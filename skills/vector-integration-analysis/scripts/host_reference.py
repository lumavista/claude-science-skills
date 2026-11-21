#!/usr/bin/env python3
"""Stage 2b -- settle the host genome, fetch its reference, annotate sites.

Every number this skill reports about the host -- the read-count baseline,
the genome overview, which gene a site lands in -- depends on knowing which
organism and which assembly the provider mapped against. Paperwork often
omits it, and a human default silently applied to a CHO line is wrong by a
fifth of the genome. So this stage decides the host from the evidence, says
how it decided, and stops when the evidence does not decide.

How the host is decided
-----------------------
1. **Declared.** The manifest's `host` block (or `--species`, `--assembly`,
   `--cell-line`) names it. A declared assembly is still checked against the
   junction tables: every host contig must exist in it and every junction
   position must fall inside its contig.
2. **Inferred.** With nothing declared, the host contig names and the largest
   junction coordinate on each (from the Stage 2 classified-junction tables,
   or `--contig NAME[:POS]`) are matched against each assembly in the
   catalogue below, under every naming convention NCBI records for it
   (RefSeq `NC_000001.11`, GenBank `CM000663.2`, UCSC `chr1`, bare `1`). A
   provider `.fai` (`host.fai`) makes this decisive, since contig lengths are
   compared too. A cell-line or species hint narrows the candidates.
3. **Asked.** If more than one assembly fits (`chr1` alone fits human, mouse
   and hamster) or none does, the `host_species` check fails and names the
   candidates. Ask the user which organism and build the line was mapped
   against and put it in the manifest's `host` block; never pick one.

The catalogue holds only the current, most-used assembly of the common
expression and research hosts. Any other NCBI assembly is reachable by
declaring its accession (`"assembly": "GCF_000001405.25"`), and an assembly
NCBI does not hold by declaring local `fai` and `annotation` files.

What it fetches and writes (beside the manifest, under `host_reference/`)
------------------------------------------------------------------------
* the NCBI assembly report (cached; source of names, lengths, genome size)
* `host.fai` -- every host sequence, named the way the junction tables name
  them, for the genome-overview figure
* `host_genes.tsv.gz` -- one row per annotated gene with its merged exons,
  built from the NCBI RefSeq GTF (or `host.annotation`), for site annotation
* `genome_size` -- the haploid primary-assembly length, for Stage 3

and records all of it under `stages.host` in the manifest. The declared
`genome_size` is never overwritten; a declared value that disagrees with the
assembly by more than 10% fails the `host_genome_size` check.

Usage
-----
    python host_reference.py --run run.json --strict
    python host_reference.py --run run.json --cell-line CHO-K1
    python host_reference.py --contig chr7:55019017 --contig chr12 --offline \\
        --cache host_reference --out-dir host_reference

Checks: `host_species`, `host_contigs`, `host_genome_size`, `host_annotation`.

Standard library only. Network access goes to ftp.ncbi.nlm.nih.gov only, and
`--offline` forbids it (everything must already be in the cache).
"""
from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import json
import os
import re
import sys
import time
import urllib.error
import urllib.request
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from via_run import (FAIL, PASS, WARN, Checks, RunError,  # noqa: E402
                     add_check_arguments, load_run)

NCBI = "https://ftp.ncbi.nlm.nih.gov/genomes/all"

# --------------------------------------------------------------------------
# the catalogue: the current, most-used assembly of each common host
# --------------------------------------------------------------------------

CATALOGUE: List[Dict] = [
    {"key": "human", "species": "Homo sapiens", "common": "human",
     "assembly": "GRCh38.p14", "accession": "GCF_000001405.40",
     "aliases": ["grch38", "grch38.p14", "hg38"],
     "cell_lines": ["hek293", "293", "293t", "293f", "293ft", "expi293",
                    "freestyle293", "hela", "jurkat", "k562", "hepg2",
                    "ht1080", "hct116", "a549", "u2os", "mcf7", "perc6",
                    "cap", "ipsc", "hipsc", "hesc", "thp1", "huvec", "rpe1",
                    "nk92"]},
    {"key": "mouse", "species": "Mus musculus", "common": "mouse",
     "assembly": "GRCm39", "accession": "GCF_000001635.27",
     "aliases": ["grcm39", "mm39"],
     "cell_lines": ["ns0", "sp20", "nih3t3", "3t3", "c2c12", "raw2647",
                    "mesc", "baf3", "l929", "hybridoma", "el4", "b16",
                    "b16f10", "ct26", "4t1", "neuro2a", "n2a"]},
    {"key": "chinese_hamster", "species": "Cricetulus griseus",
     "common": "Chinese hamster", "assembly": "CriGri-PICRH-1.0",
     "accession": "GCF_003668045.3",
     "aliases": ["crigri-picrh-1.0", "picrh", "picrh-1.0"],
     "cell_lines": ["cho", "chok1", "chos", "chodg44", "dg44", "dxb11",
                    "chodxb11", "expicho", "chozn", "chok1sv", "chok1gs",
                    "cho-gs"]},
    {"key": "rat", "species": "Rattus norvegicus", "common": "rat",
     "assembly": "mRatBN7.2", "accession": "GCF_015227675.2",
     "aliases": ["mratbn7.2", "rn7"],
     "cell_lines": ["pc12", "h9c2", "ins1", "ins1e", "rbl2h3", "yb2/0",
                    "yb20"]},
    {"key": "african_green_monkey", "species": "Chlorocebus sabaeus",
     "common": "African green monkey", "assembly": "Chlorocebus_sabeus_1.1",
     "accession": "GCF_000409795.2",
     "aliases": ["chlorocebus_sabeus_1.1", "chlsab2"],
     "cell_lines": ["vero", "veroe6", "verocc81", "cos7", "cos1", "cv1"],
     "note": "COS and CV-1 derive from Chlorocebus aethiops; the C. sabaeus "
             "assembly is the closest reference NCBI annotates."},
    {"key": "dog", "species": "Canis lupus familiaris", "common": "dog",
     "assembly": "ROS_Cfam_1.0", "accession": "GCF_014441545.1",
     "aliases": ["ros_cfam_1.0"],
     "cell_lines": ["mdck", "mdckii", "mdck2", "mdcksiat1"]},
]

#: host-side names that are legitimately not in any assembly report
EXTRA_RE = re.compile(r"(decoy|ebv|^hla[-_]|phix|lambda|^chrun_jtfh|"
                      r"_alt$|^kn707|^jtfh)", re.I)


def _norm(s: str) -> str:
    return re.sub(r"[^a-z0-9.]", "", (s or "").lower())


def _tok(s: str) -> str:
    """Cell-line / species token: case, spaces, dashes and underscores off."""
    return re.sub(r"[\s_\-]", "", (s or "").lower())


def catalogue_entry(text: str) -> Optional[Dict]:
    """The catalogue assembly named by an accession, name or synonym."""
    t = _norm(text)
    for e in CATALOGUE:
        if t in (_norm(e["accession"]), _norm(e["assembly"]),
                 *(_norm(a) for a in e["aliases"])):
            return e
    return None


def entries_for_hint(species: str = "", cell_line: str = "") -> List[Dict]:
    """Catalogue assemblies a species or cell-line hint allows (all if none)."""
    out = list(CATALOGUE)
    if species:
        s = _tok(species)
        out = [e for e in out if s and (
            s in (_tok(e["species"]), _tok(e["common"]), _tok(e["key"]))
            or s in _tok(e["species"]) + "|" + _tok(e["common"]))]
    if cell_line:
        c = _tok(cell_line)
        hit = [e for e in out if c in e["cell_lines"]]
        if not hit:  # 'CHO-K1 GS knockout' -> longest known prefix
            hit = [e for e in out
                   if any(c.startswith(x) for x in e["cell_lines"]
                          if len(x) >= 3)]
        out = hit
    return out


# --------------------------------------------------------------------------
# fetching
# --------------------------------------------------------------------------

class FetchError(RuntimeError):
    pass


def _get(url: str, dest: str, offline: bool, tries: int = 3) -> str:
    """Download `url` to `dest` unless it is already there."""
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return dest
    if offline:
        raise FetchError("%s is not in the cache and --offline forbids the "
                         "download. Fetch it by hand from %s, or rerun "
                         "without --offline." % (os.path.basename(dest), url))
    os.makedirs(os.path.dirname(dest) or ".", exist_ok=True)
    tmp = dest + ".part"
    last = None
    for attempt in range(tries):
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": "vector-integration-analysis"})
            with urllib.request.urlopen(req, timeout=120) as r, \
                    open(tmp, "wb") as fh:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    fh.write(chunk)
            os.replace(tmp, dest)
            return dest
        except (urllib.error.URLError, OSError) as exc:
            last = exc
            time.sleep(2 * (attempt + 1))
    if os.path.exists(tmp):
        os.remove(tmp)
    raise FetchError("could not download %s (%s)" % (url, last))


def ncbi_dir(accession: str, cache: str, offline: bool,
             known_name: str = "") -> Tuple[str, str]:
    """(directory URL, directory name) of an NCBI assembly accession."""
    m = re.match(r"^(GC[AF])_(\d{3})(\d{3})(\d{3})\.\d+$", accession)
    if not m:
        raise FetchError("%r is not an NCBI assembly accession (GCF_/GCA_ "
                         "followed by nine digits and a version)" % accession)
    parent = "%s/%s/%s/%s/%s" % (NCBI, m.group(1), m.group(2), m.group(3),
                                 m.group(4))
    if known_name:
        name = "%s_%s" % (accession, known_name)
        return "%s/%s" % (parent, name), name
    listing = _get(parent + "/", os.path.join(cache, accession + ".listing"),
                   offline)
    with open(listing, encoding="utf-8", errors="replace") as fh:
        names = re.findall(r'href="(%s_[^"/]+)/?"' % re.escape(accession),
                           fh.read())
    if not names:
        raise FetchError("NCBI has no directory for %s under %s"
                         % (accession, parent))
    return "%s/%s" % (parent, names[0]), names[0]


def ncbi_files(accession: str, cache: str, offline: bool,
               known_name: str = "") -> Dict[str, str]:
    base, name = ncbi_dir(accession, cache, offline, known_name)
    return {"report_url": "%s/%s_assembly_report.txt" % (base, name),
            "report": os.path.join(cache, name + "_assembly_report.txt"),
            "gtf_url": "%s/%s_genomic.gtf.gz" % (base, name),
            "gtf": os.path.join(cache, name + "_genomic.gtf.gz"),
            "name": name}


# --------------------------------------------------------------------------
# the assembly report: names, lengths, genome size
# --------------------------------------------------------------------------

STYLES = ("refseq", "genbank", "ucsc", "name", "chr", "bare")


class Assembly:
    """One NCBI assembly report: every sequence under every name."""

    def __init__(self, path: str):
        self.meta: Dict[str, str] = {}
        self.seqs: List[Dict[str, object]] = []
        with open(path, encoding="utf-8", errors="replace") as fh:
            for ln in fh:
                if ln.startswith("#"):
                    m = re.match(r"#\s*([^:]+):\s*(.*)", ln)
                    if m:
                        self.meta[m.group(1).strip()] = m.group(2).strip()
                    continue
                f = ln.rstrip("\n").split("\t")
                if len(f) < 10:
                    continue
                try:
                    length = int(f[8])
                except ValueError:
                    continue
                role, mol, loc = f[1], f[2], f[3]
                names = {"name": f[0], "genbank": f[4], "refseq": f[6],
                         "ucsc": f[9]}
                if role == "assembled-molecule" and mol not in ("na", ""):
                    names["bare"] = mol
                    names["chr"] = "chr" + mol
                    if loc.lower().startswith("mito"):
                        names["chr"] = "chrM"
                        names["bare"] = "MT"
                self.seqs.append({
                    "names": {k: v for k, v in names.items()
                              if v and v.lower() != "na"},
                    "length": length, "role": role, "unit": f[7],
                    "mito": loc.lower().startswith("mito")})
        # The header lists the units in a layout that varies by assembly, so
        # take the primary unit from the rows: the one the chromosomes are in.
        units = [str(s["unit"]) for s in self.seqs
                 if s["role"] == "assembled-molecule" and not s["mito"]] \
            or [str(s["unit"]) for s in self.seqs]
        self.primary_unit = max(set(units), key=units.count) if units else ""
        self.alias: Dict[str, int] = {}
        for i, s in enumerate(self.seqs):
            for v in s["names"].values():  # type: ignore[union-attr]
                self.alias.setdefault(_norm(v), i)

    def find(self, name: str) -> Optional[Dict]:
        i = self.alias.get(_norm(name))
        return None if i is None else self.seqs[i]

    def style_of(self, name: str) -> Optional[str]:
        s = self.find(name)
        if not s:
            return None
        for st in STYLES:
            if _norm(s["names"].get(st, "")) == _norm(name):  # type: ignore
                return st
        return None

    def genome_size(self) -> int:
        """Haploid primary-assembly length: placed, unlocalized and unplaced
        sequence of the primary unit, without alt loci, patches or the
        mitochondrion."""
        return sum(int(s["length"]) for s in self.seqs
                   if s["unit"] == self.primary_unit and not s["mito"]
                   and s["role"] in ("assembled-molecule",
                                     "unlocalized-scaffold",
                                     "unplaced-scaffold"))

    def organism(self) -> str:
        return self.meta.get("Organism name", "")

    def name_in(self, seq: Dict, style: str) -> str:
        n = seq["names"]  # type: ignore[index]
        return n.get(style) or n.get("refseq") or n.get("genbank") \
            or n.get("name")


# --------------------------------------------------------------------------
# evidence: host contigs from the junction tables
# --------------------------------------------------------------------------

def contigs_from_tables(paths: Iterable[str]) -> Dict[str, int]:
    """{host contig: largest host position} from classified-junction tables."""
    out: Dict[str, int] = {}
    for p in paths:
        with open(p, encoding="utf-8", newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                c = (row.get("host_contig") or "").strip()
                if not c or row.get("kind", "host-vector") not in (
                        "host-vector", "host-host"):
                    continue
                try:
                    pos = int(float(row.get("host_pos") or 0))
                except ValueError:
                    pos = 0
                out[c] = max(out.get(c, 0), pos)
    return out


def parse_contig_args(items: Sequence[str]) -> Dict[str, int]:
    out: Dict[str, int] = {}
    for it in items:
        name, _, pos = it.rpartition(":") if ":" in it else (it, "", "0")
        try:
            p = int(pos.replace(",", "")) if pos else 0
        except ValueError:
            name, p = it, 0
        out[name] = max(out.get(name, 0), p)
    return out


def read_fai_lengths(path: str) -> Dict[str, int]:
    out: Dict[str, int] = {}
    with open(path, encoding="utf-8") as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) >= 2 and f[0]:
                try:
                    out[f[0]] = int(f[1])
                except ValueError:
                    continue
    return out


def score(asm: Assembly, contigs: Dict[str, int],
          fai: Dict[str, int]) -> Dict:
    """How well the evidence fits one assembly."""
    unknown, outside, extras, styles = [], [], [], {}
    for c, pos in contigs.items():
        s = asm.find(c)
        if s is None:
            (extras if EXTRA_RE.search(c) else unknown).append(c)
            continue
        st = asm.style_of(c) or "refseq"
        styles[st] = styles.get(st, 0) + 1
        if pos > int(s["length"]):
            outside.append("%s:%d > %d" % (c, pos, int(s["length"])))
    len_match = len_mismatch = 0
    for c, L in fai.items():
        s = asm.find(c)
        if s is None:
            continue
        if int(s["length"]) == L:
            len_match += 1
        else:
            len_mismatch += 1
    n = len(contigs) - len(extras)
    found = n - len(unknown)
    fits = (n > 0 or len_match > 0) and not unknown and not outside \
        and len_mismatch == 0
    return {"found": found, "n": n, "unknown": unknown, "outside": outside,
            "extras": extras, "styles": styles, "len_match": len_match,
            "len_mismatch": len_mismatch, "fits": fits}


# --------------------------------------------------------------------------
# genes: GTF -> compact gene table -> site annotation
# --------------------------------------------------------------------------

_ATTR = re.compile(r'(\S+)\s+"([^"]*)"')
GENE_COLS = ["contig", "start", "end", "strand", "gene_id", "gene", "biotype",
             "exons"]


def _merge(iv: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    out: List[Tuple[int, int]] = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1] + 1:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def build_gene_table(gtf: str, out: str, rename) -> int:
    """GTF (gz or plain) -> one row per gene with merged exons. `rename` maps
    a GTF sequence name to the name the junction tables use."""
    genes: Dict[Tuple[str, str], Dict] = {}
    opener = gzip.open if gtf.endswith(".gz") else open
    with opener(gtf, "rt", encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("gene", "exon"):
                continue
            a = dict(_ATTR.findall(f[8]))
            gid = a.get("gene_id", "")
            if not gid:
                continue
            key = (f[0], gid)
            g = genes.get(key)
            if g is None:
                g = genes[key] = {"contig": rename(f[0]), "start": 0, "end": 0,
                                  "strand": f[6], "gene_id": gid, "gene": "",
                                  "biotype": "", "exons": []}
            s, e = int(f[3]), int(f[4])
            if f[2] == "gene":
                g["start"], g["end"] = s, e
                g["gene"] = a.get("gene_name") or a.get("gene") or gid
                g["biotype"] = a.get("gene_biotype") or a.get("gene_type", "")
            else:
                g["exons"].append((s, e))
                g["gene"] = g["gene"] or a.get("gene_name") or a.get("gene", "")
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    n = 0
    with gzip.open(out, "wt", encoding="utf-8") as fh:
        fh.write("\t".join(GENE_COLS) + "\n")
        for g in sorted(genes.values(), key=lambda x: (x["contig"], x["start"])):
            ex = _merge(g["exons"])
            if not g["start"] and ex:
                g["start"], g["end"] = ex[0][0], ex[-1][1]
            if not g["start"]:
                continue
            fh.write("\t".join([
                g["contig"], str(g["start"]), str(g["end"]), g["strand"],
                g["gene_id"], g["gene"] or g["gene_id"], g["biotype"],
                ";".join("%d-%d" % x for x in ex)]) + "\n")
            n += 1
    return n


ANNOT_COLS = ["host_region", "host_gene", "host_gene_biotype", "host_strand",
              "nearest_gene", "nearest_gene_distance", "nearest_gene_side"]

_BIOTYPE_RANK = {"protein_coding": 0, "lncRNA": 1, "lnc_RNA": 1}


class GeneIndex:
    """The gene table, queryable by host interval."""

    def __init__(self, path: str):
        self.by: Dict[str, List[Dict]] = {}
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt", encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                row["start"], row["end"] = int(row["start"]), int(row["end"])
                row["exons"] = [tuple(int(x) for x in iv.split("-"))
                                for iv in row["exons"].split(";") if iv]
                self.by.setdefault(row["contig"], []).append(row)
        self.starts: Dict[str, List[int]] = {}
        self.maxlen: Dict[str, int] = {}
        for c, gs in self.by.items():
            gs.sort(key=lambda g: g["start"])
            self.starts[c] = [g["start"] for g in gs]
            self.maxlen[c] = max(g["end"] - g["start"] + 1 for g in gs)

    def contigs(self) -> List[str]:
        return list(self.by)

    def annotate(self, contig: str, start: int, end: int) -> Dict[str, object]:
        """What the host interval [start, end] (1-based) falls in."""
        out: Dict[str, object] = {c: "" for c in ANNOT_COLS}
        gs, st = self.by.get(contig), self.starts.get(contig)
        if not gs:
            out["host_region"] = "no-annotation"
            return out
        lo = bisect.bisect_left(st, start - self.maxlen[contig])
        hi = bisect.bisect_right(st, end)
        hits = [g for g in gs[lo:hi] if g["start"] <= end and g["end"] >= start]
        if hits:
            def exonic(g):
                return any(s <= end and e >= start for s, e in g["exons"])
            hits.sort(key=lambda g: (not exonic(g),
                                     _BIOTYPE_RANK.get(g["biotype"], 2),
                                     g["start"]))
            g = hits[0]
            out.update({"host_region": "exonic" if exonic(g) else "intronic",
                        "host_gene": ",".join(dict.fromkeys(
                            h["gene"] for h in hits)),
                        "host_gene_biotype": g["biotype"],
                        "host_strand": g["strand"]})
            return out
        out["host_region"] = "intergenic"
        best = None
        i = bisect.bisect_left(st, start)
        cand = gs[max(0, i - 50):i + 50]
        for g in cand:
            d = g["start"] - end if g["start"] > end else start - g["end"]
            if d > 0 and (best is None or d < best[0]):
                best = (d, g)
        if best:
            d, g = best
            before = g["start"] > end
            side = "upstream" if before == (g["strand"] != "-") else \
                "downstream"
            out.update({"nearest_gene": g["gene"], "nearest_gene_distance": d,
                        "nearest_gene_side": side})
        return out


def annotate_rows(rows: Sequence[Dict], genes_path: str) -> List[Dict]:
    """Add ANNOT_COLS to each called-site row (keys contig, start, end)."""
    idx = GeneIndex(genes_path)
    out = []
    for r in rows:
        a = idx.annotate(str(r["contig"]), int(r["start"]), int(r["end"]))
        out.append(dict(r, **a))
    return out


# --------------------------------------------------------------------------
# deciding the host
# --------------------------------------------------------------------------

def decide(declared: Dict, contigs: Dict[str, int], fai: Dict[str, int],
           cache: str, offline: bool, checks: Checks
           ) -> Tuple[Optional[Dict], Optional[Assembly], Dict, str]:
    """(entry, assembly, score, how) -- or entry None when undecided.

    `entry` is a catalogue-shaped dict (accession, assembly, species ...).
    """
    asm_text = declared.get("assembly", "")
    species = declared.get("species", "")
    cell = declared.get("cell_line", "")

    # 1. an explicitly named assembly
    if asm_text:
        e = catalogue_entry(asm_text)
        if e is None and re.match(r"^GC[AF]_\d{9}\.\d+$", asm_text):
            e = {"key": asm_text, "accession": asm_text, "assembly": "",
                 "species": species, "common": species}
        if e is None:
            raise RunError(
                "host assembly %r is not in the catalogue (%s). Give its NCBI "
                "accession instead (e.g. GCF_000001405.25 for GRCh37), or "
                "supply host.fai and host.annotation files."
                % (asm_text, ", ".join(x["assembly"] for x in CATALOGUE)))
        asm = load_assembly(e, cache, offline)
        if not e.get("assembly"):
            e = dict(e, assembly=asm.meta.get("Assembly name", ""),
                     species=asm.organism() or species)
        return e, asm, score(asm, contigs, fai), "declared"

    # 2. inference over the catalogue, narrowed by any hint
    pool = entries_for_hint(species, cell)
    if not pool:
        checks.add("host_species", FAIL,
                   "the hint (species %r, cell line %r) matches no catalogue "
                   "host (%s). Ask the user which organism and assembly the "
                   "reads were mapped against, and set host.assembly to its "
                   "NCBI accession." % (species, cell, ", ".join(
                       "%s %s" % (x["common"], x["assembly"])
                       for x in CATALOGUE)))
        return None, None, {}, ""
    if not contigs and not fai:
        if len(pool) == 1:
            e = pool[0]
            asm = load_assembly(e, cache, offline)
            return e, asm, score(asm, contigs, fai), "hint"
        checks.add("host_species", FAIL,
                   "no host evidence: no host contigs in the junction tables "
                   "and no host.fai, and no species or cell line declared. "
                   "Ask the user which organism the line is and which "
                   "assembly the reads were mapped against.")
        return None, None, {}, ""

    scored, unreachable = [], []
    for e in pool:
        try:
            asm = load_assembly(e, cache, offline)
        except FetchError as exc:
            unreachable.append("%s %s (%s)" % (e["common"], e["assembly"], exc))
            continue
        scored.append((e, asm, score(asm, contigs, fai)))
    if unreachable:
        if not scored:
            raise FetchError("no catalogue assembly report could be read: "
                             + "; ".join(unreachable))
        checks.add("host_catalogue", WARN,
                   "left out of the inference, report unavailable: "
                   + "; ".join(unreachable))
    fits = [x for x in scored if x[2]["fits"]]
    if fai and len(fits) > 1:  # lengths decide between same-named assemblies
        top = max(x[2]["len_match"] for x in fits)
        fits = [x for x in fits if x[2]["len_match"] == top]
    if len(fits) == 1:
        e, asm, sc = fits[0]
        return e, asm, sc, "hint" if (species or cell) else "inferred"

    def line(x):
        e, _, sc = x
        why = ("fits" if sc["fits"] else
               "%d/%d contigs found%s" % (
                   sc["found"], sc["n"],
                   ", positions past contig end" if sc["outside"] else ""))
        return "%s %s (%s): %s" % (e["common"], e["assembly"],
                                   e["accession"], why)
    if fits:
        msg = ("the host contigs fit %d catalogue assemblies equally well, so "
               "the host cannot be inferred: %s. Ask the user which organism "
               "the line is (or its cell line) and set host.species or "
               "host.cell_line; a provider .fai in host.fai would also "
               "decide it." % (len(fits), "; ".join(line(x) for x in fits)))
    else:
        msg = ("no catalogue assembly fits the host contigs (%s): %s. The "
               "line may be another organism or the reads another build. Ask "
               "the user, then set host.assembly to the NCBI accession, or "
               "supply host.fai and host.annotation."
               % (", ".join(sorted(contigs)[:8]) or "from host.fai",
                  "; ".join(line(x) for x in scored)))
    checks.add("host_species", FAIL, msg)
    return None, None, {"candidates": [
        {"species": x[0]["species"], "assembly": x[0]["assembly"],
         "accession": x[0]["accession"], "fits": x[2]["fits"]}
        for x in (fits or scored)]}, ""


def load_assembly(e: Dict, cache: str, offline: bool) -> Assembly:
    known = e["assembly"] if e.get("assembly") else ""
    files = ncbi_files(e["accession"], cache, offline, known)
    _get(files["report_url"], files["report"], offline)
    return Assembly(files["report"])


def naming_style(sc: Dict) -> str:
    st = sc.get("styles") or {}
    return max(st, key=lambda k: st[k]) if st else "ucsc"


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Settle the host organism and assembly, fetch its "
                    "sequence names, lengths and gene annotation, and record "
                    "them for Stages 3 and 4.",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--species", default="",
                    help="host species, e.g. 'Cricetulus griseus' or 'mouse'")
    ap.add_argument("--assembly", default="",
                    help="assembly name, synonym or NCBI accession")
    ap.add_argument("--cell-line", default="", help="e.g. CHO-K1, HEK293")
    ap.add_argument("--fai", default="",
                    help="the provider's host .fai: names AND lengths, the "
                         "strongest evidence for which build was used")
    ap.add_argument("--annotation", default="",
                    help="local GTF to use instead of the NCBI RefSeq GTF")
    ap.add_argument("--table", action="append", default=[],
                    help="classified-junction table (Stage 2 output); "
                         "default: every table Stage 2 recorded in the run")
    ap.add_argument("--contig", action="append", default=[],
                    metavar="NAME[:POS]", help="a host contig as evidence")
    ap.add_argument("--cache", default="",
                    help="download cache (default: <out-dir>)")
    ap.add_argument("--out-dir", default="",
                    help="output directory (default: host_reference/ beside "
                         "the run manifest)")
    ap.add_argument("--offline", action="store_true",
                    help="never download; use only what is in the cache")
    ap.add_argument("--no-annotation", action="store_true",
                    help="skip the gene annotation (and its download)")
    ap.add_argument("--json", default="", help="write the result here too")
    ap.add_argument("--list", action="store_true",
                    help="print the catalogue and exit")
    add_check_arguments(ap)
    return ap


def main(argv: Optional[Sequence[str]] = None) -> int:
    a = build_parser().parse_args(argv)
    if a.list:
        for e in CATALOGUE:
            print("%-22s %-24s %-18s %s" % (e["common"], e["assembly"],
                                            e["accession"],
                                            ", ".join(e["cell_lines"][:6])))
        return 0
    checks = Checks("host_reference", a.accept)
    try:
        run = load_run(a.run)
        declared: Dict = dict((run.fact("host") or {}) if run else {})
        for k in ("species", "assembly", "cell_line", "fai", "annotation"):
            v = getattr(a, k)
            if v:
                declared[k] = v
        res = (lambda p: run.resolve(p)) if run else (lambda p: p)
        out_dir = a.out_dir or (os.path.join(run.base, "host_reference")
                                if run else "host_reference")
        cache = a.cache or out_dir
        tables = list(a.table)
        if run and not tables:
            tables = [run.resolve(r["junctions"])
                      for r in run.stage("classify").values()
                      if r.get("junctions")]
        contigs = contigs_from_tables(tables)
        contigs.update(parse_contig_args(a.contig))
        fai_path = res(declared["fai"]) if declared.get("fai") else ""
        fai = read_fai_lengths(fai_path) if fai_path else {}
        if run is None and not (contigs or fai or declared):
            raise RunError("need --run, or host evidence (--contig/--table/"
                           "--fai) or a declaration (--species/--assembly/"
                           "--cell-line)")
    except (RunError, OSError) as exc:
        print("error: %s" % exc, file=sys.stderr)
        return 2

    os.makedirs(out_dir, exist_ok=True)
    local_only = bool(fai_path and declared.get("annotation")
                      and not declared.get("assembly"))
    try:
        if local_only:
            entry = {"species": declared.get("species", ""),
                     "assembly": "local files",
                     "accession": "", "common": declared.get("species", "")}
            asm, sc, how = None, {"fits": True, "styles": {}}, "local files"
        else:
            entry, asm, sc, how = decide(declared, contigs, fai, cache,
                                         a.offline, checks)
    except (RunError, FetchError) as exc:
        print("error: %s" % exc, file=sys.stderr)
        return 2

    result: Dict = {"declared": {k: v for k, v in declared.items()},
                    "evidence": {"tables": [os.path.relpath(t, out_dir)
                                            for t in tables],
                                 "contigs": len(contigs),
                                 "fai": bool(fai)}}
    if entry is None:
        result.update(sc)
        return finish(a, run, checks, result, out_dir)

    style = naming_style(sc)
    result.update({"species": entry.get("species", ""),
                   "assembly": entry.get("assembly", ""),
                   "accession": entry.get("accession", ""),
                   "decided_by": how, "naming": style})
    if entry.get("note"):
        result["note"] = entry["note"]
    checks.add("host_species", PASS,
               "%s %s (%s), %s" % (result["species"], result["assembly"],
                                   result["accession"] or "local files", how))

    # ---- contigs against the assembly
    if asm is not None and (contigs or fai):
        bad = sc["unknown"] + sc["outside"]
        if sc["len_mismatch"]:
            bad.append("%d .fai contig length(s) differ from %s"
                       % (sc["len_mismatch"], result["assembly"]))
        if bad:
            checks.add("host_contigs", FAIL,
                       "the junction tables do not fit %s: %s. The reads were "
                       "probably mapped to another build; confirm the build "
                       "with the provider and set host.assembly."
                       % (result["assembly"], "; ".join(bad[:10])))
        elif sc["extras"]:
            checks.add("host_contigs", WARN,
                       "%d host contig(s) are decoy/viral/HLA sequences not "
                       "in %s: %s" % (len(sc["extras"]), result["assembly"],
                                      ", ".join(sc["extras"][:6])))
        else:
            checks.add("host_contigs", PASS,
                       "all %d host contigs found in %s, named %s-style"
                       % (sc["n"], result["assembly"], style))

    # ---- the .fai and the genome size
    fai_out = os.path.join(out_dir, "host.fai")
    if asm is not None:
        with open(fai_out, "w", encoding="utf-8") as fh:
            for s in asm.seqs:
                if s["role"] in ("alt-scaffold", "fix-patch", "novel-patch"):
                    continue
                fh.write("%s\t%d\t0\t0\t0\n" % (asm.name_in(s, style),
                                                int(s["length"])))
        gsize = asm.genome_size()
    else:
        with open(fai_out, "w", encoding="utf-8") as fh:
            for n, L in fai.items():
                fh.write("%s\t%d\t0\t0\t0\n" % (n, L))
        gsize = sum(L for n, L in fai.items() if not EXTRA_RE.search(n))
    result["fai"] = fai_out
    result["genome_size"] = gsize
    declared_size = run.fact("genome_size") if run else None
    if declared_size:
        dev = abs(float(declared_size) - gsize) / gsize
        status = FAIL if dev > 0.10 else PASS
        checks.add("host_genome_size", status,
                   "declared genome_size %s vs %s primary assembly %d (%.0f%% "
                   "apart)%s" % (declared_size, result["assembly"], gsize,
                                 dev * 100,
                                 ": correct genome_size in the manifest or "
                                 "delete it so Stage 3 uses the assembly's"
                                 if status == FAIL else ""))
    else:
        checks.add("host_genome_size", PASS,
                   "genome_size %d from the %s primary assembly"
                   % (gsize, result["assembly"]))

    # ---- annotation
    if a.no_annotation:
        checks.add("host_annotation", WARN,
                   "annotation skipped (--no-annotation): sites will carry "
                   "no gene context")
    else:
        try:
            if declared.get("annotation"):
                gtf = res(declared["annotation"])
            else:
                files = ncbi_files(result["accession"], cache, a.offline,
                                   result["assembly"])
                gtf = _get(files["gtf_url"], files["gtf"], a.offline)
            rename = ((lambda n: asm.name_in(asm.find(n), style)
                       if asm.find(n) else n) if asm is not None
                      else (lambda n: n))
            genes_out = os.path.join(out_dir, "host_genes.tsv.gz")
            n = build_gene_table(gtf, genes_out, rename)
            result["genes"] = genes_out
            result["annotation_source"] = gtf
            result["n_genes"] = n
            checks.add("host_annotation", PASS if n else WARN,
                       "%d genes from %s" % (n, os.path.basename(gtf)))
        except (FetchError, OSError) as exc:
            checks.add("host_annotation", WARN,
                       "no gene annotation (%s): sites will carry no gene "
                       "context" % exc)
    return finish(a, run, checks, result, out_dir)


def finish(a, run, checks: Checks, result: Dict, out_dir: str) -> int:
    result["checks"] = checks.as_list()
    if run is not None:
        rec = {k: v for k, v in result.items() if k != "evidence"}
        for k in ("fai", "genes", "annotation_source"):
            if rec.get(k):
                rec[k] = run.relative(rec[k])
        run.record("host", rec)
    path = a.json or os.path.join(out_dir, "host.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(result, fh, indent=2)
        fh.write("\n")
    if result.get("assembly"):
        print("host: %s %s (%s), decided by %s; genome_size %s; %s genes"
              % (result.get("species"), result.get("assembly"),
                 result.get("accession") or "local", result.get("decided_by"),
                 result.get("genome_size"), result.get("n_genes", "no")),
              file=sys.stderr)
    return checks.finish(a.strict)


if __name__ == "__main__":
    raise SystemExit(main())
