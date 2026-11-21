#!/usr/bin/env python3
"""Assemble the final report from the pipeline outputs and the skill's figures.

This is the last stage, and the authoritative report producer. It fills
assets/report_template.md (or a copy into which the interpretive prose has
been written) deterministically:

  [[FIG:<id>]]     the SVG figures make_figures.py drew -- the coverage
                   tracks, genome_overview, copy_number, architecture
  [[TABLE:<id>]]   tables built from the stage outputs: the called sites,
                   copy numbers, site bounds, site accounting, checks, ...
  [[VAL:<id>]]     single values: the site range, the unresolved count, the
                   shared-block closure deficit, the host, the modality, ...

`<construct A>`, `<construct B>`, ... inside a placeholder name the constructs
in report order: A is the primary mapping, the rest follow the identity map.
`--list-placeholders` prints every id.

Outputs, all beside the manifest by default:

  report.html        the deliverable: one self-contained file, the SVGs
                     embedded inline, no external file, script or font;
                     prints to PDF from any browser
  report_filled.md   the same report as Markdown, linking the SVG files
  report.docx        OPTIONAL, only when an SVG rasteriser (cairosvg or
                     rsvg-convert) is present: the same report with the
                     skill's figures rasterised. Without one it is skipped,
                     with the reason, and the HTML stands
  report.json        what was built, what was embedded, and the checks

Every figure in every output is the skill's own: an SVG listed in
figures/figures.json and carrying make_figures.py's figure id. This script
draws nothing. An image in the prose that the skill did not draw fails the
`foreign_figure` check and is removed; a skill figure no placeholder asked for
is appended, never dropped.

    python scripts/build_report.py --run run.json --report report.md --strict

Python 3.8+. Standard library only; the optional docx path additionally uses a
rasteriser if one is installed.
"""
from __future__ import annotations

import argparse
import csv
import html
import json
import os
import re
import shlex
import shutil
import struct
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
import zipfile
from typing import Dict, List, Optional, Sequence, Tuple

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from via_run import (ACCEPTED, FAIL, PASS, WARN, Checks,  # noqa: E402
                     InputError, RunError, add_check_arguments, load_run)

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_TEMPLATE = os.path.join(HERE, "..", "assets", "report_template.md")

PLACEHOLDER = re.compile(r"\[\[(FIG|TABLE|VAL):([^\]]+)\]\]")
ANY_PLACEHOLDER = re.compile(r"\[\[[^\[\]\n]*\]\]")
LETTER = re.compile(r"<\s*construct\s+([A-Za-z])\s*>", re.I)
COMMENT = re.compile(r"<!--.*?-->", re.S)
FIG_TOKEN = re.compile(r"\x00FIG(\d+)\x00")
# an image the prose brought with it, in any form a renderer would show
FOREIGN_IMAGE = re.compile(
    r"!\[[^\]]*\]\([^)]*\)|<\s*img\b[^>]*>|<\s*svg\b.*?</\s*svg\s*>|<\s*svg\b[^>]*>"
    r"|<\s*(?:object|embed|iframe|picture)\b[^>]*>", re.I | re.S)
# template text the prose still has to replace: <line name>, <N>, <...>
TEMPLATE_MARK = re.compile(r"<(?!!--)[^<>]{1,800}>")

STAGE_ORDER = ("identity", "classify", "host", "copy_number", "reconcile")


class NotAvailable(Exception):
    """A placeholder whose input is missing; the reason goes in its slot."""


# --------------------------------------------------------------------------
# formatting
# --------------------------------------------------------------------------

def num(x, nd: int = 3) -> str:
    """A script output as printed: integers bare, floats to `nd` places with
    trailing zeros dropped. Presentation only; never a recomputation."""
    if x is None or x == "":
        return "—"
    if isinstance(x, bool):
        return str(x).lower()
    if isinstance(x, (int, float)):
        if float(x).is_integer():
            return str(int(x))
        return ("%.*f" % (nd, x)).rstrip("0").rstrip(".")
    return str(x)


def cell(s) -> str:
    """One Markdown table cell: no newline, no unescaped pipe."""
    return str(s).replace("\n", " ").replace("|", "\\|").strip() or "—"


def md_table(header: Sequence[str], rows: Sequence[Sequence]) -> str:
    out = ["| " + " | ".join(cell(h) if h else " " for h in header) + " |",
           "|" + "|".join("---" for _ in header) + "|"]
    for r in rows:
        out.append("| " + " | ".join(cell(c) for c in r) + " |")
    return "\n".join(out)


def plural(n: int, word: str, many: Optional[str] = None) -> str:
    return "%d %s" % (n, word if n == 1 else (many or word + "s"))


def locus_name(r: Dict) -> str:
    s, e = r.get("start", ""), r.get("end", "")
    return "%s:%s" % (r.get("contig", ""), s if not e or e == s else "%s-%s" % (s, e))


# --------------------------------------------------------------------------
# inputs
# --------------------------------------------------------------------------

def read_json(path: Optional[str], what: str) -> Optional[Dict]:
    if not path:
        return None
    if not os.path.exists(path):
        return None
    try:
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, ValueError) as e:
        raise InputError("%s %s unreadable (%s)" % (what, path, e))


def read_tsv(path: Optional[str], what: str) -> Optional[List[Dict[str, str]]]:
    if not path or not os.path.exists(path):
        return None
    try:
        with open(path, encoding="utf-8", newline="") as fh:
            return list(csv.DictReader(fh, delimiter="\t"))
    except (OSError, csv.Error) as e:
        raise InputError("%s %s unreadable (%s)" % (what, path, e))


class Inputs:
    """Everything the report draws on, read once."""

    def __init__(self, a: argparse.Namespace):
        self.run = None
        if a.run:
            try:
                self.run = load_run(a.run)
            except RunError as e:
                raise InputError("run manifest: %s" % e)
        run = self.run
        self.paths = self._input_paths(a)
        self.base = run.base if run else os.getcwd()
        self.fig_dir = os.path.abspath(a.figures or os.path.join(self.base, "figures"))

        self.ident = read_json(self.paths["identity"], "identity map")
        self.cn = read_json(self.paths["copy_number"], "copy-number JSON")
        self.summary = read_json(self.paths["site_summary"], "site summary")
        self.called = read_tsv(self.paths["called"], "called-site table")
        if a.host:
            self.host = read_json(a.host, "host record")
        else:
            self.host = dict(run.stage("host")) if run and run.stage("host") else None
        self.figures = read_json(os.path.join(self.fig_dir, "figures.json"),
                                 "figures manifest")
        self.manifest_facts = dict(run.data) if run else {}
        #: report order: construct A, B, C, ...
        self.constructs = self._construct_order()

    def _input_paths(self, a: argparse.Namespace) -> Dict[str, Optional[str]]:
        """Each input: the explicit flag, else what its stage recorded."""
        run = self.run

        def st(stage: str, key: str) -> Optional[str]:
            if run is None or not run.stage(stage).get(key):
                return None
            return run.resolve(run.stage(stage)[key])

        return {"identity": a.identity or (run.identity_path() if run else None),
                "copy_number": a.copy_number or st("copy_number", "json"),
                "called": a.called or st("reconcile", "called"),
                "site_summary": a.site_summary or st("reconcile", "summary"),
                "host": a.host}

    def _construct_order(self) -> List[str]:
        """Identity-map order, then any construct only a later stage names;
        the primary mapping first, so it is construct A."""
        order: List[str] = list((self.ident or {}).get("construct_order")
                                or (self.ident or {}).get("constructs", {}).keys())
        later = list((self.summary or {}).get("mappings") or [])
        for s in ((self.cn or {}).get("samples") or {}).values():
            later += list(s.get("copies_per_haploid") or {})
        order += [c for c in dict.fromkeys(later) if c not in order]
        primary = (self.manifest_facts.get("primary_mapping")
                   or (self.cn or {}).get("primary_mapping"))
        if primary in order:
            order.remove(primary)
            order.insert(0, primary)
        return order

    def ordered(self, names) -> List[str]:
        """`names` in report order; any the report order lacks go last."""
        names = list(names)
        return [c for c in self.constructs if c in names] + \
            [c for c in names if c not in self.constructs]

    # ---- each stage's output, or why it is not there
    def site_summary(self) -> Dict:
        if not self.summary:
            raise NotAvailable("no site summary (reconcile_sites.py has not run)")
        return self.summary

    def identity(self) -> Dict:
        if not self.ident:
            raise NotAvailable("no identity map (construct_identity.py has not run)")
        return self.ident

    def classify(self) -> Dict[str, Dict]:
        """Stage 2's per-mapping records from the manifest; {} without one."""
        return (self.run.stage("classify") if self.run else {}) or {}

    # ---- derived lookups
    def samples(self) -> List[str]:
        names = list((self.manifest_facts.get("samples") or {}).keys())
        for s in ((self.cn or {}).get("samples") or {}):
            if s not in names:
                names.append(s)
        for k in (self.called[0].keys() if self.called else []):
            if k.startswith("support_") and k[8:] not in names:
                names.append(k[8:])
        return names

    def per_construct(self) -> Dict[str, Dict]:
        return self.site_summary().get("per_construct") or {}

    def called_rows(self) -> List[Dict[str, str]]:
        if self.called is None:
            raise NotAvailable("no called_sites.tsv (reconcile_sites.py has not run)")
        return self.called

    def cn_samples(self) -> Dict[str, Dict]:
        if not self.cn:
            raise NotAvailable("no copy-number JSON (copy_number.py has not run)")
        return self.cn.get("samples") or {}

    def stage_checks(self) -> List[Tuple[str, Dict]]:
        """(stage label, check) for every check every stage recorded."""
        out: List[Tuple[str, Dict]] = []
        if self.run is not None:
            stages = self.run.data.get("stages") or {}
            for name in STAGE_ORDER:
                rec = stages.get(name) or {}
                if name == "classify":
                    for m, r in rec.items():
                        out += [("classify (%s)" % m, c) for c in r.get("checks") or []]
                else:
                    out += [(name, c) for c in rec.get("checks") or []]
            return out
        for name, d in (("identity", self.ident), ("host", self.host),
                        ("copy_number", self.cn), ("reconcile", self.summary)):
            out += [(name, c) for c in (d or {}).get("checks") or []]
        return out


# --------------------------------------------------------------------------
# construct letters
# --------------------------------------------------------------------------

def resolve_letters(key: str, constructs: Sequence[str]) -> str:
    def sub(m):
        i = ord(m.group(1).upper()) - ord("A")
        if i >= len(constructs):
            raise NotAvailable(
                "the line carries %s (%s), so there is no construct %s"
                % (plural(len(constructs), "construct"), ", ".join(constructs) or "none",
                   m.group(1).upper()))
        return constructs[i]
    return LETTER.sub(sub, key)


# --------------------------------------------------------------------------
# [[TABLE:...]]
# --------------------------------------------------------------------------

def t_summary(I: Inputs, arg: str) -> str:
    pc = I.per_construct()
    cons = I.ordered(pc)
    rows: List[List[str]] = []
    rows.append(["Integration sites (bounded range)"] + [pc[c].get("site_range", "—") for c in cons])
    rows.append(["certain / vector-end / copy-ceiling"] + [
        "%s / %s / %s" % (pc[c].get("certain", 0), pc[c].get("vector_end", 0),
                          pc[c].get("copy_ceiling", 0)) for c in cons])
    rows.append(["Unresolved loci naming it a candidate"] + [
        num(pc[c].get("unresolved_candidate_for")) for c in cons])
    if I.called is not None:
        both = []
        for c in cons:
            mine = [r for r in I.called if r.get("construct") == c]
            n2 = sum(1 for r in mine if (r.get("flanks") or "").strip() in ("2", "both"))
            both.append("%d of %d" % (n2, len(mine)))
        rows.append(["Both flanks recovered"] + both)
    for s, sd in ((I.cn or {}).get("samples") or {}).items():
        cph = sd.get("copies_per_haploid") or {}
        rows.append(["Copies per haploid genome, %s" % s] + [num(cph.get(c)) for c in cons])
    rows.append(["Copy-number ceiling (sites)"] + [num(pc[c].get("ceiling")) for c in cons])
    acc = {r.get("construct"): r for r in (I.summary or {}).get("site_accounting") or []}
    if acc:
        rows.append(["Do the called sites account for the DNA?"] + [
            (acc.get(c) or {}).get("verdict", "—") for c in cons])
    return md_table([""] + cons, rows)


def _called_table(I: Inputs, rows: List[Dict[str, str]], construct: Optional[str]) -> str:
    samples = [k[8:] for k in (rows[0].keys() if rows else
                               (I.called[0].keys() if I.called else []))
               if k.startswith("support_")]
    host = bool(rows) and "host_region" in rows[0]
    header = ["#", "Locus", "Construct", "Tier", "Evidence", "Vector position",
              "Footprint", "Flanks"] + ["Support %s" % s for s in samples] + ["Total"]
    if host:
        header += ["Host region", "Gene / nearest gene"]
    header += ["Reason"]
    if not rows:
        return "*No called locus%s.*" % (" is attributed to or names %s" % construct
                                         if construct else "")
    return md_table(header, [_called_row(i, r, samples, host) for i, r in enumerate(rows, 1)])


def _called_row(i: int, r: Dict[str, str], samples: List[str], host: bool) -> List[str]:
    cons = r.get("construct") or ("— (candidates: %s)" % r.get("candidates", "")
                                 if r.get("confidence") == "unresolved" else "—")
    line = [str(i), locus_name(r), cons, r.get("confidence", ""),
            r.get("evidence", ""), (r.get("vector_positions") or "").replace(";", "; "),
            r.get("footprint") or "—", r.get("flanks", "")]
    line += [r.get("support_" + s, "") for s in samples] + [r.get("total_support", "")]
    if host:
        gene = r.get("host_gene") or ""
        if not gene and r.get("nearest_gene"):
            gene = "near %s (%s bp %s)" % (r["nearest_gene"], r.get("nearest_gene_distance", "?"),
                                           r.get("nearest_gene_side", ""))
        line += [r.get("host_region") or "—", gene or "—"]
    return line + [r.get("note", "")]


def t_called_sites(I: Inputs, arg: str) -> str:
    rows = I.called_rows()
    if arg:
        rows = [r for r in rows if r.get("construct") == arg
                or (r.get("confidence") == "unresolved"
                    and arg in (r.get("candidates") or "").split(","))]
        return _called_table(I, rows, arg)
    return _called_table(I, rows, None)


def t_unresolved(I: Inputs, arg: str) -> str:
    un = I.site_summary().get("unresolved") or []
    if not un:
        return "*No locus is unresolved: every called locus is attributed.*"
    return md_table(["#", "Locus", "Candidates", "Total support", "Samples", "Why unresolved"],
                    [[i, locus_name(u), ", ".join(u.get("candidates") or []),
                      u.get("total_support", ""), u.get("n_samples", ""), u.get("note", "")]
                     for i, u in enumerate(un, 1)])


def t_copy_number(I: Inputs, arg: str) -> str:
    samples = I.cn_samples()
    rows = []
    for s, sd in samples.items():
        solve = sd.get("solve") or {}
        ub = solve.get("unique_blocks") or {}
        src = solve.get("source") or {}
        base = sd.get("baseline_depth") or sd.get("genome_coverage")
        for c in I.ordered(sd.get("copies_per_haploid") or {}):
            u = ub.get(c) or {}
            spread = "—"
            if u.get("n_blocks"):
                spread = "%s block%s, %s–%s" % (u["n_blocks"], "" if u["n_blocks"] == 1 else "s",
                                               num(u.get("min")), num(u.get("max")))
                if u.get("relative_range") is not None:
                    spread += " (range %s%%)" % num(100 * u["relative_range"], 1)
            rows.append([s, c, num((sd.get("copies_per_haploid") or {}).get(c)),
                         num((sd.get("copies_per_cell") or {}).get(c)),
                         src.get(c) or solve.get("method", "—"),
                         num((solve.get("unique_bp") or {}).get(c)), spread,
                         "%s× (%s)" % (num(base), sd.get("baseline_source") or "read count")
                         if base else "—"])
    if not rows:
        raise NotAvailable("the copy-number JSON holds no per-construct copy number")
    return md_table(["Sample", "Construct", "Copies / haploid", "Copies / cell", "Source",
                     "Unique bp", "Per-block estimates", "Baseline"], rows)


def t_site_bounds(I: Inputs, arg: str) -> str:
    pc = I.per_construct()
    p = (I.summary or {}).get("parameters") or {}
    rows = [[c, num(d.get("copies_per_haploid")), num(d.get("ceiling")),
             num(d.get("site_floor")), num(d.get("attributed")),
             "%s / %s / %s" % (d.get("certain", 0), d.get("vector_end", 0),
                               d.get("copy_ceiling", 0)),
             num(d.get("unresolved_candidate_for")), d.get("site_range", "—")]
            for c, d in ((c, pc[c]) for c in I.ordered(pc))]
    t = md_table(["Construct", "Copies / haploid", "Ceiling", "Floor", "Attributed",
                  "certain / vector-end / copy-ceiling", "Unresolved candidate for",
                  "Site range"], rows)
    if p:
        t += ("\n\nCeiling = ceil(copies × ploidy %s ÷ copies per site %s); "
              "copy-number margin %s." % (num(p.get("ploidy")), num(p.get("min_copies_per_site")),
                                          num(p.get("cn_margin"))))
    return t


def t_site_accounting(I: Inputs, arg: str) -> str:
    acc = I.site_summary().get("site_accounting") or []
    if not acc:
        raise NotAvailable("the site summary carries no site accounting")
    return md_table(["Construct", "Copies / haploid", "Sites the DNA implies", "Attributed",
                     "Unaccounted", "Verdict"],
                    [[r.get("construct"), num(r.get("copies_per_haploid")),
                      num(r.get("implied_sites")), num(r.get("attributed")),
                      num(r.get("unaccounted")), r.get("verdict", "")] for r in acc])


def t_checks(I: Inputs, arg: str) -> str:
    allc = I.stage_checks()
    if not allc:
        raise NotAvailable("no stage recorded its checks (no --run manifest and no stage JSON)")
    flagged = [(s, c) for s, c in allc if c.get("status") != PASS]
    if not flagged:
        return "*Every check of every stage passed (%d checks).*" % len(allc)
    rank = {FAIL: 0, ACCEPTED: 1, WARN: 2}
    flagged.sort(key=lambda sc: rank.get(sc[1].get("status"), 3))
    return md_table(["Stage", "Check", "Status", "Message"],
                    [[s, c.get("id"), c.get("status", "").upper(), c.get("message", "")]
                     for s, c in flagged])


def t_caveats(I: Inputs, arg: str) -> str:
    rows = []
    for c in (I.cn or {}).get("caveats") or []:
        rows.append(["copy_number.py", c])
    for s, sd in ((I.cn or {}).get("samples") or {}).items():
        for n in (sd.get("architecture") or {}).get("notes") or []:
            rows.append(["copy_number.py, %s architecture" % s, n])
    for c in (I.summary or {}).get("caveats") or []:
        rows.append(["reconcile_sites.py", c])
    if I.cn is None and I.summary is None:
        raise NotAvailable("neither the copy-number JSON nor the site summary is available")
    if not rows:
        return "*Neither copy_number.py nor reconcile_sites.py raised a caveat.*"
    return md_table(["From", "Caveat"], rows)


def t_constructs(I: Inputs, arg: str) -> str:
    ends = I.manifest_facts.get("vector_ends") or {}
    cons = I.identity().get("constructs") or {}
    rows = []
    for c in I.ordered(cons):
        d = cons[c]
        e = ends.get(c)
        rows.append([c, num(d.get("length")), "circular" if d.get("circular") else "linear",
                     "%s (%s%%)" % (num(d.get("unique_bp")), num(d.get("attributable_pct"))),
                     ", ".join(str(x) for x in e) if isinstance(e, list) else (e or "none declared")])
    return md_table(["Construct", "Length (bp)", "Reference", "Unique bp (attributable)",
                     "Declared vector ends"], rows)


def t_identity(I: Inputs, arg: str) -> str:
    cons = I.identity().get("constructs") or {}
    rows = []
    for c in I.ordered(cons):
        d = cons[c]
        uniq = ["%s-%s" % (b["start"], b["end"]) for b in d.get("blocks") or [] if b.get("unique")]
        shared = ["%s: %s bp" % (" + ".join(m.get("membership") or []), num(m.get("bp")))
                  for m in d.get("membership_classes") or [] if len(m.get("membership") or []) > 1]
        bnd = ["%s (%s clip)" % (b.get("pos"), ", ".join(b.get("clipping_constructs") or []))
               for b in d.get("boundaries") or []]
        rows.append([c, ", ".join(uniq) or "none", "; ".join(shared) or "none",
                     ", ".join(bnd) or "none"])
    return md_table(["Construct", "Construct-unique blocks", "Shared sequence",
                     "Predicted artefact boundaries (who clips)"], rows)


def t_fingerprint(I: Inputs, arg: str) -> str:
    rec = I.classify()
    if not rec:
        raise NotAvailable("no Stage 2 record in the run manifest (classify_junctions.py --run)")
    rows = []
    for m, r in rec.items():
        chk = {c.get("id"): c for c in r.get("checks") or []}
        fp = chk.get("modality_fingerprint") or {}
        rows.append([m, r.get("modality", "—"), r.get("fingerprint") or fp.get("message", "—"),
                     (fp.get("status") or "—").upper(), num(r.get("detected_loci")),
                     num(r.get("loci_both_flanks"))])
    return md_table(["Mapping", "Declared modality", "Junction fingerprint", "Check",
                     "Loci (Stage 2)", "Both flanks"], rows)


def t_validation_loci(I: Inputs, arg: str) -> str:
    rows = I.called_rows()
    un = [r for r in rows if r.get("confidence") == "unresolved"]
    inf = sorted([r for r in rows if r.get("confidence") == "inferred"],
                 key=lambda r: int(r.get("total_support") or 0))
    cert = sorted([r for r in rows if r.get("confidence") == "certain"],
                  key=lambda r: -int(r.get("total_support") or 0))[:1]
    out = []
    for r in un:
        out.append([locus_name(r), "candidates: " + (r.get("candidates") or ""), "unresolved",
                    r.get("evidence", ""), "decides which construct, if any, sits here"])
    for r in inf:
        out.append([locus_name(r), r.get("construct", ""), "inferred", r.get("evidence", ""),
                    "confirms an attribution made by argument (%s)" % (r.get("note") or "")])
    for r in cert:
        out.append([locus_name(r), r.get("construct", ""), "certain", r.get("evidence", ""),
                    "positive control: highest-support certain locus"])
    if not out:
        return "*No called locus to validate.*"
    return md_table(["#", "Locus", "Construct", "Tier", "Evidence", "Junction PCR tests"],
                    [[i] + r for i, r in enumerate(out, 1)])


TABLES = {
    "summary": (t_summary, "per construct: site range, tiers, flanks, copies, ceiling, accounting"),
    "called_sites": (t_called_sites, "called_sites.tsv; `called_sites:<construct A>` for one construct"),
    "unresolved": (t_unresolved, "unresolved loci, candidates named"),
    "copy_number": (t_copy_number, "copies per haploid/cell per sample, source, per-block spread, baseline"),
    "site_bounds": (t_site_bounds, "floor, attributed, ceiling, unresolved per construct"),
    "site_accounting": (t_site_accounting, "do the called sites account for the DNA?"),
    "checks": (t_checks, "every non-passing check of every stage"),
    "caveats": (t_caveats, "every caveat from copy_number.py and reconcile_sites.py"),
    "constructs": (t_constructs, "length, unique bp, declared ends per construct"),
    "identity": (t_identity, "unique blocks, shared sequence, predicted artefact boundaries"),
    "fingerprint": (t_fingerprint, "Stage 2 junction fingerprint per mapping"),
    "validation_loci": (t_validation_loci, "junction-PCR priorities: unresolved, inferred, a control"),
}


# --------------------------------------------------------------------------
# [[VAL:...]]
# --------------------------------------------------------------------------

def v_name(I: Inputs, arg: str) -> str:
    if not arg:
        raise NotAvailable("name needs a construct, e.g. name:<construct A>")
    return arg


def v_samples(I: Inputs, arg: str) -> str:
    s = I.samples()
    if not s:
        raise NotAvailable("no sample is named in the manifest or the stage outputs")
    return "%s (reference sample: %s)" % (", ".join(s), s[0])


def v_host(I: Inputs, arg: str) -> str:
    h = I.host
    if not h:
        raise NotAvailable("no host stage was run (host_reference.py)")
    parts = ", ".join(str(h[k]) for k in ("species", "assembly")
                      if h.get(k) and h[k] != "local files")
    acc = h.get("accession")
    s = "%s%s" % (parts or "host from local files", " (%s)" % acc if acc else "")
    if h.get("decided_by"):
        s += ", decided by %s" % h["decided_by"]
    if h.get("cell_line"):
        s += ", cell line %s" % h["cell_line"]
    return s


def v_ploidy(I: Inputs, arg: str) -> str:
    p = ((I.summary or {}).get("parameters") or {})
    if p.get("ploidy") is not None:
        return "%s (%s)" % (num(p["ploidy"]), p.get("ploidy_source") or "reconcile_sites.py")
    if I.manifest_facts.get("ploidy") is not None:
        return "%s (run manifest)" % num(I.manifest_facts["ploidy"])
    raise NotAvailable("no ploidy in the site summary or the manifest")


def v_modality(I: Inputs, arg: str) -> str:
    mod = I.manifest_facts.get("modality") or (I.cn or {}).get("modality")
    if not mod:
        raise NotAvailable("no modality in the manifest or the copy-number JSON")
    rec = I.classify()
    parts = []
    for m, r in rec.items():
        st = {c.get("id"): c.get("status") for c in r.get("checks") or []}
        parts.append("%s %s" % (m, st.get("modality_fingerprint", "not checked")))
    return "%s (declared)%s" % (mod, "; junction fingerprint check: " + ", ".join(parts)
                                if parts else "")


def v_site_range(I: Inputs, arg: str) -> str:
    pc = I.per_construct()
    cons = I.ordered(pc)
    parts = []
    for c in cons:
        d = pc[c]
        rng = d.get("site_range", "?")
        n = d.get("site_high") or 0
        parts.append("%s %s site%s (%s certain, %s inferred; ceiling %s)" % (
            c, rng, "" if n == 1 and "-" not in str(rng) else "s", d.get("certain", 0),
            num((d.get("vector_end") or 0) + (d.get("copy_ceiling") or 0)), num(d.get("ceiling"))))
    nu = len((I.summary or {}).get("unresolved") or [])
    return "%s; %s, attributed to no construct" % (
        "; ".join(parts), plural(nu, "locus unresolved", "loci unresolved"))


def v_sites(I: Inputs, arg: str) -> str:
    pc = I.per_construct()
    if arg not in pc:
        raise NotAvailable("the site summary has no construct %s" % arg)
    return str(pc[arg].get("site_range", "—"))


def _count(I: Inputs, tier: str) -> str:
    return str(sum(1 for r in I.called_rows() if r.get("confidence") == tier))


def v_site_derivation(I: Inputs, arg: str) -> str:
    s = I.site_summary()
    p = s.get("parameters") or {}
    rec = I.classify()
    per = ", ".join("%s %s" % (m, num(r.get("detected_loci"))) for m, r in rec.items()
                    if r.get("detected_loci") is not None)
    rows = I.called or []
    tiers = {t: sum(1 for r in rows if r.get("confidence") == t)
             for t in ("certain", "inferred", "unresolved")}
    return ("Stage 2 reported %s locus entries%s; collapsing one event seen by several "
            "mappings left %s distinct loci; %s met retention (≥ %s reads in ≥ %s samples) "
            "and %s went to the sub-threshold tail; of the %s called, %s are certain, "
            "%s inferred and %s unresolved" % (
                num(s.get("n_entries")), " (%s)" % per if per else "", num(s.get("n_unified")),
                num(s.get("n_called")), num(p.get("min_reads")), num(p.get("min_samples")),
                num(s.get("n_sub")), num(s.get("n_called")), tiers["certain"],
                tiers["inferred"], tiers["unresolved"]))


def v_copy_number_route(I: Inputs, arg: str) -> str:
    out = []
    for s, sd in I.cn_samples().items():
        solve = sd.get("solve") or {}
        base = sd.get("baseline_depth") or sd.get("genome_coverage")
        out.append("%s: %s, baseline %s× (%s)" % (
            s, solve.get("method", "route not recorded"), num(base),
            sd.get("baseline_source") or "read-count estimate"))
    if not out:
        raise NotAvailable("the copy-number JSON holds no sample")
    return "; ".join(out)


def v_copy_number(I: Inputs, arg: str) -> str:
    out = []
    for s, sd in I.cn_samples().items():
        cph = (sd.get("copies_per_haploid") or {}).get(arg)
        if cph is None:
            continue
        solve = sd.get("solve") or {}
        u = (solve.get("unique_blocks") or {}).get(arg) or {}
        txt = "%s copies per haploid genome in %s (%s per cell)" % (
            num(cph), s, num((sd.get("copies_per_cell") or {}).get(arg)))
        if u.get("n_blocks"):
            txt += ", from %s of %s retained unique bp; per-block estimates %s–%s" % (
                plural(u["n_blocks"], "block"), num(u.get("retained_bp")),
                num(u.get("min")), num(u.get("max")))
        else:
            txt += ", source: %s" % ((solve.get("source") or {}).get(arg)
                                     or solve.get("method", "not recorded"))
        out.append(txt)
    if not out:
        raise NotAvailable("the copy-number JSON has no copy number for %s" % arg)
    return "; ".join(out)


def v_deficit(I: Inputs, arg: str) -> str:
    out = []
    for s, sd in I.cn_samples().items():
        cl = (sd.get("solve") or {}).get("closure") or {}
        if cl.get("deficit_pct") is not None:
            out.append("%s %+.1f%% (observed %s× vs %s× predicted)" % (
                s, cl["deficit_pct"], num(cl.get("observed"), 1), num(cl.get("predicted"), 1)))
        else:
            out.append("%s n/a (%s)" % (
                s, "no shared class with every carrier measured" if sd.get("solve")
                else "copy number not measured from depth"))
    if not out:
        raise NotAvailable("the copy-number JSON holds no sample")
    return "; ".join(out)


def v_ceiling_divisor(I: Inputs, arg: str) -> str:
    p = ((I.summary or {}).get("parameters") or {})
    if p.get("min_copies_per_site") is None:
        raise NotAvailable("no site summary parameters (reconcile_sites.py has not run)")
    return "%s copies per site (%s)" % (num(p["min_copies_per_site"]),
                                        p.get("min_copies_per_site_source") or "default")


def v_flank_detection(I: Inputs, arg: str) -> str:
    out = []
    for s, sd in I.cn_samples().items():
        fd = sd.get("flank_detection") or []
        rates = [f.get("flank_detection_rate") for f in fd if f.get("flank_detection_rate") is not None]
        if rates:
            out.append("%s: %s–%s across anchors of %s–%s reads" % (
                s, num(min(rates)), num(max(rates)), num(min(f.get("anchor") for f in fd)),
                num(max(f.get("anchor") for f in fd))))
    if not out:
        raise NotAvailable("the copy-number JSON reports no flank detection rate")
    return "flank detection rate " + "; ".join(out)


def v_multi_locus_contigs(I: Inputs, arg: str) -> str:
    counts: Dict[str, int] = {}
    for r in I.called_rows():
        counts[r.get("contig", "")] = counts.get(r.get("contig", ""), 0) + 1
    many = [(c, n) for c, n in counts.items() if n > 1]
    if not many:
        return "no contig carries more than one called locus"
    return ", ".join("%s (%d loci)" % cn for cn in many)


VALUES = {
    "name": (v_name, "a construct's name: name:<construct A>"),
    "samples": (v_samples, "sample names and the reference sample"),
    "host": (v_host, "species, assembly, accession and how decided (Stage 2b)"),
    "ploidy": (v_ploidy, "ploidy and where it came from"),
    "modality": (v_modality, "declared modality and the Stage 2 fingerprint check"),
    "site_range": (v_site_range, "every construct's bounded site range with the unresolved count"),
    "sites": (v_sites, "one construct's site range: sites:<construct A>"),
    "n_certain": (lambda I, a: _count(I, "certain"), "called loci tiered certain"),
    "n_inferred": (lambda I, a: _count(I, "inferred"), "called loci tiered inferred"),
    "n_unresolved": (lambda I, a: _count(I, "unresolved"), "called loci left unresolved"),
    "site_derivation": (v_site_derivation, "Stage 2 entries -> collapsed -> retained -> tiers"),
    "copy_number_route": (v_copy_number_route, "how copy number was obtained, and the baseline"),
    "copy_number": (v_copy_number, "one construct's copy number: copy_number:<construct A>"),
    "deficit": (v_deficit, "the shared-block closure deficit per sample"),
    "ceiling_divisor": (v_ceiling_divisor, "copies per site the ceiling divides by, and its source"),
    "flank_detection": (v_flank_detection, "flank detection rate"),
    "multi_locus_contigs": (v_multi_locus_contigs, "contigs carrying more than one called locus"),
}


# --------------------------------------------------------------------------
# figures
# --------------------------------------------------------------------------

def load_figure(fig_dir: str, rec: Dict) -> Tuple[Optional[str], str]:
    """The SVG markup ready to inline, or (None, why it is not the skill's)."""
    path = os.path.join(fig_dir, rec.get("file") or "")
    try:
        with open(path, encoding="utf-8") as fh:
            text = fh.read()
    except OSError as e:
        return None, "%s unreadable (%s)" % (path, e.strerror or e)
    try:
        root = ET.fromstring(text.encode("utf-8"))
    except ET.ParseError as e:
        return None, "%s is not well-formed SVG (%s)" % (path, e)
    want = re.sub(r"[^A-Za-z0-9_-]", "_", rec.get("id", ""))
    if root.tag != "{http://www.w3.org/2000/svg}svg" or root.get("data-figure") != want:
        return None, ("%s does not carry make_figures.py's figure id %r: it is not the "
                      "skill's figure" % (path, want))
    svg = re.sub(r"^\s*<\?xml[^>]*\?>\s*", "", text)
    return svg.strip(), ""


# --------------------------------------------------------------------------
# filling
# --------------------------------------------------------------------------

class Filled:
    def __init__(self):
        self.text = ""
        self.figs: List[Dict] = []        # in order of appearance, with "svg"
        self.placed: List[str] = []
        self.unavailable: List[Tuple[str, str]] = []
        self.unknown: List[str] = []
        self.foreign: List[str] = []
        self.unplaced: List[str] = []
        self.template_marks: List[str] = []


class Filler:
    """Replaces each placeholder in the prose with its script output,
    recording on `F` everything the checks need to know."""

    def __init__(self, I: Inputs):
        self.I, self.F = I, Filled()
        self.figs = {f["id"]: f for f in (I.figures or {}).get("figures", [])}
        self.skipped = {s.get("id"): s.get("reason")
                        for s in (I.figures or {}).get("skipped", [])}
        self._loaded: Dict[str, Tuple[Optional[str], str]] = {}

    # ---- slot text
    def fig_token(self, f: Dict) -> str:
        if f["id"] not in self._loaded:
            self._loaded[f["id"]] = load_figure(self.I.fig_dir, f)
        svg, why = self._loaded[f["id"]]
        if svg is None:
            return "\n\n> **Figure %s not embedded:** %s\n\n" % (f["id"], why)
        self.F.figs.append(dict(f, svg=svg))
        self.F.placed.append(f["id"])
        return "\n\n\x00FIG%d\x00\n\n" % (len(self.F.figs) - 1)

    def note(self, kind: str, key: str, why: str, inline: str = "not available") -> str:
        """A visible note in the slot: a block for figures and tables, inline
        for values."""
        self.F.unavailable.append(("%s:%s" % (kind, key), why))
        if kind == "VAL":
            return "*(%s: %s)*" % (inline, why)
        label = {"FIG": "Figure", "TABLE": "Table"}[kind]
        return "\n\n> **%s %s not generated:** %s\n\n" % (label, key, why)

    def drop_foreign(self, m) -> str:
        self.F.foreign.append(re.sub(r"\s+", " ", m.group(0))[:80])
        return "\n\n> **Image removed:** the report embeds only the figures the skill " \
               "drew (make_figures.py). Add a figure there, not here.\n\n"

    def leftover(self, m) -> str:
        """Anything still bracketed is a placeholder this script does not know."""
        self.F.unknown.append(m.group(0))
        return "*(unknown placeholder %s)*" % m.group(0)[2:-2]

    # ---- placeholders
    def sub(self, m) -> str:
        kind, raw = m.group(1), m.group(2).strip()
        try:
            key = resolve_letters(raw, self.I.constructs)
        except NotAvailable as e:
            return self.note(kind, re.sub(r"[<>]", "", raw), str(e), "not applicable")
        if kind == "FIG":
            return self.figure(key)
        return self.data(kind, key, m.group(0))

    def figure(self, key: str) -> str:
        if self.I.figures is None:
            return self.note("FIG", key, "no figures/figures.json: run make_figures.py "
                             "--run run.json before build_report.py")
        if key in self.figs:
            return self.fig_token(self.figs[key])
        if key.startswith("coverage_"):
            hits = [f for f in self.figs.values()
                    if f.get("kind") == "coverage" and f.get("construct") == key[9:]]
            if hits:
                return "".join(self.fig_token(f) for f in hits)
        return self.note("FIG", key, self.skipped.get(key)
                         or self.skipped.get(key.split("_")[0])
                         or "make_figures.py drew no such figure")

    def data(self, kind: str, key: str, whole: str) -> str:
        name, _, arg = key.partition(":")
        table = TABLES if kind == "TABLE" else VALUES
        if name not in table:
            self.F.unknown.append(whole)
            if kind == "VAL":
                return "*(unknown value %s)*" % key
            return self.note(kind, key, "unknown placeholder id (build_report.py "
                             "--list-placeholders lists them)")
        try:
            out = table[name][0](self.I, arg)
        except NotAvailable as e:
            return self.note(kind, key, str(e))
        if kind == "TABLE":
            return "\n\n%s\n\n" % out
        return out.replace("|", "/").replace("\n", " ")

    # ---- after the placeholders
    def appendix(self, text: str) -> str:
        """Every skill figure appears: what no placeholder asked for is appended."""
        unplaced = [f for f in self.figs.values() if f["id"] not in self.F.placed]
        if not unplaced:
            return text
        self.F.unplaced = [f["id"] for f in unplaced]
        return text.rstrip() + "\n\n## Appendix: figures no template slot asked for\n\n" \
            "These were drawn by make_figures.py for this run but no placeholder in the " \
            "report placed them.\n\n" + "".join(self.fig_token(f) for f in unplaced)


def template_marks(text: str) -> List[str]:
    """The <...> template text still waiting for prose, each once, shortened."""
    marks: List[str] = []
    for m in TEMPLATE_MARK.finditer(FIG_TOKEN.sub("", text)):
        s = re.sub(r"\s+", " ", m.group(0))
        s = s if len(s) <= 60 else s[:56] + "...>"
        if s not in marks:
            marks.append(s)
    return marks


def fill(prose: str, I: Inputs) -> Filled:
    fl = Filler(I)
    # the template's own instructions are documentation, not report text
    text = COMMENT.sub("", prose)
    text = FOREIGN_IMAGE.sub(fl.drop_foreign, text)
    text = PLACEHOLDER.sub(fl.sub, text)
    text = ANY_PLACEHOLDER.sub(fl.leftover, text)
    text = fl.appendix(text)
    fl.F.template_marks = template_marks(text)
    fl.F.text = re.sub(r"\n{3,}", "\n\n", text).strip() + "\n"
    return fl.F


# --------------------------------------------------------------------------
# Markdown subset -> blocks (shared by the HTML and the docx renderers)
# --------------------------------------------------------------------------

INLINE = re.compile(r"`([^`]+)`|\*\*(.+?)\*\*|(?<![\w*])\*(?!\s)([^*\n]+?)(?<!\s)\*(?![\w*])"
                    r"|\[([^\]]+)\]\(([^)\s]+)\)")
TABLE_SEP = re.compile(r"^\|?\s*:?-{3,}:?\s*(\|\s*:?-{3,}:?\s*)*\|?$")
LIST_ITEM = re.compile(r"^(\d+\.|[-*])\s+(.*)$")
LABEL_LINE = re.compile(r"^\*\*[^*]+:\*\*")


def segments(text: str, b: bool = False, i: bool = False, href: Optional[str] = None):
    """(text, bold, italic, code, href) runs of one line of inline Markdown."""
    out, pos = [], 0
    for m in INLINE.finditer(text):
        if m.start() > pos:
            out.append((text[pos:m.start()], b, i, False, href))
        if m.group(1) is not None:
            out.append((m.group(1), b, i, True, href))
        elif m.group(2) is not None:
            out += segments(m.group(2), True, i, href)
        elif m.group(3) is not None:
            out += segments(m.group(3), b, True, href)
        else:
            out += segments(m.group(4), b, i, m.group(5))
        pos = m.end()
    if pos < len(text):
        out.append((text[pos:], b, i, False, href))
    return out


def split_row(line: str) -> List[str]:
    s = line.strip()
    if s.startswith("|"):
        s = s[1:]
    if s.endswith("|") and not s.endswith("\\|"):
        s = s[:-1]
    return [c.strip().replace("\\|", "|") for c in re.split(r"(?<!\\)\|", s)]


HEADING = re.compile(r"(#{1,6})\s+(.*)$")
RULE = re.compile(r"-{3,}|\*{3,}")


def _is_table(lines: List[str], j: int) -> bool:
    return (lines[j].strip().startswith("|") and j + 1 < len(lines)
            and bool(TABLE_SEP.match(lines[j + 1].strip())))


def _starts_block(lines: List[str], j: int) -> bool:
    s = lines[j].strip()
    return bool(HEADING.match(s) or FIG_TOKEN.fullmatch(s) or s.startswith(">")
                or LIST_ITEM.match(s) or RULE.fullmatch(s) or _is_table(lines, j))


def _take_table(lines: List[str], i: int) -> Tuple[Tuple, int]:
    header, rows = split_row(lines[i]), []
    i += 2
    while i < len(lines) and lines[i].strip().startswith("|"):
        rows.append(split_row(lines[i]))
        i += 1
    return ("table", header, rows), i


def _take_quote(lines: List[str], i: int) -> Tuple[Tuple, int]:
    parts = []
    while i < len(lines) and lines[i].strip().startswith(">"):
        parts.append(re.sub(r"^>\s?", "", lines[i].strip()))
        i += 1
    return ("quote", " ".join(p for p in parts if p)), i


def _take_list(lines: List[str], i: int) -> Tuple[Tuple, int]:
    """Items of one kind (ordered or not); indented lines continue an item,
    and one blank line between items of the same kind does not end the list."""
    ordered = LIST_ITEM.match(lines[i].strip()).group(1)[0].isdigit()

    def same_kind(s: str) -> Optional[re.Match]:
        m = LIST_ITEM.match(s)
        return m if m and m.group(1)[0].isdigit() == ordered else None

    items: List[str] = []
    while i < len(lines):
        s = lines[i].strip()
        m = same_kind(s)
        if m:
            items.append(m.group(2))
        elif s and lines[i][:1] in (" ", "\t") and items:
            items[-1] += " " + s
        elif not (not s and i + 1 < len(lines) and same_kind(lines[i + 1].strip())):
            break
        i += 1
    return ("ol" if ordered else "ul", items), i


def _take_paragraph(lines: List[str], i: int) -> Tuple[Tuple, int]:
    """A paragraph; a `**Label:**` line starts a new line within it."""
    pieces: List[str] = []
    while i < len(lines) and lines[i].strip() and (not pieces or not _starts_block(lines, i)):
        s = lines[i].strip()
        if not pieces or LABEL_LINE.match(s):
            pieces.append(s)
        else:
            pieces[-1] += " " + s
        i += 1
    return ("p", pieces), i


def parse_blocks(md: str) -> List[Tuple]:
    lines = md.split("\n")
    blocks: List[Tuple] = []
    i = 0
    while i < len(lines):
        s = lines[i].strip()
        fig, head = FIG_TOKEN.fullmatch(s), HEADING.match(s)
        if not s:
            i += 1
            continue
        if fig or head or RULE.fullmatch(s):
            blocks.append(("fig", int(fig.group(1))) if fig else
                          ("h", len(head.group(1)), head.group(2).strip()) if head else ("hr",))
            i += 1
            continue
        take = (_take_table if _is_table(lines, i) else _take_quote if s.startswith(">")
                else _take_list if LIST_ITEM.match(s) else _take_paragraph)
        block, i = take(lines, i)
        blocks.append(block)
    return blocks


# --------------------------------------------------------------------------
# HTML
# --------------------------------------------------------------------------

def inline_html(text: str) -> str:
    out = []
    for t, b, i, code, href in segments(text):
        h = html.escape(t, quote=False)
        if code:
            h = "<code>%s</code>" % h
        if i:
            h = "<em>%s</em>" % h
        if b:
            h = "<strong>%s</strong>" % h
        if href:
            h = '<a href="%s">%s</a>' % (html.escape(href), h)
        out.append(h)
    return "".join(out)


def slug(text: str, used: Dict[str, int]) -> str:
    s = re.sub(r"[^a-z0-9]+", "-", re.sub(r"<[^>]*>", "", text.lower())).strip("-") or "s"
    used[s] = used.get(s, 0) + 1
    return s if used[s] == 1 else "%s-%d" % (s, used[s])


PAGE_CSS = """
:root{--bg:#fcfcfb;--fg:#0b0b0b;--fg2:#52514e;--rule:#e4e3df;--accent:#2a78d6;
--note-bg:#f4f3ef;--warn:#b54708;--table-head:#f4f3ef}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--bg:#1a1a19;
--fg:#f0efec;--fg2:#c3c2b7;--rule:#33332f;--accent:#3987e5;--note-bg:#242422;
--warn:#f79009;--table-head:#242422}}
:root[data-theme="dark"]{--bg:#1a1a19;--fg:#f0efec;--fg2:#c3c2b7;--rule:#33332f;
--accent:#3987e5;--note-bg:#242422;--warn:#f79009;--table-head:#242422}
*{box-sizing:border-box}
html{-webkit-text-size-adjust:100%}
body{margin:0;background:var(--bg);color:var(--fg);
font:15px/1.55 system-ui,-apple-system,"Segoe UI",Helvetica,Arial,sans-serif}
main{max-width:980px;margin:0 auto;padding:32px 16px 64px}
h1{font-size:1.75rem;line-height:1.25;margin:0 0 .75rem}
h2{font-size:1.3rem;margin:2.2rem 0 .6rem;padding-top:.6rem;border-top:1px solid var(--rule)}
h3{font-size:1.08rem;margin:1.6rem 0 .5rem}
p,ul,ol{margin:.55rem 0}
a{color:var(--accent)}
code{font:.88em ui-monospace,SFMono-Regular,Menlo,Consolas,monospace}
hr{border:0;border-top:1px solid var(--rule);margin:1.5rem 0}
blockquote{margin:.8rem 0;padding:.55rem .9rem;background:var(--note-bg);
border-left:3px solid var(--warn);border-radius:2px}
blockquote p{margin:0}
.table-wrap{overflow-x:auto;margin:.8rem 0;max-width:100%}
table{border-collapse:collapse;font-size:.84rem;min-width:50%}
th,td{border:1px solid var(--rule);padding:.32rem .5rem;text-align:left;vertical-align:top}
th{background:var(--table-head);font-weight:600}
td{font-variant-numeric:tabular-nums}
figure{margin:1.2rem 0;break-inside:avoid;page-break-inside:avoid}
figure svg{display:block;max-width:100%;height:auto}
figcaption{color:var(--fg2);font-size:.86rem;margin-top:.4rem}
.provenance{color:var(--fg2);font-size:.8rem;margin-top:3rem;border-top:1px solid var(--rule);
padding-top:.8rem}
@media print{:root{--bg:#fff;--fg:#000;--fg2:#333;--rule:#bbb;--note-bg:#f4f4f4;
--table-head:#eee}main{max-width:none;padding:0}h2{break-after:avoid}
table{font-size:8.5pt}.table-wrap{overflow:visible}}
"""


def print_light_css(figs: List[Dict]) -> str:
    """make_figures.py's SVGs switch palette with prefers-color-scheme; on
    paper they must stay light. Its own light rules, re-applied for print."""
    for f in figs:
        m = re.search(r"<style>(.*?)@media", f["svg"], re.S)
        if m:
            rules = re.sub(r"\b(fill|stroke):(#[0-9a-fA-F]{3,8})", r"\1:\2 !important",
                           re.sub(r"^text\{[^}]*\}", "", m.group(1)))
            return "@media print{%s}" % rules
    return ""


def _html_figure(f: Dict, n: int) -> str:
    return ('<figure id="fig-%s" data-figure-id="%s" data-source="make_figures.py">\n'
            "%s\n<figcaption><strong>Figure %d.</strong> %s. %s</figcaption>\n</figure>"
            % (html.escape(f["id"]), html.escape(f["id"]), f["svg"], n,
               html.escape(f.get("title") or f["id"]), html.escape(f.get("caption") or "")))


def _html_table(header: List[str], rows: List[List[str]]) -> str:
    head = "".join("<th>%s</th>" % inline_html(c) for c in header)
    body = "".join("<tr>%s</tr>" % "".join("<td>%s</td>" % inline_html(c) for c in r)
                   for r in rows)
    return ('<div class="table-wrap"><table><thead><tr>%s</tr></thead>'
            "<tbody>%s</tbody></table></div>" % (head, body))


def _html_block(b: Tuple, used: Dict[str, int]) -> str:
    kind = b[0]
    if kind == "h":
        return '<h%d id="%s">%s</h%d>' % (b[1], slug(b[2], used), inline_html(b[2]), b[1])
    if kind == "p":
        return "<p>%s</p>" % "<br>\n".join(inline_html(p) for p in b[1])
    if kind == "hr":
        return "<hr>"
    if kind == "quote":
        return "<blockquote><p>%s</p></blockquote>" % inline_html(b[1])
    if kind in ("ul", "ol"):
        return "<%s>%s</%s>" % (kind, "".join("<li>%s</li>" % inline_html(x) for x in b[1]),
                                kind)
    return _html_table(b[1], b[2])


def render_html(blocks: List[Tuple], figs: List[Dict], title: str, footer: str) -> str:
    out: List[str] = []
    used: Dict[str, int] = {}
    nfig = 0
    for b in blocks:
        if b[0] == "fig":
            nfig += 1
            out.append(_html_figure(figs[b[1]], nfig))
        else:
            out.append(_html_block(b, used))
    return ("<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n"
            "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
            "<meta name=\"generator\" content=\"vector-integration-analysis build_report.py\">\n"
            "<title>%s</title>\n<style>%s%s</style>\n</head>\n<body>\n<main>\n%s\n"
            "<p class=\"provenance\">%s</p>\n</main>\n</body>\n</html>\n"
            % (html.escape(title), PAGE_CSS, print_light_css(figs), "\n".join(out),
               html.escape(footer)))


def render_md(text: str, figs: List[Dict], fig_dir: str, md_path: str) -> str:
    base = os.path.dirname(os.path.abspath(md_path))

    def sub(m):
        f = figs[int(m.group(1))]
        rel = os.path.relpath(os.path.join(fig_dir, f["file"]), base).replace(os.sep, "/")
        return "![%s](%s)\n\n*Figure: %s*" % (f.get("title", f["id"]), rel, f.get("caption", ""))
    return FIG_TOKEN.sub(sub, text)


# --------------------------------------------------------------------------
# optional docx: the same blocks, the skill's SVGs rasterised
# --------------------------------------------------------------------------

def find_rasteriser(choice: str, command: str) -> Tuple[Optional[Tuple[str, str]], str]:
    if command:
        return ("command", command), ""
    if choice in ("auto", "cairosvg"):
        try:
            import cairosvg  # noqa: F401  (optional, docx only)
            return ("cairosvg", ""), ""
        except Exception:  # noqa: BLE001 -- ImportError, or a missing libcairo
            if choice == "cairosvg":
                return None, "cairosvg is not importable"
    if choice in ("auto", "rsvg-convert"):
        exe = shutil.which("rsvg-convert")
        if exe:
            return ("rsvg-convert", exe), ""
        if choice == "rsvg-convert":
            return None, "rsvg-convert is not on PATH"
    if choice == "none":
        return None, "--rasteriser none"
    return None, "no SVG rasteriser found: neither the cairosvg package nor rsvg-convert"


def rasterise(r: Tuple[str, str], svg_path: str, png_path: str) -> None:
    kind, arg = r
    if kind == "cairosvg":
        import cairosvg
        cairosvg.svg2png(url=svg_path, write_to=png_path, scale=2)
    elif kind == "rsvg-convert":
        subprocess.run([arg, "--format=png", "--zoom=2", "--output", png_path, svg_path],
                       check=True, capture_output=True)
    else:
        cmd = [part.format(svg=svg_path, png=png_path) for part in shlex.split(arg)]
        subprocess.run(cmd, check=True, capture_output=True)


def png_size(data: bytes) -> Tuple[int, int]:
    if data[:8] != b"\x89PNG\r\n\x1a\n" or data[12:16] != b"IHDR":
        raise ValueError("not a PNG")
    return struct.unpack(">II", data[16:24])


def _x(s: str) -> str:
    return html.escape(s, quote=False)


def _runs(text: str, size: Optional[int] = None, italic_all: bool = False) -> str:
    out = []
    for t, b, i, code, _href in segments(text):
        rpr = ""
        if b:
            rpr += "<w:b/>"
        if i or italic_all:
            rpr += "<w:i/>"
        if code:
            rpr += '<w:rFonts w:ascii="Consolas" w:hAnsi="Consolas"/>'
        if size:
            rpr += '<w:sz w:val="%d"/><w:szCs w:val="%d"/>' % (size, size)
        out.append('<w:r>%s<w:t xml:space="preserve">%s</w:t></w:r>'
                   % ("<w:rPr>%s</w:rPr>" % rpr if rpr else "", _x(t)))
    return "".join(out)


def _p(body: str, ppr: str = "") -> str:
    return "<w:p>%s%s</w:p>" % ("<w:pPr>%s</w:pPr>" % ppr if ppr else "", body)


TEXT_W_TWIPS = 9638          # A4, 2 cm margins
TEXT_W_EMU = 6000000


def _docx_table(header: List[str], rows: List[List[str]]) -> str:
    ncol = max([len(header)] + [len(r) for r in rows])
    w = TEXT_W_TWIPS // max(ncol, 1)

    def row(cells, head=False):
        cells = list(cells) + [""] * (ncol - len(cells))
        shade = '<w:shd w:val="clear" w:color="auto" w:fill="EEEEEA"/>' if head else ""
        tcs = "".join('<w:tc><w:tcPr><w:tcW w:w="%d" w:type="dxa"/>%s</w:tcPr>%s</w:tc>'
                      % (w, shade, _p(_runs("**%s**" % c if head and c else c, size=16),
                                      '<w:spacing w:before="0" w:after="0"/>'))
                      for c in cells)
        return "<w:tr>%s%s</w:tr>" % ("<w:trPr><w:tblHeader/></w:trPr>" if head else "", tcs)

    bd = "".join('<w:%s w:val="single" w:sz="4" w:space="0" w:color="BFBEB8"/>' % e
                 for e in ("top", "left", "bottom", "right", "insideH", "insideV"))
    grid = "".join('<w:gridCol w:w="%d"/>' % w for _ in range(ncol))
    return ('<w:tbl><w:tblPr><w:tblW w:w="%d" w:type="dxa"/><w:tblBorders>%s</w:tblBorders>'
            '<w:tblLayout w:type="autofit"/></w:tblPr><w:tblGrid>%s</w:tblGrid>%s%s</w:tbl>'
            % (TEXT_W_TWIPS, bd, grid, row(header, True), "".join(row(r) for r in rows))
            + _p(""))


def _docx_figure(k: int, f: Dict, png: bytes, n: int) -> str:
    """The rasterised skill figure, inline at text width, then its caption."""
    pw, ph = png_size(png)
    cx = min(TEXT_W_EMU, int(pw / 2 * 9525))
    cy = int(cx * ph / pw)
    pic = ('<w:r><w:drawing><wp:inline distT="0" distB="0" distL="0" distR="0">'
           '<wp:extent cx="%d" cy="%d"/><wp:docPr id="%d" name="%s"/>'
           '<wp:cNvGraphicFramePr><a:graphicFrameLocks noChangeAspect="1"/>'
           '</wp:cNvGraphicFramePr><a:graphic><a:graphicData '
           'uri="http://schemas.openxmlformats.org/drawingml/2006/picture"><pic:pic>'
           '<pic:nvPicPr><pic:cNvPr id="%d" name="%s.png"/><pic:cNvPicPr/></pic:nvPicPr>'
           '<pic:blipFill><a:blip r:embed="rIdImg%d"/><a:stretch><a:fillRect/></a:stretch>'
           '</pic:blipFill><pic:spPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="%d" cy="%d"/>'
           '</a:xfrm><a:prstGeom prst="rect"><a:avLst/></a:prstGeom></pic:spPr></pic:pic>'
           "</a:graphicData></a:graphic></wp:inline></w:drawing></w:r>"
           % (cx, cy, n, _x(f["id"]), n, _x(f["id"]), k, cx, cy))
    caption = _runs("**Figure %d.** %s. %s" % (n, f.get("title") or f["id"],
                                                f.get("caption") or ""), size=18)
    return _p(pic, "<w:keepNext/>") + _p(caption, '<w:pStyle w:val="Caption"/>')


def _docx_block(b: Tuple) -> str:
    kind = b[0]
    if kind == "h":
        return _p(_runs(b[2]), '<w:pStyle w:val="Heading%d"/>' % min(b[1], 3))
    if kind == "p":
        return _p('<w:r><w:br/></w:r>'.join(_runs(x) for x in b[1]))
    if kind == "hr":
        return _p("", '<w:pBdr><w:bottom w:val="single" w:sz="4" w:space="1" '
                  'w:color="A3A29C"/></w:pBdr>')
    if kind == "quote":
        return _p(_runs(b[1]), '<w:pBdr><w:left w:val="single" w:sz="18" w:space="8" '
                  'w:color="B54708"/></w:pBdr><w:ind w:left="284"/>')
    if kind in ("ul", "ol"):
        return "".join(_p('<w:r><w:t xml:space="preserve">%s\t</w:t></w:r>'
                          % ("%d." % k if kind == "ol" else "•") + _runs(item),
                          '<w:ind w:left="425" w:hanging="425"/>')
                       for k, item in enumerate(b[1], 1))
    return _docx_table(b[1], b[2])


def docx_document(blocks: List[Tuple], figs: List[Dict], pngs: Dict[int, bytes]) -> str:
    body: List[str] = []
    nfig = 0
    for b in blocks:
        if b[0] == "fig":
            nfig += 1
            body.append(_docx_figure(b[1], figs[b[1]], pngs[b[1]], nfig))
        else:
            body.append(_docx_block(b))
    return ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
            '<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" '
            'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" '
            'xmlns:wp="http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing" '
            'xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" '
            'xmlns:pic="http://schemas.openxmlformats.org/drawingml/2006/picture">'
            "<w:body>%s<w:sectPr><w:pgSz w:w=\"11906\" w:h=\"16838\"/><w:pgMar w:top=\"1134\" "
            "w:right=\"1134\" w:bottom=\"1134\" w:left=\"1134\" w:header=\"708\" "
            "w:footer=\"708\" w:gutter=\"0\"/></w:sectPr></w:body></w:document>"
            % "".join(body))


DOCX_STYLES = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:styles xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">
<w:docDefaults><w:rPrDefault><w:rPr><w:rFonts w:ascii="Calibri" w:hAnsi="Calibri" w:cs="Calibri"/>
<w:sz w:val="21"/><w:szCs w:val="21"/></w:rPr></w:rPrDefault>
<w:pPrDefault><w:pPr><w:spacing w:after="120" w:line="264" w:lineRule="auto"/></w:pPr></w:pPrDefault>
</w:docDefaults>
<w:style w:type="paragraph" w:default="1" w:styleId="Normal"><w:name w:val="Normal"/></w:style>
<w:style w:type="paragraph" w:styleId="Heading1"><w:name w:val="heading 1"/><w:basedOn w:val="Normal"/>
<w:next w:val="Normal"/><w:qFormat/><w:pPr><w:keepNext/><w:spacing w:before="240" w:after="120"/>
<w:outlineLvl w:val="0"/></w:pPr><w:rPr><w:b/><w:sz w:val="34"/></w:rPr></w:style>
<w:style w:type="paragraph" w:styleId="Heading2"><w:name w:val="heading 2"/><w:basedOn w:val="Normal"/>
<w:next w:val="Normal"/><w:qFormat/><w:pPr><w:keepNext/><w:spacing w:before="320" w:after="100"/>
<w:outlineLvl w:val="1"/></w:pPr><w:rPr><w:b/><w:sz w:val="27"/></w:rPr></w:style>
<w:style w:type="paragraph" w:styleId="Heading3"><w:name w:val="heading 3"/><w:basedOn w:val="Normal"/>
<w:next w:val="Normal"/><w:qFormat/><w:pPr><w:keepNext/><w:spacing w:before="240" w:after="80"/>
<w:outlineLvl w:val="2"/></w:pPr><w:rPr><w:b/><w:sz w:val="23"/></w:rPr></w:style>
<w:style w:type="paragraph" w:styleId="Caption"><w:name w:val="caption"/><w:basedOn w:val="Normal"/>
<w:qFormat/><w:pPr><w:spacing w:after="240"/></w:pPr><w:rPr><w:color w:val="52514E"/></w:rPr></w:style>
</w:styles>
"""


def _docx_parts(doc: str, pngs: Dict[int, bytes]) -> Dict[str, str]:
    """The package's XML parts: content types, relationships, document, styles."""
    rels = ['<Relationship Id="rIdStyles" Type="http://schemas.openxmlformats.org/'
            'officeDocument/2006/relationships/styles" Target="styles.xml"/>']
    for k in sorted(pngs):
        rels.append('<Relationship Id="rIdImg%d" Type="http://schemas.openxmlformats.org/'
                    'officeDocument/2006/relationships/image" Target="media/figure%d.png"/>'
                    % (k, k))
    return {
        "[Content_Types].xml":
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n<Types xmlns="http://'
            'schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" '
            'ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
            '<Default Extension="xml" ContentType="application/xml"/><Default Extension="png" '
            'ContentType="image/png"/><Override PartName="/word/document.xml" ContentType="'
            'application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>'
            '<Override PartName="/word/styles.xml" ContentType="application/vnd.openxmlformats-'
            'officedocument.wordprocessingml.styles+xml"/></Types>',
        "_rels/.rels":
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n<Relationships xmlns="'
            'http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="'
            'rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/'
            'officeDocument" Target="word/document.xml"/></Relationships>',
        "word/_rels/document.xml.rels":
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n<Relationships xmlns="'
            'http://schemas.openxmlformats.org/package/2006/relationships">%s</Relationships>'
            % "".join(rels),
        "word/document.xml": doc,
        "word/styles.xml": DOCX_STYLES,
    }


def write_docx(path: str, blocks: List[Tuple], figs: List[Dict], pngs: Dict[int, bytes]) -> None:
    parts = _docx_parts(docx_document(blocks, figs, pngs), pngs)
    tmp = path + ".tmp"
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as z:
        for name, data in list(parts.items()) + [("word/media/figure%d.png" % k, pngs[k])
                                                 for k in sorted(pngs)]:
            zi = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
            zi.compress_type = zipfile.ZIP_DEFLATED
            z.writestr(zi, data)
    os.replace(tmp, path)


def build_docx(path: str, blocks, figs, fig_dir: str, r) -> Tuple[bool, str]:
    """(written, why not). Never a partial docx, never a substitute figure."""
    pngs: Dict[int, bytes] = {}
    with tempfile.TemporaryDirectory() as td:
        for k, f in enumerate(figs):
            out = os.path.join(td, "f%d.png" % k)
            try:
                rasterise(r, os.path.join(fig_dir, f["file"]), out)
                with open(out, "rb") as fh:
                    pngs[k] = fh.read()
                png_size(pngs[k])
            except Exception as e:  # noqa: BLE001 -- any rasteriser failure skips the docx
                return False, "rasterising %s with %s failed (%s)" % (
                    f["id"], r[0], str(e).strip().splitlines()[-1] if str(e).strip() else
                    type(e).__name__)
    write_docx(path, blocks, figs, pngs)
    return True, ""


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def list_placeholders() -> str:
    L = ["[[FIG:<id>]]      a figure from make_figures.py: coverage_<construct A>, "
         "coverage_<sample>_<construct>, genome_overview, copy_number, architecture", ""]
    L += ["[[TABLE:%s]]%s%s" % (k, " " * max(1, 22 - len(k)), v[1]) for k, v in TABLES.items()]
    L.append("")
    L += ["[[VAL:%s]]%s%s" % (k, " " * max(1, 24 - len(k)), v[1]) for k, v in VALUES.items()]
    L += ["", "<construct A>, <construct B>, ... name constructs in report order: A is the "
          "primary mapping, the rest follow the identity map."]
    return "\n".join(L)


def parse_args(argv: Optional[Sequence[str]]) -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Assemble the final report (self-contained HTML; optional docx) from "
                    "the pipeline outputs and make_figures.py's figures.")
    add_check_arguments(ap)
    ap.add_argument("--report", metavar="REPORT.md",
                    help="the template with the interpretive prose written in, placeholders "
                         "left in place (default: assets/report_template.md as it stands)")
    ap.add_argument("--figures", metavar="DIR",
                    help="make_figures.py's output directory, holding figures.json "
                         "(default: figures beside the manifest)")
    ap.add_argument("--identity", help="identity.json (default: from --run)")
    ap.add_argument("--copy-number", help="copy_number.py JSON (default: from --run)")
    ap.add_argument("--called", help="called_sites.tsv (default: from --run)")
    ap.add_argument("--site-summary", help="site_summary.json (default: from --run)")
    ap.add_argument("--host", metavar="HOST.json", help="host_reference.py's host.json "
                    "(default: the host stage of --run, if any)")
    _add_output_arguments(ap)
    ap.add_argument("--list-placeholders", action="store_true",
                    help="print every placeholder id and exit")
    return ap.parse_args(argv)


def _add_output_arguments(ap: argparse.ArgumentParser) -> None:
    ap.add_argument("--out-html", help="the HTML report (default: report.html beside the "
                    "manifest, else in the current directory)")
    ap.add_argument("--out-md", help="the filled Markdown (default: report_filled.md "
                    "beside the HTML)")
    ap.add_argument("--json", help="build record (default: report.json beside the HTML)")
    ap.add_argument("--docx", help="the optional Word copy (default: report.docx beside "
                    "the HTML)")
    ap.add_argument("--no-docx", action="store_true", help="do not attempt the docx")
    ap.add_argument("--rasteriser", choices=("auto", "cairosvg", "rsvg-convert", "none"),
                    default="auto", help="SVG rasteriser for the docx (default: auto)")
    ap.add_argument("--rasteriser-cmd", default="", metavar='"CMD {svg} {png}"',
                    help="any other rasteriser, as a command with {svg} and {png}")


def output_paths(a: argparse.Namespace, I: Inputs) -> Dict[str, str]:
    """Every output path, and the input report; refuses an output that is it."""
    html_path = os.path.abspath(a.out_html or os.path.join(I.base, "report.html"))
    here = os.path.dirname(html_path)
    p = {"html": html_path,
         "md": os.path.abspath(a.out_md or os.path.join(here, "report_filled.md")),
         "json": os.path.abspath(a.json or os.path.join(here, "report.json")),
         "docx": os.path.abspath(a.docx or os.path.join(here, "report.docx"))}
    src = os.path.abspath(a.report or DEFAULT_TEMPLATE)
    if src in p.values():
        raise InputError("refusing to overwrite the input report %s with an output; "
                         "name the output differently" % src)
    p["source"] = src
    return p


def read_prose(src: str) -> str:
    try:
        with open(src, encoding="utf-8") as fh:
            return fh.read()
    except OSError as e:
        raise InputError("cannot read the report %s (%s)" % (src, e.strerror or e))


def check_inputs(I: Inputs, checks: Checks) -> None:
    """The stage results, their checks, and the figures, before filling."""
    missing = [k for k, v in (("identity map", I.ident), ("copy-number JSON", I.cn),
                              ("called_sites.tsv", I.called), ("site summary", I.summary))
               if v is None]
    checks.add("report_inputs", FAIL if missing else PASS,
               ("missing %s: the report would be built without the stage results. Run the "
                "stages with --run, or name the files explicitly" % ", ".join(missing))
               if missing else "identity map, copy number, called sites and site summary read")
    flagged = [(s, c) for s, c in I.stage_checks() if c.get("status") == FAIL]
    checks.add("upstream_checks", FAIL if flagged else PASS,
               ("%s failed upstream and was not accepted: %s. That output must not be "
                "reported as it stands; resolve it, or accept it at that stage, and rerun"
                % (plural(len(flagged), "check"),
                   "; ".join("%s [%s]" % (s, c.get("id")) for s, c in flagged)))
               if flagged else "no stage carries a failing check", acceptable=False)
    if I.figures is None:
        checks.add("figures_drawn", FAIL,
                   "no %s: run `make_figures.py --run run.json` first (or point --figures at "
                   "its --out). The report's figures are the skill's, never substitutes"
                   % os.path.join(I.fig_dir, "figures.json"))
        return
    sk = I.figures.get("skipped") or []
    checks.add("figures_drawn", WARN if sk else PASS,
               "%d figure(s) drawn%s" % (len(I.figures.get("figures") or []),
                                         "; skipped by make_figures.py: " + "; ".join(
                                             "%s (%s)" % (s.get("id"), s.get("reason"))
                                             for s in sk) if sk else ""))


def check_fill(I: Inputs, F: Filled, checks: Checks) -> None:
    """What filling found: provenance, foreign images, unplaced figures,
    unfilled placeholders and unwritten prose."""
    bad_prov = [f["id"] for f in (I.figures or {}).get("figures", [])
                if load_figure(I.fig_dir, f)[0] is None]
    checks.add("figure_provenance", FAIL if bad_prov else PASS,
               ("not embedded, because the file is missing or is not the SVG make_figures.py "
                "wrote: %s. Rerun make_figures.py" % ", ".join(bad_prov)) if bad_prov
               else "every embedded figure is make_figures.py's SVG", acceptable=False)
    checks.add("foreign_figure", FAIL if F.foreign else PASS,
               ("the prose carried %s the skill did not draw, removed from the report: %s. "
                "Figures come only from make_figures.py"
                % (plural(len(F.foreign), "image"),
                   "; ".join("`%s`" % f.replace("`", "'") for f in F.foreign)))
               if F.foreign else "no image other than the skill's figures", acceptable=False)
    checks.add("figures_placed", WARN if F.unplaced else PASS,
               ("no placeholder asked for %s; appended under 'Appendix' so none is lost"
                % ", ".join(F.unplaced)) if F.unplaced else "every figure has a template slot")
    na = F.unavailable + [(u, "unknown placeholder") for u in F.unknown]
    checks.add("placeholders", WARN if na else PASS,
               ("%s could not be filled, and say why in place: %s"
                % (plural(len(na), "placeholder"),
                   "; ".join("%s (%s)" % (k, v) for k, v in na)))
               if na else "every placeholder filled from a script output")
    marks = F.template_marks
    checks.add("prose_unfilled", WARN if marks else PASS,
               ("%s still hold template text to replace with prose: %s"
                % (plural(len(marks), "marker"),
                   ", ".join(marks[:12]) + (", ..." if len(marks) > 12 else "")))
               if marks else "no template marker left in the prose")


def with_banner(text: str, checks: Checks) -> str:
    """The failing/accepted-check banner, under the title, before anything is read."""
    banner = checks.banner_md()
    if not banner:
        return text
    head, _, rest = text.partition("\n")
    if head.startswith("# "):
        return head + "\n\n" + "\n".join(banner) + "\n" + rest
    return "\n".join(banner) + "\n" + text


def write_reports(text: str, F: Filled, I: Inputs, p: Dict[str, str]) -> List[Tuple]:
    """report.html and report_filled.md; returns the parsed blocks."""
    blocks = parse_blocks(text)
    title = next((re.sub(r"[*`]", "", b[2]) for b in blocks if b[0] == "h" and b[1] == 1),
                 "Vector integration analysis")
    footer = ("Built by build_report.py from the stage outputs; figures drawn by "
              "make_figures.py and embedded unchanged (%s)."
              % (", ".join(f["id"] for f in F.figs) or "none"))
    os.makedirs(os.path.dirname(p["html"]), exist_ok=True)
    with open(p["html"], "w", encoding="utf-8") as fh:
        fh.write(render_html(blocks, F.figs, title, footer))
    with open(p["md"], "w", encoding="utf-8") as fh:
        fh.write(render_md(text, F.figs, I.fig_dir, p["md"]))
    embedded = list(dict.fromkeys(f["id"] for f in F.figs))
    print("report: %s (self-contained HTML, %d figure(s) inline: %s)"
          % (p["html"], len(embedded), ", ".join(embedded) or "none"))
    print("        %s (Markdown, linking the SVG files)" % p["md"])
    return blocks


def maybe_docx(a: argparse.Namespace, blocks, F: Filled, I: Inputs,
               p: Dict[str, str]) -> Tuple[bool, str]:
    """The optional docx: (written, why not)."""
    written, why = False, ""
    if a.no_docx:
        why = "--no-docx"
    elif not F.figs:
        why = "no skill figure to embed"
    else:
        r, why = find_rasteriser(a.rasteriser, a.rasteriser_cmd)
        if r:
            written, why = build_docx(p["docx"], blocks, F.figs, I.fig_dir, r)
    if written:
        print("        %s (optional docx, the skill's figures rasterised)" % p["docx"])
        return written, why
    if os.path.exists(p["docx"]) and not a.no_docx:
        os.remove(p["docx"])        # never leave a stale docx beside a new report
    print("docx: skipped (%s). The HTML report is the deliverable: %s; open it in "
          "any browser and print to PDF." % (why, p["html"]))
    return written, why


def write_record(I: Inputs, F: Filled, p: Dict[str, str], docx: Tuple[bool, str],
                 checks: Checks) -> None:
    """report.json, and the `report` stage in the manifest."""
    written, why = docx
    rec = {"html": p["html"], "md": p["md"], "docx": p["docx"] if written else None,
           "docx_skipped": None if written else why, "source": p["source"],
           "figures_embedded": list(dict.fromkeys(f["id"] for f in F.figs)),
           "figures_appended": F.unplaced, "checks": checks.as_list()}
    with open(p["json"], "w", encoding="utf-8") as fh:
        json.dump(rec, fh, indent=2)
        fh.write("\n")
    if I.run is not None:
        rel = {k: (I.run.relative(v) if isinstance(v, str) and os.path.isabs(v) else v)
               for k, v in rec.items()}
        rel["json"] = I.run.relative(p["json"])
        I.run.record("report", rel)


def build(argv: Optional[Sequence[str]] = None) -> int:
    a = parse_args(argv)
    if a.list_placeholders:
        print(list_placeholders())
        return 0
    I = Inputs(a)
    p = output_paths(a, I)
    prose = read_prose(p["source"])
    checks = Checks("build_report", a.accept)
    check_inputs(I, checks)
    F = fill(prose, I)
    check_fill(I, F, checks)
    blocks = write_reports(with_banner(F.text, checks), F, I, p)
    docx = maybe_docx(a, blocks, F, I, p)
    write_record(I, F, p, docx, checks)
    return checks.finish(a.strict)


def main(argv: Optional[Sequence[str]] = None) -> int:
    try:
        return build(argv)
    except InputError as e:
        print("error: %s" % e, file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
