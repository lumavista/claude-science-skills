"""Shared helpers for the vector-integration-analysis skill.

Standard library only: parses GenBank/FASTA and delimited junction tables with
no third-party packages, so it runs anywhere Python 3.8+ does.

Three design rules, each of which exists because its absence has produced a
wrong answer in a real analysis:

1. **Nothing is hardcoded to a caller's vocabulary.** Clip/orientation values
   are resolved from the data against a table of known conventions, or declared
   explicitly. When they cannot be resolved the code degrades to single-flank
   mode with a loud warning; it never crashes and never guesses.
2. **Construct names match exactly, or by a declared alias.** Loose prefix
   matching silently reclassifies host contigs as vector (a construct called
   `A` otherwise swallows the real GRCh38 contig `AC233755.1`), which deletes
   genuine integration sites from the locus table without any error.
3. **Any number of constructs and samples.** Two of each is a common case, not
   a structural assumption.

See references/input-formats.md.
"""
from __future__ import annotations

import csv
import json
import re
from dataclasses import asdict, dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

#: tokens a caller may use for "the other side did not map anywhere"
UNMAPPED_TOKENS = {
    "flankingdna", "flanking", "unknown", "unmapped", "unaligned", "none",
    "na", "n/a", "nan", "null", ".", "-", "",
}


def revcomp(s: str) -> str:
    return s.translate(COMPLEMENT)[::-1]


def natural_key(s: str):
    """Sort contig names the way a human would: 1, 2, 10, X, GL000220.1.

    Each chunk carries a type rank so numeric and alphabetic chunks are never
    compared directly; mixed assemblies (`1`, `MT`, `GL000220.1`,
    `scaffold_3`) would otherwise raise TypeError.
    """
    parts = re.split(r"(\d+)", str(s))
    return [(0, int(p), "") if p.isdigit() else (1, 0, p.lower())
            for p in parts if p != ""]


_NORM = re.compile(r"[^a-z0-9]+")


def norm(s: str) -> str:
    return _NORM.sub("_", str(s).strip().lower()).strip("_")


# --------------------------------------------------------------------------
# Clip conventions
# --------------------------------------------------------------------------
#
# Canonical internal values, and the only two that matter:
#
#   "R"  the reference continues UP TO this base, then novel sequence begins.
#        On the host side this is the LEFT flank of an insert.
#   "L"  novel sequence, then the reference RESUMES at this base.
#        On the host side this is the RIGHT flank of an insert.
#
# Callers disagree about spelling and, worse, some emit a strand rather than a
# clip side. A strand column cannot be converted to a clip side, so it is
# rejected rather than reinterpreted.

CLIP_R, CLIP_L = "R", "L"

#: raw value -> canonical, for conventions seen in the wild
CLIP_VOCABULARIES: Sequence[Dict[str, str]] = (
    {"right_clip": CLIP_R, "left_clip": CLIP_L},
    {"rightclip": CLIP_R, "leftclip": CLIP_L},
    {"right": CLIP_R, "left": CLIP_L},
    {"rightmost": CLIP_R, "leftmost": CLIP_L},
    {"r": CLIP_R, "l": CLIP_L},
    {"3": CLIP_R, "5": CLIP_L},
    {"3prime": CLIP_R, "5prime": CLIP_L},
    {"downstream": CLIP_R, "upstream": CLIP_L},
    {"soft_clip_right": CLIP_R, "soft_clip_left": CLIP_L},
    {"sc_right": CLIP_R, "sc_left": CLIP_L},
    {"end": CLIP_R, "start": CLIP_L},
    {"tail": CLIP_R, "head": CLIP_L},
)

#: values that are a STRAND, not a clip side. Present so they can be named in
#: the error rather than silently treated as a clip convention.
STRAND_VALUES = {"+", "-", "plus", "minus", "fwd", "rev", "forward", "reverse",
                 "1", "-1", "f", "w", "c"}


def norm_clip(s: str) -> str:
    """Normalise a clip value without destroying it.

    `norm()` would reduce `+` and `-` to the empty string, which turns a strand
    column into an apparently absent one and hides the real diagnosis. Clip
    values are therefore only lowercased and stripped of whitespace and
    quoting; `'` is dropped so `3'`/`5'` collapse onto `3`/`5`.
    """
    v = str(s).strip().strip("\"'").lower().replace("'", "").replace(" ", "_")
    # `-` is a strand value on its own but a separator inside a word
    return v.replace("-", "_") if len(v) > 1 else v


@dataclass
class ClipResolution:
    """How the clip column was interpreted, and whether it can be trusted."""
    mapping: Dict[str, str] = field(default_factory=dict)
    ok: bool = False
    reason: str = ""
    observed: Tuple[str, ...] = ()

    def convert(self, raw: str) -> str:
        return self.mapping.get(norm_clip(raw), "")


def resolve_clip_vocabulary(values: Iterable[str],
                            override: Optional[Dict[str, str]] = None
                            ) -> ClipResolution:
    """Decide what the clip column means, from the values actually present.

    `override` maps raw values onto `R`/`L` and always wins. Without it, the
    observed value set is matched against `CLIP_VOCABULARIES`; a set that looks
    like a strand, or that matches nothing, is reported as unusable so the
    caller can degrade rather than invent flank pairings.
    """
    obs = tuple(sorted({norm_clip(v) for v in values if norm_clip(v)}))

    if override:
        m = {norm_clip(k): v.upper()[:1] for k, v in override.items()}
        bad = {k: v for k, v in m.items() if v not in (CLIP_R, CLIP_L)}
        if bad:
            return ClipResolution({}, False,
                                  f"--clip-map values must be R or L, got {sorted(bad.values())}",
                                  obs)
        unmapped = [v for v in obs if v not in m]
        return ClipResolution(
            m, True,
            f"declared via --clip-map"
            + (f"; {len(unmapped)} value(s) left unmapped: {unmapped[:4]}" if unmapped else ""),
            obs)

    if not obs:
        return ClipResolution({}, False, "clip column absent or empty", obs)

    if set(obs) <= STRAND_VALUES:
        return ClipResolution(
            {}, False,
            "values look like a STRAND (%s), not a clip side. A strand cannot be "
            "converted to a clip side: it says which way the read aligned, not "
            "which side of the breakpoint the reference continues on. Supply the "
            "caller's clip column, or declare the meaning with --clip-map."
            % ",".join(obs[:4]),
            obs)

    for vocab in CLIP_VOCABULARIES:
        if set(obs) <= set(vocab):
            return ClipResolution(dict(vocab), True,
                                  "recognised convention %s" % "/".join(sorted(vocab)),
                                  obs)

    return ClipResolution(
        {}, False,
        "unrecognised clip values %s. Declare them with "
        "--clip-map '%s=R,%s=L' once you have confirmed the convention against "
        "a junction of known orientation." % (list(obs[:4]), obs[0],
                                              obs[1] if len(obs) > 1 else "<other>"),
        obs)


# --------------------------------------------------------------------------
# Reference sequences
# --------------------------------------------------------------------------

@dataclass
class Feature:
    label: str
    spans: List[Tuple[int, int]]          # 1-based inclusive
    complement: bool = False

    def covers(self, pos: int) -> bool:
        return any(a <= pos <= b for a, b in self.spans)


@dataclass
class Construct:
    name: str
    seq: str
    features: List[Feature] = field(default_factory=list)
    circular: bool = False

    @property
    def length(self) -> int:
        return len(self.seq)

    def span_features(self, a: int, b: int, default: str = "unannotated") -> str:
        """Feature labels overlapping [a, b], excluding whole-construct spans."""
        labs: List[str] = []
        for f in self.features:
            if any(not (b < s or a > e) for s, e in f.spans):
                if sum(e - s + 1 for s, e in f.spans) < 0.75 * self.length:
                    labs.append(f.label)
        return ", ".join(dict.fromkeys(labs)) or default


def _parse_genbank(text: str, name_hint: str) -> Construct:
    text = text.replace("\r", "")
    m = re.search(r"^LOCUS\s+(\S+)", text, re.M)
    name = m.group(1) if m else name_hint

    parts = text.split("ORIGIN", 1)
    if len(parts) < 2:
        raise ValueError("%s: no ORIGIN block found; is this GenBank?" % name_hint)
    seq = "".join(re.findall(r"^\s*\d+\s+([acgtnACGTN ]+)$", parts[1], re.M))
    seq = seq.replace(" ", "").upper()

    circular = bool(re.search(r"^LOCUS\s.*\bcircular\b", text, re.M | re.I))

    feats: List[Feature] = []
    for fm in re.finditer(
        r"^     (\S+)\s+(\S+)\n((?:^                     .*\n)*)",
        parts[0], re.M,
    ):
        kind, loc, quals = fm.group(1), fm.group(2), fm.group(3)
        if kind.lower() == "source":
            continue
        lab = (re.search(r'/label="([^"]+)"', quals)
               or re.search(r'/gene="([^"]+)"', quals)
               or re.search(r'/product="([^"]+)"', quals)
               or re.search(r'/note="([^"]+)"', quals))
        label = lab.group(1) if lab else kind
        nums = [int(x) for x in re.findall(r"\d+", loc)]
        spans = [(nums[i], nums[i + 1]) for i in range(0, len(nums) - 1, 2)]
        if spans:
            feats.append(Feature(label, spans, "complement" in loc))
    return Construct(name, seq, feats, circular)


def _parse_fasta(text: str, name_hint: str) -> Construct:
    blocks = [b for b in text.split(">") if b.strip()]
    if not blocks:
        raise ValueError("%s: empty FASTA" % name_hint)
    if len(blocks) > 1:
        raise ValueError(
            "%s: FASTA holds %d records. Supply one file per construct, so each "
            "one can be named and compared separately."
            % (name_hint, len(blocks)))
    header, *rest = blocks[0].splitlines()
    seq = "".join(rest).replace(" ", "").upper()
    return Construct(header.split()[0] if header.split() else name_hint, seq)


def load_construct(path: str, name: Optional[str] = None) -> Construct:
    with open(path, encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    hint = name or path.rsplit("/", 1)[-1].rsplit(".", 1)[0]
    c = (_parse_fasta(text, hint) if text.lstrip().startswith(">")
         else _parse_genbank(text, hint))
    if name:
        c.name = name
    if not c.seq:
        raise ValueError("%s: parsed zero bases" % path)
    if len(c.name) < 3:
        raise ValueError(
            "%s: construct name %r is too short to match safely against contig "
            "names. Use a distinctive name (the plasmid ID, e.g. pLV-EF1a-GFP)."
            % (path, c.name))
    return c


# --------------------------------------------------------------------------
# Junction tables — schema driven
# --------------------------------------------------------------------------

#: canonical field names used throughout the skill
CANONICAL = ["ref1", "pos1", "clip1", "filler", "ref2", "pos2", "clip2",
             "support", "flag"]

#: header spellings seen in the wild, mapped onto the canonical names
ALIASES: Dict[str, Sequence[str]] = {
    "ref1":    ("ref1", "ref_a", "chrom", "chromosome", "chr", "seqnames",
                "contig", "contig1", "reference1", "seq1", "genome", "host"),
    "pos1":    ("pos1", "position1", "pos_a", "start", "start_position",
                "breakpoint1", "bp1", "integration_site", "position",
                "genome_pos", "host_pos"),
    "clip1":   ("clip1", "clip_type_1", "cliptype1", "clip_a", "clip_side_1",
                "clip_side", "softclip1", "soft_clip_1", "side1"),
    "filler":  ("filler", "filler_length", "flank_length",
                "filler/flank_length", "microhomology", "mh",
                "microhomology_length", "insert_len", "untemplated"),
    "ref2":    ("ref2", "ref_b", "plasmid", "vector", "construct", "contig2",
                "reference2", "seq2", "partner", "transgene"),
    "pos2":    ("pos2", "position2", "pos_b", "end", "end_position",
                "breakpoint2", "bp2", "vector_pos", "plasmid_pos"),
    "clip2":   ("clip2", "clip_type_2", "cliptype2", "clip_b", "clip_side_2",
                "softclip2", "soft_clip_2", "side2"),
    "support": ("support", "supporting_read", "supporting_reads", "reads",
                "read_count", "readcount", "umi_count", "count", "n_reads",
                "split_reads", "n_split"),
    "flag":    ("flag", "status", "filter", "note", "comment", "qc"),
}

#: deliberately NOT aliased onto clip1/clip2 -- see resolve_clip_vocabulary
STRAND_HEADERS = {"strand", "strand1", "strand2", "orientation", "orientation1",
                  "orientation2", "sense", "direction"}


def _alias_map() -> Dict[str, str]:
    out: Dict[str, str] = {}
    for canon, names in ALIASES.items():
        for n in names:
            out.setdefault(norm(n), canon)
    return out


def parse_column_spec(spec: str) -> Dict[str, object]:
    """Parse `--columns ref1=chrom,pos1=3,support=reads`.

    Values are either a header name or a 1-based column index.
    """
    mapping: Dict[str, object] = {}
    for item in spec.split(","):
        item = item.strip()
        if not item:
            continue
        k, _, v = item.partition("=")
        k = k.strip()
        if k not in CANONICAL:
            raise ValueError("unknown field %r; expected one of %s" % (k, CANONICAL))
        v = v.strip()
        mapping[k] = int(v) - 1 if v.isdigit() else norm(v)
    return mapping


def parse_clip_map(spec: str) -> Dict[str, str]:
    """Parse `--clip-map right_clip=R,left_clip=L`."""
    out: Dict[str, str] = {}
    for item in spec.split(","):
        item = item.strip()
        if not item:
            continue
        k, _, v = item.partition("=")
        if not v:
            raise ValueError("--clip-map needs RAW=R or RAW=L, got %r" % item)
        out[k.strip()] = v.strip()
    return out


def _sniff(path: str) -> Tuple[str, List[List[str]]]:
    with open(path, encoding="utf-8", errors="replace", newline="") as fh:
        sample = fh.read(8192)
        fh.seek(0)
        delim = "\t"
        first = sample.split("\n")[0] if sample else ""
        if first and "\t" not in first:
            try:
                delim = csv.Sniffer().sniff(sample, delimiters="\t,;|").delimiter
            except csv.Error:
                delim = ","
        rows = [r for r in csv.reader(fh, delimiter=delim)
                if r and not (r[0].startswith("#") and len(r) < 4)]
    return delim, rows


def _resolve_columns(rows: List[List[str]], override: Dict[str, object]
                     ) -> Tuple[Dict[str, int], int, List[str]]:
    """Return {canonical: column index}, index of the first data row, warnings."""
    amap = _alias_map()
    header = rows[0]
    hnorm = [norm(h) for h in header]
    warnings: List[str] = []

    by_name: Dict[str, int] = {}
    for i, h in enumerate(hnorm):
        if h in amap and amap[h] not in by_name:
            by_name[amap[h]] = i

    strand_cols = [header[i] for i, h in enumerate(hnorm) if h in STRAND_HEADERS]
    if strand_cols and "clip1" not in by_name:
        warnings.append(
            "column(s) %s look like a strand, not a clip side, so they were NOT "
            "used for flank pairing. If your caller really encodes the clip side "
            "there, map it with --columns clip1=%s and declare the values with "
            "--clip-map." % (strand_cols, strand_cols[0]))

    # A header row names at least two canonical fields and its second cell is
    # not itself a number.
    second = (header[1].strip() if len(header) > 1 else "x") or "x"
    has_header = len(by_name) >= 2 and not re.fullmatch(r"-?\d+(\.\d+)?", second)

    cols: Dict[str, int] = {}
    if has_header:
        cols.update(by_name)
        first_data = 1
    else:
        for i, name in enumerate(CANONICAL):
            if i < len(header):
                cols[name] = i
        first_data = 0
        warnings.append(
            "no header recognised; columns read positionally as %s"
            % ", ".join(CANONICAL[:len(header)]))

    for k, v in override.items():
        if isinstance(v, int):
            cols[k] = v
        else:
            if v not in hnorm:
                raise ValueError("column %r not found in header %s" % (v, header))
            cols[k] = hnorm.index(v)
    return cols, first_data, warnings


@dataclass
class Junction:
    sample: str
    mapping: str
    ref_a: str = ""
    pos_a: int = 0
    clip_a: str = ""
    ref_b: str = ""
    pos_b: int = 0
    clip_b: str = ""
    filler: int = 0
    support: int = 0
    flag: str = ""
    kind: str = ""               # host-vector | vector-vector | vector-unmapped
    host_contig: str = ""
    host_pos: int = 0
    host_clip: str = ""
    vector: str = ""             # which construct the vector side is
    vector_pos: int = 0
    vector_clip: str = ""
    region: str = ""             # unique | shared | vector-end
    membership: str = ""         # which constructs carry this vector position
    feature: str = ""
    artefact: str = ""

    def as_dict(self) -> Dict:
        return asdict(self)


class NameMatcher:
    """Decide whether a reference name in a table denotes one of our constructs.

    Matching is exact on normalised names, plus declared aliases. Prefix
    matching is offered only as a last resort, only for names long enough to be
    distinctive, only when the match is unique, and always with a warning --
    because a construct called `A` would otherwise claim the GRCh38 contig
    `AC233755.1`, and the resulting misclassification is silent.
    """

    MIN_PREFIX_LEN = 5

    def __init__(self, names: Iterable[str],
                 aliases: Optional[Dict[str, str]] = None,
                 allow_prefix: bool = True):
        self.names = {norm(n): n for n in names}
        self.allow_prefix = allow_prefix
        self.warnings: List[str] = []
        self._cache: Dict[str, Optional[str]] = {}

        # Documented direction is CONSTRUCT=TABLENAME. A user who writes it the
        # other way round is told so instead of getting silence.
        self.alias: Dict[str, str] = {}
        for k, v in (aliases or {}).items():
            nk, nv = norm(k), norm(v)
            if nk in self.names:
                self.alias[nv] = self.names[nk]
            elif nv in self.names:
                self.alias[nk] = self.names[nv]
                self.warnings.append(
                    "alias %r reads TABLENAME=CONSTRUCT; the documented order is "
                    "CONSTRUCT=TABLENAME. Accepted as written." % (k + "=" + v))
            else:
                self.warnings.append(
                    "alias %r names no construct in the identity map (%s); ignored."
                    % (k + "=" + v, ", ".join(sorted(self.names.values()))))

    def match(self, ref: str) -> Optional[str]:
        r = norm(ref)
        if r in self._cache:
            return self._cache[r]
        hit = self._match(r, ref)
        self._cache[r] = hit
        return hit

    def _match(self, r: str, raw: str) -> Optional[str]:
        if not r:
            return None
        if r in self.alias:
            return self.alias[r]
        if r in self.names:
            return self.names[r]
        if not self.allow_prefix:
            return None
        hits = {orig for n, orig in self.names.items()
                if len(n) >= self.MIN_PREFIX_LEN
                and (r.startswith(n) or n.startswith(r))
                and len(r) >= self.MIN_PREFIX_LEN}
        if len(hits) == 1:
            hit = hits.pop()
            self.warnings.append(
                "reference %r matched construct %r by prefix, not exactly. "
                "Confirm this is right; declare it with --construct-alias "
                "%s=%s to silence this." % (raw, hit, hit, raw))
            return hit
        if len(hits) > 1:
            self.warnings.append(
                "reference %r is an ambiguous prefix of %s; treated as host "
                "sequence. Declare it with --construct-alias."
                % (raw, sorted(hits)))
        return None


def read_junction_table(
    path: str,
    sample: str,
    mapping: str,
    matcher: NameMatcher,
    min_support: int = 1,
    columns: Optional[Dict[str, object]] = None,
    unmapped_tokens: Iterable[str] = (),
    clip_map: Optional[Dict[str, str]] = None,
) -> Tuple[List[Junction], Dict[str, object]]:
    """Read any delimited junction table.

    Which side is vector and which is host is decided by `matcher`, so there is
    no assumption that host contigs are called `chr...`. A side matching no
    construct and no unmapped token is host sequence.

    Returns (junctions, info). `info` records how the columns and the clip
    vocabulary were resolved, plus every warning, so the caller can show the
    user exactly what was assumed.
    """
    delim, rows = _sniff(path)
    if not rows:
        return [], {"path": path, "rows_read": 0, "delimiter": delim,
                    "warnings": ["file is empty"], "clip": ClipResolution(),
                    "vector_side_found": False, "columns": {}}

    cols, first, warnings = _resolve_columns(rows, columns or {})

    missing = [f for f in ("ref1", "pos1", "ref2", "pos2", "support")
               if f not in cols]
    if missing:
        raise ValueError(
            "%s: could not identify column(s) %s. Header was %r. "
            "Map them explicitly with --columns." % (path, missing, rows[0]))

    unmapped = UNMAPPED_TOKENS | {norm(t) for t in unmapped_tokens}

    def cell(r: List[str], fieldname: str, default: str = "") -> str:
        i = cols.get(fieldname)
        return r[i].strip() if i is not None and i < len(r) else default

    data = rows[first:]

    clip_raw = [cell(r, "clip1") for r in data] + [cell(r, "clip2") for r in data]
    clip = resolve_clip_vocabulary(clip_raw, clip_map)
    if not clip.ok:
        warnings.append(
            "clip side unusable (%s). Flanks cannot be paired, so every site "
            "will be reported as single-flank and the site count loses its main "
            "consistency check. The analysis continues; the loss is recorded."
            % clip.reason)

    def as_int(v: str) -> int:
        v = v.strip()
        if not v or norm(v) in {"na", "n_a", "nan", "none", "null", "."}:
            return 0
        return int(round(float(v)))

    out: List[Junction] = []
    unparsable = 0
    below = 0
    for r in data:
        if not r or len(r) < 2:
            continue
        try:
            support = as_int(cell(r, "support", "0"))
            if support < min_support:
                below += 1
                continue
            j = Junction(
                sample=sample, mapping=mapping,
                ref_a=cell(r, "ref1"), pos_a=as_int(cell(r, "pos1", "0")),
                clip_a=clip.convert(cell(r, "clip1")),
                ref_b=cell(r, "ref2"), pos_b=as_int(cell(r, "pos2", "0")),
                clip_b=clip.convert(cell(r, "clip2")),
                filler=as_int(cell(r, "filler", "0")),
                support=support,
                flag=cell(r, "flag"),
            )
        except (ValueError, TypeError):
            unparsable += 1
            continue
        out.append(j)

    oriented = _orient(out, matcher, unmapped)
    vectors = sorted({j.vector for j in oriented if j.vector})

    info = {
        "path": path,
        "delimiter": delim,
        "header": first == 1,
        "columns": {k: (rows[0][v] if first == 1 and v < len(rows[0])
                        else "col%d" % (v + 1))
                    for k, v in sorted(cols.items(), key=lambda kv: kv[1])},
        "rows_read": len(data),
        "rows_unparsable": unparsable,
        "rows_below_threshold": below,
        "junctions": len(oriented),
        "clip": clip,
        "vector_names_seen": vectors,
        "vector_side_found": bool(vectors),
        "warnings": warnings + list(matcher.warnings),
    }
    return oriented, info


def _orient(recs: List[Junction], matcher: NameMatcher,
            unmapped: set) -> List[Junction]:
    """Put the host side first, label the junction kind, collapse duplicate
    rows keeping maximum support."""
    best: Dict[Tuple, Junction] = {}
    for r in recs:
        a_con, b_con = matcher.match(r.ref_a), matcher.match(r.ref_b)
        a_unm, b_unm = norm(r.ref_a) in unmapped, norm(r.ref_b) in unmapped
        a_host = not a_con and not a_unm
        b_host = not b_con and not b_unm

        # host side first
        if b_host and not a_host:
            r = Junction(sample=r.sample, mapping=r.mapping,
                         ref_a=r.ref_b, pos_a=r.pos_b, clip_a=r.clip_b,
                         ref_b=r.ref_a, pos_b=r.pos_a, clip_b=r.clip_a,
                         filler=r.filler, support=r.support, flag=r.flag)
            a_con, b_con = b_con, a_con
            a_unm, b_unm = b_unm, a_unm
            a_host, b_host = True, False

        if a_host and b_host:
            r.kind = "host-host"
            r.host_contig, r.host_pos, r.host_clip = r.ref_a, r.pos_a, r.clip_a
        elif a_host:
            r.kind = "host-vector"
            r.host_contig, r.host_pos, r.host_clip = r.ref_a, r.pos_a, r.clip_a
            r.vector, r.vector_pos, r.vector_clip = b_con or "", r.pos_b, r.clip_b
        elif a_con and b_con:
            r.kind = "vector-vector"
            r.vector, r.vector_pos, r.vector_clip = a_con, r.pos_a, r.clip_a
        else:
            r.kind = "vector-unmapped"
            if a_con:
                r.vector, r.vector_pos, r.vector_clip = a_con, r.pos_a, r.clip_a
            else:
                r.vector, r.vector_pos, r.vector_clip = b_con or "", r.pos_b, r.clip_b

        k = (r.sample, r.mapping, r.kind, r.ref_a, r.pos_a, r.ref_b, r.pos_b)
        if k not in best or r.support > best[k].support:
            best[k] = r
    return sorted(best.values(), key=lambda r: -r.support)


# --------------------------------------------------------------------------
# Tiny linear algebra: weighted least squares by normal equations
# --------------------------------------------------------------------------

def solve_wls(rows: Sequence[Sequence[float]], y: Sequence[float],
              w: Optional[Sequence[float]] = None,
              non_negative: bool = True) -> Optional[List[float]]:
    """Solve `rows · x ≈ y` in the weighted least-squares sense.

    Systems here have as many unknowns as constructs -- a handful -- so normal
    equations with Gaussian elimination are ample and keep the skill
    dependency-free. `non_negative` clamps negative solutions to zero and
    re-solves on the remaining free variables, since a copy number below zero
    is not a physical answer.

    Returns None when the system is singular (the block classes do not
    distinguish the constructs).
    """
    n = len(rows[0]) if rows else 0
    if n == 0 or len(rows) < n:
        return None
    w = list(w) if w is not None else [1.0] * len(rows)

    free = list(range(n))
    for _ in range(n + 1):
        k = len(free)
        ata = [[0.0] * k for _ in range(k)]
        atb = [0.0] * k
        for row, yi, wi in zip(rows, y, w):
            for i in range(k):
                ri = row[free[i]]
                if ri == 0.0:
                    continue
                atb[i] += wi * ri * yi
                for jj in range(k):
                    ata[i][jj] += wi * ri * row[free[jj]]
        sol = _gauss(ata, atb)
        if sol is None:
            return None
        if not non_negative or all(v >= -1e-9 for v in sol):
            out = [0.0] * n
            for i, idx in enumerate(free):
                out[idx] = max(0.0, sol[i])
            return out
        drop = min(range(k), key=lambda i: sol[i])
        free.pop(drop)
        if not free:
            return [0.0] * n
    return None


def _gauss(a: List[List[float]], b: List[float]) -> Optional[List[float]]:
    n = len(b)
    m = [row[:] + [b[i]] for i, row in enumerate(a)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(m[r][c]))
        if abs(m[p][c]) < 1e-12:
            return None
        m[c], m[p] = m[p], m[c]
        pv = m[c][c]
        for r in range(c + 1, n):
            f = m[r][c] / pv
            if f:
                for cc in range(c, n + 1):
                    m[r][cc] -= f * m[c][cc]
    x = [0.0] * n
    for r in range(n - 1, -1, -1):
        s = m[r][n] - sum(m[r][c] * x[c] for c in range(r + 1, n))
        x[r] = s / m[r][r]
    return x


# --------------------------------------------------------------------------
# Small helpers
# --------------------------------------------------------------------------

def read_depth_file(path: str, contig: Optional[str] = None) -> Dict[int, float]:
    """Read `samtools depth` output (contig, pos, depth) for one contig.

    When `contig` is None the single most frequent contig in the file is used,
    which is the common case: depth extracted over one vector reference.
    """
    per: Dict[str, Dict[int, float]] = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 3:
                continue
            try:
                pos, dep = int(f[1]), float(f[2])
            except ValueError:
                continue
            per.setdefault(f[0], {})[pos] = dep
    if not per:
        raise ValueError("%s: no usable rows; expected `samtools depth` output "
                         "(contig, position, depth)" % path)
    if contig is not None:
        for name, d in per.items():
            if norm(name) == norm(contig) or norm(contig) in norm(name) \
                    or norm(name) in norm(contig):
                return d
        raise ValueError("%s: contig %r not present; file holds %s"
                         % (path, contig, sorted(per)))
    return max(per.values(), key=len)


def dump_json(obj, path: str) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=2, default=str)


def load_json(path: str):
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def fmt_warnings(warnings: Sequence[str], prefix: str = "  warning: ") -> str:
    return "\n".join(prefix + w for w in warnings)
