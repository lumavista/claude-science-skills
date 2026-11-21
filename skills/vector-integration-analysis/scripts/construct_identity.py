#!/usr/bin/env python3
"""Stage 1 -- build the construct identity map.

Compares any number of vector references and reports, for each one:

  * a block map labelled by **which constructs carry that sequence**, not
    merely whether it is shared. A block present in {A,B} and a block present
    in {A,B,C} behave differently, and the copy-number solve needs to tell
    them apart.
  * the boundaries between membership classes, which predict exactly where
    cross-mapping artefact junctions will appear
  * the attributable fraction: how much of each construct can be assigned on
    sequence alone
  * optionally, a coordinate correspondence, so a junction reported against
    one reference can be matched to the same event reported against another

Run this before looking at a single junction. Everything downstream depends on
knowing which positions can be attributed and which cannot.

Usage
-----
    python construct_identity.py A.gb B.gb C.gb --out identity.json --md map.md
    python construct_identity.py A.gb B.gb --correspond pA:4809
    python construct_identity.py A.gb B.gb C.gb --run run.json

Options
-------
    --k N          seed length for exact matching (default 40). Lower it for
                   diverged constructs; raise it to ignore short incidental
                   matches such as restriction and primer sites.
    --min-block N  blocks shorter than this are reported but flagged as too
                   short to cross-map (default 150, about one read length).
                   Set it to your read length.
"""
from __future__ import annotations

import argparse
import os
import sys
from typing import Dict, List, Optional, Sequence, Tuple  # noqa: F401

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from via_common import (Construct, dump_json, load_construct,  # noqa: E402
                        revcomp)
from via_run import (FAIL, PASS, WARN, Checks, RunError,  # noqa: E402
                     add_check_arguments, load_run)


def kmer_index(c: Construct, k: int) -> set:
    """Every k-mer of `c`, both strands. Circular constructs wrap."""
    s = c.seq + c.seq[:k - 1] if c.circular else c.seq
    idx = set()
    for strand in (s, revcomp(s)):
        for i in range(len(strand) - k + 1):
            idx.add(strand[i:i + k])
    return idx


def membership_mask(target: Construct, others: Sequence[Construct],
                    k: int) -> List[frozenset]:
    """For each base of `target`, the set of construct names carrying it.

    `target` is always a member of its own bases. A base belongs to another
    construct when it falls inside a k-mer of `target` that construct also
    contains -- the same criterion that decides whether a read spanning that
    base can cross-map.
    """
    n = target.length
    sets: List[set] = [{target.name} for _ in range(n)]
    s = target.seq + target.seq[:k - 1] if target.circular else target.seq
    for o in others:
        idx = kmer_index(o, k)
        for i in range(len(s) - k + 1):
            if s[i:i + k] in idx:
                for j in range(i, min(i + k, n)):
                    sets[j].add(o.name)
                if target.circular:
                    for j in range(i, i + k):
                        sets[j % n].add(o.name)
    return [frozenset(x) for x in sets]


def blocks_from_membership(mask: Sequence[frozenset]
                           ) -> List[Tuple[int, int, frozenset]]:
    """Maximal runs of constant membership, as 1-based inclusive spans."""
    out: List[Tuple[int, int, frozenset]] = []
    if not mask:
        return out
    start, cur = 0, mask[0]
    for i, m in enumerate(mask[1:], start=1):
        if m != cur:
            out.append((start + 1, i, cur))
            start, cur = i, m
    out.append((start + 1, len(mask), cur))
    return out


def correspond(pos: int, src: Construct, dst: Construct, window: int = 40
               ) -> Optional[Dict[str, object]]:
    """Map a coordinate on `src` onto `dst`.

    Returns every placement found, not just the first: duplicated elements (a
    poly(A) signal present twice in one vector, say) produce genuinely
    ambiguous answers, and silently taking the first is a real source of error.
    """
    lo = max(0, pos - window)
    probe = src.seq[lo:pos + window]
    if len(probe) < 20:
        return None
    offset = (pos - 1) - lo
    hay = dst.seq + dst.seq[:window * 2] if dst.circular else dst.seq
    for length in range(len(probe), 19, -1):
        for i in range(0, len(probe) - length + 1):
            frag = probe[i:i + length]
            hits: List[int] = []
            for rc, seq in ((False, frag), (True, revcomp(frag))):
                start = 0
                while True:
                    j = hay.find(seq, start)
                    if j < 0:
                        break
                    mapped = (j + length - (offset - i) if rc
                              else j + 1 + (offset - i))
                    if dst.circular:
                        mapped = ((mapped - 1) % dst.length) + 1
                    hits.append(-mapped if rc else mapped)
                    start = j + 1
            if hits:
                uniq = sorted(set(hits), key=abs)
                return {"matched_bp": length,
                        "fraction": round(length / len(probe), 3),
                        "positions": uniq,
                        "ambiguous": len(uniq) > 1}
    return None


def build(constructs: Sequence[Construct], k: int, min_block: int) -> Dict:
    names = [c.name for c in constructs]
    out: Dict = {"k": k, "min_block": min_block, "construct_order": names,
                 "constructs": {}}
    for c in constructs:
        others = [o for o in constructs if o.name != c.name]
        mask = membership_mask(c, others, k)
        raw = blocks_from_membership(mask)

        blocks = []
        for a, b, mem in raw:
            bp = b - a + 1
            blocks.append({
                "start": a, "end": b, "bp": bp,
                "membership": sorted(mem),
                "unique": len(mem) == 1,
                "cross_mappable": bp >= min_block and len(mem) > 1,
                "features": c.span_features(a, b),
            })

        boundaries = []
        for i in range(len(raw) - 1):
            a1, b1, m1 = raw[i]
            a2, b2, m2 = raw[i + 1]
            # Only an edge where cross-mappable sequence meets something else
            # predicts a cross-mapping artefact. An edge between two blocks
            # both too short to carry a read predicts nothing.
            if max(b1 - a1 + 1, b2 - a2 + 1) < min_block:
                continue
            leaving = sorted(m1 - m2)
            entering = sorted(m2 - m1)
            if not leaving and not entering:
                continue
            boundaries.append({
                "pos": b1,
                "next_pos": a2,
                "from": sorted(m1),
                "to": sorted(m2),
                "constructs_lost": leaving,
                "constructs_gained": entering,
                # Reads of these constructs clip here when mapped to `c`.
                "clipping_constructs": leaving or entering,
            })

        unique_bp = sum(b["bp"] for b in blocks if b["unique"])
        out["constructs"][c.name] = {
            "length": c.length,
            "circular": c.circular,
            "unique_bp": unique_bp,
            "attributable_pct": round(100 * unique_bp / c.length, 1),
            "blocks": blocks,
            "boundaries": boundaries,
            # class -> total bp, the design matrix for the copy-number solve
            "membership_classes": _classes(blocks),
        }
    return out


def _classes(blocks: Sequence[Dict]) -> List[Dict]:
    agg: Dict[Tuple[str, ...], Dict] = {}
    for b in blocks:
        key = tuple(b["membership"])
        d = agg.setdefault(key, {"membership": list(key), "bp": 0, "blocks": []})
        d["bp"] += b["bp"]
        d["blocks"].append([b["start"], b["end"]])
    return sorted(agg.values(), key=lambda d: (-d["bp"], d["membership"]))


def render_md(ident: Dict) -> str:
    names = ident["construct_order"]
    L = ["# Construct identity map", "",
         "Seed length %d bp, exact matching, both strands. Blocks shorter than "
         "%d bp are flagged as too short for a read to cross-map."
         % (ident["k"], ident["min_block"]), ""]
    for name in names:
        d = ident["constructs"][name]
        L += ["## %s (%s bp%s)" % (name, format(d["length"], ","),
                                   ", circular" if d["circular"] else ""),
              "",
              "**%s bp (%s%%) is unique to %s and therefore attributable on "
              "sequence alone.**" % (format(d["unique_bp"], ","),
                                     d["attributable_pct"], name),
              "",
              "| Range | bp | Carried by | Content | Attributable? |",
              "|---|---|---|---|---|"]
        for b in d["blocks"]:
            if b["bp"] < 30 and not b["unique"]:
                continue
            if b["unique"]:
                tag, rng = "**unique**", "**%s-%s**" % (format(b["start"], ","),
                                                        format(b["end"], ","))
            else:
                tag = "shared" if b["cross_mappable"] else "shared (too short to cross-map)"
                rng = "%s-%s" % (format(b["start"], ","), format(b["end"], ","))
            L.append("| %s | %s | %s | %s | %s |"
                     % (rng, format(b["bp"], ","), " + ".join(b["membership"]),
                        b["features"], tag))
        L.append("")
        if d["boundaries"]:
            L += ["Predicted cross-mapping artefact positions:", ""]
            for bd in d["boundaries"]:
                who = ", ".join(bd["clipping_constructs"]) or "?"
                L.append("- **%s** - reads of %s clip here (%s -> %s)"
                         % (format(bd["pos"], ","), who,
                            " + ".join(bd["from"]), " + ".join(bd["to"])))
            L.append("")
        else:
            L += ["No cross-mapping boundaries predicted: no shared block is "
                  "long enough for a read to cross out of.", ""]
        L += ["Membership classes (the design matrix for the copy-number solve):",
              "", "| Carried by | total bp |", "|---|---|"]
        for cl in d["membership_classes"]:
            L.append("| %s | %s |" % (" + ".join(cl["membership"]),
                                      format(cl["bp"], ",")))
        L.append("")

    if ident.get("correspondence"):
        L += ["## Coordinate correspondence", "",
              "| From | To | Position(s) | Matched bp | |", "|---|---|---|---|---|"]
        for c in ident["correspondence"]:
            r = c["result"]
            if isinstance(r, dict):
                L.append("| %s:%s | %s | %s | %d | %s |"
                         % (c["from"], c["pos"], c["to"],
                            ", ".join(str(p) for p in r["positions"]),
                            r["matched_bp"],
                            "**ambiguous**" if r["ambiguous"] else ""))
            else:
                L.append("| %s:%s | %s | no match | - | |"
                         % (c["from"], c["pos"], c["to"]))
        L.append("")
    return "\n".join(L)


#: a construct with less unique sequence than this is unattributable
MIN_ATTRIBUTABLE_PCT = 5.0


def stage_checks(checks: Checks, ident: Dict, n_refs: int) -> Checks:
    """What the map says about whether the line can be analysed at all."""
    thin = {n: d["attributable_pct"] for n, d in ident["constructs"].items()
            if d["attributable_pct"] < MIN_ATTRIBUTABLE_PCT}
    if thin:
        checks.add("attributable", FAIL,
                   "under %g%% of %s is unique sequence (%s). Junctions in "
                   "those constructs cannot be attributed on sequence alone "
                   "and the copy-number solve will be weakly determined: "
                   "check that every construct in the line has its own "
                   "reference and that --k suits how diverged they are, and "
                   "say so in the report if this is simply what the line is."
                   % (MIN_ATTRIBUTABLE_PCT, ", ".join(sorted(thin)),
                      ", ".join("%s %.1f%%" % kv for kv in sorted(thin.items()))))
    else:
        checks.add("attributable", PASS,
                   "every construct carries unique sequence: %s"
                   % ", ".join("%s %.1f%%" % (n, d["attributable_pct"])
                               for n, d in ident["constructs"].items()))

    if n_refs == 1:
        checks.add("constructs", WARN,
                   "one construct in the map. If the line carries more, add "
                   "their references -- including any the provider did not map "
                   "against -- before Stage 2, because their sequence is what "
                   "explains the artefacts their absence created.")
    else:
        checks.add("constructs", PASS, "%d constructs compared" % n_refs)
    return checks


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("refs", nargs="+", help="GenBank or FASTA references, one per construct")
    ap.add_argument("--name", action="append", default=[],
                    help="override construct name, in the same order as refs. "
                         "Use the real plasmid ID; one- or two-letter names are "
                         "rejected because they cannot be matched safely "
                         "against contig names.")
    ap.add_argument("--k", type=int, default=40)
    ap.add_argument("--min-block", type=int, default=150,
                    help="set this to your read length (default 150)")
    ap.add_argument("--circular", action="store_true",
                    help="treat every reference as circular")
    ap.add_argument("--out", help="write the identity map as JSON")
    ap.add_argument("--md", help="write a human-readable map as Markdown")
    ap.add_argument("--correspond", default="",
                    help="comma-separated NAME:POS to map onto every other construct")
    add_check_arguments(ap)
    a = ap.parse_args(argv)

    try:
        run = load_run(a.run)
    except RunError as exc:
        print("error: %s" % exc, file=sys.stderr)
        return 2
    if run is not None:
        a.out = run.output(a.out, "identity.json")
        a.md = run.output(a.md, "identity_map.md")

    if a.k < 15:
        print("error: --k below 15 produces chance matches in random sequence",
              file=sys.stderr)
        return 2

    constructs: List[Construct] = []
    for i, p in enumerate(a.refs):
        nm = a.name[i] if i < len(a.name) else None
        try:
            c = load_construct(p, nm)
        except (OSError, ValueError) as exc:
            print("error: %s" % exc, file=sys.stderr)
            return 2
        if a.circular:
            c.circular = True
        constructs.append(c)

    if len({c.name for c in constructs}) != len(constructs):
        print("error: construct names are not unique; set them with --name",
              file=sys.stderr)
        return 2

    total = sum(c.length for c in constructs)
    if total > 2_000_000:
        print("warning: %s bp of reference across %d constructs. The k-mer "
              "index is held in memory; this is sized for vectors, not for "
              "genomes." % (format(total, ","), len(constructs)),
              file=sys.stderr)

    ident = build(constructs, a.k, a.min_block)

    if a.correspond:
        ident["correspondence"] = []
        by = {c.name: c for c in constructs}
        for item in a.correspond.split(","):
            nm, _, ps = item.strip().partition(":")
            if nm not in by:
                print("error: unknown construct %r; have %s"
                      % (nm, ", ".join(by)), file=sys.stderr)
                return 2
            if not ps.strip().isdigit():
                print("error: --correspond needs NAME:POSITION, got %r" % item,
                      file=sys.stderr)
                return 2
            for other in constructs:
                if other.name == nm:
                    continue
                r = correspond(int(ps), by[nm], other)
                ident["correspondence"].append(
                    {"from": nm, "pos": int(ps), "to": other.name,
                     "result": r or "no match"})

    md = render_md(ident)
    if a.out:
        dump_json(ident, a.out)
    if a.md:
        with open(a.md, "w", encoding="utf-8") as fh:
            fh.write(md)
    if not a.out and not a.md:
        print(md)
    if run is not None:
        run.record("identity", {
            "path": run.relative(a.out), "md": run.relative(a.md),
            "constructs": {n: ident["constructs"][n]["length"]
                           for n in ident["construct_order"]}})

    checks = stage_checks(Checks("construct_identity", a.accept), ident,
                          len(constructs))
    if a.out:
        ident["checks"] = checks.as_list()
        dump_json(ident, a.out)
    if run is not None:
        run.record("identity", dict(run.stage("identity"),
                                    checks=checks.as_list()))

    for name, d in ident["constructs"].items():
        print("%s: %s bp | attributable %s bp (%s%%) | boundaries %s"
              % (name, format(d["length"], ","), format(d["unique_bp"], ","),
                 d["attributable_pct"],
                 ", ".join(str(b["pos"]) for b in d["boundaries"]) or "none"),
              file=sys.stderr)
        if d["attributable_pct"] < 5:
            print("  warning: under 5%% of %s is unique. Junctions in this "
                  "construct are essentially unattributable on sequence alone; "
                  "the copy-number solve will be weakly determined and the "
                  "report must say so." % name, file=sys.stderr)
    for c in ident.get("correspondence", []):  # noqa: E501
        r = c["result"]
        if isinstance(r, dict):
            print("  %s:%s -> %s:%s (%d bp)%s"
                  % (c["from"], c["pos"], c["to"],
                     ",".join(str(p) for p in r["positions"]), r["matched_bp"],
                     " AMBIGUOUS" if r["ambiguous"] else ""), file=sys.stderr)
        else:
            print("  %s:%s -> %s: no match" % (c["from"], c["pos"], c["to"]),
                  file=sys.stderr)
    return checks.finish(a.strict)


if __name__ == "__main__":
    raise SystemExit(main())
