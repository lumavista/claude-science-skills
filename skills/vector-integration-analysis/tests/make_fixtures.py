#!/usr/bin/env python3
"""Generate synthetic fixtures for the test suite.

Every sequence here is pseudo-random with a fixed seed. **No real vector,
plasmid or genome sequence is committed to this repository** — fixtures are
generated, not stored, so the tests can be public without distributing anyone's
construct sequences.

The three constructs deliberately share a backbone and a selection marker, so
the shared-sequence problem the skill exists to solve is present in the tests.

    python3 make_fixtures.py <output_dir>
"""
from __future__ import annotations

import os
import random
import sys

SEED = 20260101
TRUE_DEPTHS = {"pAlpha": 290.0, "pBeta": 116.0, "pGamma": 58.0}
GENOME_COVERAGE = 29.0           # -> true copies/haploid 10, 4, 2
DEDUP_PAIRS = 299_666_666        # gives C = 29.0 at read_length 150


def randseq(rng: random.Random, n: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def genbank(name: str, seq: str, feats) -> str:
    out = ["LOCUS       %-16s %d bp    DNA     linear   SYN" % (name, len(seq)),
           "FEATURES             Location/Qualifiers"]
    for label, (start, end) in feats:
        out.append("     misc_feature    %d..%d" % (start, end))
        out.append('                     /label="%s"' % label)
    out.append("ORIGIN")
    for i in range(0, len(seq), 60):
        chunk = seq[i:i + 60].lower()
        groups = " ".join(chunk[j:j + 10] for j in range(0, len(chunk), 10))
        out.append("%9d %s" % (i + 1, groups))
    out.append("//")
    return "\n".join(out) + "\n"


HEADER = ("chromosome\tposition1\tclip_type_1\tfiller_length\tref2\t"
          "position2\tclip_type_2\tsupporting_reads\tflag\n")


def write_table(path: str, rows, right: str, left: str) -> None:
    with open(path, "w") as fh:
        fh.write(HEADER)
        for contig, pos, clip, filler, ref2, pos2, clip2, support in rows:
            fh.write("\t".join(str(x) for x in (
                contig, pos, {"R": right, "L": left}.get(clip, clip), filler,
                ref2, pos2, {"R": right, "L": left}.get(clip2, clip2),
                support, "PASS")) + "\n")


def main() -> int:
    out = sys.argv[1] if len(sys.argv) > 1 else "fixtures"
    os.makedirs(out, exist_ok=True)
    rng = random.Random(SEED)

    backbone = randseq(rng, 3000)
    marker = randseq(rng, 800)
    cargo = {"pAlpha": randseq(rng, 1500),
             "pBeta": randseq(rng, 1600),
             "pGamma": randseq(rng, 1400)}

    seqs = {}
    for name, c in cargo.items():
        seqs[name] = backbone + c + marker
        feats = [("backbone", (1, 3000)),
                 ("cargo-%s" % name[1:], (3001, 3000 + len(c))),
                 ("selection marker", (3001 + len(c), len(seqs[name])))]
        with open(os.path.join(out, "%s.gb" % name), "w") as fh:
            fh.write(genbank(name, seqs[name], feats))

    alpha_len = len(seqs["pAlpha"])

    # ---- end-joining line: scattered vector breakpoints, a tandem array,
    #      one host contig that is an alt-contig accession, plus artefacts
    nhej = [
        ("chr1", 1_000_000, "R", 0, "pAlpha", 3200, "L", 30),
        ("chr1", 1_000_014, "L", 2, "pAlpha", 4100, "R", 26),
        ("chr7", 55_500_000, "R", 1, "pAlpha", 1500, "L", 18),
        ("chrX", 12_000_000, "R", 0, "pAlpha", 3900, "L", 9),
        ("AC233755.1", 41_000, "R", 0, "pAlpha", 3450, "L", 12),
        ("pAlpha", alpha_len - 1, "R", 0, "pAlpha", 2, "L", 400),
        ("pAlpha", 3001, "R", 110, "FlankingDNA", 0, "", 250),
        ("pAlpha", 2200, "R", -8, "pAlpha", 2280, "R", 4),
        ("pAlpha", 4000, "R", 0, "pAlpha", 4000, "R", 40),
    ]

    # ---- integrase line: every host junction at a vector terminus,
    #      5 bp target-site duplication, one off-end junction
    lenti = [
        ("chr1", 1_000_000, "R", 0, "pAlpha", 1, "L", 30),
        ("chr1", 999_996, "L", 0, "pAlpha", alpha_len, "R", 28),
        ("chr7", 55_500_000, "R", 0, "pAlpha", 1, "L", 22),
        ("chr7", 55_499_996, "L", 0, "pAlpha", alpha_len, "R", 20),
        ("chrX", 12_000_000, "R", 0, "pAlpha", 1, "L", 15),
        ("chr3", 7_000_000, "R", 0, "pAlpha", 1, "L", 11),
        ("chr3", 6_999_996, "L", 0, "pAlpha", alpha_len, "R", 10),
        # a host junction away from either terminus: a rearranged provirus,
        # backbone carryover, or an artefact -- never a canonical junction
        ("chr9", 3_000_000, "R", 0, "pAlpha", 2600, "L", 6),
        # a vector-vector junction, which single-copy delivery should not
        # produce at all; the accounting must query it rather than ignore it
        ("pAlpha", 1800, "R", -2, "pAlpha", 3600, "L", 7),
    ]

    # clip vocabularies the reader must cope with
    write_table(os.path.join(out, "junctions_nhej.tsv"), nhej,
                "right_clip", "left_clip")
    write_table(os.path.join(out, "junctions_nhej_altvocab.tsv"), nhej,
                "rightmost", "leftmost")
    write_table(os.path.join(out, "junctions_nhej_strand.tsv"), nhej, "+", "-")
    write_table(os.path.join(out, "junctions_nhej_noclip.tsv"), nhej, "", "")
    write_table(os.path.join(out, "junctions_nhej_exotic.tsv"), nhej,
                "clipped_at_end", "clipped_at_start")
    write_table(os.path.join(out, "junctions_integrase.tsv"), lenti,
                "right_clip", "left_clip")

    # second sample, support scaled, for the cross-sample path
    write_table(os.path.join(out, "junctions_nhej_s2.tsv"),
                [r[:-1] + (max(1, int(r[-1] * 0.6)),) for r in nhej],
                "right_clip", "left_clip")

    # ---- per-base depth consistent with TRUE_DEPTHS
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "..", "scripts"))
    import subprocess
    ident_path = os.path.join(out, "identity.json")
    script = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "..", "scripts", "construct_identity.py")
    subprocess.run([sys.executable, script,
                    os.path.join(out, "pAlpha.gb"),
                    os.path.join(out, "pBeta.gb"),
                    os.path.join(out, "pGamma.gb"),
                    "--out", ident_path], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    import json
    ident = json.load(open(ident_path))
    drng = random.Random(SEED + 1)
    for mapping, d in ident["constructs"].items():
        with open(os.path.join(out, "depth_%s.txt" % mapping), "w") as fh:
            for b in d["blocks"]:
                expected = sum(TRUE_DEPTHS[n] for n in b["membership"])
                for p in range(b["start"], b["end"] + 1):
                    val = max(0, round(drng.gauss(expected, expected ** 0.5)))
                    fh.write("%s\t%d\t%d\n" % (mapping, p, val))

    print("fixtures written to %s" % out)
    print("true copies per haploid: %s"
          % {k: round(v / GENOME_COVERAGE, 2) for k, v in TRUE_DEPTHS.items()})
    print("dedup_pairs for C=%.1fx at read_length 150: %d"
          % (GENOME_COVERAGE, DEDUP_PAIRS))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
