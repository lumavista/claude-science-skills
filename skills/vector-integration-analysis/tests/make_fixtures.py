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

# Concatemer-deficit case: two constructs sharing ~4.6 kb of backbone, one at
# low copy and one at high copy, with the shared blocks carrying 15% less depth
# than additivity predicts (backbone loss in concatemers). An all-block solve
# drives the low-copy construct negative; clamping that reported it as 0.
CONCAT_COPIES = {"pDelta": 1.5, "pEpsilon": 28.0}
CONCAT_SHARED_FACTOR = 0.85      # observed / (A + B) over shared blocks
CONCAT_DEDUP_PAIRS = 323_640_000  # read-count C = 31.32x, ~8% above the 29.0x
                                  # genome median: the read-count bias


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

    write_concatemer_fixture(out, ident_script=script)
    write_reconcile_fixture(out)
    write_spread_fixture(out, ident_script=script)

    print("fixtures written to %s" % out)
    print("true copies per haploid: %s"
          % {k: round(v / GENOME_COVERAGE, 2) for k, v in TRUE_DEPTHS.items()})
    print("dedup_pairs for C=%.1fx at read_length 150: %d"
          % (GENOME_COVERAGE, DEDUP_PAIRS))
    return 0


def write_concatemer_fixture(out: str, ident_script: str) -> None:
    """Two constructs, shared backbone under-represented by 15%."""
    import json
    import subprocess
    rng = random.Random(SEED + 2)
    backbone = randseq(rng, 3800)
    marker = randseq(rng, 800)                       # 4.6 kb shared in total
    cargo = {"pDelta": randseq(rng, 1200), "pEpsilon": randseq(rng, 2000)}
    for name, c in cargo.items():
        seq = backbone + c + marker
        feats = [("backbone", (1, 3800)),
                 ("cargo-%s" % name[1:], (3801, 3800 + len(c))),
                 ("selection marker", (3801 + len(c), len(seq)))]
        with open(os.path.join(out, "%s.gb" % name), "w") as fh:
            fh.write(genbank(name, seq, feats))
    ident_path = os.path.join(out, "identity_concat.json")
    subprocess.run([sys.executable, ident_script]
                   + [os.path.join(out, "%s.gb" % n) for n in cargo]
                   + ["--out", ident_path], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    ident = json.load(open(ident_path))
    true_depth = {n: c * GENOME_COVERAGE for n, c in CONCAT_COPIES.items()}
    drng = random.Random(SEED + 3)
    for mapping, d in ident["constructs"].items():
        with open(os.path.join(out, "depth_concat_%s.txt" % mapping), "w") as fh:
            for b in d["blocks"]:
                expected = sum(true_depth[n] for n in b["membership"])
                if len(b["membership"]) > 1:
                    expected *= CONCAT_SHARED_FACTOR
                for p in range(b["start"], b["end"] + 1):
                    val = max(0, round(drng.gauss(expected, expected ** 0.5)))
                    fh.write("%s\t%d\t%d\n" % (mapping, p, val))


# --------------------------------------------------------------------------
# reconcile_sites fixtures
# --------------------------------------------------------------------------
#
# Junction tables for the two constructs built by write_concatemer_fixture
# (pDelta 1.5 copies, pEpsilon 28 copies, sharing a 3,800 bp backbone and an
# 800 bp marker), mapped against each separately. Every feature the
# reconciliation has to get right is present:
#
#   * two loci in pEpsilon-unique cargo              -> certain, pEpsilon
#   * one locus at pDelta's declared end, in SHARED  -> vector-end, pDelta
#     backbone, reported by BOTH mappings
#   * five high-support shared-backbone loci, each   -> the ceiling decides
#     reported by BOTH mappings (so collapse must
#     halve them)
#   * a single-read, single-sample locus             -> sub-threshold tail
#
# The declared ends differ between the constructs because they were
# linearized at different sites -- which is what makes the vector-end tier
# able to tell them apart inside shared sequence.

RECON_ENDS = {"pDelta": 1200, "pEpsilon": 3000}

#: (contig, host_pos, vector_pos, support_s1, support_s2, in_pDelta_table,
#:  in_pEpsilon_table)
RECON_LOCI = [
    # pEpsilon-unique cargo (3801..5800 of pEpsilon): attributable on sequence
    ("chr2",  5_000_000, 4200, 40, 36, False, True),
    ("chr4",  8_000_000, 5200, 33, 30, False, True),
    # pDelta's linearization site, inside shared backbone: both mappings see
    # it, but only pDelta has a declared end there
    ("chr5",  2_000_000, RECON_ENDS["pDelta"], 25, 22, True, True),
    # shared backbone, away from either declared end: nothing but the copy
    # number can place these. Descending support, so which one falls off the
    # ceiling is deterministic.
    ("chr6",  1_000_000, 2000, 60, 55, True, True),
    ("chr8",  1_000_000, 2200, 50, 45, True, True),
    ("chr10", 1_000_000, 2400, 42, 38, True, True),
    ("chr11", 1_000_000, 2600, 31, 28, True, True),
    ("chr12", 1_000_000, 2800, 12, 10, True, True),
    # the noise tail: one read, one sample
    ("chr13", 1_000_000, 2500,  1,  0, False, True),
]


def write_reconcile_fixture(out: str) -> None:
    for which, mapping in enumerate(("pDelta", "pEpsilon")):
        for si, sample in enumerate(("s1", "s2")):
            rows = []
            for contig, hpos, vpos, s1, s2, in_delta, in_eps in RECON_LOCI:
                sup = (s1, s2)[si]
                if sup > 0 and (in_delta, in_eps)[which]:
                    rows.append((contig, hpos, "R", 0, mapping, vpos, "L", sup))
            write_table(os.path.join(out, "junctions_recon_%s_%s.tsv"
                                     % (mapping, sample)),
                        rows, "right_clip", "left_clip")


# --------------------------------------------------------------------------
# per-unique-block spread fixture (copy_number change C)
# --------------------------------------------------------------------------
#
# pTheta carries TWO separate unique blocks at DIFFERENT depths, so the
# bp-weighted mean hides a real disagreement unless the per-block estimates
# are reported. pIota carries ONE small unique block that survives trimming
# with 150 bp -- under the old `min_class_bp` of 200 it was silently dropped
# and the construct fell back to the flagged solve.

SPREAD_DEPTHS = {"pTheta_cargoA": 20.0, "pTheta_cargoB": 30.0,
                 "pIota_unique": 50.0, "shared": 75.0}


def write_spread_fixture(out: str, ident_script: str) -> None:
    import json
    import subprocess
    rng = random.Random(SEED + 4)
    backbone = randseq(rng, 2500)
    marker = randseq(rng, 800)
    spacer = randseq(rng, 600)
    cargo_a = randseq(rng, 1200)
    cargo_b = randseq(rng, 1200)
    iota_u = randseq(rng, 450)

    theta = backbone + cargo_a + marker + cargo_b + spacer
    iota = backbone + marker + spacer + iota_u
    for name, seq, feats in (
            ("pTheta", theta, [("backbone", (1, 2500)),
                               ("cargo-A", (2501, 3700)),
                               ("selection marker", (3701, 4500)),
                               ("cargo-B", (4501, 5700)),
                               ("spacer", (5701, 6300))]),
            ("pIota", iota, [("backbone", (1, 2500)),
                             ("selection marker", (2501, 3300)),
                             ("spacer", (3301, 3900)),
                             ("cargo-iota", (3901, 4350))])):
        with open(os.path.join(out, "%s.gb" % name), "w") as fh:
            fh.write(genbank(name, seq, feats))

    ident_path = os.path.join(out, "identity_spread.json")
    subprocess.run([sys.executable, ident_script,
                    os.path.join(out, "pTheta.gb"),
                    os.path.join(out, "pIota.gb"),
                    "--out", ident_path], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    ident = json.load(open(ident_path))
    drng = random.Random(SEED + 5)

    def expected(mapping: str, block) -> float:
        if len(block["membership"]) > 1:
            return SPREAD_DEPTHS["shared"]
        if mapping == "pIota":
            return SPREAD_DEPTHS["pIota_unique"]
        # pTheta: the first unique block is cargo-A, the second cargo-B
        return SPREAD_DEPTHS["pTheta_cargoA" if block["start"] < 4000
                             else "pTheta_cargoB"]

    for mapping, d in ident["constructs"].items():
        with open(os.path.join(out, "depth_spread_%s.txt" % mapping), "w") as fh:
            for b in d["blocks"]:
                e = expected(mapping, b)
                for p in range(b["start"], b["end"] + 1):
                    val = max(0, round(drng.gauss(e, e ** 0.5)))
                    fh.write("%s\t%d\t%d\n" % (mapping, p, val))


if __name__ == "__main__":
    raise SystemExit(main())
