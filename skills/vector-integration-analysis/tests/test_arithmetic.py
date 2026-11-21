#!/usr/bin/env python3
"""Unit regressions for the Stage 3/4 arithmetic the shell suite does not reach.

Every case here is a wrong biological number that passed a fully green
353-assertion run, which is the point: the shell suite exercises the manifest
plumbing end to end, so an error in a formula survives it as long as the
pipeline still runs. These assert on the formulas directly.

Run standalone (`python3 tests/test_arithmetic.py`) or via run_tests.sh.
Prints one `PASS`/`FAIL` line per case and exits non-zero on any failure.
"""

from __future__ import annotations

import math
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "scripts"))

import copy_number  # noqa: E402
import reconcile_sites  # noqa: E402
import via_common  # noqa: E402

PASS = FAIL = 0


def check(label: str, condition: bool, detail: str = "") -> None:
    global PASS, FAIL
    if condition:
        PASS += 1
        print("  \033[32mPASS\033[0m %s" % label)
    else:
        FAIL += 1
        print("  \033[31mFAIL\033[0m %s%s" % (label, " -- " + detail if detail else ""))


def close(a, b, tol=1e-6) -> bool:
    return a is not None and abs(float(a) - float(b)) <= tol


# ---------------------------------------------------------------- ploidy units

def test_implied_sites_carries_ploidy() -> None:
    """implied_sites must be per CELL, like Stage 4's ceiling.

    `copies` is per haploid genome. Stage 4 caps sites at
    ceil(copies * ploidy / divisor). When Stage 3 computed copies / divisor it
    was a factor of `ploidy` low, so a correct, self-consistent diploid line --
    4 single-copy loci, 2.0 copies per haploid, 1 copy per site -- was
    reported as "EXCESS: 4 loci detected but the measured DNA only supports
    2.0 (0.5x). Loci are probably being double-counted", i.e. a correct result
    flagged as a 2x artefact with remediation pointing at flank pairing.
    """
    print("\n== implied_sites is in copies-per-cell units ==")
    copies_per_haploid, divisor, ploidy = 2.0, 1.0, 2

    # Stage 3's number, via the real code path.
    cfg = {"ploidy": ploidy, "mappings": {"pA": {}}}
    jrows = {"pA": [
        # Flank reads, no internal reads; integrase => copies_per_site fixed 1.
        {"sample": "s1", "kind": "host-vector", "support": 10,
         "vector_pos": 1, "artefact": ""},
        {"sample": "s1", "kind": "host-vector", "support": 10,
         "vector_pos": 4000, "artefact": ""},
    ]}
    arch = copy_number.architecture(cfg, jrows, "pA", "s1", copies_per_haploid,
                                    "integrase", False)
    implied = arch["implied_sites"]
    check("copies_per_site is 1 with no internal reads",
          close(arch["copies_per_site"], 1.0), "got %r" % arch["copies_per_site"])
    check("implied_sites = copies x ploidy / divisor (4, not 2)",
          close(implied, 4.0), "got %r" % implied)

    # Stage 4's ceiling, via the real code path. The two must agree.
    caps = reconcile_sites.ceilings({"pA": copies_per_haploid}, ["pA"],
                                    ploidy, divisor)
    check("Stage 3 implied_sites agrees with Stage 4's ceiling",
          caps["pA"] == math.ceil(implied), "cap %r vs implied %r" % (caps["pA"], implied))

    # And the verdict on a self-consistent line must not cry EXCESS.
    verdict = copy_number.reconcile(implied, 4)
    v = verdict["verdict"] if isinstance(verdict, dict) else str(verdict)
    check("4 real loci at 4 implied sites is not called an excess",
          "EXCESS" not in v.upper(), v)


# -------------------------------------------------------- missing depth bases

def test_missing_depth_positions_count_as_zero() -> None:
    """A position absent from a depth file is at zero depth, not absent from n.

    `samtools depth` without -a omits zero-coverage bases. Averaging only the
    present positions deletes the evidence of a deletion: 700 uncovered bp of a
    3000 bp class read 2.00 copies/haploid where the truth was 1.53.
    """
    print("\n== missing depth positions are zero, not excluded ==")
    ident = {"constructs": {"pA": {
        "length": 1000, "circular": False,
        "blocks": [{"membership": ["pA"], "start": 1, "end": 1000}],
    }}}
    # 600 of 1000 bases present at 10x; 400 genuinely uncovered and omitted.
    depth = {p: 10.0 for p in range(1, 601)}
    blocks = copy_number.block_depths(ident, "pA", depth, trim=0)
    check("one block returned", len(blocks) == 1, "got %d" % len(blocks))
    if blocks:
        b = blocks[0]
        check("mean averages over the full span (6.0, not 10.0)",
              close(b["mean"], 6.0, 1e-9), "got %r" % b["mean"])
        check("bases is the span, not the count of present positions",
              b["bases"] == 1000, "got %r" % b["bases"])
        check("observed_bases records how many were actually in the file",
              b.get("observed_bases") == 600, "got %r" % b.get("observed_bases"))


# ------------------------------------------------------- depth contig lookup

def test_depth_contig_lookup_is_not_a_loose_substring() -> None:
    """An exact contig name must win, and an ambiguous one must refuse.

    Bidirectional substring matching on dict order answered a request for
    `pAlphaX` with `pAlphaX_v2`'s 99x instead of its own 5x -- a 20x wrong
    copy number, silently. Near-miss names are routine here.
    """
    print("\n== depth-file contig selection ==")
    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, "depth.txt")
        with open(path, "w") as fh:
            for p in range(1, 4):
                fh.write("pAlphaX_v2\t%d\t99\n" % p)
            for p in range(1, 4):
                fh.write("pAlphaX\t%d\t5\n" % p)
        d = via_common.read_depth_file(path, "pAlphaX")
        got = sum(d.values()) / len(d)
        check("exact name wins over a longer contig listed first",
              close(got, 5.0), "got %r" % got)

        d2 = via_common.read_depth_file(path, "pAlphaX_v2")
        got2 = sum(d2.values()) / len(d2)
        check("the longer name still resolves to itself", close(got2, 99.0),
              "got %r" % got2)

        # Ambiguous: matches both by substring, neither exactly.
        try:
            via_common.read_depth_file(path, "pAlpha")
            check("an ambiguous contig request is refused", False,
                  "returned a value instead of raising")
        except ValueError as exc:
            check("an ambiguous contig request is refused",
                  "matches 2 contigs" in str(exc), str(exc))


def main() -> int:
    test_implied_sites_carries_ploidy()
    test_missing_depth_positions_count_as_zero()
    test_depth_contig_lookup_is_not_a_loose_substring()
    print("\n  test_arithmetic.py: %d passed, %d failed" % (PASS, FAIL))
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
