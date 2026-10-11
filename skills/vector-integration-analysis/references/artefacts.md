# Artefact catalogue

Junction callers report everything that looks like a breakpoint. Most of what
they report in a multi-construct stable line is not a rearrangement. Classify
before counting.

`classify_junctions.py` flags the classes below that are marked with a flag
name. **The flags are starting points.** Confirm each one by the stated test
before accepting it, and override when the test fails.

Thresholds quoted are the script's defaults, scaled from `--read-length`. They
are starting points too.

---

## 1. Boundary cross-mapping — flag `boundary-crossmap`

**What it is.** The mapping reference contained one construct but the cells
carry two or more. A read from an absent construct aligns across shared
sequence, then runs into sequence unique to its true source, which matches
nothing. The aligner clips it. Many such reads clip at the same base, and the
caller reports a junction to unknown DNA (`FlankingDNA`, `unknown`, `NA`, `.` —
callers differ).

**Signature.**
- vector position sits at a predicted membership boundary (`identity.json`
  names which construct's reads clip there)
- filler ≈ read length − aligner anchor
- support can be very large, in proportion to the absent construct's copy
  number

**Confirming test.** Take the clipped consensus sequence and the vector-side
segment, and ask whether they are **contiguous in the other construct**. If
yes, the read is an ordinary unrearranged molecule of that construct and the
junction is pure artefact. If no, something genuinely joined two distant parts
of the other construct — a real rearrangement, but of a different vector. Both
outcomes remove it from the current construct's evidence; only the second is
biology.

Most callers emit the junction consensus or the clipped sequence alongside the
table. This test is only possible with the sequence in hand, so ask for it if
the deliverable omits it.

**Why this matters quantitatively.** These junctions are the basis of the
fallback copy-number estimator in `quantification.md §3`. Do not simply discard
them.

---

## 2. Long filler without a boundary — flag `long-filler`

Same mechanism, but the clip point is inside a shared block rather than at its
edge — usually because the two constructs diverge internally, or a third
sequence is present.

**Signature.** Filler at or above the cross-mapping threshold, position not
near a predicted boundary.

**Action.** Recover the consensus and identify it. If it matches a construct
you know about, update the identity map and rerun. If it matches nothing you
expected, you may have an unlisted construct, a helper plasmid or bacterial
backbone carried over — flag it in the report rather than ignoring it. A
handful of unexplained junctions with partial identity to a selection marker is
a real finding and worth a database search.

---

## 3. Short-range inverted junctions — flag `short-range-inverted`

**What it is.** Fold-back structures or library chimeras.

**Signature.** Two vector positions tens to a couple of hundred bases apart,
**the same clip side on both sides**, several bases of apparent microhomology
(negative filler), support at a small fraction of local depth.

**Confirming test.** They appear in the deeper sample and not the shallower
one, at the same positions, because doubling depth pushes them over the
reporting threshold. If the same profile appears when the reads are mapped
against a *different* reference, that confirms the diagnosis.

**Do not confuse with a self-junction.** A vector position joined to itself
(span ~0, same clip side) is a palindromic head-to-head or tail-to-tail
concatemer junction and is **real**. The artefact class has a short but
non-zero span, and the script enforces a lower span bound for exactly this
reason.

---

## 4. Vector-end religation — flag `vector-end-religation`

**What it is.** A junction joining the two declared ends of the delivered
molecule, or the two ends of a linear reference.

**Why it is ambiguous — three sources, indistinguishable by short reads:**
1. a genuine head-to-tail concatemer formed by end joining before integration;
2. uncut or religated input molecule still present in the prep;
3. **any co-resident construct whose equivalent site is intact** — if construct
   A was linearized at one site and construct B at a different one, every B
   molecule carries A's cut site contiguously and generates this junction when
   mapped to a linear A reference.

Source 3 is the one people miss, and it can dominate: a junction carrying
hundreds of reads can be entirely explained by a co-resident construct.

**Confirming test.** Reconstruct the junction consensus. If it equals the
*native* circular sequence, it is uninformative. A junction that has **lost
bases** relative to native — ends resected — cannot come from an intact
molecule and is a genuine end-joining concatemer junction. That distinction is
the whole game.

**Second test.** Map the reads against the other construct. If the junction
disappears, it was native sequence there, confirming it carried no information.

**Default:** excluded from copies-per-site estimates. State which you did.

**Modality note.** For integrase- and transposase-mediated delivery this class
should not appear at all, and the script marks it
`vector-end-religation(unexpected-for-modality)`. Its presence means
unintegrated vector, carried-over producer or donor plasmid, or a rearranged
insert.

---

## 5. Reference-edge clipping — flag `reference-edge`

Reads pile up clipping at position 1 and at the last base of any **linear**
reference, simply because the reference ends. Breakpoint tables often show
these as the highest-support entries in the dataset. They are not evidence of
structure.

**Important exception.** For a provirus or transposon, the LTR or ITR termini
*are* the genuine integration junctions, and on a linear reference they sit at
exactly these coordinates. This is why `--vector-ends` must be supplied: the
script checks declared ends before applying this rule, so real integrase
junctions are not thrown away as edge artefacts. Without `--vector-ends`,
every real junction in a lentiviral line will be flagged here.

---

## 6. Off-end host junctions — flag `off-end-host-junction`

**Only raised for fixed-end modalities.** For integrase- and
transposase-mediated delivery the vector breakpoint is not free: it is at a
terminus. A host–vector junction elsewhere in the vector is therefore not a
canonical integration junction.

**What it can be, in order of likelihood:**
- a rearranged or internally deleted insert;
- plasmid backbone that integrated by end joining alongside the intended event
  — the carryover check in `modalities.md`;
- a library chimera;
- the reference being the transfer plasmid rather than the proviral form, which
  shifts the expected junction a couple of bases and makes every real junction
  look off-end. Check this first: if *all* junctions are off-end by the same
  small offset, the reference is wrong, not the line.

---

## 7. Depth-driven threshold effects — not a junction class

Not an artefact in a single junction, but an artefact in the *comparison*.

A sample sequenced twice as deep reports far more junctions at a fixed support
cutoff, because low-frequency events cross it. A report saying one sample has
roughly twice the junctions of another is usually reporting sensitivity, not
biology.

**Test.** Compare **host–vector** junction counts between samples, not total
counts, and check how many of the extra calls fall below threshold in the
shallower sample. If most do, it is sensitivity.

---

## 8. Library chimeras at genuine positions — the one to worry about

With vector at many times the concentration of a single-copy locus, accidental
ligation of vector to genomic fragments during library prep is a real
mechanism. This is the artefact that would produce *false integration sites*,
and no flag can catch it, because such a junction looks structurally normal.

**Why genuine junctions survive the challenge — apply all five:**

- **Support after deduplication.** Support in the tens means that many
  independently captured fragments. A ligation chimera is a one-off accident,
  and dozens of accidents landing on the same host base *and* the same vector
  base is not credible.
- **Reproduction at the identical base pair** across independently prepared
  libraries.
- **Reproduction at the coordinates predicted by sequence correspondence** when
  reads are mapped against a different reference.
- **Junction microhomology** consistent with end joining, rather than the blunt
  profile of blunt-end ligation.
- **Target-site footprint matching the modality.** The strongest of the five
  where it applies: a chimera has no reason to carry a 4 or 5 bp target-site
  duplication of exactly the enzyme's signature. See `modalities.md`.

A junction that passes these is real. One supported by a handful of reads in a
single library is not, whatever its position.
