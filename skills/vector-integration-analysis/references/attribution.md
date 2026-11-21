# Site attribution — the priority order, and why the ceiling is hard

`scripts/reconcile_sites.py` is the authoritative site producer. It collapses
the per-mapping loci, applies retention, and attributes every called locus to a
construct or to nobody. **Report its output. Do not re-collapse, re-threshold
or re-attribute by hand**, and do not assemble a site list from the per-mapping
`loci_<mapping>.tsv` tables — those are per-mapping intermediates, and a site
in shared sequence appears in every one of them.

This file is the argument behind the procedure. Read it before changing a
parameter or defending an attribution.

## Why this is a script and not a judgement

On a real two-construct line, two different analysts working from the same
per-mapping tables produced two different wrong site results:

- one over-assigned shared-backbone loci to the **low-copy** construct, by
  reasoning from support ratios, giving it several times more sites than its
  measured copy number could source;
- the other never collapsed the two mapping runs, so every shared-backbone
  locus was listed **twice**, once under each construct, roughly doubling the
  site count.

Both failures have one cause: locus collapse, retention and attribution were
left to whoever was writing the report. They are now three deterministic steps
in one script, with every threshold a documented parameter. The same inputs
always give the same sites.

---

## 1. Collapse across mappings

One physical host–vector junction is reported once per mapping run that
detected it — same host coordinate, different vector coordinate. A locus seen
in both runs is **one event, never two**.

Entries from *different* mappings whose host intervals lie within `--window`
(default 500, the `classify_junctions --pair-window`) are linked, and each
linked group becomes one unified locus.

Entries from the *same* mapping are never linked to each other directly.
`classify_junctions.py` has already decided whether two host junctions are one
locus (a paired `R`+`L` flank) or two (two flanks of the same side), and
re-merging them here would silently deflate the site count. A group can still
end up holding two entries of one mapping by being linked through a third
mapping's entry; that is recorded on the locus (`same_mapping_merge`) and
printed in the table's note column, so it is visible rather than silent.

**Support is combined per sample as the maximum over entries, never the sum.**
The mappings report the same reads. Summing them inflates every shared-sequence
site in proportion to the number of mappings — which is exactly the kind of
quiet multiplication this stage exists to stop.

## 2. Retention

A unified locus is called only if

- total support ≥ `--min-reads` (default 5), **and**
- it is supported in ≥ `--min-samples` independently prepared samples
  (default 2).

Everything else goes to `subthreshold_sites.tsv`, which is the single-read and
noise tail. It is reported as a tail — never folded into the called set, and
never silently discarded either.

Two conditions, because they fail for different reasons. Read support separates
signal from noise within one library. The sample count asks whether the event
survived an independent library preparation, which is what a library chimera
does not do. A locus seen once, in one sample, at one read has the shape of a
chimera at a genuine position (`artefacts.md` §8).

With a single sample, `--min-samples 2` calls nothing. The script says so
loudly rather than reporting an empty site list. Lower it deliberately and
state that you did.

## 3. Attribution, in strict priority

| Priority | Tier | Evidence | Confidence | Capped by the ceiling? |
|---|---|---|---|---|
| a | `unique` | vector side in sequence only one construct carries | `certain` | no |
| b | `vector-end` | within `--end-tol` of exactly one construct's declared end | `inferred` | no |
| c | `copy-ceiling` | the measured copy number is what can source it | `inferred` | **yes** |
| d | — | two or more constructs remain admissible | `unresolved` | — |

Each locus is tested in that order and takes the first tier that fires.

### (a) unique — sequence evidence

If, in any mapping, the vector side of the junction sits in a singleton
membership class of the identity map, the locus belongs to that construct and
is `certain`. This is sequence evidence and nothing downstream overrules it.

### (b) vector-end — structural evidence

Otherwise, if the junction sits within `--end-tol` of exactly one construct's
declared end (`--vector-ends CONSTRUCT=a,b`), it is attributed to that
construct as `inferred`.

A construct's end can only be tested where that construct was mapped against,
because `vector_pos` is in each mapping's own coordinates. Constructs that
could not be tested are named in the locus's note rather than assumed away.
If the junction sits at the declared end of **two** constructs, the tier does
not fire and the locus falls through.

This tier is what lets a junction inside *shared* sequence still be attributed:
constructs linearized at different sites have different declared ends, so a
junction at one of them names its construct even though the sequence does not.

### (c) copy-ceiling — the hard cap

A locus with no unique and no vector-end signal can only be placed by asking
which construct's measured copy number can source it.

```
cap_X = ceil(copies_per_haploid_X × ploidy ÷ min_copies_per_site)
```

A candidate is **admissible** when it has a reported copy number and its cap is
not yet full. Among the admissible candidates, the locus goes to the
highest-copy one **only if it leads the next by at least `--cn-margin`**
(default 2.0, i.e. twice the copy number). Otherwise nobody gets it.

Loci are processed best-supported first, so cap consumption is deterministic.
Support decides only the *order* in which interchangeable loci consume
capacity; it never decides *which construct* a locus goes to.

Three properties of this step matter:

- **The cap is hard.** A construct whose cap is full is not admissible at all.
  A construct measured at N copies cannot be assigned more copy-ceiling sites
  than its ceiling allows, **whatever the support ratios look like**.
- **No support ratio and no cross-sample scaling enters it.** The script
  computes neither here, so neither can override the ceiling. If you use a
  scaling-type heuristic in your write-up, it may only break ties among
  constructs that are *all* admissible under the ceiling. It may never place a
  locus on a construct whose cap is full.
- **An unknown copy number is not an unlimited one.** A construct with no
  reported copy number gets no ceiling and is attributed nothing by this tier.
  The whole premise of the tier is that the copy number sources the site.

`--min-copies-per-site` is the fewest copies a real site is assumed to carry,
and it is the ceiling's divisor. The default 1.0 is the most permissive honest
value and the only one valid for integrase or transposase delivery. On a
concatemer line, where each site holds an array, pass the measured copies per
site from `copy_number.py`; that is what makes the ceiling binding rather than
nominal.

The script never takes that value from the copy-number stage by itself, because
changing the divisor changes every site range and that is a decision, not a
derivation. What it does instead, when the run manifest carries a measured
copies-per-site above 1 and the divisor is still 1, is raise the
`ceiling_divisor` check naming the measured value. Rerun with it, or keep 1 and
say in the report that the ceiling is nominal.

### (d) unresolved

Everything left is listed **once**, with its candidate constructs named, and
attributed to no construct. A locus never appears under two constructs.

`unresolved` is a result, not a failure. When two constructs have similar copy
numbers and a locus sits in sequence they share, away from either declared end,
the data genuinely cannot place it. Forcing a choice there is how a confident
wrong answer gets made. These loci are the validation priorities — one junction
PCR settles each one (`validation.md`).

### When the tiers contradict the ceiling

Tiers (a) and (b) are not capped. If sequence-level evidence gives a construct
more sites than its ceiling admits, the sites stand and the **contradiction is
reported** as a caveat. Suppressing a locus that sequence evidence places would
be hiding the disagreement inside a site count. Something else is then wrong:
the copy number is underestimated, the ploidy is wrong, or those sites carry
fewer copies each than `--min-copies-per-site` assumes.

---

## 4. Site counts

Per construct, one bounded range — not per sample, not per mapping:

```
floor = certain + vector-end loci
high  = cap_X  (the copy-number ceiling)
```

The floor counts only the tiers that do not depend on the copy number. The
upper bound is what the measured DNA can source. Copy-ceiling attributions sit
between the two: inside the range, but not part of the floor, because the
argument for them is arithmetic rather than sequence. The unresolved count is
stated separately, and never silently distributed over the constructs.

Report the range and the unresolved count together. A site count without its
uncertainty, or without saying which loci are certain and which inferred, is
not a result.

---

## 5. Parameters

Every threshold is a parameter with a documented default. None is tuned to any
dataset, and none encodes a construct, copy number, coordinate, chromosome or
sample name.

| Parameter | Default | What it changes |
|---|---|---|
| `--window` | 500 | host distance within which two mappings' entries are the same junction |
| `--min-reads` | 5 | support a called locus needs |
| `--min-samples` | 2 | independent samples a called locus needs |
| `--end-tol` | 3 | tolerance for "at a declared vector end" |
| `--ploidy` | from the copy-number JSON, else 2 with a caveat | scales every ceiling |
| `--min-copies-per-site` | 1.0 | the ceiling's divisor |
| `--cn-margin` | 2.0 | how far the top candidate must lead before the ceiling tier attributes |

Raising `--cn-margin` moves loci from `copy-ceiling` to `unresolved`: more
conservative, more bench follow-up. Lowering it does the opposite, and below
about 1.5 it starts attributing on copy-number differences that are inside the
copy-number estimate's own uncertainty. Check the per-construct spread in the
`copy_number.py` report (`quantification.md` §2) before trusting a small
margin.
