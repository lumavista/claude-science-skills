# Sensitivity for donor-level pseudobulk DE

The sample-size question for single-cell and bulk transcriptomics of vaccinated
or infected subjects. It is here rather than in `scanpy` or `pydeseq2` because
it is a power calculation, and because every skill in that chain needs to cite
the same numbers.

**Donors are the replicates, not cells.** A per-cell test across conditions
treats thousands of cells from one donor as independent samples. It does not
produce a weak p-value, it produces a wrong one: in a simulated null with 8
donors, a per-cell Wilcoxon called roughly **1,390 of 2,000 genes** falsely
significant against **0–3** by donor pseudobulk — and donor pseudobulk is
itself mildly anti-conservative at 4 donors per arm. So aggregate to one
pseudobulk profile per donor per cell type, then test across donors.

## The arithmetic

With `k` donors per arm and a donor-level standard deviation `sigma_donor` on
the natural-log scale, the standard error of the arm difference in **log2**
units is

```
SE = sigma_donor * sqrt(2 / k) / ln(2)
```

A gene is detectable when `|log2FC| / SE` clears the multiple-testing
threshold. That threshold is a **range, not a number**: Benjamini–Hochberg
rejects at `p <= alpha * R / m` for `R` discoveries among `m` genes, so it sits
near Bonferroni when almost nothing is truly changed and relaxes as the number
of real signals grows. Quoting one power figure implies a precision the design
does not have.

Run it rather than re-deriving it:

```python
import sys; sys.path.insert(0, "scripts")
from power import donor_de_sensitivity

r = donor_de_sensitivity(sigma_donor=0.35, donors_per_arm=4, n_genes=2000)
r["se"]                                  # 0.357 log2 units
r["min_detectable_log2fc"]               # {'optimistic': 1.07, 'conservative': 1.50}
r["detectable"][1.0]["power_optimistic"] # 0.42
r["verdict"]
```

## What it gives, at sigma_donor = 0.35 and m ≈ 2000

Power envelope (conservative = BH≈Bonferroni, optimistic = BH with ~100
discoveries). The simulated values sit inside each envelope, which is why the
envelope rather than a point estimate:

| donors/arm | SE (log2) | log2FC 0.5 | log2FC 1.0 | log2FC 2.0 | simulated |
|---|---|---|---|---|---|
| 2 | 0.505 | 0.1–3% | 1.3–15% | 40–82% | ~2% at log2FC 1.0 |
| 4 | 0.357 | 0.2–6% | 8–42% | 92–99% | ~34% at 1.0, ~0% at 0.5 |
| 6 | 0.292 | 0.6–12% | 22–66% | ~100% | — |
| 8 | 0.252 | 1.3–20% | 40–82% | ~100% | — |
| 12 | 0.206 | 3.7–36% | 74–97% | ~100% | — |

**Read it as two rules.**

- **Fewer than 4 donors per arm cannot support gene-level conclusions at all.**
  Not "has low power" — cannot support them. Say so instead of reporting a
  gene list.
- **At 4 per arm, a gene that was not called is uninformative, not unchanged.**
  State which of the two you mean, every time. Never read a null result as
  evidence of no effect without this number in hand.

## Three things that make it worse than the table

1. **Pseudobulking per cell type splits the counts further.** The table is an
   optimistic bound for the most abundant population; sensitivity for any
   non-abundant cell type is worse. Report the donor count *per cell type*,
   not for the experiment.
2. **`sigma_donor` is an assumption.** 0.35 is a reasonable default for human
   PBMC pseudobulk, but it is the input that most changes the answer. Estimate
   it from pilot or public data of the same tissue and platform and say where
   it came from, as `references/effect_sizes.md` requires of any effect size.
3. **A responder cutoff derived from the same data it is then tested on is
   circular.** Pre-specify response groups. If you choose a threshold, say that
   you chose it and show the sensitivity to it.

## Where this is enforced

The two scripts that would otherwise emit a confident gene list from an
underpowered design refuse or warn, quoting these figures:

- `scanpy/scripts/pseudobulk.py` logs replicates per arm per stratum, refuses
  below 2 per arm (`--allow-underpowered` to override) and warns below 4.
- `pydeseq2/scripts/run_deseq2_analysis.py` prints samples per tested level,
  raises below 2, warns below 4, and records `samples_per_level` in the run
  manifest.

Related: `references/effect_sizes.md` for choosing a SESOI,
`references/simulation_based_power.md` when the design needs a mixed model
rather than this closed form.
