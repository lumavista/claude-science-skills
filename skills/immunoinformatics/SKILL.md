---
name: immunoinformatics
description: Predicts and triages T-cell and B-cell epitopes for vaccine design — MHC class I and II binding (NetMHCpan, NetMHCIIpan, MHCflurry, MixMHCpred), linear and conformational B-cell epitopes (BepiPred, DiscoTope, SEMA), IEDB validated epitopes, HLA population coverage, human-proteome self-similarity screening, and mosaic or cocktail immunogen design with junctional neo-epitope checks.
license: MIT
compatibility: Python 3.10+. MHCflurry 2.1 and the IEDB web APIs need network access and run anywhere. NetMHCpan-4.1, NetMHCIIpan-4.3, BepiPred-3.0 and DiscoTope-3.0 are DTU academic-licence tools installed separately; scripts use them if present, else fall back. Conformational prediction needs a structure.
metadata:
  version: "1.1"
  skill-author: Drafted for vaccine-design subset
  last-reviewed: "2026-10-16"
---

# Immunoinformatics

Takes antigen sequences to a ranked, filtered, population-aware epitope
shortlist, and assesses how well a proposed immunogen covers a target
population and the variants actually circulating.

Keep five kinds of result separate and never substitute one for another:
**predicted binder** (a model score; most are not presented), **presented
peptide** (seen by immunopeptidomics; presentation ≠ recognition), **epitope**
(experimentally recognised by a TCR or antibody — IEDB has these, predictors
cannot produce them), **immunogenic epitope** (elicits a response in a given
construct, dose, adjuvant and host), **protective epitope** (correlates with
protection in a challenge or efficacy study). Binding prediction is the weakest
and where almost all computational output sits, so the ranked list is a
screening shortlist, not a design. Class II prediction is substantially less
accurate than class I (open groove, shifting registers) — do not present class
II rankings with equal confidence.

## Install

```bash
pip install "mhcflurry>=2.1" "biopython>=1.83" "pandas>=2.0" "requests>=2.31"
mhcflurry-downloads fetch models_class1_presentation
```

MixMHCpred and the DTU Health Tech tools (NetMHCpan, NetMHCIIpan, BepiPred-3.0,
DiscoTope-3.0, from services.healthtech.dtu.dk) are optional, installed
separately, and detected by executable name on PATH. **They are academic-only,
and MHCflurry (Apache 2.0) is the permissive class I substitute with no class II
equivalent**, so read `references/tool-licensing.md` — which carries every
tool's exact version and licence — before wiring any into a published or
productised pipeline.

## MHC prediction backends and fallback routing

[scripts/predict_mhc_with_fallback.py](scripts/predict_mhc_with_fallback.py)
provides intelligent routing for MHC Class I prediction. **Never simulate MHC
predictions** — if no backend is available, fail with clear instructions instead.

Route order (first available wins):
1. **Local MHCflurry** (if model successfully loads) — works offline, no dependencies
2. **NetMHCpan** (if installed and on PATH) — higher accuracy, requires license
3. **Fail with clear instructions** (no synthetic predictions)

Example usage:

```python
import sys; sys.path.insert(0, "scripts")
from predict_mhc_with_fallback import predict_mhc_class1_with_fallback, MHCPredictionError

peptides = ["FVFSTSFYL", "YVFSTSFYL", "MTMPLSCTK"]
alleles = ["HLA-A*02:01", "HLA-A*11:01"]

try:
    results = predict_mhc_class1_with_fallback(
        peptides,
        alleles,
        include_percentile=True,
    )
    
    for r in results:
        print(f"{r['peptide']} → {r['allele']}: "
              f"{r['percentile_rank']:.3f}% ({r['backend']})")

except MHCPredictionError as e:
    print(f"MHC prediction not available:\n{e}")
    sys.exit(1)
```

Returns `{"peptide": "...", "allele": "...", "percentile_rank": ..., "backend": "mhcflurry"|"netmhcpan"}`.

**Why no simulation:** Synthetic MHC predictions waste experimental resources.
If the model cannot load or download (common in sandboxes), fail loudly and
report what is needed. This prevents silent errors that appear correct until
validation.

## Workflow

Run every step through `scripts/epitope_triage.py`, which encodes the ordering,
column conventions and safety filter and records tool versions in a
`.meta.json` sidecar. If a subcommand fails, report the error rather than
reimplementing it ad hoc.

### 1. Establish the antigen set

Start from actual protein sequences with provenance — strain or isolate,
accession, date — and resolve signal peptides, cleavage products and polyprotein
boundaries first, since a peptide spanning a cleavage site cannot be presented
from the mature protein.

**Never write a sequence from memory.** Every residue must come from a retrieved
record with an accession, because every downstream step runs cleanly on a wrong
sequence and the error is invisible in the output; if retrieval fails, stop,
report which database was queried and what came back, and ask for an accession
or FASTA rather than writing out a "canonical" sequence. Before accepting one,
state that it is the protein asked for (record name and organism), that the
length is right, and that the quoted accession resolves.

```bash
python scripts/epitope_triage.py prepare \
  --fasta antigens.fasta --out work/antigens.clean.fasta
```

For a drifting pathogen, pull the currently circulating variants first
(`pathogen-variant-surveillance`) and carry them through to step 6 — designing
against a reference that no longer circulates is the most common avoidable error
in this workflow.

### 2. MHC class I binding

```bash
python scripts/epitope_triage.py mhc1 \
  --fasta work/antigens.clean.fasta \
  --alleles HLA-A*01:01,HLA-A*02:01,HLA-A*03:01,HLA-A*24:02,HLA-B*07:02,HLA-B*08:01,HLA-B*44:02 \
  --lengths 8,9,10,11 \
  --out work/mhc1.tsv
```

That supertype panel covers a large fraction of most populations; for a specific
region, ancestry or cohort use its allele frequencies instead (`afnd.net`).

Rank on **percentile rank**, not affinity — raw nM values are not comparable
across alleles because some bind more promiscuously — and treat the ≤0.5%
strong / ≤2% weak cutoffs as convention balancing sensitivity against list
length, not biology. Prefer a presentation score (MHCflurry
`presentation_score`, NetMHCpan `%Rank_EL`), since eluted-ligand models are
trained on immunopeptidomics and reflect processing.

### 3. MHC class II binding

```bash
python scripts/epitope_triage.py mhc2 \
  --fasta work/antigens.clean.fasta \
  --alleles HLA-DRB1*01:01,HLA-DRB1*03:01,HLA-DRB1*04:01,HLA-DRB1*07:01,HLA-DRB1*15:01 \
  --out work/mhc2.tsv
```

Requires NetMHCIIpan; the script reports what is missing rather than silently
substituting. Class II drives B-cell help and affinity maturation, so a
construct with good class I coverage and no class II content often disappoints.
The permissive fallback is MHCnuggets (Apache-2.0, class I and II, unmaintained
since 2023) — name it and state its limits rather than claiming no open option
exists.

### 4. Cross-check IEDB

```bash
python scripts/epitope_triage.py iedb \
  --organism "<taxon id or name>" --out work/iedb_known.tsv
```

Do this **before ranking anything**. A validated epitope outranks any
prediction, and a high-ranking prediction in a well-studied antigen with no IEDB
support deserves scepticism because it has probably already been looked for.

### 5. B-cell epitopes

Antibody epitopes are overwhelmingly conformational, so use linear predictors
only to triage and prefer structure-based prediction wherever a structure
exists or can be predicted.

- **Linear** — BepiPred-3.0, read as a smoothed propensity track, not epitopes.
- **Conformational** — DiscoTope-3.0 or SEMA on an experimental or
  high-confidence predicted structure; check per-residue pLDDT, since a
  low-confidence loop gives a meaningless surface.
- **Accessibility** — compute RSA and exclude buried residues regardless of
  score.

For a glycosylated viral surface protein run `glycoengineering` in parallel: a
site under a dense glycan shield can score well and be antibody-inaccessible,
the dominant failure mode here.

### 6. Conservation across circulating variants

**Use `conservation-mapping` whenever you have or can build an alignment** — it
applies Henikoff weighting, reports effective sample size and reconciles
alignment positions to PDB numbering, none of which the fallback can do. Align
with `phylogenetics`, then run its `conservation_score.py`, join candidates on
`ref_pos` for a per-residue profile across each peptide, and report the
effective sample size it prints rather than the raw sequence count.

The alignment-free fallback, where no alignment is practical:

```bash
python scripts/epitope_triage.py conservation \
  --candidates work/mhc1.tsv --variants work/variant_sequences.fasta \
  --out work/conserved.tsv
```

It counts exact peptide occurrences in unaligned sequences without correcting
for redundant sampling, so a clade-dominated panel reports clade-specific
peptides as highly conserved; if you use it, say so and report panel
composition. Either way, check *why* a region is conserved — functional
constraint is a good reason, shallow sampling is not, and the fraction alone
cannot distinguish them.

### 7. Self-similarity and allergenicity

A safety filter, not optional for anything intended to enter a human.

```bash
python scripts/epitope_triage.py selfcheck \
  --candidates work/conserved.tsv \
  --human-proteome UP000005640_9606.fasta \
  --out work/flagged.tsv
```

Flag and review any candidate with close identity to a human protein: exact or
near-exact self-matches carry a theoretical autoimmunity risk and must be
dropped or escalated to someone qualified, not silently ranked. Screen
separately against known allergens (AllergenOnline, or the FAO/WHO
35%-identity-over-80-residues criterion).

### 8. Population coverage

```bash
python scripts/epitope_triage.py coverage \
  --candidates work/flagged.tsv --population "World" --out work/coverage.tsv
```

Calls IEDB's population coverage tool. Report the fraction expected to present
at least one candidate **per region** as well as globally — a construct at 95%
world coverage can sit far lower in a specific population, which is exactly the
gap that matters for deployment.

### 9. Multivalent construct, if breadth demands one

When no single antigen covers the panel, use a cocktail of natural sequences or
a synthetic mosaic set; both are k-mer coverage optimisations in
`mosaic_design.py`.

```bash
python scripts/mosaic_design.py cocktail --fasta work/panel.fasta \
  -n 4 -k 9 --out work/cocktail
python scripts/mosaic_design.py mosaic --fasta work/panel.fasta \
  -n 2 -k 9 --beam 200 --out work/mosaic
python scripts/mosaic_design.py evaluate --fasta work/panel.fasta \
  --designs work/my_construct.fasta -k 9 --out work/construct_eval
```

**Mosaic mode cannot create junctional neo-epitopes** — it walks a de Bruijn
graph built from the panel, so every output k-mer occurs in some real strain by
construction, verified by the script, which fails loudly otherwise. **Check any
other construct with `evaluate`**, since naive concatenation, string-of-beads
and fusion designs all create junction-spanning peptides that exist in no real
strain; run it at k=9 and again at k=15 for class I and class II registers.

**Read the coverage output honestly.** The aggregate is the misleading number —
a design can reach 90% weighted coverage while an entire clade sits near zero,
because weighted optimisation correctly ignores a sparsely sampled clade. Report
the minimum and 10th percentile, name the uncovered strains, and say whether the
gap was deliberate; `--unweighted` trades common-variant depth for rare-variant
breadth.

Cocktails are foldable and processable, being real proteins; mosaics usually are
not — fine for mRNA or peptide delivery, not for anything that must present a
conformational surface. Carry that into `construct-assembly`, whose QC differs
depending on whether a construct must fold.

### 10. Report

Per candidate: sequence, source antigen and position, restricting alleles with
percentile ranks, IEDB support, conservation fraction, self-similarity flag,
accessibility or structural context. Record tool versions and the allele panel —
predictor versions change rankings, so a result without its version is not
reproducible. Use `scientific-writing` for the write-up, stating in methods that
candidates are computational predictions pending experimental validation.

## Interpretation

- **Rank ordering is not calibrated** — rank 1 and rank 20 are not reliably
  distinguishable, so carry tens of candidates into screening.
- **Negative predictions are weak**; absence from a shortlist is not evidence a
  peptide is not an epitope, particularly for class II.
- **Allele coverage is uneven** — predictors are trained disproportionately on
  common European-ancestry alleles, and under-reported performance on
  under-represented alleles propagates into coverage inequity.
- **Processing is modelled crudely** — proteasomal cleavage and TAP predictors
  add little over a good eluted-ligand model, so do not stack them as
  independent filters.
- **A construct is not a list.** Epitope order, linkers, flanking residues and
  delivery format all change what is processed and presented; check
  junction-spanning peptides with `mosaic_design.py evaluate` (step 9).

## Reference files

- `references/tool-licensing.md` — licence and current version for every
  predictor, folding model and design tool named here, including
  code-vs-weights splits. It records NetMHCpan as 4.2e, superseding 4.1 above.
- `references/compute-routes.md` — which steps run on local CPU, local GPU, or
  Tamarind/a hosted route. Prediction and triage are CPU work and the GPU tier
  is folding, embedding, inverse folding and MD, so on an 8 GiB card most of
  this workflow is local.

## Official references

- IEDB — https://www.iedb.org/ ; tools and API at https://nextgen-tools.iedb.org/docs/
- MHCflurry 2.0 — O'Donnell et al., *Cell Systems* 2020; https://openvax.github.io/mhcflurry/
- NetMHCpan-4.1 / NetMHCIIpan-4.0 — Reynisson et al., *NAR* 2020
- BepiPred-3.0 — Clifford et al., *Protein Science* 2022
- DiscoTope-3.0 — Høie et al., *Front. Immunol.* 2024
- Allele Frequency Net Database — http://www.allelefrequencies.net/
- Human reference proteome — https://www.uniprot.org/proteomes/UP000005640
