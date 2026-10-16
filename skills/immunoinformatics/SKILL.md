---
name: immunoinformatics
description: Predicts and triages T-cell and B-cell epitopes for protective vaccine design. Covers MHC class I and II binding prediction (NetMHCpan, NetMHCIIpan, MHCflurry, MixMHCpred), linear and conformational B-cell epitope prediction (BepiPred, DiscoTope, SEMA), retrieval of experimentally validated epitopes from IEDB, HLA population coverage analysis, conservation across circulating variants, self-similarity screening against the human proteome, and mosaic or cocktail construct design with junctional neo-epitope checking. Use for immunogen candidate selection, epitope-based construct design, and coverage assessment across target populations. For per-position conservation scoring and mapping onto structure, this skill defers to conservation-mapping; for linkers, tags and construct developability it defers to construct-assembly.
license: MIT
compatibility: Python 3.10+. MHCflurry 2.1 and the IEDB web APIs need network access and run anywhere. NetMHCpan-4.1, NetMHCIIpan-4.3, BepiPred-3.0 and DiscoTope-3.0 are DTU Health Tech tools under an academic-only license and must be installed separately; the scripts here call them if present and otherwise fall back. Conformational B-cell prediction needs a structure (experimental or predicted).
metadata:
  version: "1.1"
  skill-author: Drafted for vaccine-design subset
  last-reviewed: "2026-10-16"
---

# Immunoinformatics

## When to use

Use this workflow to go from a set of antigen sequences to a ranked, filtered,
population-aware set of candidate epitopes, and to assess how well a proposed
immunogen covers a target population and the variants actually circulating.

Keep five kinds of result separate, and never let one stand in for another:

1. **A predicted binder** is a sequence a model scores as likely to bind an
   allele. Most predicted binders are not presented.
2. **A presented peptide** is one observed on the cell surface by
   immunopeptidomics. Presentation does not imply recognition.
3. **An epitope** is a peptide experimentally shown to be recognised by a T-cell
   receptor or antibody. IEDB records these; predictors do not produce them.
4. **An immunogenic epitope** elicits a response when delivered in a given
   construct, dose, adjuvant and host. This is an experimental result.
5. **A protective epitope** correlates with protection in a challenge or
   efficacy study. Nothing upstream of a trial establishes this.

Binding prediction is the weakest of the five and is where almost all
computational output sits. Treat the ranked list this skill produces as a
shortlist for experimental screening, not as a design.

Class II prediction is substantially less accurate than class I — the binding
groove is open at both ends, registers shift, and benchmark AUCs run well below
class I. Do not present class II rankings with the same confidence.

## Versions and installation

```bash
# Permissively licensed, works everywhere.
pip install "mhcflurry>=2.1" "biopython>=1.83" "pandas>=2.0" "requests>=2.31"
mhcflurry-downloads fetch models_class1_presentation

# Optional: MixMHCpred (academic license, separate download)
#   https://github.com/GfellerLab/MixMHCpred

# Optional DTU Health Tech tools — academic use only, register and download:
#   NetMHCpan-4.1, NetMHCIIpan-4.3, BepiPred-3.0, DiscoTope-3.0
#   https://services.healthtech.dtu.dk/
# Install to a directory on PATH. The scripts detect them by executable name.
```

**Licensing matters here.** The DTU predictors are free for academic use and
require a commercial license otherwise. MHCflurry (Apache 2.0) is the permissive
substitute for class I; there is no equally permissive class II substitute, so
commercial class II work needs a DTU license. Check before wiring any of these
into a pipeline you intend to publish or productise.

## Workflow

Run the steps through `scripts/epitope_triage.py`. It encodes the ordering, the
column conventions and the safety filter, and it records tool versions in a
sidecar `.meta.json` for reproducibility. Ad-hoc prediction code skips all of
that. If a subcommand fails, report the error rather than reimplementing it.


### 1. Establish the antigen set

Start from the actual protein sequences, with provenance. Record the strain or
isolate, the accession, and the date. Resolve signal peptides, cleavage products
and polyprotein boundaries before predicting — a peptide spanning a cleavage site
cannot be presented from the mature protein.

**Never write a sequence from memory.** Every residue must come from a retrieved
record with an accession attached. A recalled sequence will be subtly wrong, and
because every downstream step still runs cleanly on a wrong sequence, the error
is invisible in the output and poisons the entire result set.

If retrieval fails, stop and report the failure. Do not substitute a
"canonical", "reference", "representative" or "well-known" sequence written
out directly — those words are the tell that no record was retrieved. Say which
database was queried, what came back, and ask the user for an accession or a
FASTA file. Stopping costs minutes; a fabricated antigen costs the whole study.

Before accepting a retrieved sequence, verify three things and state them:

1. **It is the protein you asked for.** Check the record's name and organism.
   Hits returning a different protein from the same organism, or the same
   protein family from a different organism, are common when an ID is guessed
   rather than looked up.
2. **The length is right.** Confirm against what the protein should be.
3. **The accession resolves.** Quote it. If the ID was guessed rather than
   obtained from a search result, that is not provenance.

Prefer the UniProt connector over constructing REST URLs by hand. If a
connector returns something unexpected, that is a signal to search properly,
not to try another guessed identifier.

```bash
python scripts/epitope_triage.py prepare \
  --fasta antigens.fasta \
  --out work/antigens.clean.fasta
```

For a drifting pathogen, pull the currently circulating variants first
(`pathogen-variant-surveillance` in this collection) and carry the variant set
through to step 6. Designing against a reference sequence that no longer
circulates is the most common avoidable error in this workflow.

### 2. Predict MHC class I binding

Choose the allele set deliberately. The default supertype panel below covers a
large fraction of most populations, but if your target population is specific —
a region, an ancestry group, a cohort — use its allele frequencies instead
(Allele Frequency Net Database, `afnd.net`).

```bash
python scripts/epitope_triage.py mhc1 \
  --fasta work/antigens.clean.fasta \
  --alleles HLA-A*01:01,HLA-A*02:01,HLA-A*03:01,HLA-A*24:02,HLA-B*07:02,HLA-B*08:01,HLA-B*44:02 \
  --lengths 8,9,10,11 \
  --out work/mhc1.tsv
```

Rank on **percentile rank**, not affinity. Percentile rank normalises for the
fact that some alleles bind more promiscuously than others; raw nM values are not
comparable across alleles. Conventional cutoffs: strong binder ≤ 0.5%, weak
binder ≤ 2%. These are convention, not biology — they were chosen to balance
sensitivity and list length.

Prefer a presentation score (MHCflurry's `presentation_score`, or NetMHCpan's
`%Rank_EL` eluted-ligand rank) over a pure binding-affinity rank. Eluted-ligand
models are trained on mass-spec immunopeptidomics and reflect antigen processing,
not just groove affinity.

### 3. Predict MHC class II binding

```bash
python scripts/epitope_triage.py mhc2 \
  --fasta work/antigens.clean.fasta \
  --alleles HLA-DRB1*01:01,HLA-DRB1*03:01,HLA-DRB1*04:01,HLA-DRB1*07:01,HLA-DRB1*15:01 \
  --out work/mhc2.tsv
```

This requires NetMHCIIpan; the script reports what is missing rather than
silently substituting. Class II matters for vaccine design well beyond CD4
help — it drives B-cell help and affinity maturation, so a construct with good
class I coverage and no class II content often disappoints.

### 4. Check against experimentally validated epitopes

Before ranking anything, find out what is already known. IEDB holds curated
experimental data for most studied pathogens, and a validated epitope outranks
any prediction.

```bash
python scripts/epitope_triage.py iedb \
  --organism "<taxon id or name>" \
  --out work/iedb_known.tsv
```

Cross-reference the prediction output against this. Overlap with a known
epitope is the strongest evidence available at this stage; a high-ranking
prediction in a well-studied antigen with no IEDB support deserves scepticism,
because it has probably already been looked for.

### 5. Predict B-cell epitopes

Antibody epitopes are overwhelmingly conformational, so linear predictors have
limited value on their own. Use them only to triage, and prefer structure-based
prediction whenever a structure exists or can be predicted.

- **Linear:** BepiPred-3.0. Treat the output as a smoothed propensity track, not
  a set of epitopes.
- **Conformational:** DiscoTope-3.0 or SEMA, on an experimental structure or a
  high-confidence predicted model. Check pLDDT per residue if using a predicted
  structure — a low-confidence loop produces a meaningless surface.
- **Accessibility:** compute relative solvent accessibility and exclude buried
  residues regardless of score.

For a glycosylated viral surface protein, run the glycan analysis
(`glycoengineering` in this collection) in parallel. A site under a dense glycan
shield can score well and be inaccessible to antibody; this is the dominant
failure mode for envelope glycoprotein epitope prediction.

### 6. Assess conservation across circulating variants

**Use the `conservation-mapping` skill when you have, or can build, an
alignment.** It applies Henikoff sequence weighting, reports effective sample
size, and reconciles alignment positions to PDB residue numbering — none of
which the fallback below can do.

```bash
# align (phylogenetics skill), then score
python ../conservation-mapping/scripts/conservation_score.py \
  --aln work/panel.afa --ref <REFERENCE_ACCESSION> \
  --out work/conservation.tsv \
  --pdb input/antigen.pdb --pdb-chain A --pdb-out work/conservation.pdb
```

Join candidates to that table on `ref_pos` to get a per-residue conservation
profile across each peptide, rather than a single intact/not-intact call.

The alignment-free fallback, for when no alignment is practical:

```bash
python scripts/epitope_triage.py conservation \
  --candidates work/mhc1.tsv \
  --variants work/variant_sequences.fasta \
  --out work/conserved.tsv
```

This counts exact peptide occurrences in unaligned sequences. It does not
correct for redundant sampling, so a panel dominated by one clade will report
clade-specific peptides as highly conserved. If you use it, say so, and report
the panel composition alongside the fraction.

Either way: weight selection toward conserved candidates, but check *why* a
region is conserved. Functional constraint is a good reason; shallow or skewed
sampling is not, and the two are indistinguishable from the fraction alone.

### 7. Screen for self-similarity and allergenicity

This step is a safety filter and is not optional for anything intended to enter
a human.

```bash
python scripts/epitope_triage.py selfcheck \
  --candidates work/conserved.tsv \
  --human-proteome UP000005640_9606.fasta \
  --out work/flagged.tsv
```

Flag and review any candidate with close identity to a human protein. Exact or
near-exact matches to self-peptides carry a theoretical autoimmunity risk and
should be dropped or escalated to someone qualified to assess them, not
silently ranked. Molecular mimicry is a hypothesis-generating signal, not a
verdict — but an unexamined self-match in a vaccine construct is indefensible.

Screen separately against known allergens (AllergenOnline, or the FAO/WHO
criterion of 35% identity over an 80-residue window).

### 8. Population coverage

```bash
python scripts/epitope_triage.py coverage \
  --candidates work/flagged.tsv \
  --population "World" \
  --out work/coverage.tsv
```

This calls IEDB's population coverage tool. Report coverage as the fraction of
the target population expected to present at least one candidate, and give the
figure per region rather than only globally — a construct at 95% world coverage
can sit far lower in a specific population, which is exactly the kind of gap
that matters for deployment.

### 9. Design a multivalent construct, if breadth demands one

When no single antigen covers the strain panel, the options are a cocktail of
natural sequences or a synthetic mosaic set. Both are coverage-optimisation
problems over k-mers, and both are handled by `mosaic_design.py`.

```bash
# natural-sequence cocktail: every component is a real strain
python scripts/mosaic_design.py cocktail --fasta work/panel.fasta \
  -n 4 -k 9 --out work/cocktail

# synthetic mosaics: more coverage per component
python scripts/mosaic_design.py mosaic --fasta work/panel.fasta \
  -n 2 -k 9 --beam 200 --out work/mosaic
```

**Mosaic mode cannot create junctional neo-epitopes.** It walks a de Bruijn
graph built from the panel, so every k-mer in the output occurs in some real
strain by construction. The script verifies this and fails loudly if it ever
does not.

**Check any other construct with `evaluate`.** Naive concatenation,
string-of-beads designs and fusion constructs all create junction-spanning
peptides that exist in no real strain:

```bash
python scripts/mosaic_design.py evaluate --fasta work/panel.fasta \
  --designs work/my_construct.fasta -k 9 --out work/construct_eval
```

Run this at k=9 and again at k=15 to cover class I and class II registers.

**Read the coverage output honestly.** The script reports aggregate coverage
and the per-strain distribution. The aggregate is the misleading one — a design
can reach 90% weighted coverage while an entire clade sits near zero, because
weighted optimisation correctly ignores a sparsely sampled clade. Report the
minimum and 10th percentile, name the uncovered strains, and state whether the
gap was accepted deliberately. `--unweighted` trades common-variant depth for
rare-variant breadth if that is the call you want to make.

Cocktails are guaranteed foldable and processable because the components are
real proteins. Mosaics usually are not foldable — fine for mRNA or peptide
delivery, not for anything that must present a conformational surface. Carry
that distinction into `construct-assembly`, which applies different QC
depending on whether a construct must fold.

### 10. Report

Produce a table with, per candidate: sequence, source antigen and position,
restricting alleles with percentile ranks, IEDB support, conservation fraction,
self-similarity flag, and accessibility or structural context. Record the tool
versions and the allele panel — predictor versions change rankings, and a result
without its version is not reproducible.

Use `scientific-writing` in this collection for the write-up, and state
explicitly in any methods section that candidates are computational predictions
pending experimental validation.

## Interpretation

- **Rank ordering is not calibrated.** A candidate at rank 1 and one at rank 20
  are not reliably distinguishable. Carry enough candidates into screening to
  absorb that, typically tens rather than a handful.
- **Negative predictions are weak.** Absence from a shortlist is not evidence a
  peptide is not an epitope, particularly for class II.
- **Allele coverage is uneven.** Predictors are trained disproportionately on
  common European-ancestry alleles. Performance on under-represented alleles is
  worse and under-reported, which propagates directly into coverage inequity.
- **Processing is modelled crudely.** Proteasomal cleavage and TAP transport
  predictors add little over a good eluted-ligand model and should not be
  stacked as independent filters.
- **A construct is not a list.** Epitope order, linkers, flanking residues and
  delivery format all change what is processed and presented. String-of-beads
  designs create junctional neo-epitopes at every join; check those with
  `mosaic_design.py evaluate` (step 9) rather than by eye.

## Reference files

- `references/tool-licensing.md` — licensing and current versions for every
  predictor, folding model and design tool named here, including the code-vs-
  weights splits that a single licence field would misrepresent. Read it before
  wiring any tool into a pipeline intended for publication or product. Note it
  records NetMHCpan as 4.2e, superseding the 4.1 named below.
- `references/compute-routes.md` — which steps run on CPU locally, which go to
  Tamarind, and which free hosted routes exist when Tamarind does not fit.

## Official references

- IEDB Analysis Resource — https://nextgen-tools.iedb.org/ ; API docs at
  https://nextgen-tools.iedb.org/docs/
- IEDB database — https://www.iedb.org/
- MHCflurry 2.0 — O'Donnell et al., *Cell Systems* 2020; docs at
  https://openvax.github.io/mhcflurry/
- NetMHCpan-4.1 / NetMHCIIpan-4.0 — Reynisson et al., *NAR* 2020
- BepiPred-3.0 — Clifford et al., *Protein Science* 2022
- DiscoTope-3.0 — Høie et al., *Front. Immunol.* 2024
- Allele Frequency Net Database — http://www.allelefrequencies.net/
- Human reference proteome UP000005640 — https://www.uniprot.org/proteomes/UP000005640
