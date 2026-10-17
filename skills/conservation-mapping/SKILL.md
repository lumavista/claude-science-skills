---
name: conservation-mapping
description: Scores per-position sequence conservation across a panel of strains, variants or homologues, and maps it onto protein structure. Covers Henikoff-weighted Shannon entropy with effective sample size, reconciliation of alignment positions to PDB residue numbering, conservation written into B-factors for visualisation, relative solvent accessibility, and clustering of conserved exposed residues into surface patches sized against an antibody footprint. Use when asked which residues or regions of a protein are conserved, where the conserved surface is, whether an epitope or binding site is variable across circulating strains, or to produce a per-residue conservation track for downstream epitope or immunogen work. This is the conservation engine for the immunoinformatics workflow; use it in place of per-candidate substring counting whenever an alignment is available.
license: MIT
compatibility: Python 3.10+, numpy, biopython (for PDB numbering reconciliation). Surface patches additionally need freesasa (conda-forge) or an externally supplied RSA table, e.g. from DSSP. Alignment is not performed here — use the phylogenetics skill. No network access required. Scoring stays CPU-bound — Henikoff weighting and per-column entropy are pure-Python loops over a few MB, so a realistic 2000 × 1500 panel is ~10^6–10^7 character operations, too small for a GPU to pay off. The available speedup is numpy vectorisation, not CUDA.
metadata:
  version: "1.0"
  skill-author: Drafted for vaccine-design subset
  last-reviewed: "2026-10-16"
---

# Conservation Mapping

Turns a strain panel into a per-residue conservation track in the numbering of
a reference sequence, and projects that track onto structure so conserved
surface can be seen rather than inferred.

This skill owns conservation scoring for the collection. `immunoinformatics`
step 6 delegates here whenever an alignment is available; its
`epitope_triage.py conservation` subcommand remains only as an alignment-free
fallback for counting intact peptides, and it cannot see per-position
structure or correct for sampling skew.

## The failure mode this exists to prevent

**Strain panels reflect surveillance effort, not pathogen diversity.** Sequence
databases carry thousands of entries from clades someone funded a study on and
a handful from clades circulating where no one sequences. Compute unweighted
conservation on that and every clade-dominant residue reads as conserved —
which sends the entire downstream design after the wrong residues, and nothing
later in the pipeline catches it.

Two defences, both on by default:

1. **Henikoff position-based sequence weighting**, which downweights redundant
   sequences.
2. **Effective sample size**, printed at the top of every run and recorded in
   the `.meta.json` sidecar. An effective N of 12 from 400 sequences is the
   single most useful number this skill produces.

State the effective N, not the raw sequence count, whenever you report
conservation. The specialist contract for Immuno-Designer already requires
stating how many sequences a conservation figure was computed over; effective
N is the honest version of that number.

## Workflow

### 1. Assemble the panel

Sequences come from a retrieved record with an accession, never from memory —
the same rule as `immunoinformatics` step 1, and for the same reason: every
downstream step runs cleanly on a wrong sequence, so the error is invisible.

For a drifting pathogen, pull the currently circulating set with
`pathogen-variant-surveillance` first. Designing against a reference that no
longer circulates is the most common avoidable error in this workflow.

Before aligning, check: are these all the same protein (fragments and
misannotated paralogs are common in keyword-assembled panels); are there exact
duplicates (collapse them, record the counts); is the clade structure what the
user believes it to be; and is the antigen you care about in the panel, so its
numbering anchors everything downstream.

### 2. Align

Alignment belongs to `phylogenetics` — do not align here.

```bash
python ../phylogenetics/scripts/phylogenetic_analysis.py panel.fasta \
  --aligner famsa --type aa --output-dir work --no-visualization
```

For panels spanning high divergence or with large indels, prefer MAFFT
`--mafft-method einsi` and inspect the result. Automatic alignment of divergent
envelope glycoproteins across variable loops is genuinely unreliable, and an
alignment error in a variable loop propagates into a false conserved call in
the flanking region. Note which regions you do not trust and carry that note
forward — a design built on a badly aligned loop fails for reasons no later
step will attribute correctly.

### 3. Score conservation

```bash
python scripts/conservation_score.py \
  --aln work/panel.afa \
  --ref <REFERENCE_ACCESSION> \
  --out work/conservation.tsv \
  --pdb input/antigen.pdb --pdb-chain A \
  --pdb-out work/conservation.pdb
```

| Column | Use |
|---|---|
| `ref_pos` | reference numbering |
| `pdb_resnum` | the structure's own residue ID, reconciled by alignment — the key for every structural join |
| `entropy` | normalised Shannon entropy, 0 invariant, 1 maximally variable |
| `class_entropy` | entropy over biochemical classes; low class entropy with high residue entropy means only conservative swaps are tolerated |
| `gap_frac` | weighted gap fraction — a high-gap column is an indel region, not a conserved one |
| `consensus_aa`, `consensus_frac` | majority residue and weighted support |
| `n_variants` | distinct residues observed |
| `ref_differs_from_consensus` | flags an atypical reference |
| `nglyc_sequon_N` | N-X-S/T, X≠P, annotation only |

**On PDB numbering.** FASTA offsets are not PDB residue IDs — construct tags,
disordered gaps and author numbering schemes all shift them, and
`glycoengineering/references/glycoshield.md` flags this as a hazard without
providing a fix. This script resolves it by aligning the structure's own
CA-derived sequence to the reference and reporting the identity achieved. If
that identity drops below 0.8 it warns; treat that as a signal the PDB is the
wrong protein or the wrong chain, not as a number to work around.

**On B-factors.** The written PDB carries `100*(1-entropy)`, inverted so that
standard blue-to-red B-factor colouring puts conserved surface in the hot
colour. GlycoSHIELD also encodes its shielding into B-factors — keep the two
files separate and label them, because the scales are not comparable.

**On sequons.** The `nglyc_sequon_N` column is annotation for convenience only.
`glycoengineering` owns sequon logic, occupancy and shield geometry; send any
actual glycan work there rather than reasoning from this column.

### 4. Find conserved surface patches

```bash
python scripts/surface_patches.py \
  --pdb work/conservation.pdb --chain A \
  --conservation work/conservation.tsv \
  --out work/patches.tsv
```

Conservation alone is not a target — buried residues are conserved because the
fold requires it. This intersects conservation with relative solvent
accessibility and clusters the survivors by CA distance into spatially
contiguous patches.

**Area is the number that matters**, not residue count. A typical Fab buries
roughly 600–900 Å²; the script flags patches reaching that threshold. A
conserved surface that fragments into several small patches is telling you
something real about the antigen, and the honest response is to report it
rather than loosening thresholds until a large patch appears. If you do change
`--entropy-max` or `--linkage`, say which and why.

**For envelope glycoproteins this result is incomplete on its own.** A patch
that is conserved, exposed and well-sized can still be antibody-inaccessible
under a dense glycan shield — the dominant failure mode for envelope site
selection. Run `glycoengineering`'s GlycoSHIELD workflow and read shielding
alongside these patches, not instead of them.

### 5. Consensus and ancestral sequences

If a single sequence representing the panel is wanted:

- **Consensus** — the `consensus_aa` column. Cheapest, usually adequate above
  roughly 85% identity.
- **Ancestral reconstruction** — `phylogenetics`, IQ-TREE `--ancestral`
  (see its `references/iqtree_inference.md`). Better across deeper divergence.
- **Center-of-tree** — the panel member minimising distance to all others.
  Guaranteed natural, therefore guaranteed foldable, at some cost in breadth.

**The consensus trap:** a per-position majority can produce a residue
combination that never co-occurs in nature and does not fold. Check consensus
designs for co-variation violations — positions whose consensus residues are
never observed together in any panel member — and read a consensus that fails
to refold as evidence of epistasis rather than a tool failure. For envelope
glycoproteins this is common enough that center-of-tree is often the safer
start.

## Interpretation

- **Low entropy is not evidence of functional constraint.** It can equally mean
  the panel is shallow or clade-skewed. Check effective N before attributing
  conservation to selection.
- **High entropy is not evidence of irrelevance.** Immunodominant variable
  loops are variable *because* they are under antibody pressure, which makes
  them functionally important and a poor target at the same time.
- **A single structure is one conformer.** Accessibility and loop position
  move. Patches computed on one model are a hypothesis about shape.
- **Conservation across circulating strains is not conservation under
  vaccine pressure.** A site conserved today may be free to vary once a
  response targets it; escape potential is a separate question this skill does
  not answer.

## Official references

- Henikoff & Henikoff, position-based sequence weights, *J Mol Biol* 1994
- Tien et al., maximum allowed solvent accessibility, *PLoS ONE* 2013
- Valdar, scoring residue conservation, *Proteins* 2002
- `phylogenetics` (this collection) — alignment and ancestral reconstruction
- `glycoengineering` (this collection) — sequons, glycan shield, GlycoSHIELD
- `immunoinformatics` (this collection) — the epitope workflow this feeds
