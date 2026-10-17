---
name: conservation-mapping
description: Scores per-position conservation across a strain, variant or homologue panel and maps it onto protein structure — Henikoff-weighted Shannon entropy with effective sample size, alignment-to-PDB numbering reconciliation, conservation in B-factors, and clustering of conserved exposed residues into antibody-sized surface patches. Use for conserved residues, conserved surface or a conservation track.
license: MIT
compatibility: Python 3.10+, numpy, biopython for PDB numbering reconciliation. Surface patches also need freesasa (conda-forge) or an externally supplied RSA table. Alignment is not performed here — use the phylogenetics skill. No network access required.
metadata:
  version: "1.0"
  skill-author: Drafted for vaccine-design subset
  last-reviewed: "2026-10-16"
---

# Conservation Mapping

Turns a strain panel into a per-residue conservation track in a reference
sequence's numbering, and projects it onto structure so conserved surface can be
seen rather than inferred. This skill owns conservation scoring for the
collection: `immunoinformatics` step 6 delegates here whenever an alignment is
available, and its `epitope_triage.py conservation` subcommand remains only as
an alignment-free fallback blind to per-position structure and sampling skew.

Scoring is CPU-bound — Henikoff weighting and per-column entropy are
pure-Python loops over a few MB, so a 2000 × 1500 panel is ~10^6–10^7 character
operations, too small for a GPU to pay off; the speedup available is numpy
vectorisation, not CUDA.

## The failure mode this exists to prevent

**Strain panels reflect surveillance effort, not pathogen diversity**, so
unweighted conservation makes every clade-dominant residue read as conserved and
sends the whole downstream design after the wrong residues, with nothing later
to catch it.

Two defences, both on by default: **Henikoff position-based sequence weighting**,
which downweights redundant sequences, and **effective sample size**, printed
each run and recorded in the `.meta.json` sidecar. An effective N of 12 from 400
sequences is the single most useful number this skill produces, so **report the
effective N, not the raw sequence count** — the Immuno-Designer contract already
requires saying how many sequences a figure covers, and this is the honest
version of that.

## Workflow

### 1. Assemble the panel

Sequences come from a retrieved record with an accession, never from memory —
the same rule as `immunoinformatics` step 1, because every downstream step runs
cleanly on a wrong sequence. For a drifting pathogen, pull the currently
circulating set with `pathogen-variant-surveillance` first; designing against a
reference that no longer circulates is the most common avoidable error here.

Before aligning: confirm one protein (keyword-assembled panels collect fragments
and misannotated paralogs), collapse exact duplicates recording counts, check
the clade structure against what the user believes, and include the antigen of
interest so its numbering anchors everything downstream.

### 2. Align

Alignment belongs to `phylogenetics` — do not align here.

```bash
python ../phylogenetics/scripts/phylogenetic_analysis.py panel.fasta \
  --aligner famsa --type aa --output-dir work --no-visualization
```

For high divergence or large indels prefer MAFFT `--mafft-method einsi` and
inspect the result: **automatic alignment of divergent envelope glycoproteins
across variable loops is unreliable, and an alignment error in a variable loop
becomes a false conserved call in the flanking region**, so note which regions
you do not trust and carry that note forward.

### 3. Score conservation

```bash
python scripts/conservation_score.py \
  --aln work/panel.afa --ref <REFERENCE_ACCESSION> \
  --out work/conservation.tsv \
  --pdb input/antigen.pdb --pdb-chain A \
  --pdb-out work/conservation.pdb
```

Columns: `ref_pos`, `consensus_aa`, `consensus_frac`, `n_variants` and
`ref_differs_from_consensus` are self-describing. The rest need reading:

| Column | Use |
|---|---|
| `pdb_resnum` | the structure's own residue ID, reconciled by alignment — the key for every structural join |
| `entropy` | normalised Shannon entropy, 0 invariant, 1 maximally variable |
| `class_entropy` | entropy over biochemical classes; low here with high `entropy` means only conservative swaps are tolerated |
| `gap_frac` | weighted gap fraction — a high-gap column is an indel region, not a conserved one |
| `nglyc_sequon_N` | N-X-S/T, X≠P; annotation only, since `glycoengineering` owns sequon logic |

**PDB numbering.** FASTA offsets are not PDB residue IDs — construct tags,
disordered gaps and author numbering schemes shift them — so the script aligns
the structure's CA-derived sequence to the reference and reports the identity
achieved; identity below 0.8 warns, and that means the wrong protein or chain,
not a number to work around.

**B-factors.** The written PDB carries `100*(1-entropy)`, so standard
blue-to-red colouring puts conserved surface in the hot colour. GlycoSHIELD also
writes shielding into B-factors on an incomparable scale — keep the files
separate and labelled.

### 4. Find conserved surface patches

```bash
python scripts/surface_patches.py \
  --pdb work/conservation.pdb --chain A \
  --conservation work/conservation.tsv \
  --out work/patches.tsv
```

Conservation alone is not a target — buried residues are conserved because the
fold requires it — so this intersects conservation with relative solvent
accessibility and clusters survivors by CA distance.

**Area is the number that matters**, not residue count: a typical Fab buries
roughly 600–900 Å² and the script flags patches reaching that. A conserved
surface fragmenting into small patches is a real property of the antigen, so
report it rather than loosening thresholds until a large patch appears, and say
so if you change `--entropy-max` or `--linkage`.

**For envelope glycoproteins this is incomplete on its own**: a conserved,
exposed, well-sized patch can still be antibody-inaccessible under a dense
glycan shield, the dominant failure mode for envelope site selection. Read
`glycoengineering`'s GlycoSHIELD shielding alongside these patches.

### 5. Consensus and ancestral sequences

For one representative sequence: **consensus** (`consensus_aa`, cheapest and
usually adequate above ~85% identity), **ancestral reconstruction**
(`phylogenetics`, IQ-TREE `--ancestral` — better across deeper divergence), or
**center-of-tree** (the member minimising distance to all others — natural and
therefore foldable, at some cost in breadth).

**The consensus trap:** a per-position majority can produce a residue
combination that never co-occurs in nature and does not fold, so check consensus
designs for co-variation violations and read a failure to refold as epistasis
rather than tool failure. Center-of-tree is often safer for envelope
glycoproteins.

## Interpretation

- **Low entropy is not evidence of functional constraint** — check effective N
  before ruling out a shallow or clade-skewed panel.
- **High entropy is not evidence of irrelevance** — immunodominant loops are
  variable *because* of antibody pressure, which makes them important and a poor
  target at once.
- **A single structure is one conformer**, so patches are a hypothesis about
  shape.
- **Conservation across circulating strains is not conservation under vaccine
  pressure**; escape potential is a separate question this skill does not answer.

## Official references

- Henikoff & Henikoff, position-based sequence weights, *J Mol Biol* 1994
- Tien et al., maximum allowed solvent accessibility, *PLoS ONE* 2013
- Valdar, scoring residue conservation, *Proteins* 2002
- `phylogenetics` — alignment and ancestral reconstruction
- `glycoengineering` — sequons, glycan shield, GlycoSHIELD
- `immunoinformatics` — the epitope workflow this feeds
