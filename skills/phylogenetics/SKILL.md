---
name: phylogenetics
description: Builds and analyzes phylogenetic trees using MAFFT or pyfamsa multiple sequence alignment, IQ-TREE maximum likelihood with ModelFinder and branch support, and FastTree approximate inference. Uses ETE3 for tree summaries and visualization. Applies to homologous nucleotide or protein sequences, microbial gene trees, protein families, and cautiously interpreted dated phylogenies.
license: Unknown
compatibility: Requires MAFFT and IQ-TREE 3 or FastTree executables on Linux/macOS. On Windows MAFFT, IQ-TREE and trimAl are unavailable in practice — use pyfamsa and FastTree. Optional tree summaries need Python 3.12 and ete3 3.1.3; rendering also needs PyQt5. Network only for installation.
metadata:
  version: "1.6"
  skill-author: Kuan-lin Huang, modified by Lumavista.BIO
  last-reviewed: "2026-10-11"
---

# Phylogenetics

## When to use

For an alignment of homologous loci or proteins, a gene tree, branch support,
and an interpretable tree figure. Multi-locus partitioning, concordance factors,
ancestral states and dating are in
[IQ-TREE reference](references/iqtree_inference.md).

A gene tree does not by itself establish a species tree, transmission direction
or horizontal transfer, and sampling dates alone do not justify a molecular
clock.

## Versions and installation

Reviewed against official documentation on 2026-10-01. The smoke workflow used
MAFFT **7.526**, IQ-TREE **3.1.4**, FastTree **2.2.0** and ETE3 **3.1.3** on
Python 3.12. ETE3 is deliberate — ETE4 is a separate API, not a drop-in
replacement — and ETE3 imports `cgi`, removed in Python 3.13.

Inference is CPU-bound: threads are the only parallelism, since neither FastTree
nor IQ-TREE has a CUDA path and BEAGLE's CUDA kernels serve Bayesian MCMC in
BEAST/MrBayes, a different method rather than an accelerator for this one.

```bash
conda create -n phylo -c conda-forge -c bioconda mafft iqtree fasttree trimal
uv pip install ete3==3.1.3 numpy six PyQt5   # PyQt5 only for rendering
```

**Verify that each executable launches, not merely that the environment built.**
A solver that cannot satisfy a spec may silently drop it and report the
environment as partially available, surfacing much later as a missing-command
error attributed to the wrong cause. Print every version before using it.

The helper defaults to `iqtree3` and `FastTree`; use `--iqtree-bin iqtree2` or
`--fasttree-bin fasttree` for differently named installs (version 2
compatibility is source-checked, not runtime-tested). All analyses are local.

### Platform availability

bioconda publishes no win-64 builds of MAFFT, IQ-TREE, FastTree or trimAl, so on
Windows all four specs are dropped, and vendor Windows binaries are not a
substitute: MAFFT, MUSCLE, IQ-TREE and trimAl ship non-relocatable PE images (no
`.reloc` section) that fail with `WinError 623` under forced ASLR, including
sandboxed environments, and cannot be repaired with `editbin /DYNAMICBASE`.
Check a candidate binary's PE header for `DYNAMICBASE` plus a base-relocation
directory — FastTree 2.2.0 passes, which is why it runs where the others do not.
Python extension modules are built with ASLR enabled, so
`pip install pyfamsa pyhmmer parasail biopython` (all win_amd64 wheels) loads
where a standalone executable does not.

| Role | Linux/macOS | Windows |
|---|---|---|
| MSA | MAFFT | `pyfamsa` (FAMSA) |
| Profile alignment | HMMER | `pyhmmer` |
| ML inference | IQ-TREE | FastTree only |
| Pairwise | — | `parasail`, Biopython `PairwiseAligner` |
| Trimming | trimAl | no wheel; filter columns in Biopython |

If a conda `activate.d` directory is read-only, invoke executables by absolute
path.

**Losing IQ-TREE on Windows costs ModelFinder, ultrafast bootstrap and
SH-aLRT**, leaving FastTree's approximate search and SH-like local supports as
the ceiling. State that limitation in any write-up rather than letting the
reader assume a full ML analysis.

## Workflow

### 1. Establish the biological input

- Select homologs with appropriate taxon coverage, and investigate paralogs,
  contamination and recombination before assuming one history describes all
  sites.
- Use unique, whitespace-free FASTA IDs with a separate metadata table; the
  helper rejects empty/all-missing records and Newick punctuation in IDs.
- Specify `--type nt` or `--type aa`.
- Coding DNA needs codon-aware alignment if using a codon model — ordinary MAFFT
  nucleotide alignment does not preserve reading frames.
- Keep the original sequences, accession versions, filters, alignment and exact
  commands; inspect alignment length, gaps, composition and taxon identity.

### 2. Align and inspect

```bash
mafft --auto --thread 1 --inputorder sequences.fasta > aligned.fasta
```

`--auto` delegates the strategy choice to MAFFT, and sequence count alone is not
an accuracy criterion; choose L-INS-i, E-INS-i or an FFT-NS variant deliberately
per the MAFFT manual. Single-thread runs simplify reproducibility, since
iterative refinement with multiple threads can differ between runs.

Where MAFFT is unavailable, align with `pyfamsa` (FAMSA, progressive, designed
for large families, in-process with no file round-trip):

```python
from Bio import SeqIO
from pyfamsa import Aligner, Sequence

seqs = [Sequence(rec.id.encode(), str(rec.seq).encode())
        for rec in SeqIO.parse("sequences.fasta", "fasta")]
msa = Aligner(guide_tree="upgma", threads=1).align(seqs)
# Write FASTA for FastTree, then confirm the count matches the input.
```

Name FAMSA and its version in the methods — a wheel-based aligner is still an
aligner choice that changes the alignment, and "aligned in Python" is not a
method.

**Reference-anchored projection is not a multiple sequence alignment**, and is
acceptable only for sequences verified indel-free and length-uniform; once
indels are present, independently placed gaps are mutually inconsistent. Report
the indel-freedom evidence, and prefer a real aligner.

Trimming is optional and must be justified by alignment quality; the trimAl
command below is documentation-checked, not executed in this review, and with no
`pytrimal` win-64 wheel you filter columns by gap fraction and occupancy in
Biopython instead, recording the rule.

```bash
trimal -in aligned.fasta -out trimmed.fasta -automated1 -fasta
```

Record retained columns/taxa and compare sensitivity to trimming. **A failed
trim must stop the trimming step** — never copy the untrimmed input under a
"trimmed" filename. The helper does not trim.

### 3. Infer a tree

Use the maintained helper rather than recreating subprocess wrappers, from the
skill directory with the executables on `PATH`:

```bash
python scripts/phylogenetic_analysis.py sequences.fasta --type nt --threads 1 \
  --bootstrap 1000 --seed 42 --output-dir results --no-visualization
# --aligned skips realignment; --fasttree selects FastTree (which ignores
# --threads). Without MAFFT and IQ-TREE, --aligned --fasttree is the only path.
python scripts/phylogenetic_analysis.py aligned.fasta --aligned --type aa \
  --fasttree --output-dir reviewed --no-visualization
```

**Account for every input sequence through alignment and inference**, and report
any dropped: a silently shortened alignment produces a clean-looking tree of the
wrong taxon set, and the count is the cheapest place to catch it.

IQ-TREE defaults are `-m MFP`, `-B 1000`, `--alrt 1000`; use `--model` only for
a justified fixed model, and a new output prefix when data or model assumptions
change. FastTree uses `-nt -gtr -gamma` for DNA and `-lg -gamma` for proteins,
searches with the CAT approximation before Gamma20 rescaling (so it is **not** a
full Gamma-model ML search), and reports **SH-like local support in the 0–1
range**, not UFBoot percentages. Choose by workload benchmarking, not a taxon
threshold. Flags, output files and checkpoint resumption:
[IQ-TREE reference](references/iqtree_inference.md).

### 4. Interpret support and rooting

Read `.iqtree` and `.log`, not only the image, and verify the order of IQ-TREE's
slash-delimited SH-aLRT/UFBoot labels in the report. UFBoot ≥95 and SH-aLRT ≥80
are screening thresholds, not proof of correctness or posterior probabilities,
and the standard-bootstrap 70% rule does not apply to UFBoot. Model
misspecification, alignment error and sampling can all produce strongly
supported wrong branches; consider `-bnni` when assessing model violations and
UFBoot overestimation.

**Reversible substitution models do not infer the evolutionary root.** The
helper preserves the supplied orientation unless `--outgroup` (one taxon,
existence verified) or `--midpoint` is given; midpoint rooting assumes clock-like
path lengths and is only a heuristic, and display rerooting preserves support
labels on their splits without creating a dated phylogeny or overwriting the
inferred tree. Keep the tree file, alignment and IQ-TREE artifacts together.

### 5. Summarize and visualize

ETE must read combined IQ-TREE labels as internal **names** with `format=1` —
its default numeric-support parser cannot parse `95.2/99`. The helper renders
the original text labels, and branch statistics use all original non-root edges,
including zero-length ones, without rerooting.

```python
# Run from the skill directory in the ETE3 environment.
import sys; sys.path.insert(0, "scripts")
from phylogenetic_analysis import load_tree, tree_summary
print(tree_summary("results/sequences.treefile"))
```

Omit `--no-visualization` to render. ETE3 needs PyQt5 and, headless, a Qt
platform plugin; missing rendering support leaves the inferred tree intact and
is reported explicitly, and Qt rendering was not exercised in this review.
FigTree and iTOL also display Newick — verify the viewer retains support labels
and the branch-length scale.

## Official references

- [MAFFT manual](https://mafft.cbrc.jp/alignment/software/manual/manual.html), [threading](https://mafft.cbrc.jp/alignment/software/multithreading.html)
- [IQ-TREE quickstart](https://iqtree.github.io/doc/Quickstart), [releases](https://github.com/iqtree/iqtree3/releases), [model assumptions](https://iqtree.github.io/doc/Assessing-Phylogenetic-Assumptions)
- [FastTree manual](https://morgannprice.github.io/fasttree/)
- [ETE3 parsing](https://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html), [drawing](https://etetoolkit.org/docs/latest/tutorial/tutorial_drawing.html)
- [trimAl manual](https://vicfero.github.io/trimal/), [pyhmmer docs](https://pyhmmer.readthedocs.io/en/stable/)
- [pyfamsa docs](https://pyfamsa.readthedocs.io/en/stable/); FAMSA — Deorowicz, Debudaj-Grabysz & Gudyś, *Scientific Reports* 2016
