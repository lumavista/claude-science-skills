---
name: phylogenetics
description: Builds and analyzes phylogenetic trees using MAFFT or pyfamsa multiple sequence alignment, IQ-TREE maximum likelihood with ModelFinder and branch support, and FastTree approximate inference. Uses ETE3 for tree summaries and visualization. Applies to homologous nucleotide or protein sequences, microbial gene trees, protein families, and cautiously interpreted dated phylogenies.
license: Unknown
compatibility: Requires MAFFT and IQ-TREE 3 or FastTree executables on Linux/macOS. On Windows, MAFFT, IQ-TREE and trimAl are unavailable in practice (no bioconda win-64 builds; vendor binaries fail to load under forced ASLR) — use pyfamsa for alignment and FastTree for inference. Optional tree summaries use Python 3.12 and ete3 3.1.3; rendering also needs PyQt5. Network access is required for installation, not local inference.
metadata:
  version: "1.6"
  skill-author: Kuan-lin Huang, modified by Lumavista.BIO
  last-reviewed: "2026-10-11"
---

# Phylogenetics

## When to use

Use this workflow for an alignment of homologous loci or proteins, a gene tree,
branch support, and an interpretable tree figure. Multi-locus partitioning,
concordance factors, ancestral states and dating are covered in
[IQ-TREE reference](references/iqtree_inference.md).

A gene tree does not establish a species tree, transmission direction, or horizontal
transfer by itself. Sampling dates alone do not justify a molecular clock.

## Versions and installation

Reviewed against current official documentation on 2026-10-01. The executed
smoke workflow used MAFFT **7.526**, IQ-TREE **3.1.4**, FastTree **2.2.0** and
ETE3 **3.1.3** on Python 3.12. ETE3 remains deliberate for this helper; ETE4 is a
separate API and is not a drop-in replacement. ETE3 imports `cgi`, removed in
Python 3.13. The native executables are not Python packages.

Inference is CPU-bound: neither FastTree nor IQ-TREE has a CUDA path, so threads
are the only parallelism available. BEAGLE's CUDA kernels do exist, but they serve
Bayesian MCMC in BEAST or MrBayes — a different inference method than this skill
performs, so using them would mean replacing the method rather than accelerating it.

Illustrative installation commands; verify executable versions after installation:

```bash
conda create -n phylo -c conda-forge -c bioconda mafft iqtree fasttree trimal
conda activate phylo
uv venv --python 3.12 .venv-phylo
uv pip install --python .venv-phylo/bin/python ete3==3.1.3 numpy six
# Optional rendering backend; it can require a working Qt platform plugin.
uv pip install --python .venv-phylo/bin/python PyQt5
```

**Verify that each executable launches, not merely that the environment built.**
A solver that cannot satisfy a spec may report the environment as partially
available and drop it silently; the failure then surfaces much later as a
missing-command error attributed to the wrong cause. Print the version of every
executable before using it.

### Platform availability

The line above is a Linux/macOS command. bioconda publishes no win-64 builds of
MAFFT, IQ-TREE, FastTree or trimAl, so on Windows all four specs are dropped.
Vendor-supplied Windows binaries are not a general substitute either: MAFFT,
MUSCLE, IQ-TREE and trimAl ship non-relocatable PE images
(`DllCharacteristics 0x8120`, no `.reloc` section) and fail to load with
`WinError 623` wherever forced ASLR is in effect, which includes sandboxed
execution environments. Such binaries carry no relocation data and cannot be
repaired with `editbin /DYNAMICBASE`. Before investing in any Windows binary,
read its PE header and confirm `DYNAMICBASE` and a base-relocation directory;
FastTree 2.2.0 passes this check, which is why it runs where the others do not.

Compiled **Python extension modules are built with ASLR enabled by default**, so
a pip wheel loads where a standalone executable does not. On Windows, prefer a
wheel-based aligner:

```bash
# Alignment and pairwise work, all with win_amd64 wheels.
pip install pyfamsa pyhmmer parasail biopython
```

| Role | Linux/macOS | Windows |
|---|---|---|
| MSA | MAFFT | `pyfamsa` (FAMSA) |
| Profile alignment | HMMER | `pyhmmer` |
| ML inference | IQ-TREE | FastTree only |
| Pairwise | — | `parasail`, Biopython `PairwiseAligner` |
| Trimming | trimAl | no wheel; filter columns in Biopython |

On Windows the venv interpreter is `.venv-phylo\Scripts\python.exe`, not
`.venv-phylo/bin/python`. If a conda `activate.d` directory is read-only,
invoke executables by absolute path rather than attempting to install a PATH
hook.

Losing IQ-TREE on Windows costs ModelFinder, ultrafast bootstrap and SH-aLRT;
FastTree's approximate search and SH-like local supports are the ceiling there.
State that limitation in any write-up rather than letting the reader assume a
full ML analysis.

The helper defaults to `iqtree3` and `FastTree`; use `--iqtree-bin iqtree2` for an
existing IQ-TREE 2 installation, or `--fasttree-bin fasttree` if that is the
installed name. Version 2 compatibility is source-checked, not runtime-tested
in this refresh. All analyses are local; this skill invokes no HTTP API or remote
inference endpoint.

## Workflow

### 1. Establish the biological input

- Select homologs with appropriate taxon coverage; investigate paralogs,
  contamination and recombination before assuming one history describes all sites.
- Use unique, whitespace-free FASTA IDs and keep a separate sample metadata table.
  The helper rejects empty/all-missing records and Newick punctuation in IDs.
- Specify `--type nt` or `--type aa`. The helper sends `-st DNA` or `-st AA` to
  IQ-TREE and selects nucleotide/protein FastTree flags explicitly.
- Coding DNA needs codon-aware alignment if using a codon model. Ordinary MAFFT
  nucleotide alignment does not preserve reading frames automatically.
- Keep the original sequences, accession versions, filters, alignment and exact
  commands. Inspect alignment lengths, gaps, composition and taxon identities.

### 2. Align and inspect

These MAFFT commands were exercised on small synthetic FASTA data:

```bash
mafft --auto --thread 1 --inputorder sequences.fasta > aligned.fasta
# L-INS-i: locally alignable homologs, typically small datasets.
mafft --localpair --maxiterate 1000 --thread 1 --inputorder sequences.fasta > linsi.fasta
# E-INS-i: homologs with large internal gaps.
mafft --genafpair --maxiterate 1000 --thread 1 --inputorder sequences.fasta > einsi.fasta
# FFT-NS-i (two refinement cycles), then FFT-NS-2 (progressive).
mafft --retree 2 --maxiterate 2 --thread 1 --inputorder sequences.fasta > fftnsi.fasta
mafft --retree 2 --maxiterate 0 --thread 1 --inputorder sequences.fasta > fftns.fasta
```

`auto` delegates the choice to MAFFT. Sequence count alone is not an accuracy
criterion. Single-thread runs simplify reproducibility; MAFFT iterative refinement
with multiple threads can produce different results between runs.

Where MAFFT is unavailable, align with `pyfamsa`. FAMSA is a progressive
aligner designed for large families; it handles indels normally and returns
objects in process, with no subprocess or file round-trip:

```python
from Bio import SeqIO
from pyfamsa import Aligner, Sequence

seqs = [Sequence(rec.id.encode(), str(rec.seq).encode())
        for rec in SeqIO.parse("sequences.fasta", "fasta")]
aligner = Aligner(guide_tree="upgma", threads=1)
msa = aligner.align(seqs)
# Write FASTA for FastTree, then confirm the count matches the input.
```

Set `threads=1` for reproducibility, for the same reason as MAFFT. Name FAMSA
and its version in the methods — a wheel-based aligner is still an aligner
choice that changes the alignment, and "aligned in Python" is not a method.

**Reference-anchored alternatives are not multiple sequence alignments.**
Projecting each sequence onto one reference's coordinates by pairwise alignment
is acceptable only for sequences verified to be indel-free and length-uniform,
where column correspondence is unambiguous and the result is exact. The moment
indels are present, gaps placed independently per sequence are inconsistent
between sequences, and the matrix must not be described as an MSA. Verify
indel-freedom and report the evidence before relying on this shortcut; prefer a
real aligner whenever one is available.

Trimming is optional and must be justified by alignment quality, not assumed to
improve every tree. No `pytrimal` wheel exists for win-64; where trimAl is
unavailable, filter columns by gap fraction and occupancy directly in Biopython
and record the rule applied. The following trimAl command is
documentation-checked and illustrative; trimAl was not installed for this
review:

```bash
trimal -in aligned.fasta -out trimmed.fasta -automated1 -fasta
```

Record retained columns/taxa and compare sensitivity to trimming. A failed trim
must stop the trimming step; do not silently copy the untrimmed input under a
“trimmed” filename. The bundled helper does not trim automatically.

### 3. Infer a tree

Use the maintained helper instead of recreating subprocess wrappers. From the
skill directory, with the executable names on `PATH`:

```bash
python scripts/phylogenetic_analysis.py sequences.fasta --type nt --threads 1 \
  --bootstrap 1000 --seed 42 --output-dir results --no-visualization
# Use a reviewed alignment (including any explicit trimming) without realigning.
python scripts/phylogenetic_analysis.py aligned.fasta --aligned --type nt \
  --threads 1 --output-dir reviewed --no-visualization
# Protein FastTree alternative; --threads applies to MAFFT/IQ-TREE, not FastTree.
python scripts/phylogenetic_analysis.py proteins.fasta --type aa --fasttree \
  --threads 1 --output-dir fast_protein --no-visualization
```

Where IQ-TREE is unavailable, align separately (MAFFT or `pyfamsa`) and run the
helper on the reviewed alignment with `--aligned --fasttree`, which is the only
path the helper can take without the MAFFT and IQ-TREE executables.

Account for every input sequence through alignment and inference, and report any
that were dropped. A silently shortened alignment produces a clean-looking tree
of the wrong taxon set, and the count is the cheapest place to catch it.

The standard IQ-TREE command uses `-m MFP` (ModelFinder including FreeRate
models), `-B 1000` (UFBoot) and `--alrt 1000`. `-m TEST` remains supported but
searches a narrower model set. Use `--model` for a scientifically justified fixed
model. IQ-TREE checkpoints are preserved; `--redo` explicitly restarts and
replaces results. Use a new output prefix when data or model assumptions change.

FastTree uses `-nt -gtr -gamma` for DNA and `-lg -gamma` for proteins. Its search
uses the CAT approximation before Gamma20 rescaling; it is not a full Gamma-model
ML search. JTT is FastTree's protein default (no `-jtt` flag); nucleotide JC is
selected by `-nt` without `-gtr`. FastTree defaults to **SH-like local support**
in the **0–1** range, not UFBoot percentages. Use workload-based benchmarking,
not a universal taxon threshold, to decide between inference tools.

### 4. Interpret support, rooting and outputs

Read `.iqtree` and `.log`, not only the image. In the default combined test,
IQ-TREE writes slash-delimited SH-aLRT/UFBoot labels; verify their order in the
report. UFBoot >=95 and SH-aLRT >=80 are common screening thresholds, not proof
of correctness or posterior probabilities. Do not apply a standard-bootstrap
70% rule to UFBoot. Model misspecification, alignment error and sampling can
produce strongly supported wrong branches; consider `-bnni` when assessing
model violations and UFBoot overestimation.

Reversible substitution models do not infer the evolutionary root. The helper
preserves the supplied tree orientation unless an explicit `--outgroup` or
`--midpoint` display choice is made. Select outgroups biologically; midpoint
rooting assumes sufficiently clock-like path lengths and is only a heuristic.
`--outgroup` passes one taxon to IQ-TREE and verifies that it exists. Display
rerooting preserves support labels on their splits. Neither display option
creates a dated phylogeny or overwrites the original inferred tree.

Primary outputs are `<prefix>.treefile` (IQ-TREE) or `<prefix>.tree` (FastTree),
the alignment, and IQ-TREE's `.iqtree`, `.log`, `.ckp.gz`, plus optional support
files. `.model.gz` caches model-selection work; it is not an estimated model
parameter report. Keep these artifacts together.

### 5. Summarize and visualize

ETE must read combined IQ-TREE labels as internal **names** with `format=1`;
its default numeric-support parser cannot parse `95.2/99`. The helper renders
the original text labels, not ETE's default support value. Branch statistics
use all original non-root edges, including zero-length edges, without rerooting.

```python
# Run from the skill directory in the ETE3 environment.
import sys
sys.path.insert(0, "scripts")
from phylogenetic_analysis import load_tree, tree_summary

tree = load_tree("results/sequences.treefile")
print(tree.get_leaf_names())
print(tree_summary("results/sequences.treefile"))
```

To request rendering, omit `--no-visualization`. ETE3 needs PyQt5; on a headless
host the appropriate Qt platform plugin must be available. Missing rendering
support leaves the inferred tree intact and is reported explicitly. Qt rendering
was not exercised in this review. FigTree or iTOL can also display Newick; verify
that a viewer retains the support labels and branch-length scale.

## Official references

- [MAFFT manual](https://mafft.cbrc.jp/alignment/software/manual/manual.html) and [threading](https://mafft.cbrc.jp/alignment/software/multithreading.html)
- [IQ-TREE quickstart](https://iqtree.github.io/doc/Quickstart), [current releases](https://github.com/iqtree/iqtree3/releases), [model assumptions](https://iqtree.github.io/doc/Assessing-Phylogenetic-Assumptions)
- [FastTree manual](https://morgannprice.github.io/fasttree/)
- [ETE3 tree parsing](https://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html) and [drawing](https://etetoolkit.org/docs/latest/tutorial/tutorial_drawing.html)
- [trimAl manual](https://vicfero.github.io/trimal/)
- [pyfamsa documentation](https://pyfamsa.readthedocs.io/en/stable/) and [source](https://github.com/althonos/pyfamsa/); FAMSA — Deorowicz, Debudaj-Grabysz & Gudyś, *Scientific Reports* 2016
- [pyhmmer documentation](https://pyhmmer.readthedocs.io/en/stable/)
