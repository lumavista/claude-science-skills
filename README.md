# Vaccine design skill subset

A trimmed subset of K-Dense-AI/scientific-agent-skills, selected for
virology / immunology / antigen & epitope / vaccine-design workflows, plus
local fixes and one skill written from scratch.

Upstream: https://github.com/K-Dense-AI/scientific-agent-skills (MIT)
Each skill keeps its own license in its SKILL.md `license` field.

Excluded from upstream: docs/images (229 MB of screenshots), the 22 MB
demo GIF, and ~157 skills unrelated to this work. That is what pushes
the upstream repo past the 100 MB import limit.

## Skills

19 skills. Several reference each other by name — `immunoinformatics` hands
off to `pathogen-variant-surveillance` for circulating variants, to
`glycoengineering` for glycan shielding, and to `scientific-writing` for the
write-up — so removing one degrades the others' workflows.

| | |
|---|---|
| Epitopes & antigen | `immunoinformatics`, `glycoengineering`, `esm`, `tamarind`, `adaptyv` |
| Sequences & evolution | `biopython`, `phylogenetics`, `pathogen-variant-surveillance` |
| Structure & dynamics | `molecular-dynamics` |
| Transcriptomics | `scanpy`, `anndata`, `pydeseq2`, `pathway-enrichment` |
| Literature & writing | `paper-lookup`, `scientific-writing`, `scientific-visualization` |
| Study design | `experimental-design`, `statistical-power` |

`immunoinformatics` is not from upstream. It was written for this repo and
covers MHC class I/II binding prediction, B-cell epitope prediction, IEDB
cross-referencing, conservation across circulating variants, self-similarity
screening against the human proteome, and HLA population coverage. Its
`scripts/epitope_triage.py` encodes the step ordering and writes a
`.meta.json` sidecar per output so a result carries its tool versions.

`database-lookup` from upstream is deliberately not included — it overlaps the
Claude Science Featured connectors, which are better at the same job. Skills
that reference it by name resolve to those connectors.

## Removed as redundant with Claude Science built-ins

- `scvi-tools` — Claude Science ships a featured scvi-tools skill covering the
  same models (scVI, scANVI, totalVI, MultiVI). Other skills here that reference
  `scvi-tools` by name resolve to the built-in.
- `literature-review` — Claude Science ships a featured literature review skill.
  `paper-lookup` is kept: it is raw multi-API search (18 scholarly APIs with
  provenance), not review methodology, and does not collide.

## Local modifications

### `tamarind`

Description narrowed to de novo design, binder design, antibody design and
developability, with trigger keywords cut from 27 to 14 and a `## Scope`
section added. As written upstream it competed with the Claude Science
featured skills for plain structure prediction, docking and embeddings
(AlphaFold2, Boltz, Chai-1, ESMFold2, OpenFold3, DiffDock, ESM-2), which are
better choices when they apply. The skill now defers to them and states
which of three conditions justifies a cloud submission — explicit request,
no local GPU, or batch size — before submitting.

### `phylogenetics` (v1.6)

Upstream assumes MAFFT and IQ-TREE executables, which cannot run on a
Windows host: bioconda publishes no win-64 builds, and the vendor Windows
binaries of mafft, muscle, iqtree and trimal ship non-relocatable PE images
that fail to load under forced ASLR with `WinError 623`. The skill now
states platform availability, carries a per-role fallback matrix, and adds a
PE-header pre-screen (`DYNAMICBASE` plus a base-relocation directory) so a
candidate binary can be rejected in seconds instead of after a download.

Compiled Python extension modules are ASLR-safe where standalone executables
are not, so the skill prefers a wheel-based aligner: `pyfamsa` for MSA,
`pyhmmer` for profile alignment, `parasail` or Biopython for pairwise work.
`scripts/phylogenetic_analysis.py` gained an `--aligner famsa` backend and
runs end to end with FastTree where IQ-TREE is unavailable. Losing IQ-TREE
costs ModelFinder, ultrafast bootstrap and SH-aLRT, and the skill says so
rather than letting a write-up imply a full ML analysis.

Four defects were found by running it, not by reading it:

- `_run_to_file` called `os.replace()` while the `NamedTemporaryFile` handle
  was still open. POSIX allows it; Windows refuses with `WinError 32`, and
  the `finally` clause then raised a second one. **This is upstream code and
  breaks every MAFFT and IQ-TREE output path on Windows**, not just the new
  FAMSA path. A regression test now patches `os.replace` to refuse a rename
  while a handle to the source is open, reproducing the failure on any
  platform — it was invisible to testing on Linux.
- `check_dependencies` was dead code, never called, and hardcoded `mafft` as
  required even for runs that never invoke it. It is now called and gated on
  `--aligned` and `--aligner`, and accepts an absolute path in place of a
  name on `PATH`.
- The tree leaf-count check depended on ete3, which is not installable on
  some hosts, so it degraded to a warning — and a skipped check reads like a
  passed one in a log. It now prefers dendropy and says explicitly when the
  check did not run.
- The amino-acid alphabet omitted J (Leu/Ile), U (selenocysteine) and O
  (pyrrolysine), so one legitimate residue aborted a 469-sequence run.
  FAMSA is the real constraint — verified against pyfamsa 0.7.0, it accepts
  B, X, Z and `*` but rejects J, U and O — so the helper refuses by default,
  names the substitution that would fix it, and applies `J→X, U→C, O→K` only
  under `--map-unsupported`, reporting what it changed. Editing the source
  FASTA to get past this is the wrong fix: a declared substitution in a run
  log is recoverable, a quietly altered input file is not.

The first two are worth reporting upstream — they affect any Windows user of
that skill.

## Specialists

`specialists.yaml` records the Claude Science specialists built on these
skills: each one's description, instructions, selected skills and connectors.
Six specialists:

| Specialist | Covers | Skills |
|---|---|---|
| **Claude Science** (`OPERON`) | Default profile — the full skill and connector set | all |
| **Immuno-Designer** | T-/B-cell epitope prediction, conservation, population coverage, candidate triage | `immunoinformatics`, `glycoengineering`, `tamarind`, `biopython`, `phylogenetics`, `pathogen-variant-surveillance` |
| **Literature-Search** | Provenance-tracked multi-API search, evidence-traceable write-ups, study sizing | `paper-lookup`, `literature-review`*, `scientific-writing`, `scientific-visualization`, `experimental-design`, `statistical-power` |
| **Pathogen-Surveillance** | Circulating lineages, growth, mutations, phylogenies (GenSpectrum/LAPIS) | `pathogen-variant-surveillance`, `phylogenetics`, `biopython`, `scientific-visualization` |
| **Immune-Profiler** | Systems-vaccinology transcriptomics: composition, DE, enrichment | `scanpy`, `anndata`, `pydeseq2`, `pathway-enrichment`, `statistical-power`, `scvi-tools`*, `scientific-visualization` |
| **Antigen-Structure** | Structure retrieval, epitope accessibility, glycan shielding, embeddings, light MD, cloud design | `tamarind`, `esm`, `glycoengineering`, `molecular-dynamics`, `biopython`, `adaptyv`, `scientific-visualization` |

\* `literature-review` and `scvi-tools` are Claude Science built-ins, not skills
in this repo; a repo re-import clears the selection, so re-select them.

Claude Science has no specialist import or export, so this file is
documentation rather than something you can load. It is read from the app's
local store (`user_agents` and `directory_attachments` in
`orgs/<org>/the app database`), sanitized for publication (machine identifiers and
local UUIDs removed, absolute paths replaced with `%USERPROFILE%`), and then
hand-revised — so it is no longer a byte-for-byte app export. The current
revision splits each instruction field into a portable **scientific contract**
and a fenced **host notes** block (findings from one Windows host, safe to
delete when you move machines), and adds the `genomes` connector to
Immuno-Designer and the `statistical-power` skill plus the `genes-ontologies`
connector to Immune-Profiler. Re-enter a specialist by hand under
**Settings > Specialists**, and note that re-importing this repo clears each
specialist's skill selections, so re-select them afterwards.

`vector-integration-analysis` is selected in the default profile but in none of
the five focused specialists — an open scope decision, left as-is.

Two findings from the Immune-Profiler scientific contract are worth reading
even if you never use these specialists. They are portable statistical
properties — reproduced in simulation, not artifacts of any one host:

- **Pseudoreplication.** On a simulated null with 8 donors and no true
  condition effect, a per-cell Wilcoxon test called roughly 1,390 of 2,000
  genes significant at FDR 0.05; donor pseudobulk with PyDESeq2 called 0–3.
  With donor variance removed, the per-cell test called 0, which locates the
  cause in donor variance rather than in the test. Median intraclass
  correlation 0.032 gives a design effect of about 33, so ~4,000 cells per
  arm are worth ~123 effective observations.
- **Power.** The standard error of an arm difference is
  `sigma_donor * sqrt(2/k) / ln(2)` in log2 units for k donors per arm. At
  donor SD 0.35, pseudobulk recovered ~34% of true log2FC 1.0 genes with 4
  donors per arm and ~2% with 2 per arm. Below 4 per arm, gene-level
  conclusions are not supportable, and at 4 per arm a gene not called is
  uninformative rather than unchanged.

Regenerate `specialists.yaml` after any specialist change — it goes stale
silently, and nothing in the app warns you.
