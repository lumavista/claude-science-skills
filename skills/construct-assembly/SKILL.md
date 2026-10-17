---
name: construct-assembly
description: Turns antigen or protein designs into orderable constructs and screens them for sequence-level liabilities — linkers, tags, trimerisation and display domains, signal peptides, junction scanning, ThermoMPNN stability, and developability (free cysteines, designed-in sequons, deamidation motifs, hydrophobic patches, GRAVY, pI). Inverse folding, de novo design and refolding go to the tamarind skill.
license: MIT
compatibility: Python 3.10+, standard library only for the developability screen — nothing here needs a GPU or network. ThermoMPNN is optional, MIT, installed from source by git SHA (no release tags), and uses a GPU when torch finds one while staying fast on CPU for small mutant sets.
metadata:
  version: "1.0"
  skill-author: Drafted for vaccine-design subset
  last-reviewed: "2026-10-16"
---

# Construct Assembly

The step between "we have a designed sequence" and "we can order this". It is
mostly unglamorous — linkers, tags, cysteines, sequons — and it is where more
constructs fail than at the design stage.

## Scope

**This skill does:** linker and domain architecture, junction checking,
stability prediction on point mutants, local sequence developability, and the
final ranked construct table.

**This skill does not** run inverse folding (ProteinMPNN, LigandMPNN,
SolubleMPNN), de novo backbone generation (RFdiffusion), or structure
prediction for refold validation. Those are GPU jobs and belong to `tamarind`,
which discovers tool schemas live — reimplementing submission here would drift
out of sync with it. Follow `tamarind`'s own rule: prefer the dedicated
structure-prediction skills where they can run, and say which of cloud
execution, insufficient local VRAM, or batch size makes Tamarind necessary before
submitting.

## First: does this construct have to fold?

Applying the wrong QC produces a page of meaningless failures, and this is the
most common waste of effort in the step.

| Construct type | Must fold? | QC focus |
|---|---|---|
| Conformational antigen | Yes | epitope geometry, oligomeric state, glycosylation, stability |
| Mosaic or epitope string | No | processing, junctions, expression, aggregation |
| Natural-sequence cocktail | Yes, per component | same as the parent antigen |
| Combined B- and T-cell construct | Partly | the folded module must fold; the appended module must not perturb it |

`immunoinformatics` step 9 states which kind it produced — a cocktail of real
sequences or a synthetic mosaic. A mosaic is meant to be expressed and
degraded; a low pLDDT on one says nothing and should not be reported as a
finding.

## Workflow

### 1. Design the architecture

See `references/linkers-and-domains.md` for linker families and lengths,
trimerisation and display domains, tag choice and placement, and a default
assembly order.

Three decisions that are easy to get wrong:

- **Size the linker to the geometry**, measuring terminus-to-terminus distance
  where structures exist, rather than defaulting to `(G4S)3`.
- **Foldon is immunogenic.** Responses against the trimerisation domain are
  routinely observed. Plan a foldon-only control or make it cleavable.
- **Place tags away from the epitope**, and decide explicitly whether they
  cleave. A non-cleavable tag on an immunogen is part of the immunogen.

### 2. Scan the junctions

Every join creates peptides that exist in neither parent. Use the junction
checker in `immunoinformatics`:

```bash
python ../immunoinformatics/scripts/mosaic_design.py evaluate \
  --fasta work/strain_panel.fasta \
  --designs work/constructs.fasta \
  -k 9 --out work/junction_k9
```

Run again at `-k 15` for the class II register. Any k-mer reported as absent
from the panel is a peptide the construct can present that no real strain
does — it diverts the response without adding protection. Review them; do not
filter them away silently.

### 3. Predict stability for point mutants

For stabilising substitutions — disulfide staples, cavity filling, proline
substitutions at a hinge — rank with **ThermoMPNN** (MIT, same architecture
class as ProteinMPNN). SHA-pin it; there are no release tags.

ThermoMPNN uses a GPU when torch finds one and CPU otherwise; CPU is adequate
for a handful of mutants, so the GPU only matters for a full site-saturation
scan, where batching many independent forward passes is the win. Nothing in
this skill's own scripts needs a GPU — the developability screen is standard
library only. Select the device instead of hardcoding it, keeping `SCI_DEVICE`
usable for a reproducible CPU rerun:

```python
import os, torch
device = os.environ.get("SCI_DEVICE") or ("cuda" if torch.cuda.is_available() else "cpu")
```

A `False` from `torch.cuda.is_available()` on a GPU machine is usually a
CPU-only wheel (the Windows PyPI default), not a hardware fault.

`ThermoMPNN-D` handles double mutants and `ThermoMPNN-I` handles indels. Pass
through the authors' own caveat on the indel model: it has had only limited
in-silico validation, so those numbers deserve visibly wider error bars than
point-mutant ΔΔG.

ΔΔG predictions are a ranking aid for choosing what to express, not a stability
measurement. `adaptyv` is where thermostability is actually measured.

### 4. Screen developability locally

```bash
python scripts/developability.py \
  --fasta work/constructs.fasta \
  --parent work/parent.fasta \
  --ph 7.4 \
  --out work/developability.tsv
```

Checks free cysteines, N-glycosylation sequons (flagging those **added**
relative to the parent), deamidation and isomerisation motifs, Met/Trp
oxidation counts, hydrophobic patches, GRAVY, pI and net charge.

Two flags deserve particular attention:

- **Odd cysteine count.** An unpaired thiol drives disulfide scrambling,
  dimerisation and heterogeneity. Surface redesign introduces these by
  accident routinely.
- **Added sequons.** Resurfacing creates N-X-S/T by accident, and a new glycan
  on or beside the target epitope can mask the very site the immunogen exists
  to present. When the construct differs in length from the parent the script
  reports a count delta and declines to attribute positions, because an indel
  has shifted the numbering and a positional claim would be wrong.

Any actual glycan engineering — adding, removing, occupancy, shield geometry —
goes to `glycoengineering`, which owns that logic. This screen only tells you a
sequon appeared.

**Every flag is a liability to weigh, not a verdict.** A construct with three
deamidation motifs is not disqualified; it is a construct whose stability study
should look for deamidation. Report them rather than filtering.

### 5. Validate the assembled construct

Refolding goes to `tamarind`. What matters in the result:

- RMSD to the intended backbone — and, for a conformational antigen,
  **epitope-local RMSD reported separately**. A global RMSD can look fine while
  the epitope has moved.
- pLDDT globally and at the epitope. Per the Antigen-Structure contract, do not
  interpret any region below ~70.
- Oligomeric state, which requires folding the assembly rather than the
  protomer.
- **Validate after assembly, not before.** An appended trimerisation domain or
  tag can perturb the module it is attached to.

Fold a design you expect to fail alongside the real ones. Without a negative
control a column of pLDDT values cannot distinguish a good round from a
saturated metric.

### 6. Ranked construct table

The deliverable. Per construct:

| Field | Why |
|---|---|
| ID and full sequence | orderable |
| Parent and what changed | reviewable |
| Rationale per element | linker length, domain choice, tag placement |
| Predicted breadth vs the strain panel | the actual objective |
| Structural evidence | RMSD, pLDDT, epitope-local metrics, source of each structure |
| Presentation evidence | coverage and allele breadth, for T-cell designs |
| Liabilities | from step 4, unfiltered |
| **Discriminating experiment** | the assay that would falsify this fastest |

The last column is the most valuable and the most often omitted. The pipeline's
job is to propose a small set of discriminating experiments — make that
explicit per construct rather than leaving it to be inferred from the
rationale. `adaptyv` is where several of those assays can actually be run.

## Interpretation

- **Ranking among the survivors is weakly informative.** This pipeline narrows
  thousands to tens; the tools discriminate plausible from implausible far
  better than they discriminate two plausible designs. Say so, so that
  position 1 is not over-trusted against position 7.
- **Carry several candidates.** Developability metrics are predictors, not
  verdicts, and ranking on any single score discards designs that fail one
  proxy and succeed at the real thing.
- **Sequence-level screening sees no structure.** A hydrophobic patch that is
  buried is not an aggregation risk; this screen cannot tell. Read it alongside
  accessibility from `conservation-mapping` or `biopython`.
- **A design is a hypothesis until assayed.**

## Reference files

- `references/linkers-and-domains.md` — linker families and sizing,
  trimerisation and display domains, tags, signal peptides, assembly order
- `immunoinformatics/references/tool-licensing.md` — licensing and pins for
  ThermoMPNN, the MPNN family and the folding models referenced here
- `immunoinformatics/references/compute-routes.md` — what runs on CPU, what
  goes to Tamarind, and the free hosted routes

## Official references

- ThermoMPNN — Dieckhaus et al., *PNAS* 2024, doi:10.1073/pnas.2314853121
- Kyte & Doolittle hydropathy, *J Mol Biol* 1982
- Foldon / T4 fibritin trimerisation — Güthe et al., *J Mol Biol* 2004
- SpyTag/SpyCatcher — Zakeri et al., *PNAS* 2012
- `tamarind` (this collection) — MPNN design, RFdiffusion, refolding
- `glycoengineering` (this collection) — sequon engineering and glycan shield
- `adaptyv` (this collection) — experimental characterisation
