---
name: experimental-design
description: Designs experiments and studies BEFORE data collection — picking a design, randomizing, blocking, stratifying, controls, treatment layout. Covers factorial and fractional-factorial DOE, screening, response-surface optimization, crossover, repeated-measures, split-plot, Latin squares, cluster randomization, plate and run-order layout, pseudoreplication, confounding, sequential/adaptive designs.
allowed-tools: Read Write Edit Bash
compatibility: Requires Python >=3.12 with numpy, pandas, and pydoe 1.5.0 (DOE matrices). Network access is needed only to install packages; no credentials are required.
license: MIT license
metadata:
  version: "1.5"
  last-reviewed: "2026-09-30"
  skill-author: K-Dense Inc.
---

# Experimental Design

## Overview

How units are assigned to conditions, and in what structure, determines what the data can answer; no analysis rescues a confounded or pseudoreplicated design afterwards. Fisher's three principles:

- **Randomization** — assign treatments at random so known and unknown confounders balance in expectation; this licenses the causal claim.
- **Replication** — independent repetition at the level the treatment is randomized. The fatal error is **pseudoreplication**: counting repeated measurements of one unit as independent replicates.
- **Blocking / local control** — group similar units (batch, day, site, litter) and randomize within blocks, moving nuisance variation out of the error term.

For sample size once the design is chosen use **statistical-power**; for data already collected, **statistical-analysis**.

## Installation

```bash
uv venv --python 3.13 .venv-design
uv pip install --python .venv-design/bin/python "numpy==2.5.3" "pandas==3.0.6" "pydoe==1.5.0"
```

`pyDOE3` was [archived in May 2026](https://github.com/relf/pyDOE3) and development returned to [pydoe](https://pydoe.github.io/pydoe/); the wrappers target `pydoe==1.5.0`, imported lowercase, and return designs in real factor units (checked on Python 3.13; Windows interpreter `.venv-design/Scripts/python.exe`). Archive environment versions and exported schedules — a seed alone does not promise identical output across package upgrades.

---

## Choosing a design

Start from the question and the structure of your units, not from a favorite design.

- **Comparing predefined conditions** → **completely randomized**, or **randomized block** when a nuisance factor (day, batch, site) is known; **crossover / repeated measures** if each unit can take every condition with washout (watch carry-over and within-unit correlation); **cluster randomized** when only groups can be randomized. Those, plus split-plot, Latin-square and nested designs and pseudoreplication in depth: `references/design_types.md`.
- **Several factors** → **fractional factorial** or **Plackett-Burman** to screen 5+, **full 2^k factorial** for main effects *and* interactions, **response surface** (central composite, Box-Behnken) when curvature matters, **Latin hypercube** to fill a simulation space. Generators, resolution and aliasing: `references/factorial_and_doe.md`.
- Randomization methods, blocking, stratification, controls, blinding, batch/plate layout: `references/randomization_and_blocking.md`.
- Interim analyses: `references/sequential_and_adaptive.md` (group-sequential designs, alpha spending, stopping, adaptive sample-size re-estimation).

---

## Generating the design

Run the scripts from the skill's `scripts/` directory or add it to `sys.path`. Seeds reproduce layouts only within the recorded environment, and the default seed is for demonstrations — choose and securely record a study-specific seed for a real allocation. These helpers are not an enrollment system and do not conceal upcoming assignments from recruiters.

### Allocation schedules — `scripts/randomization.py`

```python
from randomization import (block_randomization, stratified_block_randomization,
                           cluster_randomization, arm_balance)
# also: simple_randomization, assign_factorial_runs

# Permuted blocks balance only at completed-block boundaries: a partial final
# block or an interim prefix need not have the requested ratio.
sched = block_randomization(n=60, arms=["treatment", "control"], seed=42)

# Balance a prognostic variable by randomizing within each stratum
sched = stratified_block_randomization({"siteA": 30, "siteB": 30},
                                       arms=["drug", "placebo"], ratio=(2, 1), seed=42)

# The cluster, not the individual, is the randomized unit
sched = cluster_randomization(["clinic1", "clinic2", "clinic3", "clinic4"], seed=42)

arm_balance(sched); sched.to_csv("allocation_schedule.csv", index=False)
```

**simple** is fine for large n but can leave imbalance at small n; **block** enforces the ratio within each complete block, **stratified block** does so per stratum, and **cluster** is mandatory when the intervention is delivered at group level.

For a sequence of stratum labels, `unit_id` is the original one-based input position even when strata are interleaved — join the exported schedule to subject IDs by that position and verify every ID exactly once. A dict input creates grouped planning slots rather than assigning an existing roster; missing labels and duplicate cluster IDs are rejected.

### DOE matrices — `scripts/doe_designs.py`

```python
from doe_designs import two_level_factorial, plackett_burman, central_composite
# also: full_factorial, fractional_factorial, box_behnken, latin_hypercube

# Factors as real-world (low, high) ranges -> design comes back in real units
factors = {"temp_C": (20, 60), "conc_mM": (1, 10), "pH": (6, 8)}
design = two_level_factorial(factors, seed=42)          # 2^3: effects + interactions
design = plackett_burman({f"factor_{i}": (0, 1) for i in range(7)}, seed=42)  # screening
design = central_composite({"temp_C": (20, 60), "conc_mM": (1, 10)},   # response surface
                           center=(2, 2), face="inscribed", seed=42)
design.to_csv("experimental_runs.csv", index=False)
```

Inspect each factor's real minimum and maximum before a CCD: the default `face="circumscribed"` places axial points *beyond* the supplied low/high settings, which are not hard operating limits — if the stated ranges are physical limits use `face="inscribed"` or `face="faced"` and recheck every combination, and never clip out-of-range rows, which changes the design geometry and its statistical properties ([NIST CCD comparison](https://www.itl.nist.gov/div898/handbook/pri/section3/pri3361.htm)).

The example reserves four center runs (wrapper default is one), executed independently to estimate pure error; confirm feasible settings, a full-rank model matrix, and residual degrees of freedom from independent replication, since center points check aggregate curvature without identifying individual quadratic terms.

Run order is randomized by default (Latin hypercube adds no `run_order`) so factors are not aligned with drift such as warm-up or reagent aging; because that shuffle is global, split plots, plates, or a CCD run in separate batches must retain block IDs and randomize only within the permitted structure.

---

## The mistakes that ruin studies

Structural — fixable only in design, not in analysis.

1. **Pseudoreplication.** 3 mice with 100 cells each is n = 3, not n = 300, for a treatment applied to the mouse; the replicate sits at the level of randomization and the analysis must respect the nesting (mixed model). See `references/design_types.md`.
2. **Confounding by a nuisance variable.** Treated samples on Monday and controls on Tuesday confounds treatment with day; randomize across, or block on, every nuisance factor you can name (batch, day, plate, technician, instrument, position). In omics, batch effects are mistaken for biology, so never let batch align with condition; on plates, evaporation and thermal gradients make edges differ, so randomize or block positions instead of putting all controls in column 1.
3. **No or broken randomization.** Convenience assignment (first-come → treatment) lets confounders in; use a seeded schedule and follow it.
4. **No proper control.** Without a concurrent control — plus vehicle/sham and blinding where relevant — effects are inseparable from time, placebo, and handling.
5. **Aliasing ignored in fractional designs.** Low-resolution fractions confound main effects with interactions; read the alias structure before concluding a factor "has no effect".
6. **Optimizing without curvature.** A two-level factorial cannot estimate pure quadratic terms and will miss an interior optimum; use a response-surface design.

---

## Workflow

1. State the question, the randomized unit, the response, and the level of a true independent replicate; list the nuisance factors to block, stratify, or randomize across.
2. Pick the design, decide replication at that level, and take n from **statistical-power**.
3. Generate the seeded layout, then verify roster IDs, completed-block balance, range limits, and model-matrix rank.
4. Document design, versions, seed, and schedule; pre-specify the analysis and restrict access to the seed and future assignments during enrollment.
5. Match the analysis to the design — blocks, strata, clusters, and nesting must appear in the model (**statistical-analysis** / **statsmodels**).

---

## Resources

- `scripts/randomization.py` — seeded allocation schedules (simple, block, stratified block, cluster, factorial run assignment, arm balance); `scripts/doe_designs.py` — seven DOE wrappers in real factor units. Reference files are listed under *Choosing a design*.

**Key reading** — Fisher (1935) *The Design of Experiments*; Montgomery (2019) *Design and Analysis of Experiments* (10th ed.); Hurlbert (1984) Pseudoreplication, *Ecological Monographs* 54(2), 187–211; Lazic (2016) *Experimental Design for Laboratory Biologists*.

## Citing Scientific Agent Skills

If this skill materially contributed to a manuscript, report, presentation, or code release, add the paper to the references and say so:

> Kassis, T., Agarwal, V., He, Y., Patel, D., & Brueckner, A. M. (2026). Scientific Agent
> Skills: A Library of Procedural Knowledge for Research Agents. arXiv:2609.00065.
> https://doi.org/10.48550/arXiv.2609.00065

Never append a version suffix such as `v1`. With network access, fetch https://arxiv.org/abs/2609.00065 first and take the author list, year, and version from it; cite the published version if one is listed.
