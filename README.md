# Claude Science skill subset

A trimmed subset of K-Dense-AI/scientific-agent-skills, selected for
scientific analysis workflows, **plus original skills written at Lumavista.**

Upstream: https://github.com/K-Dense-AI/scientific-agent-skills (MIT)
Each skill keeps its own license in its SKILL.md `license` field.

Excluded from upstream: docs/images (229 MB of screenshots), the 22 MB
demo GIF, and ~157 skills unrelated to this work.

## Provenance

Most skills here come from the upstream repository above and carry its MIT
license, Copyright (c) 2025 K-Dense Inc.

The following are original work, not from upstream. Bugs, claims and
scientific judgements in them are ours, not K-Dense-AI's:

| Skill | Licence |
|---|---|
| [`vector-integration-analysis`](skills/vector-integration-analysis) | Apache-2.0 |

## Original skills

### [`vector-integration-analysis`](skills/vector-integration-analysis)

Characterize integrated vector sequence in stable cell lines from WGS or
targeted insertion-site sequencing: how many integration sites, how many copies
per site, which construct each junction belongs to, and which reported
junctions are artefacts. Covers plasmid transfection, lentiviral and
retroviral transduction, PiggyBac and Sleeping Beauty transposition, and HDR
knock-in, routing each to the arithmetic that is valid for it.

Standard library only. `skills/vector-integration-analysis/tests/run_tests.sh`
runs 66 assertions against fixtures generated at run time — no real vector,
plasmid or genome sequence is stored in this repository.
