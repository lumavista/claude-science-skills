# Tool Licensing and Versions

Verified 2026-10-16 against primary sources (PyPI metadata, raw LICENSE files,
git tag refs, DTU service pages). Re-verify anything marked UNVERIFIED before
depending on it.

SKILL.md states the core point already — the DTU predictors are academic-only
and MHCflurry is the permissive class I substitute. This file carries the full
table, because the same trap recurs across structure and design tooling and is
not obvious from any tool's front page.

## Two distinctions that cause most of the errors

**Code licence is not weights licence.** Chroma is Apache-2.0 code with
academic-only weights. ColabFold is MIT glue over weights with separate terms.
RFdiffusion deliberately goes the other way — its BSD text explicitly covers
the weights. Recording one licence per tool will encode a falsehood.

**"Freely available" is not open source.** Across immunology tooling this
nearly always means free for academics. NetMHCpan, NetMHCIIpan, MixMHCpred,
MixMHC2pred and BigMHC are all academic-gated.

## MHC binding and presentation

| Tool | Pin | Licence | Class | CPU |
|---|---|---|---|---|
| **MHCflurry** | `mhcflurry>=2.3,<2.4` (2.3.14) | **Apache-2.0** | I | Yes, incl. ARM |
| **MHCnuggets** | `mhcnuggets==2.4.1` | **Apache-2.0** | **I and II** | Yes |
| **pVACtools** | `pvactools==7.1.5` | BSD-3-Clause-Clear | orchestrator | varies |
| NetMHCpan | `4.2e` — exact, incl. letter | Academic; commercial = paid DTU licence | I | Yes |
| NetMHCIIpan | `4.3k` — exact, incl. letter | Academic | II | Yes |
| MixMHCpred | `3.0`, git SHA | Academic (LICR) | I | Yes |
| MixMHC2pred | `2.0`, git SHA | Academic (LICR) | II | Yes |
| BigMHC | git SHA | Bespoke JHU academic licence | I | Yes, no GPU needed |

**NetMHCpan is at 4.2, not 4.1.** The current package is `4.2e`; the service
page records an update of 2026-07-06 (NetMHCpan-4.2, *Front Immunol* 2025,
PMC12367478). SKILL.md's compatibility line and reference list still name 4.1 —
treat 4.2e as current and pin the full letter-suffixed string. NetMHCIIpan is
at `4.3k`; the letter suffixes are real releases, and the `static` versus
non-`static` build distinction matters for containerisation.

**MHCnuggets is the only permissive class II option**, and it has been
unmaintained since 2023. That is worth saying out loud to anyone doing
commercial CD4 work: the permissive path exists but is thin, and class II
prediction is the weaker half of the field to begin with.

**Pin MHCflurry's weight release, not just the package.** Weights version
independently: `mhcflurry downloads fetch models_class1_presentation`, then
record the release from `mhcflurry downloads info` and pass `--model-release`.
Otherwise the same package version gives different numbers on different
machines, silently.

## B-cell epitope prediction

| Tool | Pin | Licence | Hosted |
|---|---|---|---|
| BepiPred 3.0 | `3.0b` (source, DTU request form) | Assume DTU academic — UNVERIFIED | DTU server, 50 seqs / 300k aa per job |
| DiscoTope 3.0 | `3.0` | Assume DTU academic — UNVERIFIED | DTU server |
| SEMA 2.0 | git SHA | **NO LICENCE FILE — all rights reserved** | `sema.airi.net` |

**SEMA 2.0's licence absence is verified, not assumed** — both `LICENSE` and
`LICENSE.md` 404 in the repo. No licence means no rights granted, not "open".
Use the web interface, do not vendor the code, and do not build a commercial
pipeline on it without contacting the authors. SKILL.md names SEMA as a
conformational option; this is the constraint on that recommendation.

Make "confirm a LICENSE file exists" an explicit step for any 2025–26
PLM-based epitope predictor. That cohort very often ships without one.

**DiscoTope 3.0 works on predicted structures** — it uses inverse-folding
latent representations rather than requiring experimental coordinates, which
matters when the only structure available is an AlphaFold model.

## Structure prediction

| Tool | Pin | Code | Weights | CPU |
|---|---|---|---|---|
| **Boltz** | `boltz==2.2.1` | MIT | **MIT** | Technically, but upstream calls it significantly slower |
| **OpenFold3** | `openfold3==0.5.0` | Apache-2.0 | **Apache-2.0** | No |
| **Chai-1** | `chai_lab==0.6.1` | Apache-2.0 | **Apache-2.0** | No — needs bf16 CUDA |
| ColabFold | `colabfold==1.6.3` | MIT | AF2 weights CC-BY-4.0 | No |
| ESM family | `esm==3.4.1.post1` | MIT | **MIT** | ESM-2 8M–650M only |
| ESM-2 legacy | `fair-esm==2.0.0` | MIT | MIT | Yes, small/mid |
| AlphaFold3 | — | Apache-2.0 | **Non-commercial, no redistribution, approval-gated** | No |

**AF3 weights bar commercial use entirely**, including research performed on
behalf of a commercial organisation, and bar using AF3 output to train
competing structure predictors. A university researcher on an industry-funded
contract is not covered. For commercial work the substitutes are Boltz-2 (MIT),
OpenFold3 (Apache-2.0) or Chai-1 (Apache-2.0).

**Chai-2 and Chai-3 are not open.** No public weights; `chai-lab` on PyPI is
Chai-1 and has not moved since March 2025. Do not pin a "Chai-2" package.

**ESM naming has shifted.** EvolutionaryScale now presents as Biohub; the API
moved from `forge.evolutionaryscale.ai` to `biohub.ai`, and ESM3's original
non-commercial Cambrian licence has been replaced by MIT. The `esm` skill's
`references/forge-api.md` predates this — check it against `biohub-platform.md`
before using an endpoint from it.

## Design and stability

| Tool | Pin | Code | Weights | CPU |
|---|---|---|---|---|
| **ProteinMPNN** | git SHA | MIT | MIT | **Yes, natively** |
| **LigandMPNN** | git SHA | MIT | **MIT** | **Yes, natively** |
| **ThermoMPNN** (+`-D`, `-I`) | git SHA | MIT | MIT | **Yes** |
| RFdiffusion | git SHA | BSD | **BSD covers weights** | No |
| RFdiffusion All-Atom | git SHA | BSD | **BSD covers weights** | No |
| RFantibody | git SHA | MIT | MIT | No — CUDA 11.8+, Linux |
| Genie2 | git SHA | Apache-2.0 | Apache-2.0 | No |
| Chroma | git SHA | Apache-2.0 | **Academic/non-profit, API key** | No |

**SolubleMPNN is not a separate project** — it is the `soluble_model_weights/`
directory inside the ProteinMPNN repo (`v_48_010.pt`, `v_48_020.pt`), same MIT
licence. People search for it as its own tool and conclude it is missing.

**Prefer LigandMPNN over ProteinMPNN.** It serves the ProteinMPNN models
through the same code and has materially better PDB handling: ProDy parsing
that preserves residue indices, chain letters and insertion codes, no `X`
padding across gaps, and residue-index selection syntax (`A23`, `B42D`).

**Chroma is the trap here** — "Apache-2.0" is true of the code and false of the
weights. Also, the PyPI package named `chroma` is an unrelated 2013 colour
library.

**Most of this section has no release tags.** ProteinMPNN, LigandMPNN, all
RFdiffusion variants, Genie2 and ThermoMPNN must be SHA-pinned.

## Glycan shield

| Tool | Licence | CPU |
|---|---|---|
| GlycoSHIELD | UNVERIFIED | **Yes** |
| GLYCO-2 | UNVERIFIED | Yes |
| Glycosylator (`glycosylator==5.6.3`) | UNVERIFIED (empty PyPI field) | Yes |
| glycowork (`glycowork==1.11.0`) | UNVERIFIED (empty PyPI field) | Yes |

`glycoengineering/references/glycoshield.md` pins GlycoSHIELD at commit
`f1ca31aa` and documents its CLI option-forwarding defect — that file is the
authority for running it. What is recorded here is only that its licence could
not be verified, which matters before any commercial use.

## Do not pin these

Verified unavailable or unconfirmable on 2026-10-16: **TransPHLA** (repo 404s),
**HLAthena** (no maintained repo or confirmed server), **LANL Mosaic Vaccine
Designer** and **EpiGraph** (liveness unconfirmed, web-only, no code release).

The last two are why `scripts/mosaic_design.py` implements the k-mer coverage
algorithm directly rather than depending on a hosted tool.
