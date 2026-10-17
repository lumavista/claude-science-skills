# Compute Routes: Local GPU, Local CPU, and Hosted

Immunogen work spans steps that are trivial on a laptop and steps that need a
datacentre. Being explicit about which is which is the difference between a
workflow that completes and one that stalls silently partway through.

The reference host for this collection is a **Windows 11 workstation with an
NVIDIA RTX 4060, 8 GiB VRAM** (driver 617.42, confirmed via `nvidia-smi` in
2026-10). That changes the routing advice this file used to give: folding,
inverse folding and embedding are now local operations for proteins of ordinary
size, and the cloud is for what does not fit rather than for everything.

Verify before trusting. `nvidia-smi` reports the card and driver;
`skills/esm/scripts/_device.py` reports what torch can actually reach. They
disagree more often than people expect, and the usual cause is a CPU-only torch
wheel, not a hardware fault.

## The one rule that saves the most compute

**Check whether the structure already exists before predicting it.** For a
known antigen an experimental structure or an AlphaFold DB model usually
exists, needs no GPU at all, and beats anything predicted fresh. Fetch it
through the Structures & Interactions connector. This is the step most often
skipped, and skipping it is the single most expensive mistake in the workflow.

## Three tiers

**Tier 1 — CPU, no caveat.** More of this workflow than people expect. Putting
these on a GPU buys nothing; most are bound by IO, branching or small matrices:

- MHC class I and II prediction (MHCflurry, MHCnuggets, and the DTU tools if
  licensed — small networks, CPU by design; a GPU is not the bottleneck at
  typical peptide counts)
- IEDB retrieval, self-similarity screening, population coverage prep
- Alignment and conservation (`phylogenetics`, `conservation-mapping`)
- Mosaic and cocktail design, junction checking (`mosaic_design.py`) — graph
  search over a strain panel, not arithmetic
- Relative solvent accessibility (freesasa, DSSP)
- Glycan shield grafting (GlycoSHIELD — avoiding MD is its whole contribution)

**Tier 2 — local GPU, with CPU as a correct but slow fallback.** These are the
steps that the 4060 changes from "hand off" to "run here":

| Step | Local route | 8 GiB notes |
|---|---|---|
| ESMC embeddings | `esm` skill, 300M/600M | comfortable; 6B does not fit |
| ESM3 generation | `esm3-sm-open-v1` (1.4B) | fits in bf16 (~2.8 GiB), not fp32 |
| Inverse folding | ProteinMPNN / LigandMPNN / SolubleMPNN | small models, fast |
| Point-mutant stability | ThermoMPNN (`construct-assembly`) | small; CPU is viable too |
| MD | OpenMM CUDA platform (`molecular-dynamics`) | the clearest win on this card |

**Folding is not in that table, and that is the surprise.** An 8 GiB card is
not what blocks it on native Windows — the Python wheels are. `esmfold2`
depends on flash-attn and Triton for its fused kernels, and neither publishes a
`win_amd64` wheel (flash-attn and DeepSpeed are sdist-only; upstream Triton has
no Windows build, only a community `triton-windows` fork). ColabFold-class
AlphaFold2 runs on JAX, and `jax-cuda12-plugin` / `jax-cuda13-plugin` ship
manylinux wheels only, so there is no CUDA JAX on native Windows at all. Boltz
is torch-based and so does import, but its fast `cuequivariance` kernels are
Linux shared objects, leaving `--no_kernels` at roughly half speed, and peak
VRAM is the binding limit at 8 GiB regardless.

So folding stays a Tier 3 hand-off here, for a platform reason rather than a
memory one. WSL2 recovers the JAX and Triton paths and is the right move if
local folding matters; it does not recover VRAM, so the 24-80 GiB models stay
out of reach either way.

Reduced precision is the default worth taking: on Ampere and newer, bfloat16
halves the weight footprint while keeping fp32's exponent range, so it needs no
loss scaling. For folding specifically, memory scales with sequence length
through the activations rather than with checkpoint size — so the limit reached
first is chain length. Fold a protomer rather than a trimer and state that
quaternary and trimer-interface epitopes are not covered by that model.

**Tier 3 — exceeds 8 GiB, so hand off.** Tamarind first, because it is already
wired in and schema-discoverable; free notebooks when Tamarind does not fit
(no account or budget, a one-off check, or a model Tamarind does not expose):

| Need | Route | Limits |
|---|---|---|
| Existing structure | Structures & Interactions connector | none — check this first |
| AF3-class, commercially usable | ColabFold `AlphaFold3_of3.ipynb` (OpenFold3) | free-tier GPU varies |
| Structure + binding affinity | ColabFold `Boltz1.ipynb` | same |
| AF2 monomer or complex | ColabFold `AlphaFold2_mmseqs2.ipynb` | ~2000 residues on 16 GiB |
| Complexes, no setup | `lab.chaidiscovery.com` | web UI |
| Non-commercial only | `alphafoldserver.com` | 30 jobs/day, 5000 tokens/job |
| Single-sequence triage | ColabFold `ESMFold.ipynb` | no MSA, lower accuracy |
| ESMC-6B embeddings | Biohub hosted client | ~12 GiB of weights; local is out |
| RFdiffusion at scale | Tamarind | runs locally for modest sizes |

Say which of the three Tamarind conditions applies before submitting — cloud
execution requested, insufficient local VRAM, or batch size. "This job does not
fit in 8 GiB" is a measurable claim; "no local GPU" is no longer a true one.

**Check the hand-off actually exists before planning around it.** Tier 3 is
load-bearing for everything above 8 GiB, and it is credential-gated: Tamarind
needs an API key (Customize > Credentials), and the remote-compute skills need
a connected provider — a Modal account or an SSH/Slurm target. With neither
configured, `list_compute` returns nothing and Tier 3 collapses to the free
notebooks and the structure connector. That is a workable set, but it is not
the same plan, so verify the credential before promising a routed job rather
than discovering it at submission time.

## Constraints that bite

**A CPU-only torch wheel.** The default PyPI `torch` on Windows has no CUDA.
Everything will work and be slow, which is the hardest failure to notice.
Install from the `download.pytorch.org/whl/cu###` index and assert
`torch.cuda.is_available()` before a long run.

**Windows versus WSL2.** Much of the structural-biology stack is Linux-first.
ColabFold-style installs and the RAPIDS family in particular are far easier
under WSL2, and some have no native Windows build at all. Confirm a Windows
build exists before promising one; WSL2 still reaches the same physical GPU.

**MMseqs2 server etiquette.** ColabFold's `--use_msa_server` is a shared public
resource. Serial queries from one IP are permitted; parallelising across
machines is explicitly discouraged. Forty designs is forty serial jobs — give
the realistic wall-clock before starting, and consider single-sequence folding
for first-pass triage with only survivors going to MSA-based folding. A local
GPU does not change this: the MSA server is the shared resource, not your card.

**AlphaFold Server token cap.** 5000 tokens, where each ligand or glycan *atom*
counts. A fully glycosylated trimer is not feasible; a single domain is.

**VRAM is shared with the desktop.** A display-attached 8 GiB card has several
hundred MiB already spoken for by the compositor. Plan against free memory, not
total, when a run is close to the limit.

## Ordering the work so folding is rare

A local GPU makes folding cheaper, not free. The ordering still holds, because
the point was never purely cost — it is that cheap filters are better filters
when applied first:

1. Filter on sequence-level criteria first — conservation, predicted
   presentation, self-similarity, obvious liabilities. Free.
2. Design sequences on CPU or local GPU. Cheap.
3. Rank on cheap proxies — MPNN score and confidence, ThermoMPNN ΔΔG. Cheap.
4. Fold only the top N. Decide N explicitly and say why.
5. **Fold a negative control alongside** — a design expected to fail.

Step 5 is the one most often skipped and the one that makes the round
interpretable. Without it, a column of pLDDT values cannot distinguish a good
design round from a saturated metric, and the Antigen-Structure contract's rule
about not interpreting below ~70 pLDDT gives you a floor but no ceiling.
