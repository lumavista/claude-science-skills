# Compute Routes Without a Local GPU

The Antigen-Structure and Immuno-Designer host notes record no GPU, no remote
compute and no managed model endpoints. That does not block immunogen work,
but it changes which steps run locally and which are handed off — and being
explicit about which is which is the difference between a workflow that
completes and one that stalls silently partway through.

This collection already has one no-GPU answer: **Tamarind**, which runs
folding, design and docking as cloud jobs. Tamarind is the right default when
a job needs a GPU, because it is already wired in, schema-discoverable and
covered by an existing skill. What this file adds is the free-notebook route
for cases where Tamarind is not appropriate — no account or budget, a one-off
check, or a model Tamarind does not expose.

As SKILL.md and the specialist contracts both state: prefer the dedicated
structure-prediction skills over Tamarind where they can run, and say which of
the three conditions (cloud execution needed, no local GPU, batch size) applies
before submitting.

## Three tiers

**Tier 1 — runs locally on CPU, no caveat.** More of this workflow than people
expect:

- MHC class I and II prediction (MHCflurry, MHCnuggets, and the DTU tools if
  licensed — all small networks, CPU by design)
- IEDB retrieval, self-similarity screening, population coverage prep
- Alignment and conservation (`phylogenetics`, `conservation-mapping`)
- Mosaic and cocktail design, junction checking (`mosaic_design.py`)
- Sequence design and stability (ProteinMPNN, LigandMPNN, ThermoMPNN)
- Glycan shield grafting (GlycoSHIELD — avoiding MD is its whole contribution)
- Relative solvent accessibility (freesasa, DSSP)
- ESM-2 embeddings at the 8M–650M checkpoints

**Tier 2 — hand off to a hosted GPU.** Folding and refold validation. Tamarind
first; free notebooks when Tamarind does not fit:

| Need | Route | Limits |
|---|---|---|
| Existing structure | Structures & Interactions connector (PDB / AlphaFold DB) | none — **always check this first** |
| AF3-class, commercially usable | ColabFold `AlphaFold3_of3.ipynb` (OpenFold3) | free-tier GPU availability varies |
| Structure + binding affinity | ColabFold `Boltz1.ipynb` | same |
| AF2 monomer or complex | ColabFold `AlphaFold2_mmseqs2.ipynb` | ~2000 residues on 16GB |
| Single-sequence triage | ColabFold `ESMFold.ipynb` | no MSA, lower accuracy |
| Complexes, no setup | `lab.chaidiscovery.com` | web UI |
| Non-commercial only | `alphafoldserver.com` | 30 jobs/day, 5000 tokens/job |

The first row is the one that gets skipped. For a known antigen an
experimental structure or an AlphaFold DB model usually exists, needs no GPU,
and beats anything predicted fresh.

**Tier 3 — genuinely needs a GPU.** RFdiffusion and variants, Genie2, Chroma.
Tamarind runs RFdiffusion; the others have no hosted free route. Emit a job
script and state the requirement rather than implying they can run locally.

## Constraints that bite

**MMseqs2 server etiquette.** ColabFold's `--use_msa_server` is a shared public
resource. Serial queries from one IP are permitted; parallelising across
machines is explicitly discouraged. Forty designs is forty serial jobs — give
the realistic wall-clock before starting, and consider single-sequence ESMFold
for first-pass triage with only survivors going to MSA-based folding.

**Colab free tier.** GPU assignment is luck-of-the-draw and sessions time out.
Over roughly 2000 residues on a 16GB GPU will fail — a trimeric envelope
glycoprotein exceeds that, so fold the protomer and state that quaternary and
trimer-interface epitopes are not covered by that model.

**AlphaFold Server token cap.** 5000 tokens, where each ligand or glycan *atom*
counts. A fully glycosylated trimer is not feasible; a single glycosylated
domain is.

## Ordering the work so folding is rare

Folding is the bottleneck, so order every pipeline to fold as little as
possible:

1. Filter on sequence-level criteria first — conservation, predicted
   presentation, self-similarity, obvious liabilities. Free.
2. Design sequences on CPU. Cheap.
3. Rank on cheap proxies — MPNN score and confidence, ThermoMPNN ΔΔG. Cheap.
4. Fold only the top N. Decide N explicitly and say why.
5. **Fold a negative control alongside** — a design expected to fail.

Step 5 is the one most often skipped and the one that makes the round
interpretable. Without it, a column of pLDDT values cannot distinguish a good
design round from a saturated metric, and the Antigen-Structure contract's rule
about not interpreting below ~70 pLDDT gives you a floor but no ceiling.
