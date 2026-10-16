# Linkers, Tags and Display Domains

Reference for the parts of a construct that are not the antigen. Choices here
decide whether a design expresses, folds, trimerises and purifies — and they
are made badly more often than the antigen design itself, because they feel
like boilerplate.

## Linkers

**Flexible (GS-family).** `(G4S)n` is the default when two modules must fold
independently and relative geometry does not matter. `(GGGGS)1` ≈ 18 Å
extended, and length scales roughly linearly. Glycine-rich linkers are
protease-susceptible and can be immunogenic at junctions; keep them as short
as the geometry allows.

**Rigid.** `(EAAAK)n` forms a helix and holds modules apart at a defined
distance; `(AP)n` and proline-rich `(XP)n` are stiffer still. Use when you need
to prevent two domains interacting, or to project an epitope away from a
carrier so it is accessible.

**Cleavable.** Furin (`RXKR`), TEV (`ENLYFQ/G`), thrombin, 3C. Needed when a
tag or a trimerisation domain must come off before use. Check the site does not
occur internally in the antigen.

**Choosing length.** Too short strains the fold and can prevent either module
folding; too long invites proteolysis, aggregation and conformational
heterogeneity. If a structure of both modules exists, measure the distance
between the termini to be joined and size the linker to it rather than
defaulting to `(G4S)3`.

**The junction is not neutral.** Every join creates peptides spanning it that
exist in neither parent. For any construct containing an antigen you care
about immunologically, scan the junctions —
`immunoinformatics/scripts/mosaic_design.py evaluate` at k=9 and k=15 does
this against the strain panel.

## Trimerisation and oligomerisation domains

Relevant for class I viral fusion proteins, which are trimeric and whose
neutralising epitopes are frequently quaternary — a monomeric construct simply
does not present them.

| Domain | Origin | Notes |
|---|---|---|
| **Foldon** | T4 fibritin, ~27 aa | The default. Small, folds autonomously, very well characterised. Immunogenic in its own right — responses against foldon are routinely observed, so include a foldon-only control or plan to cleave it. |
| **GCN4 variants** (IZ, IQ) | engineered coiled coil | Defined register, can set trimer geometry more precisely than foldon |
| **Ferritin** | 24-mer particle | Display rather than trimerisation; presents 8 trimeric spikes |
| **SpyTag/SpyCatcher** | engineered isopeptide | Modular particle display; lets you make antigen and scaffold separately |

**Placement.** Almost always C-terminal, after the ectodomain, joined by a
short flexible linker. Check the added domain does not occlude the target
epitope — for a conformational antigen, recompute accessibility on the
assembled construct, not on the bare antigen.

**Validate after assembly, not before.** An appended domain can perturb the
fold of the module it is attached to. Validating the bare antigen and then
shipping the fusion is a common and avoidable error.

## Purification and detection tags

| Tag | Notes |
|---|---|
| His6/His8 | Standard. Terminal, with a short linker. Can bind metals in ways that affect some assays |
| Strep-II | Gentler elution, useful when His is unsuitable |
| Avi | Site-specific biotinylation — the right choice for making probes for B-cell sorting |
| Twin-Strep, FLAG | Alternatives when the above conflict with downstream use |

**Place tags away from the epitope**, and decide explicitly whether they are
cleavable. A His-tag adjacent to a target patch changes what the immune system
sees, and a non-cleavable tag on an immunogen means the tag is part of the
immunogen.

For B-cell probes specifically, an Avi tag plus controlled biotinylation beats
random chemical biotinylation, which modifies surface lysines — often including
the ones in the epitope you are sorting on.

## Signal peptides and expression

Match the signal peptide to the expression system rather than carrying over the
native one. Common choices are tPA and IgG-kappa leaders for mammalian
secretion. For mRNA constructs the signal peptide, codon optimisation, GC
content, 5' secondary structure and uridine content all matter and are
decided together, not separately.

## Order of assembly

A workable default for a secreted trimeric antigen:

```
[signal peptide] - [antigen ectodomain] - [short flexible linker]
  - [trimerisation domain] - [cleavage site] - [purification tag]
```

State the rationale for each element in the construct table. A reader six
weeks later cannot reconstruct why a particular linker length was chosen, and
that is exactly the parameter they will want to vary in the next round.
