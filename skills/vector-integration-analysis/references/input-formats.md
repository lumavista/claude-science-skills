# Input formats

The scripts are schema-driven. Nothing assumes a particular provider's layout,
filenames, or that host contigs are called `chr...`.

## References

GenBank (`.gb`, `.gbk`) or FASTA (`.fa`, `.fasta`), **one construct per file**.

GenBank is strongly preferred: feature annotation is what turns "position
4,809" into "inside the transgene CDS", which is what makes a locus
interpretable in the report. With FASTA you still get the membership map and
all the arithmetic, but blocks are described by coordinates only.

Supply one reference per construct in the line — **including constructs the
provider did not map against, and including helper or donor plasmids.** Their
sequence is what explains the artefacts their absence created, and their
presence in the map is how the copy-number solve stays determined.

Construct names must be at least three characters. Short names cannot be
matched safely against contig names, and a construct called `A` would claim
contigs such as `AC233755.1` — a misclassification that silently deletes
genuine integration sites from the locus table. Use the real plasmid ID.

### Which reference form?

For a provirus or transposon, build the reference from the **integrated form**
if you can: the provirus between the LTR termini, or the cassette between the
ITRs. A transfer- or donor-plasmid reference puts the expected junction a
couple of bases inside the terminus and makes every real junction look
slightly off-end.

If you only have the plasmid, use it — but also include it as a *separate*
construct from the integrated form when checking backbone carryover, since that
check is precisely "is sequence outside the termini present at all?"

## Junction tables

Any delimited text file (tab, comma, semicolon or pipe; detected
automatically) with one junction per row. Both headed and headerless files
work.

### Fields

| Canonical name | Required | Meaning |
|---|---|---|
| `ref1`, `pos1` | yes | one side of the junction: reference name and coordinate |
| `ref2`, `pos2` | yes | the other side |
| `support` | yes | supporting reads, ideally after deduplication |
| `clip1`, `clip2` | strongly recommended | clip side of each end |
| `filler` | strongly recommended | bases between the two sides; negative = microhomology |
| `flag` | no | caller's own status or filter column |

Without `clip1`/`clip2` nothing can be paired, so every site looks
single-flanked and the site count loses its main consistency check. Without
`filler` the cross-mapping filler rules are weakened.

### Column resolution

1. **By header name.** Common spellings are recognised automatically —
   `chromosome`/`chrom`/`seqnames`/`contig` for `ref1`,
   `n_reads`/`read_count`/`supporting_reads`/`split_reads` for `support`, and
   so on.
2. **Positionally.** A headerless file is read in the canonical order
   `ref1, pos1, clip1, filler, ref2, pos2, clip2, support, flag`.
3. **Explicitly**, when either of the above is wrong:

```bash
--columns 'ref1=chromosome,pos1=3,clip1=side_a,support=n_reads'
```

Values are header names or 1-based column indices, mixed freely. The scripts
print and tabulate how they resolved every column — check that rather than
assuming.

## The clip side

**This is the field most likely to be misread, and misreading it is silent.**

Two canonical values:

- **`R`** — the reference continues **up to** this base, then novel sequence
  begins. On the host side this is the **left** flank of an insert.
- **`L`** — novel sequence, then the reference **resumes** at this base. On the
  host side this is the **right** flank.

Recognised automatically: `right_clip`/`left_clip`, `rightclip`/`leftclip`,
`right`/`left`, `rightmost`/`leftmost`, `R`/`L`, `3'`/`5'`,
`downstream`/`upstream`, `end`/`start`, `tail`/`head`, and the `soft_clip_*`
spellings.

Anything else, declare it:

```bash
--clip-map 'clipped_at_end=R,clipped_at_start=L'
```

**A strand is not a clip side.** Values such as `+`/`-`, `fwd`/`rev` or `1`/`-1`
say which way the read aligned, not which side of the breakpoint the reference
continues on, and cannot be converted. The scripts detect this, name it, and
degrade to single-flank mode rather than inventing pairings. If your caller
genuinely encodes the clip side in a column headed `strand`, map it explicitly
with `--columns clip1=strand` plus a `--clip-map`.

**Confirm the convention before relying on it.** Take one junction whose
orientation you know independently — an on-target HDR junction, or any locus
you have already confirmed by PCR — and check that the clip side the caller
assigns matches the definition above. Two signs you have it inverted:

- paired-flank footprints come out as large positive numbers where the modality
  predicts a small negative duplication;
- nothing pairs at all, because every junction at a site carries the same clip
  side.

## Which side is the host genome?

Decided by matching reference names against the constructs in the identity map,
**not** by name patterns. Anything that is not a construct and not an unmapped
token is host sequence, so RefSeq accessions, Ensembl names, scaffolds, alt
contigs and non-human assemblies all work unchanged.

Matching is exact on normalised names. Table names often differ from the
identity map — a construct called `pFOO` in GenBank may be `pFOO_linearized` in
the mapping FASTA — so declare it:

```bash
--construct-alias pFOO=pFOO_linearized      # CONSTRUCT=TABLENAME
```

A limited prefix match is attempted as a fallback for names of five characters
or more, and always warns. Pass `--no-prefix-match` to require exact names or
explicit aliases, which is the right setting for an automated pipeline.

If no reference name in a table matches any construct, the script **stops**
rather than silently classifying every junction as host-to-host. That failure
mode is silent and catastrophic, so it is made fatal.

## Unmapped partners

Callers denote "the other side matched nothing" variously: `FlankingDNA`,
`unknown`, `unmapped`, `NA`, `null`, `.`, `-`, or an empty field. All are
recognised. Add others with `--unmapped-token`.

**These rows matter** — they carry the cross-mapping signal used for the
fallback copy-number estimator — so do not filter them out before running the
scripts.

## Vector ends

`--vector-ends` takes the coordinates of the delivered molecule's defined ends
**in that construct's own numbering**:

- linearized plasmid: the restriction site coordinates (for a reference built
  by cutting a circular plasmid, the first and last few bases);
- provirus: the LTR termini;
- transposon: the ITR termini;
- circular or supercoiled delivery with no defined end: omit it.

Its meaning depends on `--modality` — see `modalities.md`. For fixed-end
modalities, supplying it is effectively mandatory: without it every real
junction is flagged as a reference-edge artefact.

## Support thresholds

`--min-support` defaults to **1**, deliberately. The lowest-support junctions
are frequently the most informative — a single-copy insertion in a multi-copy
background has the fewest reads — so filtering is left to you, explicitly, and
recorded in the report.

## Per-base depth

```bash
samtools depth -a -r <vector_contig> SAMPLE.bam > depth_<construct>.txt
```

Three columns: contig, position, depth. Point `copy_number.py` at one file per
vector reference per sample. This is what turns copy number from an inference
into a measurement; always ask for it.

## Coverage and read counts

`copy_number.py` takes these as numbers in its config, not as files, because
every provider formats them differently. Read them off whatever summary the
provider supplies:

- deduplicated read **pairs** per sample (if only single-end counts are
  available, halve them or set `read_length` accordingly);
- mean depth over each construct, per sample — needed only for the fallback
  routes.

## What is not supported

Aligned BAM/CRAM as a direct input. The scripts work from called junctions plus
per-base depth. If you have BAMs and no junction calls, run a structural-variant
or insertion-site caller first.
