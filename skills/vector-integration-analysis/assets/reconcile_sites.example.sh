#!/usr/bin/env bash
# The final stage: one called-site table from the per-mapping loci tables.
#
# Replace every name and number. Nothing here is specific to a construct, a
# copy number, a coordinate or a sample -- the constructs come from the
# identity map, the ends from the delivered molecule's linearization or
# terminus coordinates, and the copy numbers from copy_number.py.
#
# Run it AFTER classify_junctions.py (once per mapping) and copy_number.py.
# Report what it writes; do not re-derive loci, thresholds or attribution.
#
# This spells out every flag, which is useful for reading. In a real run the
# manifest supplies them:
#
#     python scripts/reconcile_sites.py --run run.json --strict
#
# and only the decisions below -- the retention thresholds, the ceiling's
# divisor and the margin -- are worth passing explicitly.
set -euo pipefail

python scripts/reconcile_sites.py \
    --identity identity.json \
    \
    `# one loci table per separately mapped construct` \
    --loci pAlpha=loci_pAlpha.tsv \
    --loci pBeta=loci_pBeta.tsv \
    \
    `# per-construct copies per haploid, and the ploidy that scales the ceiling` \
    --copy-number accounting.json \
    --ploidy 2 \
    \
    `# the delivered molecule's declared ends, in each construct's numbering.` \
    `# Constructs linearized at different sites are told apart by this even` \
    `# inside shared sequence.` \
    --vector-ends pAlpha=1,8412 \
    --vector-ends pBeta=1,7650 \
    \
    `# retention: support, and independently prepared samples` \
    --min-reads 5 \
    --min-samples 2 \
    \
    `# the ceiling's divisor: the fewest copies a real site is assumed to` \
    `# carry. 1 is the only valid value for integrase/transposase delivery;` \
    `# on a concatemer line pass the copies per site copy_number.py measured.` \
    --min-copies-per-site 1 \
    \
    `# how far the top candidate must lead on copy number before a shared` \
    `# locus is attributed to it at all. Raise it to be more conservative.` \
    --cn-margin 2.0 \
    \
    --out-called called_sites.tsv \
    --out-subthreshold subthreshold_sites.tsv \
    --summary site_summary.json \
    --md site_summary.md
