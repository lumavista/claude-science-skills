#!/usr/bin/env bash
# The host stage: deciding the host organism and assembly, fetching its
# reference, and annotating called sites -- offline, against a synthetic
# catalogue cache, so no network and no real genome data are needed.
#
# Each case is a way the host goes wrong in practice: a human default applied
# to a hamster line, a build guessed from `chr1`, junction coordinates past
# the end of the contig they claim, a genome size off by a fifth.
#
#   ./test_host.sh            run everything
#   KEEP=1 ./test_host.sh     keep the working directory for inspection
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="$HERE/../scripts"
WORK="$(mktemp -d)"
PY="${PYTHON:-python3}"
PASS=0; FAIL=0

cleanup() { [ -n "${KEEP:-}" ] && echo "kept: $WORK" || rm -rf "$WORK"; }
trap cleanup EXIT

ok()   { PASS=$((PASS+1)); printf '  \033[32mPASS\033[0m %s\n' "$1"; }
bad()  { FAIL=$((FAIL+1)); printf '  \033[31mFAIL\033[0m %s\n' "$1"; }
assert_has() {
  if grep -qiE -- "$3" "$2" 2>/dev/null; then ok "$1"; else
    bad "$1 (expected to match: $3)"; fi
}
assert_lacks() {
  if grep -qiE -- "$3" "$2" 2>/dev/null; then bad "$1 (should not match: $3)";
  else ok "$1"; fi
}
# expect_exit <label> <code> <log> <cmd...>
expect_exit() {
  local lbl="$1" want="$2" log="$3"; shift 3
  "$@" >"$log" 2>&1; local got=$?
  if [ "$got" -eq "$want" ]; then ok "$lbl"; else
    bad "$lbl (exit $got, expected $want)"; tail -5 "$log"; fi
}
HR() { "$PY" "$SCRIPTS/host_reference.py" --offline --cache "$WORK/cache" "$@"; }

echo "== fixtures =="
"$PY" -I "$HERE/make_fixtures.py" "$WORK" >"$WORK/fixtures.log" 2>&1 \
  && ok "pipeline fixtures generated" || { bad "fixtures"; cat "$WORK/fixtures.log"; exit 1; }
cd "$WORK" || exit 1

# A synthetic catalogue cache: one assembly report per catalogue entry, under
# the file name the script would download it to, plus a GTF for the human
# entry. Names follow NCBI's conventions; lengths and genes are invented.
"$PY" - "$SCRIPTS" <<'EOF' >cache.log 2>&1 && ok "synthetic catalogue cache" || { bad "cache"; cat cache.log; }
import gzip, os, sys
sys.path.insert(0, sys.argv[1])
from host_reference import CATALOGUE
os.makedirs("cache", exist_ok=True)
# chromosomes each fake assembly carries (hamster has no chr1, as in RefSeq)
CHROMS = {"human": [str(i) for i in range(1, 23)] + ["X", "Y"],
          "mouse": [str(i) for i in range(1, 20)] + ["X", "Y"],
          "chinese_hamster": [str(i) for i in range(2, 11)] + ["X"],
          "rat": [str(i) for i in range(1, 21)] + ["X", "Y"],
          "african_green_monkey": [str(i) for i in range(1, 30)] + ["X", "Y"],
          "dog": [str(i) for i in range(1, 39)] + ["X"]}
for k, e in enumerate(CATALOGUE):
    chroms = CHROMS[e["key"]]
    rows = []
    for i, c in enumerate(chroms):
        rows.append(["chr" + c if e["key"] == "chinese_hamster" else c,
                     "assembled-molecule", c, "Chromosome",
                     "CM9%02d%03d.1" % (k, i), "=",
                     "NC_9%02d%03d.1" % (k, i), "Primary Assembly",
                     str(100_000_000), "chr" + c if e["key"] == "human" else "na"])
    if e["key"] == "human":   # primary assembly sums to exactly 3.1 Gb
        rows[0][8] = str(3_100_000_000 - 100_000_000 * (len(rows) - 1))
    rows.append(["MT", "assembled-molecule", "MT", "Mitochondrion",
                 "J0%d.1" % k, "=", "NC_9%02d999.1" % k, "non-nuclear",
                 "16569", "chrM" if e["key"] == "human" else "na"])
    rows.append(["HSCHR6_ALT", "alt-scaffold", "6", "Chromosome",
                 "KI9%02d.1" % k, "=", "NT_9%02d.1" % k, "ALT_REF_LOCI_1",
                 "4000000", "na"])
    name = "%s_%s" % (e["accession"], e["assembly"])
    with open("cache/%s_assembly_report.txt" % name, "w") as fh:
        fh.write("# Assembly name:  %s\n# Organism name:  %s (%s)\n"
                 "#\n## Assembly-Unit names:\n## Primary Assembly\n#\n"
                 "# Sequence-Name\tSequence-Role\tAssigned-Molecule\t"
                 "Assigned-Molecule-Location/Type\tGenBank-Accn\tRelationship\t"
                 "RefSeq-Accn\tAssembly-Unit\tSequence-Length\tUCSC-style-name\n"
                 % (e["assembly"], e["species"], e["common"]))
        for r in rows:
            fh.write("\t".join(r) + "\n")
    if e["key"] == "human":
        nc = lambda c: "NC_9%02d%03d.1" % (k, chroms.index(c))
        g = [  # seqid in RefSeq naming: the script must rename to chr*
            (nc("2"), "gene", 4_990_000, 5_010_000, "+", 'gene_id "EXA"; gene "EXA"; gene_biotype "protein_coding";'),
            (nc("2"), "exon", 4_999_000, 5_001_000, "+", 'gene_id "EXA"; gene "EXA";'),
            (nc("6"), "gene", 900_000, 1_100_000, "-", 'gene_id "INB"; gene "INB"; gene_biotype "protein_coding";'),
            (nc("6"), "exon", 900_000, 901_000, "-", 'gene_id "INB"; gene "INB";'),
            (nc("6"), "exon", 1_099_000, 1_100_000, "-", 'gene_id "INB"; gene "INB";'),
            (nc("12"), "gene", 1_050_000, 1_060_000, "-", 'gene_id "NEC"; gene "NEC"; gene_biotype "lncRNA";'),
            (nc("12"), "exon", 1_050_000, 1_060_000, "-", 'gene_id "NEC"; gene "NEC";'),
        ]
        with gzip.open("cache/%s_genomic.gtf.gz" % name, "wt") as fh:
            fh.write("#gtf-version 2.2\n")
            for s, t, a, b, st, attr in g:
                fh.write("\t".join([s, "RefSeq", t, str(a), str(b), ".", st, ".", attr]) + "\n")
print("ok")
EOF

echo
echo "== deciding the host from contig names alone =="
expect_exit "a RefSeq accession names one assembly: inferred" 0 inf1.log \
  HR --contig NC_900001.1:5000 --no-annotation --out-dir o1
assert_has "inferred as human, said so" inf1.log "Homo sapiens .*decided by inferred"
expect_exit "chr-style names fit several assemblies: stop, never guess" 3 inf2.log \
  HR --contig chr6:1000000 --contig chr12:1000000 --no-annotation --out-dir o2 --strict
assert_has "the stop names the candidates" inf2.log "host_species.*human .*mouse"
assert_has "the stop says to ask the user" inf2.log "ask the user"
expect_exit "--accept host_species is refused: still a stop" 3 inf2a.log \
  HR --contig chr6:1000000 --contig chr12:1000000 --no-annotation --out-dir o2a \
  --strict --accept host_species
assert_has "the refusal is said out loud" inf2a.log "accept host_species refused"
assert_lacks "never marked accepted" inf2a.log "CHECK ACCEPTED \[host_species\]"
expect_exit "a cell-line hint settles it" 0 inf3.log \
  HR --contig chr6:1000000 --contig chr12:1000000 --cell-line HEK-293T \
  --no-annotation --out-dir o3 --strict
assert_has "hint recorded as how it was decided" inf3.log "decided by hint"
expect_exit "a contig absent from the hinted host stops the run" 3 inf4.log \
  HR --contig chr1:1000 --cell-line CHO-K1 --no-annotation --out-dir o4 --strict
assert_has "nothing fits: the build or organism is questioned" inf4.log \
  "no catalogue assembly fits"
expect_exit "positions past every contig end fit nothing" 3 inf5.log \
  HR --contig chr2:150000000 --no-annotation --out-dir o5 --strict
assert_has "the overrun is named" inf5.log "positions past contig end"
expect_exit "an unknown cell line is refused, not guessed" 3 inf6.log \
  HR --contig chr2:1000 --cell-line Sf9 --no-annotation --out-dir o6 --strict
expect_exit "--offline with an empty cache is an input error" 2 inf7.log \
  "$PY" "$SCRIPTS/host_reference.py" --offline --cache empty --contig chr2:1 \
  --no-annotation --out-dir o7
assert_has "it names what to fetch" inf7.log "ftp.ncbi.nlm.nih.gov"
expect_exit "an unknown assembly name is refused with the fix" 2 inf8.log \
  HR --assembly hg19 --contig chr2:1 --no-annotation --out-dir o8
assert_has "it asks for the accession" inf8.log "accession"

echo
echo "== a declared assembly is still checked =="
expect_exit "declared build that the junctions overrun" 3 dec1.log \
  HR --assembly GRCh38 --contig chr2:150000000 --no-annotation --out-dir d1 --strict
assert_has "host_contigs fails, naming the overrun" dec1.log "host_contigs.*chr2:150000000"
expect_exit "decoy contigs warn, they do not stop" 0 dec2.log \
  HR --assembly hg38 --contig chr2:5 --contig chrUn_JTFH01000001v1_decoy:10 \
  --no-annotation --out-dir d2 --strict
assert_has "decoy named in the warning" dec2.log "WARNING \[host_contigs\].*decoy"

echo
echo "== the full run: no genome_size declared, annotation on =="
mkdir -p mf
cp pDelta.gb pEpsilon.gb junctions_recon_*.tsv depth_concat_*.txt mf/
write_manifest() { # write_manifest <file> <extra json fields>
cat > "$1" <<JSON
{$2 "read_length":150,"ploidy":2,"modality":"nhej",
 "primary_mapping":"pEpsilon",
 "vector_ends":{"pDelta":[1200],"pEpsilon":[3000]},
 "samples":{
  "s1":{"junctions":{"pDelta":"junctions_recon_pDelta_s1.tsv",
                     "pEpsilon":"junctions_recon_pEpsilon_s1.tsv"},
        "dedup_pairs":323640000,"baseline_depth":29.0,
        "depth_files":{"pDelta":"depth_concat_pDelta.txt",
                       "pEpsilon":"depth_concat_pEpsilon.txt"}},
  "s2":{"junctions":{"pDelta":"junctions_recon_pDelta_s2.tsv",
                     "pEpsilon":"junctions_recon_pEpsilon_s2.tsv"}}}}
JSON
}
cp -r mf nohost
write_manifest mf/run.json '"host":{"cell_line":"HEK293"},'
write_manifest nohost/run.json '"genome_size":3100000000,'
for d in mf nohost; do
  "$PY" "$SCRIPTS/construct_identity.py" $d/pDelta.gb $d/pEpsilon.gb --run $d/run.json >$d/s1.log 2>&1
  for m in pDelta pEpsilon; do
    "$PY" "$SCRIPTS/classify_junctions.py" --run $d/run.json --mapping $m --strict >$d/s2_$m.log 2>&1
  done
done
expect_exit "Stage 3 without the host stage and without genome_size stops" 2 mf/s3_early.log \
  "$PY" "$SCRIPTS/copy_number.py" --run mf/run.json --strict
assert_has "and says to run the host stage" mf/s3_early.log "host_reference.py"
expect_exit "host stage on the manifest" 0 mf/host.log \
  HR --run mf/run.json --strict
got=$("$PY" -c "import json;d=json.load(open('mf/run.json'));h=d['stages']['host'];print(h['assembly'],h['genome_size'],h['naming'],h['decided_by'],'genome_size' in d)")
[ "$got" = "GRCh38.p14 3100000000 ucsc hint False" ] \
  && ok "host recorded under stages; declared facts untouched" \
  || bad "host record was '$got'"
assert_has "host .fai named as the junctions name contigs" mf/host_reference/host.fai "^chr12	100000000"
assert_lacks "alt loci left out of the .fai" mf/host_reference/host.fai "ALT|^NT_"
for d in mf nohost; do
  "$PY" "$SCRIPTS/copy_number.py" --run $d/run.json --strict >$d/s3.log 2>&1 \
    && ok "Stage 3 ($d)" || { bad "Stage 3 ($d)"; tail -5 $d/s3.log; }
  "$PY" "$SCRIPTS/reconcile_sites.py" --run $d/run.json --strict \
    --min-copies-per-site 10 >$d/s4.log 2>&1 \
    && ok "Stage 4 ($d)" || { bad "Stage 4 ($d)"; tail -5 $d/s4.log; }
done
"$PY" - <<'EOF' >same.log 2>&1
import json
a = json.load(open("mf/accounting.json")); b = json.load(open("nohost/accounting.json"))
for k in ("copies_per_haploid", "copies_per_cell", "solve"):
    assert a["samples"]["s1"][k] == b["samples"]["s1"][k], k
print("OK")
EOF
grep -q "^OK$" same.log && ok "assembly genome_size gives the declared one's copy numbers" \
  || { bad "genome_size from the host changed a number"; cat same.log; }
n=$(head -1 nohost/called_sites.tsv | awk -F'\t' '{print NF}')
if cmp -s <(cut -f1-"$n" mf/called_sites.tsv) nohost/called_sites.tsv; then
  ok "annotation only appends columns: the rest is byte-identical"
else bad "annotation changed the called table's existing columns"; fi
assert_lacks "no host stage, no annotation columns" nohost/called_sites.tsv "host_region"
assert_has "header ends with the annotation" mf/called_sites.tsv \
  "note	host_region	host_gene	host_gene_biotype	host_strand	nearest_gene	nearest_gene_distance	nearest_gene_side$"
assert_has "exon hit is exonic" mf/called_sites.tsv "^chr2	5000000	.*exonic	EXA	protein_coding	\+"
assert_has "between exons is intronic" mf/called_sites.tsv "^chr6	1000000	.*intronic	INB	protein_coding	-"
assert_has "intergenic site names nearest gene, distance and side" \
  mf/called_sites.tsv "^chr12	1000000	.*intergenic	+NEC	50000	downstream$"
assert_has "contig without genes is marked as such" mf/called_sites.tsv "^chr5	.*no-annotation"
assert_has "region counts in the site summary" mf/site_summary.json '"host_annotation"'

echo
echo "== figures pick the host .fai up from the manifest =="
"$PY" "$SCRIPTS/make_figures.py" --run mf/run.json >mf/figs.log 2>&1 \
  && ok "make_figures --run" || { bad "make_figures --run"; cat mf/figs.log; }
assert_has "genome overview drawn from the host stage" mf/figs.log "wrote .*genome_overview.svg"
assert_has "coverage drawn from the manifest's depth files" mf/figs.log "wrote .*coverage_s1_pEpsilon.svg"

echo
echo "== a declared genome_size that contradicts the assembly =="
cp -r nohost wrong
"$PY" - <<'EOF'
import json
d = json.load(open("wrong/run.json")); d["genome_size"] = 2400000000
d["host"] = {"assembly": "GRCh38"}; json.dump(d, open("wrong/run.json", "w"))
EOF
expect_exit "a hamster-sized genome on a human assembly stops the run" 3 wrong/host.log \
  HR --run wrong/run.json --no-annotation --strict
assert_has "host_genome_size names both values" wrong/host.log \
  "host_genome_size.*2400000000.*3100000000"

echo
echo "=============================="
echo " passed: $PASS   failed: $FAIL"
echo "=============================="
[ "$FAIL" -eq 0 ]
