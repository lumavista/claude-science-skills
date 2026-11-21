#!/usr/bin/env bash
# Exercise every code path in the skill against synthetic fixtures.
#
# Each case below corresponds to a way this analysis has gone wrong in
# practice. A regression here is a wrong biological answer, not a cosmetic bug,
# so the suite asserts on content rather than just on exit status.
#
#   ./run_tests.sh            run everything
#   KEEP=1 ./run_tests.sh     keep the working directory for inspection
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

# assert_has <label> <file-or-"-"> <pattern>
assert_has() {
  if grep -qiE -- "$3" "$2" 2>/dev/null; then ok "$1"; else
    bad "$1 (expected to match: $3)"; fi
}
assert_lacks() {
  if grep -qiE -- "$3" "$2" 2>/dev/null; then bad "$1 (should not match: $3)";
  else ok "$1"; fi
}
# assert_close <label> <actual> <expected> <tolerance>
assert_close() {
  "$PY" - "$2" "$3" "$4" <<'EOF' >/dev/null 2>&1
import sys
a, e, t = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
sys.exit(0 if abs(a - e) <= t else 1)
EOF
  if [ $? -eq 0 ]; then ok "$1 ($2 ~ $3)"; else bad "$1 (got $2, expected $3+-$4)"; fi
}

echo "== fixtures =="
"$PY" -I "$HERE/make_fixtures.py" "$WORK" >"$WORK/fixtures.log" 2>&1 \
  && ok "fixtures generated" || { bad "fixture generation"; cat "$WORK/fixtures.log"; exit 1; }

cd "$WORK" || exit 1

echo
echo "== Stage 1: identity map =="
for n in 1 2 3; do
  case $n in
    1) refs="pAlpha.gb" ;;
    2) refs="pAlpha.gb pBeta.gb" ;;
    3) refs="pAlpha.gb pBeta.gb pGamma.gb" ;;
  esac
  if "$PY" "$SCRIPTS/construct_identity.py" $refs --out "id$n.json" --md "map$n.md" \
       >"id$n.log" 2>&1; then ok "$n construct(s)"; else bad "$n construct(s)"; cat "id$n.log"; fi
done
assert_has "membership sets, not just 'shared'" map3.md "pAlpha \+ pBeta \+ pGamma"
assert_has "boundary names the clipping construct" map3.md "reads of pGamma clip here"
assert_has "one-construct map has no boundaries" id1.log "boundaries none"

echo
echo "== Stage 1: short construct names are rejected =="
if "$PY" "$SCRIPTS/construct_identity.py" pAlpha.gb pBeta.gb --name A --name B \
     --out bad.json >short.log 2>&1; then
  bad "short name should be refused"
else
  assert_has "short name refused with a reason" short.log "too short to match safely"
fi

echo
echo "== Stage 2: clip vocabularies =="
run2() { # run2 <label> <table> <modality> <ends> <extra...>
  local lbl="$1" tbl="$2" mod="$3" ends="$4"; shift 4
  "$PY" "$SCRIPTS/classify_junctions.py" --identity id3.json --mapping pAlpha \
    --modality "$mod" --vector-ends "$ends" --table parental="$tbl" \
    --out-loci "loci_$lbl.tsv" --out-junctions "junc_$lbl.tsv" \
    --md "md_$lbl.md" "$@" >"s2_$lbl.log" 2>&1
}

run2 canon junctions_nhej.tsv nhej 1,5300 && ok "canonical right_clip/left_clip" \
  || { bad "canonical clip vocabulary"; cat s2_canon.log; }
assert_has "canonical: a flank pair was formed" s2_canon.log "1 with both flanks"

run2 alt junctions_nhej_altvocab.tsv nhej 1,5300 \
  && ok "rightmost/leftmost auto-recognised" || bad "rightmost/leftmost"
assert_has "alt vocabulary still pairs flanks" s2_alt.log "1 with both flanks"

run2 strand junctions_nhej_strand.tsv nhej 1,5300 \
  && ok "strand values do not crash" || bad "strand values crashed"
assert_has "strand values diagnosed as a strand" s2_strand.log "look like a STRAND"
assert_has "strand: degrades to single-flank" s2_strand.log "clip side unusable"

run2 noclip junctions_nhej_noclip.tsv nhej 1,5300 \
  && ok "absent clip column does not crash" || bad "absent clip column crashed"
assert_has "absent clip reported, not hidden" s2_noclip.log "NOTHING could be paired"

run2 exotic junctions_nhej_exotic.tsv nhej 1,5300 \
  && ok "unrecognised vocabulary does not crash" || bad "unrecognised vocabulary crashed"
assert_has "unrecognised vocabulary names --clip-map" s2_exotic.log "clip.map"

run2 exotic2 junctions_nhej_exotic.tsv nhej 1,5300 \
  --clip-map 'clipped_at_end=R,clipped_at_start=L' \
  && ok "--clip-map rescues an unknown vocabulary" || bad "--clip-map"
assert_has "clip-map restores flank pairing" s2_exotic2.log "1 with both flanks"

echo
echo "== Stage 2: host contigs are not swallowed by construct names =="
assert_has "alt-contig accession kept as host" loci_canon.tsv "^AC233755\.1"

echo
echo "== Stage 2: artefact flagging =="
assert_has "boundary cross-map flagged, with construct" junc_canon.tsv "boundary-crossmap"
assert_has "vector-end religation flagged" junc_canon.tsv "vector-end-religation"
assert_has "short-range inverted flagged" junc_canon.tsv "short-range-inverted"
"$PY" - <<'EOF' >/dev/null 2>&1
import csv, sys
rows = list(csv.DictReader(open("junc_canon.tsv"), delimiter="\t"))
# the self-junction at 4000 (same clip side, zero span) is a real palindromic
# concatemer join and must NOT be flagged
hit = [r for r in rows if r["vector_pos"] == "4000"]
sys.exit(0 if hit and not hit[0]["artefact"] else 1)
EOF
[ $? -eq 0 ] && ok "self-junction left unflagged (it is real)" \
             || bad "self-junction was flagged as an artefact"

echo
echo "== Stage 2: integrase modality =="
run2 lenti junctions_integrase.tsv integrase 1,5300 \
  && ok "integrase run" || { bad "integrase run"; cat s2_lenti.log; }
assert_has "fingerprint detects fixed vector ends" s2_lenti.log "FIXED vector ends"
assert_has "LTR junctions kept, not called edge artefacts" s2_lenti.log "vector-end *7"
assert_has "5 bp target-site duplication confirmed" s2_lenti.log "expected footprint"
assert_lacks "real LTR junctions not flagged reference-edge" loci_lenti.tsv "reference-edge"
assert_has "off-end junction flagged" junc_lenti.tsv "off-end-host-junction"

echo
echo "== Stage 2: fingerprint contradicts a wrong declaration =="
run2 wrongmod junctions_integrase.tsv nhej 1,5300
assert_has "wrong modality warned about" s2_wrongmod.log "do not look like the declared modality"

echo
echo "== Stage 2: multiple samples =="
"$PY" "$SCRIPTS/classify_junctions.py" --identity id3.json --mapping pAlpha \
  --modality nhej --vector-ends 1,5300 \
  --table parental=junctions_nhej.tsv --table treated=junctions_nhej_s2.tsv \
  --norm treated=0.9 --out-loci loci_two.tsv --out-junctions junc_two.tsv \
  >s2_two.log 2>&1 && ok "two samples" || { bad "two samples"; cat s2_two.log; }
assert_has "cross-sample ratio column present" loci_two.tsv "ratio_treated"

echo
echo "== Stage 3: measured copy number (3 constructs, from depth) =="
cat > cfg_measured.json <<'JSON'
{"identity":"id3.json","genome_size":3100000000,"read_length":150,"ploidy":2,
 "primary_mapping":"pAlpha","modality":"nhej",
 "samples":{"parental":{"dedup_pairs":299666666,
   "depth_files":{"pAlpha":"depth_pAlpha.txt","pBeta":"depth_pBeta.txt",
                  "pGamma":"depth_pGamma.txt"},
   "detected_loci":4}},
 "mappings":{"pAlpha":{"junctions":"junc_canon.tsv","mean_coverage":{"parental":464.0}},
             "pBeta":{"mean_coverage":{"parental":464.0}},
             "pGamma":{"mean_coverage":{"parental":464.0}}}}
JSON
"$PY" "$SCRIPTS/copy_number.py" cfg_measured.json --md q_measured.md \
  --json q_measured.json >q_measured.log 2>&1 \
  && ok "measured solve ran" || { bad "measured solve"; cat q_measured.log; }
assert_has "route reported as measured" q_measured.md "weighted least squares"
assert_lacks "no caveat banner when measured" q_measured.md "NOT DECOMPOSED"
for c in pAlpha:10 pBeta:4 pGamma:2; do
  name="${c%%:*}"; want="${c##*:}"
  got=$("$PY" -c "import json;print(json.load(open('q_measured.json'))['samples']['parental']['copies_per_haploid']['$name'])")
  assert_close "solved copies/haploid for $name" "$got" "$want" 0.3
done

assert_has "old config states the read-count baseline" q_measured.log \
  "baseline: read-count estimate .*genome-median per-base depth .* preferred"

echo
echo "== Stage 3: concatemer backbone deficit (no silent clamped zero) =="
# pDelta 1.5 and pEpsilon 28 copies share 4.6 kb of backbone carrying 15% less
# depth than A+B. The old all-block solve with clamping reported pDelta as 0.
# The read-count genome coverage is ~8% above the genome median (29.0x).
cat > cfg_concat.json <<'JSON'
{"identity":"identity_concat.json","genome_size":3100000000,"read_length":150,
 "ploidy":2,"primary_mapping":"pDelta","modality":"nhej",
 "samples":{"s1":{"dedup_pairs":323640000,"baseline_depth":29.0,
   "depth_files":{"pDelta":"depth_concat_pDelta.txt",
                  "pEpsilon":"depth_concat_pEpsilon.txt"}}},
 "mappings":{"pDelta":{},"pEpsilon":{}}}
JSON
"$PY" "$SCRIPTS/copy_number.py" cfg_concat.json --md q_concat.md \
  --json q_concat.json >q_concat.log 2>&1 \
  && ok "concatemer-deficit case ran" || { bad "concatemer-deficit case"; cat q_concat.log; }
jq_concat() { "$PY" -c "import json;d=json.load(open('q_concat.json'))['samples']['s1'];print($1)"; }
got=$(jq_concat "d['copies_per_haploid']['pDelta']")
assert_close "low-copy construct from its unique block, not 0" "$got" 1.5 0.15
got=$(jq_concat "d['copies_per_haploid']['pEpsilon']")
assert_close "high-copy construct close to truth" "$got" 28 1.0
got=$(jq_concat "d['solve']['closure']['deficit_pct']")
assert_close "shared-block deficit reported near -15%" "$got" -15 2
"$PY" - <<'EOF' >/dev/null 2>&1
import json, sys
d = json.load(open("q_concat.json"))["samples"]["s1"]
cps = d["copies_per_haploid"]
# no reported copy number may be a clamped zero; the diagnostic stays negative
ok = all(v is None or v > 0.05 for v in cps.values()) \
     and min(d["solve"]["wls_depths"].values()) < 0
sys.exit(0 if ok else 1)
EOF
[ $? -eq 0 ] && ok "no silent clamped zero; negative solve kept as diagnostic" \
             || bad "a copy number was clamped to zero or the solve was clamped"
assert_has "negative solve raised as a caveat, not clamped" q_concat.md \
  "ADDITIVE MODEL INCOMPLETE.*backbone"
assert_has "caveat printed before the numbers" <(head -8 q_concat.md) \
  "ADDITIVE MODEL INCOMPLETE"
assert_has "genome-median baseline used and stated" q_concat.log \
  "baseline: genome-median per-base depth 29.00x"
assert_has "closure table in the report" q_concat.md "Shared-block closure: -1[0-9]\."

"$PY" -c "
import json; c=json.load(open('cfg_concat.json'))
del c['samples']['s1']['baseline_depth']; json.dump(c, open('cfg_concat_rc.json','w'))"
"$PY" "$SCRIPTS/copy_number.py" cfg_concat_rc.json --json q_concat_rc.json \
  >q_concat_rc.log 2>&1 && ok "no baseline_depth falls back to read count" \
  || bad "read-count fallback"
got=$("$PY" -c "import json;print(json.load(open('q_concat_rc.json'))['samples']['s1']['copies_per_haploid']['pEpsilon'])")
assert_close "read-count baseline biases low by ~8%" "$got" 25.9 0.6

echo
echo "== Stage 3: edge-trim mode and low unique-bp warning =="
"$PY" "$SCRIPTS/copy_number.py" cfg_concat.json --edge-trim 400 \
  --md q_trim_be.md >q_trim_be.log 2>&1 && ok "block-edge trim (default) ran" \
  || bad "block-edge trim"
assert_has "thin unique block named with its bp" q_trim_be.md \
  "pDelta's copy number rests on only 399 bp"
"$PY" "$SCRIPTS/copy_number.py" cfg_concat.json --edge-trim 400 \
  --trim-mode contig-end --md q_trim_ce.md --json q_trim_ce.json \
  >q_trim_ce.log 2>&1 && ok "contig-end trim ran" || bad "contig-end trim"
assert_lacks "contig-end keeps the unique block, no thin-block warning" \
  q_trim_ce.md "rests on only"
got=$("$PY" -c "import json;print(json.load(open('q_trim_ce.json'))['samples']['s1']['solve']['unique_bp']['pDelta'])")
assert_close "contig-end retains the whole interior unique block" "$got" 1199 1
"$PY" -c "
import json; c=json.load(open('cfg_concat.json'))
del c['samples']['s1']['depth_files']['pEpsilon']; json.dump(c, open('cfg_concat_one.json','w'))"
"$PY" "$SCRIPTS/copy_number.py" cfg_concat_one.json --md q_concat_one.md \
  >q_concat_one.log 2>&1 && ok "construct with no unique class ran" || bad "no unique class"
assert_has "solve fallback for that construct is flagged" q_concat_one.md \
  "pEpsilon .*WLS solve \(FLAGGED"

echo
echo "== Stage 3: fallback with 3 constructs must SHOUT =="
"$PY" -c "
import json; c=json.load(open('cfg_measured.json'))
del c['samples']['parental']['depth_files']
json.dump(c, open('cfg_fallback.json','w'))"
"$PY" "$SCRIPTS/copy_number.py" cfg_fallback.json --md q_fallback.md \
  >q_fallback.log 2>&1 && ok "fallback ran" || bad "fallback"
assert_has "caveat rendered in the REPORT, not just JSON" q_fallback.md "NOT DECOMPOSED"
assert_has "caveat appears before the numbers" \
  <(head -8 q_fallback.md) "UPPER BOUND"
assert_has "route named in the abundance table" q_fallback.md "upper bound from mean coverage"

echo
echo "== Stage 3: single construct =="
"$PY" -c "
import json; c=json.load(open('cfg_measured.json'))
c['identity']='id1.json'
c['mappings']={'pAlpha':c['mappings']['pAlpha']}
c['samples']['parental']['depth_files']={'pAlpha':'depth_pAlpha.txt'}
json.dump(c, open('cfg_one.json','w'))"
"$PY" "$SCRIPTS/copy_number.py" cfg_one.json --md q_one.md >q_one.log 2>&1 \
  && ok "single construct ran" || { bad "single construct"; cat q_one.log; }
assert_lacks "no phantom second construct in the header" q_one.md "A:B ratio"

echo
echo "== Stage 3: integrase must NOT apply the array estimator =="
"$PY" -c "
import json; c=json.load(open('cfg_measured.json'))
c['modality']='integrase'
c['mappings']['pAlpha']['junctions']='junc_lenti.tsv'
json.dump(c, open('cfg_lenti.json','w'))"
"$PY" "$SCRIPTS/copy_number.py" cfg_lenti.json --md q_lenti.md \
  --json q_lenti.json >q_lenti.log 2>&1 \
  && ok "integrase accounting ran" || { bad "integrase accounting"; cat q_lenti.log; }
assert_has "array estimator explicitly withheld" q_lenti.md "was NOT applied"
cps=$("$PY" -c "import json;print(json.load(open('q_lenti.json'))['samples']['parental']['architecture']['copies_per_site'])")
assert_close "copies per site fixed at 1" "$cps" 1 0.001
assert_has "vector-vector junctions queried, not ignored" q_lenti.md "worth explaining|single copies"

echo
echo "== Stage 3: undeclared modality refuses to guess =="
"$PY" -c "
import json; c=json.load(open('cfg_measured.json'))
c['modality']='auto'; json.dump(c, open('cfg_auto.json','w'))"
"$PY" "$SCRIPTS/copy_number.py" cfg_auto.json --md q_auto.md >q_auto.log 2>&1 \
  && ok "auto modality ran" || bad "auto modality"
assert_has "auto modality caveat rendered" q_auto.md "modality is undeclared"

echo
echo "== Stage 3: a sample-name mismatch must not pass silently =="
"$PY" -c "
import json; c=json.load(open('cfg_measured.json'))
c['samples']={'NAME_NOT_IN_TABLE': c['samples']['parental']}
json.dump(c, open('cfg_mismatch.json','w'))"
"$PY" "$SCRIPTS/copy_number.py" cfg_mismatch.json --md q_mismatch.md \
  >q_mismatch.log 2>&1 && ok "mismatch run completed" || bad "mismatch run"
assert_has "sample mismatch raised as a top-level caveat" q_mismatch.md \
  "does not appear in the junction table"
assert_has "mismatch says the derived numbers are meaningless" q_mismatch.md "meaningless"

echo
echo "== Stage 2: construct-name matching (regressions) =="
# A table that calls the construct something else: the documented escape hatch
# must work in the direction the docs state.
sed 's/\tpAlpha\t/\tpAlpha_linearized\t/; s/^pAlpha\t/pAlpha_linearized\t/' \
  junctions_nhej.tsv > junctions_renamed.tsv

alias_run() { # alias_run <label> <extra args...>
  local lbl="$1"; shift
  "$PY" "$SCRIPTS/classify_junctions.py" --identity id3.json --mapping pAlpha \
    --modality nhej --vector-ends 1,5300 --table parental=junctions_renamed.tsv \
    --out-loci "loci_$lbl.tsv" "$@" >"alias_$lbl.log" 2>&1
}

alias_run documented --no-prefix-match --construct-alias pAlpha=pAlpha_linearized \
  && ok "--construct-alias CONSTRUCT=TABLENAME works (documented order)" \
  || { bad "--construct-alias documented order"; cat alias_documented.log; }
assert_has "documented alias actually classified the junctions" \
  loci_documented.tsv "^chr1"

alias_run reversed --no-prefix-match --construct-alias pAlpha_linearized=pAlpha \
  && ok "reversed alias order still accepted" || bad "reversed alias order"
assert_has "reversed order warned about, not silently wrong" alias_reversed.log \
  "documented order is CONSTRUCT=TABLENAME"

alias_run prefix \
  && ok "prefix fallback matches pAlpha_linearized" || bad "prefix fallback"
assert_has "prefix match warns and names the fix" alias_prefix.log \
  "matched construct .* by prefix"

if alias_run strict --no-prefix-match; then
  bad "--no-prefix-match should refuse an unaliased renamed construct"
else
  assert_has "--no-prefix-match refuses rather than guessing" alias_strict.log \
    "silent and catastrophic"
fi

# A construct name must never claim a host contig. Checked at the matcher level
# so the assertion is about the rule, not about one fixture.
"$PY" - "$SCRIPTS" <<'EOF' >matcher.log 2>&1
import sys
sys.path.insert(0, sys.argv[1])
from via_common import NameMatcher
m = NameMatcher(["pAlpha", "pBeta"], {}, allow_prefix=True)
host = ["chr1", "chrX", "AC233755.1", "GL000220.1", "KI270733.1",
        "scaffold_7", "NC_000001.11", "MT"]
claimed = [c for c in host if m.match(c) is not None]
if claimed:
    print("FAIL: construct names claimed host contigs: %s" % claimed)
    sys.exit(1)
for name, table in (("pAlpha", "pAlpha"), ("pBeta", "pBeta_v2")):
    if m.match(table) != name:
        print("FAIL: %r did not match construct %r" % (table, name))
        sys.exit(1)
print("OK")
EOF
if grep -q "^OK$" matcher.log; then
  ok "no host contig or accession is claimed by a construct name"
else
  bad "a host contig was claimed by a construct name"; cat matcher.log
fi

echo
echo "== Error handling =="
"$PY" "$SCRIPTS/classify_junctions.py" --identity id3.json --mapping pDoesNotExist \
  --table s1=junctions_nhej.tsv >e1.log 2>&1 && bad "unknown mapping accepted" \
  || assert_has "unknown mapping refused clearly" e1.log "not in the identity map"

# a well-formed table whose reference names are all host contigs: the one case
# that must never be allowed to pass silently
printf 'chromosome\tposition1\tclip_type_1\tfiller_length\tref2\tposition2\tclip_type_2\tsupporting_reads\tflag\nchr1\t1000\tright_clip\t0\tchr5\t2000\tleft_clip\t10\tPASS\n' > unmatched.tsv
"$PY" "$SCRIPTS/classify_junctions.py" --identity id3.json --mapping pAlpha \
  --table s1=unmatched.tsv >e2.log 2>&1 && bad "table with no construct accepted" \
  || assert_has "no-construct table is fatal, not silent" e2.log "silent and catastrophic"

echo
echo "== Stage 3: per-unique-block spread and no size-threshold drops =="
# pTheta has TWO unique blocks at different depths (20x and 30x): both must
# contribute, bp-weighted, and the disagreement must be reported rather than
# averaged away. pIota's single unique block retains 150 bp, which the retired
# min_class_bp of 200 silently dropped -- it must now be used.
cat > cfg_spread.json <<'JSON'
{"identity":"identity_spread.json","genome_size":3100000000,"read_length":150,
 "ploidy":2,"primary_mapping":"pTheta","modality":"nhej",
 "samples":{"s1":{"dedup_pairs":299666666,"baseline_depth":29.0,
   "depth_files":{"pTheta":"depth_spread_pTheta.txt",
                  "pIota":"depth_spread_pIota.txt"}}},
 "mappings":{"pTheta":{},"pIota":{}}}
JSON
"$PY" "$SCRIPTS/copy_number.py" cfg_spread.json --md q_spread.md \
  --json q_spread.json >q_spread.log 2>&1 \
  && ok "per-block spread case ran" || { bad "per-block spread case"; cat q_spread.log; }
sp() { "$PY" -c "import json;d=json.load(open('q_spread.json'))['samples']['s1'];print($1)"; }
assert_close "both unique blocks contribute (2 blocks)" \
  "$(sp "d['solve']['unique_blocks']['pTheta']['n_blocks']")" 2 0
assert_close "both unique blocks' bp retained" \
  "$(sp "d['solve']['unique_blocks']['pTheta']['retained_bp']")" 1800 2
assert_close "low-depth unique block reported on its own" \
  "$(sp "d['solve']['unique_blocks']['pTheta']['min']")" 0.69 0.05
assert_close "high-depth unique block reported on its own" \
  "$(sp "d['solve']['unique_blocks']['pTheta']['max']")" 1.03 0.05
assert_close "spread reported as a relative range" \
  "$(sp "d['solve']['unique_blocks']['pTheta']['relative_range']")" 0.39 0.06
assert_close "reported copy number is their bp-weighted mean" \
  "$(sp "d['copies_per_haploid']['pTheta']")" 0.86 0.05
assert_has "disagreeing blocks flagged, not silently averaged" q_spread.md \
  "pTheta's unique blocks disagree: 2 blocks give"
assert_has "spread caveat reaches the top of the report" <(head -12 q_spread.md) \
  "unique blocks disagree"
assert_has "per-block table rendered" q_spread.md "Per-unique-block estimates"
# 150 bp of unique sequence: kept and used, with a low-confidence flag
got=$(sp "d['solve']['source']['pIota']")
if [ "$got" = "unique block" ]; then
  ok "small unique block no longer dropped on size (source: unique block)"
else
  bad "small unique block was excluded (source: $got)"
fi
assert_close "its retained bp counted in full" \
  "$(sp "d['solve']['unique_blocks']['pIota']['retained_bp']")" 150 2
assert_has "thin unique sequence still flagged" q_spread.md \
  "pIota's copy number rests on only 150 bp"

echo
echo "== Stage 3: trim_boundaries policy and deprecated keys =="
# The rule, not a fixture: an end abutting a same-membership block is not
# trimmed; an end abutting a different class, or a contig terminus, is.
"$PY" - "$SCRIPTS" <<'EOF' >trim.log 2>&1
import sys
sys.path.insert(0, sys.argv[1])
from copy_number import _trim_flags
B = [{"membership": ["A"]}, {"membership": ["A", "B"]}, {"membership": ["A"]}]
lin = _trim_flags(B, False)
assert lin == [(True, True), (True, True), (True, True)], lin
same = _trim_flags([{"membership": ["A"]}, {"membership": ["A"]}], False)
assert same == [(True, False), (False, True)], same
circ = _trim_flags([{"membership": ["A"]}, {"membership": ["A", "B"]},
                    {"membership": ["A"]}], True)
assert circ == [(False, True), (True, True), (True, False)], circ
print("OK")
EOF
if grep -q "^OK$" trim.log; then
  ok "trim_boundaries trims class edges and termini, not same-membership joins"
else
  bad "trim_boundaries policy"; cat trim.log
fi
"$PY" "$SCRIPTS/copy_number.py" cfg_concat.json --edge-trim 400 \
  --min-class-bp 200 --trim-mode block-edge --json q_depr.json \
  >q_depr.log 2>&1 && ok "deprecated keys still accepted" || bad "deprecated keys"
assert_has "--edge-trim deprecation named" q_depr.log "edge-trim is deprecated"
assert_has "--trim-mode deprecation named" q_depr.log "trim-mode is deprecated"
assert_has "min_class_bp deprecation named" q_depr.log \
  "min_class_bp is deprecated and ignored"

echo
echo "== Final stage: reconcile_sites, two constructs sharing backbone =="
# pDelta 1.5 copies and pEpsilon 28 copies. Every shared locus is reported by
# BOTH mappings and must collapse to one; attribution must come out of the
# priority order, not out of read support.
for m in pDelta pEpsilon; do
  case $m in pDelta) e=1200 ;; pEpsilon) e=3000 ;; esac
  "$PY" "$SCRIPTS/classify_junctions.py" --identity identity_concat.json \
    --mapping $m --modality nhej --vector-ends $e \
    --table s1="junctions_recon_${m}_s1.tsv" \
    --table s2="junctions_recon_${m}_s2.tsv" \
    --out-loci "loci_recon_$m.tsv" --out-junctions "junc_recon_$m.tsv" \
    >"s2_recon_$m.log" 2>&1 || bad "classify_junctions for $m"
done
cat > cfg_recon.json <<'JSON'
{"identity":"identity_concat.json","genome_size":3100000000,"read_length":150,
 "ploidy":2,"primary_mapping":"pEpsilon","modality":"nhej",
 "samples":{"s1":{"dedup_pairs":323640000,"baseline_depth":29.0,
   "depth_files":{"pDelta":"depth_concat_pDelta.txt",
                  "pEpsilon":"depth_concat_pEpsilon.txt"}}},
 "mappings":{"pDelta":{"junctions":"junc_recon_pDelta.tsv"},
             "pEpsilon":{"junctions":"junc_recon_pEpsilon.tsv"}}}
JSON
"$PY" "$SCRIPTS/copy_number.py" cfg_recon.json --json cn_recon.json \
  >cn_recon.log 2>&1 && ok "copy number for the reconcile case" \
  || { bad "copy number for the reconcile case"; cat cn_recon.log; }

# --min-copies-per-site 10 is the measured copies per site of a concatemer
# line; it is what makes the ceiling binding rather than nominal.
"$PY" "$SCRIPTS/reconcile_sites.py" --identity identity_concat.json \
  --loci pDelta=loci_recon_pDelta.tsv --loci pEpsilon=loci_recon_pEpsilon.tsv \
  --copy-number cn_recon.json --ploidy 2 --min-copies-per-site 10 \
  --vector-ends pDelta=1200 --vector-ends pEpsilon=3000 \
  --out-called called.tsv --out-subthreshold sub.tsv \
  --summary site_summary.json --md site_summary.md >recon.log 2>&1 \
  && ok "reconcile_sites ran" || { bad "reconcile_sites"; cat recon.log; }
rc() { "$PY" -c "import json;d=json.load(open('site_summary.json'));print($1)"; }

# 1. collapse: every shared locus was reported twice and must be ONE event
assert_close "entries read from both mapping tables" "$(rc "d['n_entries']")" 15 0
assert_close "cross-mapping duplicates collapsed" "$(rc "d['n_unified']")" 9 0
assert_has "collapse stated in the report" site_summary.md \
  "6 duplicate report\(s\) of one junction were collapsed"

# 2. retention: the single-read single-sample locus is the tail, not a site
assert_close "called loci" "$(rc "d['n_called']")" 8 0
assert_close "sub-threshold tail separated" "$(rc "d['n_sub']")" 1 0
assert_has "the tail locus is in the sub-threshold table" sub.tsv "^chr13"
assert_lacks "the tail locus is NOT in the called table" called.tsv "^chr13"
assert_has "tail excluded for both reasons" sub.tsv "min-reads.*min-samples"

# 3. attribution by priority
assert_has "unique-sequence locus is certain" called.tsv \
  "^chr2	5000000.*pEpsilon	certain	unique"
assert_has "vector-end locus goes to the construct whose end it is" called.tsv \
  "^chr5	2000000.*pDelta	inferred	vector-end"
assert_has "shared locus placed by the ceiling, not by support" called.tsv \
  "^chr6	1000000.*pEpsilon	inferred	copy-ceiling"
assert_has "the locus beyond every ceiling is unresolved" called.tsv \
  "^chr12	1000000	1000000	1			unresolved		pDelta,pEpsilon"
assert_has "unresolved locus names its candidates" called.tsv \
  "^chr12.*pDelta,pEpsilon"
assert_has "unresolved explained by the full ceilings" site_summary.md \
  "ceiling is already full"

# no construct may exceed its copy-number ceiling
"$PY" - <<'EOF' >cap.log 2>&1
import json
d = json.load(open("site_summary.json"))
for n, c in d["per_construct"].items():
    cap = c["ceiling"]
    assert cap is not None, n
    assert c["attributed"] <= cap, (n, c["attributed"], cap)
print("OK")
EOF
if grep -q "^OK$" cap.log; then
  ok "no construct exceeds its copy-number ceiling"
else
  bad "a construct was given more sites than its ceiling"; cat cap.log
fi
assert_close "high-copy construct holds exactly its ceiling" \
  "$(rc "d['per_construct']['pEpsilon']['attributed']")" 6 0
assert_close "low-copy construct holds exactly its ceiling" \
  "$(rc "d['per_construct']['pDelta']['attributed']")" 1 0
assert_close "ceiling of the low-copy construct is 1, not 6" \
  "$(rc "d['per_construct']['pDelta']['ceiling']")" 1 0
assert_close "exactly one unresolved locus" "$(rc "len(d['unresolved'])")" 1 0

# no locus may be listed under two constructs
"$PY" - <<'EOF' >dup.log 2>&1
import csv
from collections import Counter
rows = list(csv.DictReader(open("called.tsv"), delimiter="\t"))
seen = Counter((r["contig"], r["start"]) for r in rows)
dups = [k for k, v in seen.items() if v > 1]
assert not dups, "locus listed more than once: %s" % dups
for r in rows:
    assert "," not in r["construct"], r
    if r["confidence"] == "unresolved":
        assert r["construct"] == "", r
    else:
        assert r["construct"], r
print("OK")
EOF
if grep -q "^OK$" dup.log; then
  ok "no locus appears under two constructs, or twice under one"
else
  bad "a locus was listed under two constructs"; cat dup.log
fi
assert_has "site count is a bounded range per construct" site_summary.md \
  "\*\*2-6\*\*"
assert_has "report says it is authoritative, not to be re-derived" \
  site_summary.md "do not re-derive it"

echo
echo "== Final stage: equal copy number must NOT force attribution =="
# Same loci, two constructs of EQUAL copy number. Nothing but the copy number
# distinguishes the shared loci, so the procedure must decline to attribute
# them rather than invent a winner.
"$PY" "$SCRIPTS/reconcile_sites.py" --identity identity_concat.json \
  --loci pDelta=loci_recon_pDelta.tsv --loci pEpsilon=loci_recon_pEpsilon.tsv \
  --copies pDelta=10 --copies pEpsilon=10 --ploidy 2 \
  --vector-ends pDelta=1200 --vector-ends pEpsilon=3000 \
  --out-called called_sym.tsv --summary sym.json --md sym.md >sym.log 2>&1 \
  && ok "equal-copy case ran" || { bad "equal-copy case"; cat sym.log; }
sy() { "$PY" -c "import json;d=json.load(open('sym.json'));print($1)"; }
assert_close "every shared locus is unresolved when copies are equal" \
  "$(sy "len(d['unresolved'])")" 5 0
assert_close "no shared locus attributed by the ceiling to pDelta" \
  "$(sy "d['per_construct']['pDelta']['copy_ceiling']")" 0 0
assert_close "no shared locus attributed by the ceiling to pEpsilon" \
  "$(sy "d['per_construct']['pEpsilon']['copy_ceiling']")" 0 0
assert_has "equal copy numbers given as the reason" sym.md \
  "within 2(\.[0-9])?x"
assert_close "sequence-level evidence still attributed" \
  "$(sy "d['per_construct']['pEpsilon']['certain']")" 2 0
assert_close "vector-end evidence still attributed" \
  "$(sy "d['per_construct']['pDelta']['vector_end']")" 1 0

echo
echo "== Final stage: the ceiling outranks read support =="
# pDelta has 100x the support of pEpsilon at every shared locus, but only
# pEpsilon's copy number can source them. Support must not win.
"$PY" - <<'EOF'
import csv
rows = list(csv.DictReader(open("loci_recon_pDelta.tsv"), delimiter="\t"))
for r in rows:
    for k in r:
        if k.startswith("support_") and r[k]:
            r[k] = str(int(r[k]) * 100)
with open("loci_recon_pDelta_loud.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=rows[0].keys(), delimiter="\t")
    w.writeheader(); w.writerows(rows)
EOF
"$PY" "$SCRIPTS/reconcile_sites.py" --identity identity_concat.json \
  --loci pDelta=loci_recon_pDelta_loud.tsv \
  --loci pEpsilon=loci_recon_pEpsilon.tsv \
  --copy-number cn_recon.json --ploidy 2 --min-copies-per-site 10 \
  --vector-ends pDelta=1200 --vector-ends pEpsilon=3000 \
  --out-called called_loud.tsv --summary loud.json >loud.log 2>&1 \
  && ok "support-inflated case ran" || { bad "support-inflated case"; cat loud.log; }
assert_close "inflated support does not lift pDelta past its ceiling" \
  "$("$PY" -c "import json;print(json.load(open('loud.json'))['per_construct']['pDelta']['attributed'])")" 1 0
assert_has "shared loci still go to the construct with the copies" \
  called_loud.tsv "^chr6	1000000.*pEpsilon	inferred	copy-ceiling"

echo
echo "== Final stage: guard rails =="
"$PY" "$SCRIPTS/reconcile_sites.py" --identity identity_concat.json \
  --loci pDelta=loci_recon_pDelta.tsv --copy-number cn_recon.json \
  --min-samples 3 --summary one.json >one.log 2>&1 \
  && ok "impossible --min-samples does not crash" || bad "impossible --min-samples"
assert_has "it says no locus can be called, loudly" one.log \
  "only 2 sample\(s\) were supplied, so NO locus can be called"
"$PY" "$SCRIPTS/reconcile_sites.py" --identity identity_concat.json \
  --loci pDelta=loci_recon_pDelta.tsv --loci pEpsilon=loci_recon_pEpsilon.tsv \
  --summary nocn.json >nocn.log 2>&1 \
  && ok "missing copy number does not crash" || bad "missing copy number"
assert_has "no copy number means no ceiling attribution, and says so" nocn.log \
  "No copy number for pDelta, pEpsilon"
assert_has "assumed ploidy flagged" nocn.log "PLOIDY ASSUMED"
"$PY" "$SCRIPTS/reconcile_sites.py" --identity identity_concat.json \
  --loci pNotReal=loci_recon_pDelta.tsv >bad_recon.log 2>&1 \
  && bad "unknown construct accepted" \
  || assert_has "unknown construct refused clearly" bad_recon.log \
       "not in the identity map"
"$PY" "$SCRIPTS/reconcile_sites.py" --identity identity_concat.json \
  --loci pDelta=junc_recon_pDelta.tsv >wrongfile.log 2>&1 \
  && bad "a junction table accepted in place of a loci table" \
  || assert_has "wrong input file refused with the fix named" wrongfile.log \
       "written by --out-loci"

echo
echo "== Shared modules: one modality table, one loci-table format =="
"$PY" - "$SCRIPTS" <<'EOF' >shared.log 2>&1
import os, sys, tempfile
sys.path.insert(0, sys.argv[1])
import classify_junctions, copy_number, reconcile_sites, via_common
M = via_common.MODALITIES
assert classify_junctions.MODALITIES is M
for name, m in M.items():
    assert (name in copy_number.FIXED_COPIES_PER_SITE) == (m["concatemers"] is False), name
# what Stage 2 writes is what Stage 4 reads, through one pair of functions
assert reconcile_sites.read_loci is via_common.read_loci
path = os.path.join(tempfile.mkdtemp(), "loci.tsv")
via_common.write_loci(path, ["ref", "other"], [
    {"contig": "chr1", "start": 100, "end": 113, "flanks": 2, "footprint": 12,
     "vector_pos": 3200, "region": "unique", "membership": "pAlpha",
     "feature": "cargo", "evidence_tier": "certain", "support_ref": 56,
     "support_other": 30, "ratio_other": 0.54, "lo_other": 0.3, "hi_other": 0.9},
    {"contig": "chr2", "start": 7, "end": 7, "flanks": 1, "footprint": None,
     "vector_pos": 10, "region": "shared", "membership": "pAlpha + pBeta",
     "feature": "backbone", "evidence_tier": "inferred", "support_ref": 4}])
rows, samples = via_common.read_loci(path, "pAlpha")
assert samples == ["ref", "other"], samples
assert [r["contig"] for r in rows] == ["chr1", "chr2"]
assert rows[0]["footprint"] == 12 and rows[1]["footprint"] is None
assert rows[1]["membership"] == frozenset({"pAlpha", "pBeta"})
assert rows[0]["support"] == {"ref": 56, "other": 30}
assert rows[1]["support"] == {"ref": 4, "other": 0}
print("OK")
EOF
if grep -q "^OK$" shared.log; then
  ok "both stages read the one modality table"
  ok "a loci table written by Stage 2 round-trips through the Stage 4 reader"
else
  bad "shared modality table / loci format"; cat shared.log
fi

echo
echo "== Checks: a stop condition is an exit status, not a paragraph =="
# expect_exit <label> <wanted status> <log> <command...>
expect_exit() {
  local lbl="$1" want="$2" log="$3"; shift 3
  "$@" >"$log" 2>&1; local got=$?
  if [ "$got" -eq "$want" ]; then ok "$lbl (exit $got)"; else
    bad "$lbl (exit $got, expected $want)"; cat "$log"; fi
}
cj() { "$PY" "$SCRIPTS/classify_junctions.py" --identity id3.json --mapping pAlpha "$@"; }
# a declared modality the fingerprint contradicts
expect_exit "wrong modality still exits 0 without --strict" 0 chk_wrong.log \
  cj --modality nhej --vector-ends 1,5300 --table parental=junctions_integrase.tsv \
  --md chk_wrong.md --json chk_wrong.json
assert_has "the failing check is printed with its id" chk_wrong.log \
  "CHECK FAILED \[modality_fingerprint\]"
assert_has "the failing check is a banner at the top of the report" \
  <(head -8 chk_wrong.md) "FAILED .modality_fingerprint."
got=$("$PY" -c "import json;d=json.load(open('chk_wrong.json'));print([c['status'] for c in d['checks'] if c['id']=='modality_fingerprint'][0])")
[ "$got" = "fail" ] && ok "the failing check is machine-readable in the stage JSON" \
  || bad "stage JSON check status was '$got'"
expect_exit "--strict turns it into exit status 3" 3 chk_wrong_s.log \
  cj --modality nhej --vector-ends 1,5300 --table parental=junctions_integrase.tsv --strict
assert_has "--strict names the blocking checks" chk_wrong_s.log \
  "STOP \(classify_junctions\).*modality_fingerprint"
expect_exit "--accept acknowledges named checks and unblocks" 0 chk_wrong_a.log \
  cj --modality nhej --vector-ends 1,5300 --table parental=junctions_integrase.tsv \
  --strict --accept modality_fingerprint --accept footprint --md chk_wrong_a.md
assert_has "an accepted check stays visible in the report" chk_wrong_a.md \
  "ACCEPTED .modality_fingerprint."
expect_exit "a correct declaration passes --strict" 0 chk_ok.log \
  cj --modality integrase --vector-ends 1,5300 --table parental=junctions_integrase.tsv \
  --strict --md chk_ok.md
assert_lacks "no banner when nothing failed" chk_ok.md "Stop: resolve these"

# a fixed-end modality with no vector ends: completes as before, but says so
expect_exit "fixed-end modality without --vector-ends still exits 0" 0 chk_noends.log \
  cj --modality integrase --table parental=junctions_integrase.tsv --md chk_noends.md
assert_has "missing vector ends raised as a failing check" chk_noends.md \
  "FAILED .vector_ends."
assert_has "the message gives the construct length from the identity map" \
  chk_noends.md "pAlpha is 5300 bp.*1,5300"
expect_exit "--strict refuses a fixed-end modality without vector ends" 3 \
  chk_noends_s.log \
  cj --modality integrase --table parental=junctions_integrase.tsv --strict
expect_exit "an unusable clip side stops a --strict run" 3 chk_strand.log \
  cj --modality nhej --vector-ends 1,5300 --table parental=junctions_nhej_strand.tsv --strict
assert_has "the clip-side check names the fix" chk_strand.log \
  "CHECK FAILED \[clip_side\].*clip-map"

cn() { "$PY" "$SCRIPTS/copy_number.py" "$@"; }
expect_exit "measured copy number passes --strict" 0 chk_cn_ok.log \
  cn cfg_measured.json --json chk_cn_ok.json --strict
expect_exit "an undecomposed upper bound stops a --strict run" 3 chk_cn_fb.log \
  cn cfg_fallback.json --json chk_cn_fb.json --strict
assert_has "the route check says the number was not measured" chk_cn_fb.log \
  "CHECK FAILED \[copy_number_route\].*measured in no sample"
# a sample with depth and a sample without: the measured one still stands, so
# the check warns and names both rather than failing the whole stage
# s1 has depth, s2 has junctions and a read count but no depth: the cold-run
# case, and the one the old all-or-nothing check failed the whole stage on.
"$PY" -c "
import json; c=json.load(open('cfg_recon.json'))
s1=c['samples']['s1']
c['samples']['s2']={k:v for k,v in s1.items() if k!='depth_files'}
json.dump(c, open('cfg_mixed.json','w'))"
expect_exit "a sample without depth does not spoil the one with it" 0 \
  chk_cn_mixed.log cn cfg_mixed.json --json chk_cn_mixed.json --strict
got=$("$PY" -c "
import json
d = json.load(open('chk_cn_mixed.json'))
c = [c for c in d['checks'] if c['id'] == 'copy_number_route'][0]
print(c['status'], d['samples']['s1']['solve']['method'].split()[0],
      d['samples']['s2']['solve']['method'].split()[0])")
[ "$got" = "warn measured NOT" ] \
  && ok "the route check warns per sample and names which was measured" \
  || bad "mixed-depth route check was '$got'"
assert_has "the warning names the sample that has depth and the one that lacks it" \
  chk_cn_mixed.log "CHECK WARNING \[copy_number_route\].*depth in s1 and not in s2"
assert_close "the measured sample keeps its measured copy number" \
  "$("$PY" -c "import json;print(json.load(open('chk_cn_mixed.json'))['samples']['s1']['copies_per_haploid']['pEpsilon'])")" \
  28 1.0
expect_exit "the same run without --strict still exits 0" 0 chk_cn_fb0.log \
  cn cfg_fallback.json --json chk_cn_fb0.json
expect_exit "an undeclared modality stops a --strict run" 3 chk_cn_auto.log \
  cn cfg_auto.json --json chk_cn_auto.json --strict
expect_exit "a sample-name mismatch stops a --strict run" 3 chk_cn_mm.log \
  cn cfg_mismatch.json --json chk_cn_mm.json --strict
got=$("$PY" -c "
import json
ok = {c['id']: c['status'] for c in json.load(open('chk_cn_ok.json'))['checks']}
fb = {c['id']: c['status'] for c in json.load(open('chk_cn_fb.json'))['checks']}
print(ok['copy_number_route'], fb['copy_number_route'], ok['modality_declared'])")
[ "$got" = "pass fail pass" ] && ok "copy-number JSON carries the checks" \
  || bad "copy-number JSON checks were '$got'"

rs() { "$PY" "$SCRIPTS/reconcile_sites.py" --identity identity_concat.json "$@"; }
expect_exit "a reconciliation with nothing wrong passes --strict" 0 chk_rs_ok.log \
  rs --loci pDelta=loci_recon_pDelta.tsv --loci pEpsilon=loci_recon_pEpsilon.tsv \
  --copy-number cn_recon.json --ploidy 2 --min-copies-per-site 10 \
  --vector-ends pDelta=1200 --vector-ends pEpsilon=3000 --summary chk_rs_ok.json --strict
expect_exit "an impossible --min-samples stops a --strict run" 3 chk_rs_ms.log \
  rs --loci pDelta=loci_recon_pDelta.tsv --copy-number cn_recon.json \
  --min-samples 3 --summary chk_rs_ms.json --strict
expect_exit "a missing copy number stops a --strict run" 3 chk_rs_nocn.log \
  rs --loci pDelta=loci_recon_pDelta.tsv --loci pEpsilon=loci_recon_pEpsilon.tsv \
  --summary chk_rs_nocn.json --strict
assert_has "site summary carries the checks" chk_rs_ms.json '"id": "min_samples"'

echo
echo "== Run manifest: no value is carried between stages by hand =="
mkdir -p mf
cp pDelta.gb pEpsilon.gb junctions_recon_*.tsv depth_concat_*.txt mf/
cat > mf/run.json <<'JSON'
{"genome_size":3100000000,"read_length":150,"ploidy":2,"modality":"nhej",
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
# run from outside mf/: every path resolves against the manifest
expect_exit "Stage 1 with --run" 0 mf1.log \
  "$PY" "$SCRIPTS/construct_identity.py" mf/pDelta.gb mf/pEpsilon.gb --run mf/run.json
for m in pDelta pEpsilon; do
  expect_exit "Stage 2 with --run ($m)" 0 "mf2_$m.log" \
    "$PY" "$SCRIPTS/classify_junctions.py" --run mf/run.json --mapping $m --strict
done
expect_exit "Stage 3 with --run and no config file" 0 mf3.log \
  "$PY" "$SCRIPTS/copy_number.py" --run mf/run.json --strict
expect_exit "Stage 4 with --run" 0 mf4.log \
  "$PY" "$SCRIPTS/reconcile_sites.py" --run mf/run.json --strict --min-copies-per-site 10

# the manifest route must give exactly the flag-and-config route's answer
for pair in "loci_recon_pDelta.tsv:mf/loci_pDelta.tsv" \
            "junc_recon_pEpsilon.tsv:mf/classified_pEpsilon.tsv" \
            "called.tsv:mf/called_sites.tsv" "sub.tsv:mf/subthreshold_sites.tsv"; do
  if cmp -s "${pair%%:*}" "${pair##*:}"; then
    ok "manifest route reproduces ${pair%%:*} byte for byte"
  else
    bad "manifest route differs from the flag route for ${pair%%:*}"
  fi
done
"$PY" - <<'EOF' >mf_same.log 2>&1
import json
a = json.load(open("cn_recon.json")); b = json.load(open("mf/accounting.json"))
for k in ("copies_per_haploid", "copies_per_cell", "solve", "architecture"):
    assert a["samples"]["s1"][k] == b["samples"]["s1"][k], k
a = json.load(open("site_summary.json")); b = json.load(open("mf/site_summary.json"))
for k in ("per_construct", "unresolved", "ceilings", "caveats",
          "n_entries", "n_unified", "n_called", "n_sub"):
    assert a[k] == b[k], k
# same numbers; only the recorded provenance of the ploidy differs by route
pa, pb = dict(a["parameters"]), dict(b["parameters"])
assert pb.pop("ploidy_source") == "run manifest", b["parameters"]
pa.pop("ploidy_source")
assert pa == pb, (pa, pb)
print("OK")
EOF
grep -q "^OK$" mf_same.log && ok "manifest route reproduces the copy numbers and the site result" \
  || { bad "manifest route changed a number"; cat mf_same.log; }
mfq() { "$PY" -c "import json;d=json.load(open('mf/$1'));print($2)"; }
assert_close "detected-locus count carried from Stage 2 to Stage 3" \
  "$(mfq accounting.json "d['samples']['s1']['detected_loci']")" 9 0
assert_has "Stage 3 reconciles instead of reporting 'not reconciled'" mf3.log \
  "SHORTFALL|AGREES|EXCESS"
got=$(mfq run.json "sorted(d['stages'])")
[ "$got" = "['classify', 'copy_number', 'identity', 'reconcile']" ] \
  && ok "every stage recorded its outputs in the manifest" \
  || bad "manifest stages were $got"
got=$(mfq run.json "d['stages']['reconcile']['sites']['pEpsilon'], d['genome_size'], sorted(k for k in d if k != 'stages')")
assert_has "declared facts untouched, site range recorded" <(echo "$got") \
  "^2-6 3100000000 \['genome_size', 'modality', 'ploidy', 'primary_mapping', 'read_length', 'samples', 'vector_ends'\]"
assert_has "a sample without dedup_pairs is named, not silently dropped" \
  mf/accounting.json "s2 declare no dedup_pairs"

# an explicit flag beats the manifest
expect_exit "an explicit --modality overrides the manifest and is checked" 3 mf_over.log \
  "$PY" "$SCRIPTS/classify_junctions.py" --run mf/run.json --mapping pDelta \
  --modality integrase --strict --md mf/over.md --json mf/over.json \
  --out-loci mf/over_loci.tsv --out-junctions mf/over_junc.tsv
assert_has "override visible in the result" mf_over.log "modality=integrase"
"$PY" "$SCRIPTS/classify_junctions.py" --run mf/run.json --mapping pDelta \
  >mf_restore.log 2>&1 || bad "restoring the pDelta record"

# stopped runs name what is missing
expect_exit "copy_number with neither config nor --run is refused" 2 mf_none.log \
  "$PY" "$SCRIPTS/copy_number.py"
expect_exit "a missing manifest is refused" 2 mf_missing.log \
  "$PY" "$SCRIPTS/classify_junctions.py" --run mf/nope.json --mapping pDelta
assert_has "missing manifest names the template" mf_missing.log "run.example.json"
expect_exit "an output that would overwrite an input table is refused" 2 mf_clobber.log \
  "$PY" "$SCRIPTS/classify_junctions.py" --run mf/run.json --mapping pDelta \
  --out-junctions mf/junctions_recon_pDelta_s1.tsv
assert_has "refusal says the input would be overwritten" mf_clobber.log \
  "refusing to overwrite the input"

echo
echo "== Run manifest: a nominal ceiling is said out loud, never applied =="
# pAlpha's line has internal array junctions, so Stage 3 measures >1 copies
# per site. Stage 4 must flag that the default divisor of 1 leaves the ceiling
# nominal -- and must NOT change the divisor on its own.
mkdir -p mfa
cp pAlpha.gb pBeta.gb pGamma.gb junctions_nhej.tsv depth_pAlpha.txt \
   depth_pBeta.txt depth_pGamma.txt mfa/
cat > mfa/run.json <<'JSON'
{"genome_size":3100000000,"read_length":150,"ploidy":2,"modality":"nhej",
 "primary_mapping":"pAlpha","vector_ends":{"pAlpha":[1,5300]},
 "samples":{"parental":{"junctions":{"pAlpha":"junctions_nhej.tsv"},
   "dedup_pairs":299666666,
   "depth_files":{"pAlpha":"depth_pAlpha.txt","pBeta":"depth_pBeta.txt",
                  "pGamma":"depth_pGamma.txt"}}}}
JSON
"$PY" "$SCRIPTS/construct_identity.py" mfa/pAlpha.gb mfa/pBeta.gb mfa/pGamma.gb \
  --run mfa/run.json >mfa1.log 2>&1 || bad "manifest Stage 1 (pAlpha line)"
"$PY" "$SCRIPTS/classify_junctions.py" --run mfa/run.json --mapping pAlpha --strict \
  >mfa2.log 2>&1 || { bad "manifest Stage 2 (pAlpha line)"; cat mfa2.log; }
"$PY" "$SCRIPTS/copy_number.py" --run mfa/run.json --strict >mfa3.log 2>&1 \
  || { bad "manifest Stage 3 (pAlpha line)"; cat mfa3.log; }
"$PY" - <<'EOF' >mfa_same.log 2>&1
import json
a = json.load(open("q_measured.json")); b = json.load(open("mfa/accounting.json"))
assert a["samples"]["parental"]["copies_per_haploid"] == \
       b["samples"]["parental"]["copies_per_haploid"]
assert a["samples"]["parental"]["architecture"] == b["samples"]["parental"]["architecture"]
assert b["samples"]["parental"]["detected_loci"] == 4
print("OK")
EOF
grep -q "^OK$" mfa_same.log \
  && ok "three-construct manifest run reproduces the hand-written config's numbers" \
  || { bad "three-construct manifest run changed a number"; cat mfa_same.log; }
expect_exit "Stage 4 runs with the default divisor" 0 mfa4.log \
  "$PY" "$SCRIPTS/reconcile_sites.py" --run mfa/run.json --min-samples 1 --strict
assert_has "nominal ceiling raised as a caveat with the measured value" \
  mfa/site_summary.md "left at 1, but the copy-number stage measured"
got=$("$PY" -c "
import json, math
d = json.load(open('mfa/site_summary.json'))
st = [c['status'] for c in d['checks'] if c['id'] == 'ceiling_divisor'][0]
cn = d['per_construct']['pAlpha']['copies_per_haploid']
print(st, d['parameters']['min_copies_per_site'],
      d['per_construct']['pAlpha']['ceiling'] == math.ceil(cn * 2 / 1.0))")
[ "$got" = "warn 1.0 True" ] \
  && ok "the measured copies per site is flagged but not applied to the ceiling" \
  || bad "ceiling divisor handling was '$got'"
"$PY" "$SCRIPTS/reconcile_sites.py" --run mfa/run.json --min-samples 1 \
  --min-copies-per-site 3 --summary mfa/explicit.json --md mfa/explicit.md \
  --out-called mfa/explicit.tsv --out-subthreshold mfa/explicit_sub.tsv \
  >mfa5.log 2>&1 || bad "explicit divisor run"
assert_lacks "an explicit divisor silences the nominal-ceiling caveat" \
  mfa/explicit.md "left at 1"

echo
echo "== Figures: make_figures.py =="
# Every figure is drawn from the standard outputs above. The assertions are on
# structure (well-formed SVG, one chromosome per .fai contig, one tick per
# called locus), never on pixels.
svg_ok() { # svg_ok <label> <file>: non-empty, parses, root element is <svg>
  "$PY" - "$2" <<'EOF' >/dev/null 2>&1
import os, sys, xml.etree.ElementTree as ET
p = sys.argv[1]
assert os.path.getsize(p) > 200
assert ET.parse(p).getroot().tag == "{http://www.w3.org/2000/svg}svg"
EOF
  if [ $? -eq 0 ]; then ok "$1"; else bad "$1 ($2 missing or not well-formed SVG)"; fi
}

# multi-construct case: three constructs sharing backbone, one sample
"$PY" "$SCRIPTS/make_figures.py" --out figs3 --identity id3.json \
  --config cfg_measured.json --copy-number q_measured.json >figs3.log 2>&1 \
  && ok "make_figures ran on the three-construct case" \
  || { bad "make_figures three-construct case"; cat figs3.log; }
for c in pAlpha pBeta pGamma; do
  svg_ok "coverage track for $c is well-formed SVG" "figs3/coverage_parental_$c.svg"
done
svg_ok "three-construct copy-number figure is well-formed SVG" figs3/copy_number.svg
assert_has "absent .fai skipped with a note, not a failure" figs3.log \
  "skipped genome_overview: no host .fai"
"$PY" - <<'EOF' >figs3_check.log 2>&1
import json, xml.etree.ElementTree as ET
NS = "{http://www.w3.org/2000/svg}"
ident = json.load(open("id3.json"))
for c, d in ident["constructs"].items():
    root = ET.parse("figs3/coverage_parental_%s.svg" % c).getroot()
    g = [e for e in root.iter(NS + "g") if e.get("class") == "block"]
    assert len(g) == len(d["blocks"]), (c, len(g), len(d["blocks"]))
    mem = {e.get("data-membership") for e in g}
    assert "pAlpha,pBeta,pGamma" in mem, mem
    cls = lambda k: [e for e in root.iter() if k in (e.get("class") or "").split()]
    assert cls("baseline"), "no baseline line"
    assert cls("trimmed"), "no trimmed region"
    assert cls("depth"), "no depth trace"
root = ET.parse("figs3/copy_number.svg").getroot()
bars = [e for e in root.iter(NS + "g") if e.get("class") == "cn-bar"]
assert {b.get("data-construct") for b in bars} == {"pAlpha", "pBeta", "pGamma"}, bars
print("OK")
EOF
if grep -q "^OK$" figs3_check.log; then
  ok "coverage tracks shade every identity block, with baseline and trim"
else
  bad "coverage track content"; cat figs3_check.log
fi

# full case: two constructs, two copy-number runs (two trim policies) side by
# side, the called-site table, the site summary and a host .fai
"$PY" "$SCRIPTS/make_figures.py" --out figs --identity identity_concat.json \
  --config cfg_recon.json --copy-number cn_recon.json \
  --copy-number q_trim_ce.json --called called.tsv \
  --site-summary site_summary.json --fai contigs.fai >figs.log 2>&1 \
  && ok "make_figures ran on the full two-construct case" \
  || { bad "make_figures full case"; cat figs.log; }
for f in coverage_s1_pDelta coverage_s1_pEpsilon genome_overview copy_number \
         architecture; do
  svg_ok "$f.svg is well-formed SVG" "figs/$f.svg"
done
"$PY" - <<'EOF' >figs_check.log 2>&1
import csv, json, xml.etree.ElementTree as ET
NS = "{http://www.w3.org/2000/svg}"
def groups(path, cls):
    return [e for e in ET.parse(path).getroot().iter(NS + "g")
            if e.get("class") == cls]
fai = [l.split("\t")[0] for l in open("contigs.fai") if l.strip()]
chroms = [g.get("data-name") for g in groups("figs/genome_overview.svg", "chrom")]
assert chroms == fai, ("chromosomes drawn", chroms, "fai", fai)
called = list(csv.DictReader(open("called.tsv"), delimiter="\t"))
ticks = groups("figs/genome_overview.svg", "locus")
assert len(ticks) == len(called) == 8, (len(ticks), len(called))
conf = sorted(t.get("data-confidence") for t in ticks)
assert conf == sorted(r["confidence"] for r in called), conf
for t in ticks:
    if t.get("data-confidence") != "unresolved":
        assert t.get("data-construct"), t.attrib
bars = groups("figs/copy_number.svg", "cn-bar")
assert len(bars) == 4, len(bars)                 # 2 constructs x 2 runs
assert len({b.get("data-run") for b in bars}) == 2
cl = groups("figs/copy_number.svg", "closure")
assert cl and all(abs(float(c.get("data-deficit")) + 15) < 2 for c in cl), \
    [c.get("data-deficit") for c in cl]
rows = groups("figs/architecture.svg", "arch-row")
summ = json.load(open("site_summary.json"))["per_construct"]
assert {r.get("data-construct") for r in rows} == set(summ), rows
for r in rows:
    c = summ[r.get("data-construct")]
    assert int(r.get("data-ceiling")) == c["ceiling"]
    assert int(r.get("data-unresolved")) == c["unresolved_candidate_for"]
m = json.load(open("figs/figures.json"))
assert not m["skipped"], m["skipped"]
print("OK")
EOF
if grep -q "^OK$" figs_check.log; then
  ok "overview draws every .fai contig and one tick per called locus (8)"
  ok "copy number shows both runs and the shared-block deficit"
  ok "architecture carries each construct's ceiling and unresolved count"
else
  bad "figure content (full case)"; cat figs_check.log
fi

# no called loci: the overview still draws every chromosome
head -1 called.tsv > called_empty.tsv
"$PY" "$SCRIPTS/make_figures.py" --out figs_empty --called called_empty.tsv \
  --fai contigs.fai >figs_empty.log 2>&1 \
  && ok "overview with no called loci ran" || bad "overview with no called loci"
"$PY" - <<'EOF' >figs_empty_check.log 2>&1
import xml.etree.ElementTree as ET
NS = "{http://www.w3.org/2000/svg}"
g = list(ET.parse("figs_empty/genome_overview.svg").getroot().iter(NS + "g"))
n_fai = sum(1 for l in open("contigs.fai") if l.strip())
assert sum(e.get("class") == "chrom" for e in g) == n_fai
assert sum(e.get("class") == "locus" for e in g) == 0
print("OK")
EOF
grep -q "^OK$" figs_empty_check.log \
  && ok "no loci: all chromosomes drawn, zero ticks" \
  || { bad "no-loci overview"; cat figs_empty_check.log; }

# an input that exists but cannot be parsed is exit 2, not a skip
printf 'not json\n' >broken.json
expect_exit "make_figures: an unparsable site summary is exit 2" 2 figs_broken.log \
  "$PY" "$SCRIPTS/make_figures.py" --out figs_broken --site-summary broken.json
assert_has "the unreadable input is named" figs_broken.log "error: site summary broken.json unreadable"
expect_exit "make_figures: an unreadable --run manifest is exit 2" 2 figs_badrun.log \
  "$PY" "$SCRIPTS/make_figures.py" --out figs_badrun --run broken.json

# every input missing: nothing drawn, nothing fails, every skip explained
"$PY" "$SCRIPTS/make_figures.py" --out figs_none --fai does_not_exist.fai \
  >figs_none.log 2>&1 && ok "missing inputs degrade, exit 0" \
  || bad "missing inputs made the run fail"
"$PY" -c "
import json,sys; m=json.load(open('figs_none/figures.json'))
sys.exit(0 if not m['figures'] and len(m['skipped'])==4 else 1)" \
  && ok "each missing figure recorded as skipped with a reason" \
  || bad "skipped figures not recorded"

# the report step: placeholders filled from the manifest
cp "$HERE/../assets/report_template.md" report.md
sed -i 's/\[\[FIG:coverage_<construct A>\]\]/[[FIG:coverage_pEpsilon]]/; s/\[\[FIG:coverage_<construct B>\]\]/[[FIG:coverage_pDelta]]/' report.md
"$PY" "$SCRIPTS/make_figures.py" --out figs --identity identity_concat.json \
  --config cfg_recon.json --copy-number cn_recon.json \
  --called called.tsv --site-summary site_summary.json --fai contigs.fai \
  --fill report.md --fill-out report_filled.md >fill.log 2>&1 \
  && ok "--fill ran on the report template" || { bad "--fill"; cat fill.log; }
for f in coverage_s1_pDelta coverage_s1_pEpsilon genome_overview copy_number \
         architecture; do
  assert_has "filled report links $f" report_filled.md "\(figs/$f\.svg\)"
done
"$PY" -c "
import re; t=re.sub(r'<!--.*?-->', '', open('report_filled.md').read(), flags=re.S)
open('report_filled_body.md','w').write(t)"
assert_lacks "no figure slot left as a placeholder in the report body" \
  report_filled_body.md "\[\[FIG:"

# the template used unedited: the slot is explained in the report AND on stderr
"$PY" "$SCRIPTS/make_figures.py" --out figs --identity identity_concat.json \
  --config cfg_recon.json --copy-number cn_recon.json \
  --called called.tsv --site-summary site_summary.json --fai contigs.fai \
  --fill "$HERE/../assets/report_template.md" --fill-out report_raw.md \
  >fill_raw.out 2>fill_raw.err \
  && ok "--fill on an unedited template still exits 0" || bad "--fill unedited"
assert_has "unedited placeholder warned about on stderr, by name" fill_raw.err \
  "CHECK WARNING \[figure_placeholders\].*\[\[FIG:coverage_<construct A>\]\]"
assert_has "the warning lists the constructs that have coverage figures" \
  fill_raw.err "coverage figures exist for: pDelta, pEpsilon"
assert_lacks "an edited template draws no such warning" fill.log "WARNING"

echo
echo "== Final stage: do the called sites account for the DNA? =="
acct() { rs --loci pDelta=loci_recon_pDelta.tsv --loci pEpsilon=loci_recon_pEpsilon.tsv \
  --ploidy 2 --vector-ends pDelta=1200 --vector-ends pEpsilon=3000 "$@"; }
# 0.55 copies x ploidy 2 = 1.1 sites implied, 1 called: the ceiling rounds up
# to 2, and that rounded-up site is not missing DNA
acct --copies pDelta=0.55 --copies pEpsilon=3.5 --md acct_round.md \
  --summary acct_round.json >acct_round.log 2>&1 \
  && ok "accounting ran" || { bad "accounting"; cat acct_round.log; }
assert_has "a ceiling rounded up past the called sites is not a shortfall" \
  acct_round.md "pDelta \| 0\.55 \| 2 \| 1 \| 1 \| the called sites account"
assert_lacks "no shortfall caveat when the DNA is accounted for" acct_round.md \
  "SITE SHORTFALL"
acct_status() { "$PY" -c "
import json; d=json.load(open('$1'))
print([c['status'] for c in d['checks'] if c['id']=='site_accounting'][0])"; }
[ "$(acct_status acct_round.json)" = pass ] \
  && ok "site_accounting check passes" || bad "site_accounting should pass"
# 1.0 copies x 2 = 2 sites implied, 1 called: half the DNA has no site
acct --copies pDelta=1.0 --copies pEpsilon=3.5 --md acct_short.md \
  --summary acct_short.json >acct_short.log 2>&1
assert_has "one whole site unaccounted for is a shortfall, even at low copy" \
  acct_short.md "pDelta \| 1\.00 \| 2 \| 1 \| 1 \| SHORTFALL"
assert_has "the shortfall is on the construct it belongs to, not the other" \
  acct_short.md "pEpsilon \| 3\.50 \| 7 \| 7 \| 0 \| the called sites account"
assert_has "shortfall raised as a caveat naming the construct" acct_short.md \
  "SITE SHORTFALL \(pDelta: 1 of 2\)"
[ "$(acct_status acct_short.json)" = warn ] \
  && ok "site_accounting check warns in the JSON" || bad "site_accounting should warn"
assert_has "shortfall printed as a warning check" acct_short.log \
  "CHECK WARNING \[site_accounting\]"
# 4.6 x 2 = 9.2 implied, 7 called: within the ~25% of quantification.md s.7
acct --copies pDelta=0.55 --copies pEpsilon=4.6 --md acct_tol.md >acct_tol.log 2>&1
assert_has "agreement within ~25% accounts for the DNA" acct_tol.md \
  "pEpsilon \| 4\.60 \| 10 \| 7 \| 3 \| the called sites account"

echo
echo "== Stage 2: the locus table is a union over samples =="
# a site seen only in the second sample, and a 1-bp caller shift in another
cp junctions_nhej_s2.tsv union_s2.tsv
printf 'chr3\t40000000\tright_clip\t0\tpAlpha\t3300\tleft_clip\t50\tPASS\n' >>union_s2.tsv
awk 'BEGIN{OFS="\t"} $1=="chr7"{$2=$2+1} {print}' union_s2.tsv >union_s2_shift.tsv
for order in fwd rev; do
  if [ $order = fwd ]; then t1="parental=junctions_nhej.tsv"; t2="treated=union_s2_shift.tsv"
  else t1="treated=union_s2_shift.tsv"; t2="parental=junctions_nhej.tsv"; fi
  "$PY" "$SCRIPTS/classify_junctions.py" --identity id3.json --mapping pAlpha \
    --modality nhej --vector-ends 1,5300 --table "$t1" --table "$t2" \
    --out-loci "loci_union_$order.tsv" >"union_$order.log" 2>&1 \
    || { bad "union run ($order)"; cat "union_$order.log"; }
  assert_has "a site seen only in sample 2 is a locus ($order order)" \
    "loci_union_$order.tsv" "^chr3[[:space:]]+40000000"
done
"$PY" - <<'PYEOF' >union_check.log 2>&1
import csv
for order in ("fwd", "rev"):
    rows = list(csv.DictReader(open("loci_union_%s.tsv" % order), delimiter="\t"))
    chr1 = [r for r in rows if r["contig"] == "chr1"]
    assert len(chr1) == 1 and chr1[0]["flanks"] == "2", chr1
    # both flanks count in every sample: 30+26 and 18+15
    assert (chr1[0]["support_parental"], chr1[0]["support_treated"]) == ("56", "33"), chr1
    chr7 = [r for r in rows if r["contig"] == "chr7"]
    assert len(chr7) == 1, chr7
    assert (chr7[0]["support_parental"], chr7[0]["support_treated"]) == ("18", "10"), chr7
print("OK")
PYEOF
if grep -q "^OK$" union_check.log; then
  ok "every sample's support counts both flanks"
  ok "a 1-bp shift between samples does not split one site"
else
  bad "per-sample locus support"; cat union_check.log
fi

echo
echo "== Stage 4: argument handling =="
expect_exit "--min-copies-per-site 0 is refused, not a traceback" 2 rs_zero.log \
  rs --loci pDelta=loci_recon_pDelta.tsv --copies pDelta=1 --min-copies-per-site 0
assert_lacks "no traceback for a zero divisor" rs_zero.log "Traceback"
# one --copies under --run overrides only the construct it names
"$PY" "$SCRIPTS/reconcile_sites.py" --run mf/run.json --min-copies-per-site 10 \
  --copies pDelta=0.5 --summary mf_copies.json --md mf_copies.md \
  --out-called mf_copies_called.tsv --out-subthreshold mf_copies_sub.tsv \
  >mf_copies.log 2>&1
"$PY" - <<'PYEOF' >mf_copies_check.log 2>&1
import json
d = json.load(open("mf_copies.json"))
assert "overridden from --copies for pDelta" in d["copy_number_source"], d["copy_number_source"]
assert d["per_construct"]["pDelta"]["copies_per_haploid"] == 0.5, d["per_construct"]
assert d["per_construct"]["pEpsilon"]["copies_per_haploid"] is not None, d["per_construct"]
assert d["per_construct"]["pEpsilon"]["ceiling"] is not None, d["per_construct"]
print("OK")
PYEOF
grep -q "^OK$" mf_copies_check.log \
  && ok "--copies under --run keeps the manifest's copy numbers for other constructs" \
  || { bad "--copies under --run"; cat mf_copies_check.log mf_copies.log; }

echo
echo "== Host stage: test_host.sh =="
bash "$HERE/test_host.sh" >host_suite.log 2>&1; host_rc=$?
grep -E "FAIL" host_suite.log
hp="$(sed -n 's/.*passed: \([0-9]*\) .*/\1/p' host_suite.log | tail -1)"
hf="$(sed -n 's/.*failed: \([0-9]*\).*/\1/p' host_suite.log | tail -1)"
if [ -n "$hp" ] && [ -n "$hf" ]; then
  PASS=$((PASS+hp)); FAIL=$((FAIL+hf))
  echo "  test_host.sh: $hp passed, $hf failed (counted below)"
else
  bad "test_host.sh did not report a result (exit $host_rc)"; tail -5 host_suite.log
fi

echo
echo "=============================="
echo " passed: $PASS   failed: $FAIL"
echo "=============================="
[ "$FAIL" -eq 0 ]
