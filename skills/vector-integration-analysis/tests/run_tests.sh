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
echo "=============================="
echo " passed: $PASS   failed: $FAIL"
echo "=============================="
[ "$FAIL" -eq 0 ]
