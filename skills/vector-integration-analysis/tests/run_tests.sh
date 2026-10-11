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
echo "=============================="
echo " passed: $PASS   failed: $FAIL"
echo "=============================="
[ "$FAIL" -eq 0 ]
