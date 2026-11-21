"""The run manifest and the stage checks.

Two things every stage shares, kept in one place so no stage re-implements
them and no caller has to carry values between stages by hand.

**The run manifest** (`--run run.json`) states a run's facts once -- genome
size, read length, ploidy, modality, vector ends, the samples and their input
files -- and collects what each stage produced under `stages`. A stage given
`--run` reads its defaults from the manifest and records its outputs there, so
the next stage finds the identity map, the junction and loci tables, the
detected-locus count and the copy-number JSON without anyone retyping a path
or a number. An explicit command-line flag always wins over the manifest, and
a stage run without `--run` behaves exactly as it always has.

Paths inside a manifest are relative to the manifest's own directory. Keys
starting with `_` are comments. Stages never rewrite the facts a person
declared; they only write under `stages`.

**Checks** are the conditions under which an analysis should stop, raised by
the script that can see them instead of by a paragraph somebody has to
remember. Every stage emits the same list -- `{id, status, message}` with
status `pass`, `warn`, `fail` or `accepted` -- in its JSON output and in the
manifest, prints every one that is not a plain pass on stderr, and with
`--strict` exits with status 3 when any check fails. `--accept ID` acknowledges one failing check
by name: it stays in the output, marked `accepted`, and stops blocking.

Standard library only.
"""
from __future__ import annotations

import json
import os
import sys
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

#: exit status of a stage stopped by `--strict` on a failing check
STRICT_EXIT = 3

PASS, WARN, FAIL, ACCEPTED = "pass", "warn", "fail", "accepted"


# --------------------------------------------------------------------------
# checks
# --------------------------------------------------------------------------

class Checks:
    """The stop conditions one stage evaluated, in the order it evaluated them."""

    def __init__(self, stage: str, accept: Iterable[str] = ()):
        self.stage = stage
        self.items: List[Dict[str, str]] = []
        self._accept = {a.strip() for a in accept if a and a.strip()}

    def add(self, check_id: str, status: str, message: str) -> None:
        if status == FAIL and check_id in self._accept:
            status = ACCEPTED
        self.items.append({"id": check_id, "status": status,
                           "message": message})

    def failed(self) -> List[Dict[str, str]]:
        return [c for c in self.items if c["status"] == FAIL]

    def flagged(self) -> List[Dict[str, str]]:
        """Failing checks plus the ones a caller explicitly accepted."""
        return [c for c in self.items if c["status"] in (FAIL, ACCEPTED)]

    def unknown_accepts(self) -> List[str]:
        return sorted(self._accept - {c["id"] for c in self.items})

    def as_list(self) -> List[Dict[str, str]]:
        return [dict(c) for c in self.items]

    def banner_md(self) -> List[str]:
        """Markdown lines for the top of a report; empty when nothing failed."""
        flagged = self.flagged()
        if not flagged:
            return []
        L = ["## Stop: resolve these before using anything below", ""]
        for c in flagged:
            L.append("> **%s `%s`** -- %s"
                     % ("ACCEPTED" if c["status"] == ACCEPTED else "FAILED",
                        c["id"], c["message"]))
        return L + [""]

    def finish(self, strict: bool, stream=None) -> int:
        """Print the failing checks; return the stage's exit status.

        0 unless `strict` and at least one check is failing.
        """
        stream = stream or sys.stderr
        for name in self.unknown_accepts():
            print("  warning: --accept %s names no check of this stage (%s)"
                  % (name, ", ".join(c["id"] for c in self.items)),
                  file=stream)
        label = {FAIL: "CHECK FAILED", ACCEPTED: "CHECK ACCEPTED",
                 WARN: "CHECK WARNING"}
        for c in self.items:
            if c["status"] in label:
                print("  %s [%s]: %s"
                      % (label[c["status"]], c["id"], c["message"]),
                      file=stream)
        failed = self.failed()
        if strict and failed:
            print("STOP (%s): %d blocking check(s): %s. Exit status %d. Fix the "
                  "cause and rerun; or, having decided a check does not apply, "
                  "rerun with --accept ID."
                  % (self.stage, len(failed),
                     ", ".join(c["id"] for c in failed), STRICT_EXIT),
                  file=stream)
            return STRICT_EXIT
        return 0


def add_check_arguments(ap) -> None:
    """The three flags every stage takes for the manifest and the checks."""
    ap.add_argument("--run", default="", metavar="RUN.json",
                    help="run manifest: read this stage's defaults from it and "
                         "record this stage's outputs in it. Explicit flags "
                         "win. See assets/run.example.json.")
    ap.add_argument("--strict", action="store_true",
                    help="exit with status %d when any check fails, instead "
                         "of 0 with a warning" % STRICT_EXIT)
    ap.add_argument("--accept", action="append", default=[], metavar="CHECK_ID",
                    help="acknowledge one failing check by its id so it no "
                         "longer blocks under --strict; it stays in the "
                         "output marked `accepted`. Repeatable.")


# --------------------------------------------------------------------------
# the manifest
# --------------------------------------------------------------------------

class RunError(ValueError):
    """The manifest cannot supply what a stage asked it for."""


class Run:
    """A run manifest: declared facts, plus what each stage recorded."""

    def __init__(self, path: str):
        self.file = os.path.abspath(path)
        self.base = os.path.dirname(self.file)
        try:
            with open(self.file, encoding="utf-8") as fh:
                self.data: Dict = json.load(fh)
        except OSError as exc:
            raise RunError(
                "cannot read run manifest %s (%s). Start one by copying "
                "assets/run.example.json and filling in the facts."
                % (path, exc.strerror or exc))
        except ValueError as exc:
            raise RunError("run manifest %s is not valid JSON: %s" % (path, exc))
        if not isinstance(self.data, dict):
            raise RunError("run manifest %s must be a JSON object" % path)

    # ---- paths
    def resolve(self, p: str) -> str:
        return p if os.path.isabs(p) else os.path.join(self.base, p)

    def relative(self, p: str) -> str:
        return os.path.relpath(os.path.abspath(p), self.base)

    # ---- facts
    def fact(self, key: str, default=None):
        return self.data.get(key, default)

    def stage(self, name: str) -> Dict:
        return (self.data.get("stages") or {}).get(name) or {}

    def record(self, stage: str, data: Dict, key: Optional[str] = None) -> None:
        """Store a stage's outputs under `stages` and rewrite the manifest.

        With `key`, the stage holds one record per key (Stage 2 runs once per
        mapping). Declared facts are left exactly as they were.
        """
        stages = self.data.setdefault("stages", {})
        if key is None:
            stages[stage] = data
        else:
            stages.setdefault(stage, {})[key] = data
        tmp = self.file + ".tmp"
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(self.data, fh, indent=2)
            fh.write("\n")
        os.replace(tmp, self.file)

    # ---- what earlier stages produced
    def identity_path(self) -> Optional[str]:
        p = self.stage("identity").get("path") or self.fact("identity")
        return self.resolve(p) if p else None

    def constructs(self) -> List[str]:
        return list(self.stage("identity").get("constructs") or {})

    def vector_ends(self, mapping: str) -> List[int]:
        ends = (self.fact("vector_ends") or {}).get(mapping) or []
        try:
            return [int(x) for x in ends]
        except (TypeError, ValueError):
            raise RunError("vector_ends.%s in the run manifest must be a list "
                           "of integers, got %r" % (mapping, ends))

    def tables(self, mapping: str) -> List[Tuple[str, str]]:
        """(sample, junction table) for one mapping, in manifest order."""
        out: List[Tuple[str, str]] = []
        for sample, s in (self.fact("samples") or {}).items():
            p = (s.get("junctions") or {}).get(mapping)
            if p:
                out.append((sample, self.resolve(p)))
        return out

    def norms(self) -> Dict[str, float]:
        out: Dict[str, float] = {}
        for sample, s in (self.fact("samples") or {}).items():
            if s.get("norm") is not None:
                out[sample] = float(s["norm"])
        return out

    def output(self, explicit: Optional[str], default_name: str) -> str:
        """An output path: the explicit flag, else a named file beside the
        manifest."""
        return explicit or os.path.join(self.base, default_name)


def load_run(path: str) -> Optional[Run]:
    return Run(path) if path else None


# --------------------------------------------------------------------------
# Stage 3 and Stage 4 inputs, derived from the manifest
# --------------------------------------------------------------------------

#: per-sample manifest keys that belong to Stage 2, not to the copy-number config
_STAGE2_SAMPLE_KEYS = ("junctions", "norm")


def copy_number_config(run: Run) -> Tuple[Dict, List[str]]:
    """The copy_number.py config a careful person would have written by hand.

    Everything Stage 1 and Stage 2 already produced is taken from their
    records: the identity map, each mapping's classified-junction table, and
    the detected-locus count of the primary mapping. Returns the config and
    the names of samples left out because they declare no `dedup_pairs`.

    The constructs solved for are, in identity-map order, those that have a
    Stage 2 record, a depth file in any sample, or an entry in the manifest's
    own optional `mappings` block. A construct with none of those has no data
    to estimate it from, and naming it would only make the solve singular.
    """
    ident = run.identity_path()
    if not ident:
        raise RunError("the run manifest has no identity map yet: run "
                       "construct_identity.py --run first")
    cfg: Dict = {k: v for k, v in run.data.items()
                 if not k.startswith("_") and k not in ("stages", "vector_ends")}
    cfg["identity"] = ident
    if "genome_size" not in cfg and run.stage("host").get("genome_size"):
        # the host stage measured it from the assembly; a declared value wins
        cfg["genome_size"] = run.stage("host")["genome_size"]

    samples: Dict[str, Dict] = {}
    skipped: List[str] = []
    for name, s in (run.fact("samples") or {}).items():
        if "dedup_pairs" not in s:
            skipped.append(name)
            continue
        samples[name] = {k: v for k, v in s.items()
                         if k not in _STAGE2_SAMPLE_KEYS}
    cfg["samples"] = samples

    classified = run.stage("classify")
    declared = run.fact("mappings") or {}
    with_depth = {m for s in samples.values()
                  for m in (s.get("depth_files") or {})}
    known = run.constructs() or list(declared) or list(classified)
    names = [n for n in known
             if n in classified or n in with_depth or n in declared]
    primary = run.fact("primary_mapping") or (names[0] if names else None)
    if primary:
        cfg["primary_mapping"] = primary
    mappings: Dict[str, Dict] = {}
    for n in names:
        m = dict(declared.get(n) or {})
        rec = classified.get(n) or {}
        if rec.get("junctions"):
            m.setdefault("junctions", rec["junctions"])
        if n == primary and rec.get("detected_loci") is not None:
            m.setdefault("detected_loci", rec["detected_loci"])
        mappings[n] = m
    cfg["mappings"] = mappings
    return cfg, skipped


def reconcile_inputs(run: Run) -> Dict:
    """What reconcile_sites.py needs from the earlier stages' records."""
    classified = run.stage("classify")
    cn = run.stage("copy_number")
    return {
        "identity": run.identity_path(),
        "loci": {m: run.resolve(rec["loci"])
                 for m, rec in classified.items() if rec.get("loci")},
        "vector_ends": {m: run.vector_ends(m) for m in classified
                        if run.vector_ends(m)},
        "copy_number": run.resolve(cn["json"]) if cn.get("json") else None,
        "ploidy": run.fact("ploidy"),
        "copies_per_site": cn.get("copies_per_site"),
        "modality": run.fact("modality"),
    }


def worst(statuses: Sequence[str]) -> str:
    for s in (FAIL, WARN):
        if s in statuses:
            return s
    return PASS
