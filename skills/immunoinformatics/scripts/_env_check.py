"""Preflight for immunoinformatics prediction: report what will break first.

Run `python scripts/_env_check.py` before a first run. Every check prints PASS,
WARN or FAIL with the remedy attached; the script exits 1 if any check FAILs.

This skill is CPU work -- prediction and triage are Tier 1 in
`references/compute-routes.md`, and the GPU tier is folding, embedding and MD.
So there is deliberately no torch or CUDA check here: a CPU-only torch wheel is
not a defect for this skill, and failing on one would reject a working install.

MHCflurry does not use the Hugging Face hub. Its weights come from GitHub and
Dropbox via `mhcflurry-downloads fetch`, into the directory named by
MHCFLURRY_MODELS_DIR, so the cache checks below are about that path and not
about HF_HOME.
"""

from __future__ import annotations

import os
import pathlib
import sys

PASS, WARN, FAIL = "PASS", "WARN", "FAIL"

FETCH_HINT = ("bash skills/immunoinformatics/setup_models.sh "
              "(or: mhcflurry-downloads fetch models_class1_presentation)")


def _check_python() -> tuple[str, str, str]:
    v = sys.version_info
    got = f"{v.major}.{v.minor}.{v.micro}"
    if (v.major, v.minor) < (3, 11):
        return (FAIL, "python", f"{got}; this skill requires >=3.11.")
    return (PASS, "python", got)


def _check_mhcflurry() -> list[tuple[str, str, str]]:
    try:
        import mhcflurry
    except ImportError:
        return [(FAIL, "mhcflurry", 'not installed. pip install "mhcflurry>=2.3"')]

    version = getattr(mhcflurry, "__version__", "unknown")
    try:
        parts = tuple(int(x) for x in version.split(".")[:2])
    except ValueError:
        parts = ()
    if parts and parts < (2, 3):
        return [(FAIL, "mhcflurry", f'{version}; requires >=2.3. pip install -U "mhcflurry>=2.3"')]
    return [(PASS, "mhcflurry", version)]


def _models_dir() -> pathlib.Path | None:
    """Effective MHCflurry presentation-models directory, or None if unknown."""
    override = os.environ.get("MHCFLURRY_MODELS_DIR")
    if override:
        return pathlib.Path(override)
    try:
        from mhcflurry import downloads
        return pathlib.Path(downloads.get_default_class1_presentation_models_dir(test_exists=False))
    except Exception:
        return None


def _check_models_dir() -> list[tuple[str, str, str]]:
    results = []
    override = os.environ.get("MHCFLURRY_MODELS_DIR")
    results.append((PASS, "env.MHCFLURRY_MODELS_DIR", override) if override
                   else (WARN, "env.MHCFLURRY_MODELS_DIR",
                         "unset; using the per-user default below. Set it to keep "
                         "weights with the skill."))

    path = _models_dir()
    if path is None:
        results.append((WARN, "mhcflurry.models.dir", "could not resolve; is mhcflurry installed?"))
        return results
    results.append((PASS, "mhcflurry.models.dir", str(path)))

    # Writability is checked on the nearest existing ancestor, because the
    # models directory legitimately does not exist until the first fetch.
    probe_dir = path
    while not probe_dir.exists() and probe_dir != probe_dir.parent:
        probe_dir = probe_dir.parent
    try:
        probe = probe_dir / ".preflight-write-test"
        probe.write_text("ok", encoding="utf-8")
        probe.unlink()
        results.append((PASS, "mhcflurry.models.writable", str(probe_dir)))
    except Exception as exc:
        results.append((FAIL, "mhcflurry.models.writable",
                        f"{probe_dir} is not writable ({exc}). Set MHCFLURRY_MODELS_DIR "
                        "to a writable path."))
        return results

    if path.exists() and any(path.iterdir()):
        results.append((PASS, "mhcflurry.models.present", "weights found"))
    else:
        results.append((WARN, "mhcflurry.models.present",
                        f"no weights at {path}. Fetch them with: {FETCH_HINT}"))
    return results


def _check_predictor_loads() -> list[tuple[str, str, str]]:
    """Load the predictor for real -- the only check that proves prediction works."""
    try:
        from mhcflurry import Class1PresentationPredictor
    except ImportError:
        return [(WARN, "mhcflurry.predictor", "Class1PresentationPredictor unavailable")]
    try:
        Class1PresentationPredictor.load()
        return [(PASS, "mhcflurry.predictor", "loaded")]
    except Exception as exc:
        return [(WARN, "mhcflurry.predictor",
                 f"did not load ({type(exc).__name__}). The skill falls back to NetMHCpan "
                 f"if on PATH, else fails loudly. Fetch weights with: {FETCH_HINT}")]


def validate_env() -> list[tuple[str, str, str]]:
    """Run every check and return (status, name, detail) in report order."""
    results = [_check_python()]
    results += _check_mhcflurry()
    results += _check_models_dir()
    results += _check_predictor_loads()
    return results


def main() -> int:
    results = validate_env()
    width = max(len(name) for _, name, _ in results)
    for status, name, detail in results:
        print(f"[{status}] {name.ljust(width)}  {detail}")

    failed = sum(1 for status, _, _ in results if status == FAIL)
    warned = sum(1 for status, _, _ in results if status == WARN)
    print(f"\n{len(results)} checks: {len(results) - failed - warned} pass, {warned} warn, {failed} fail")
    if failed:
        print("Environment is not ready; fix the FAIL lines above.")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
