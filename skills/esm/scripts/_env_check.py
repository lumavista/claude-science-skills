"""Preflight for local ESM inference: report what will break before it breaks.

Run `python scripts/_env_check.py` before a first local run. Every check prints
PASS, WARN or FAIL with the remedy attached; the script exits 1 if any check
FAILs, so it is usable as a gate in a shell script.

The split between WARN and FAIL follows `_device.py`'s contract: CPU always
produces a correct result and the GPU is an accelerator, so a CPU-only torch
wheel is a WARN (slow, not wrong). FAIL is reserved for states where loading a
model raises -- a torch version the SDK rejects, or the unset-USERNAME defect
with no cache directory set to work around it.

`--prime` additionally downloads the ESMC checkpoints with `max_workers=1`.
That is the documented way around the Windows symlink race (WinError 1314), and
it is opt-in because it pulls several GiB.
"""

from __future__ import annotations

import argparse
import os
import pathlib
import sys
import tempfile

# The SDK declares torch>=2.11.0,<2.12.0, so 2.11.x is the only permitted minor.
TORCH_MIN = (2, 11)
TORCH_MAX_EXCLUSIVE = (2, 12)
ESMC_MODELS = ("biohub/ESMC-300M", "biohub/ESMC-600M")

PASS, WARN, FAIL = "PASS", "WARN", "FAIL"


def _check_python() -> tuple[str, str, str]:
    v = sys.version_info
    got = f"{v.major}.{v.minor}.{v.micro}"
    if (v.major, v.minor) < (3, 12):
        return (FAIL, "python", f"{got}; esm requires >=3.12. Build a 3.12 venv for this skill.")
    return (PASS, "python", got)


def _check_torch() -> list[tuple[str, str, str]]:
    try:
        import torch
    except ImportError:
        return [(FAIL, "torch", "not installed. pip install torch==2.11.0 "
                                "--index-url https://download.pytorch.org/whl/cu128")]

    results = []
    version = torch.__version__.split("+")[0]
    try:
        parts = tuple(int(x) for x in version.split(".")[:2])
    except ValueError:
        parts = ()
    if parts and not (TORCH_MIN <= parts < TORCH_MAX_EXCLUSIVE):
        results.append((FAIL, "torch.version",
                        f"{torch.__version__}; esm declares >=2.11.0,<2.12.0. "
                        "pip install torch==2.11.0 "
                        "--index-url https://download.pytorch.org/whl/cu128"))
    else:
        results.append((PASS, "torch.version", torch.__version__))

    # The decisive test for the silent CPU-only wheel: a CPU build reports None
    # here regardless of what hardware is present, so this distinguishes "no
    # CUDA wheel" from "no GPU" -- which is.is_available() cannot.
    if torch.version.cuda is None:
        results.append((WARN, "torch.cuda.build",
                        "CPU-only wheel (torch.version.cuda is None). Inference "
                        "will run on CPU. For GPU: pip install torch==2.11.0 "
                        "--index-url https://download.pytorch.org/whl/cu128 "
                        "(cu124 caps at 2.6.0 and cannot satisfy esm)"))
    else:
        results.append((PASS, "torch.cuda.build", f"CUDA {torch.version.cuda}"))
        try:
            if torch.cuda.is_available():
                results.append((PASS, "torch.cuda.runtime", torch.cuda.get_device_name(0)))
            else:
                results.append((WARN, "torch.cuda.runtime",
                                "CUDA build present but no device visible; check the driver"))
        except Exception as exc:  # a broken driver must report, not crash the preflight
            results.append((WARN, "torch.cuda.runtime", f"probe failed: {exc}"))
    return results


def _hf_cache_dir() -> pathlib.Path:
    try:
        from huggingface_hub import constants
        return pathlib.Path(constants.HF_HUB_CACHE)
    except Exception:
        return pathlib.Path(os.environ.get("HF_HOME", pathlib.Path.home() / ".cache" / "huggingface")) / "hub"


def _check_hf_cache() -> list[tuple[str, str, str]]:
    cache = _hf_cache_dir()
    results = []
    try:
        cache.mkdir(parents=True, exist_ok=True)
        probe = cache / ".preflight-write-test"
        probe.write_text("ok", encoding="utf-8")
        probe.unlink()
        results.append((PASS, "hf.cache.writable", str(cache)))
    except Exception as exc:
        results.append((FAIL, "hf.cache.writable",
                        f"{cache} is not writable ({exc}). Set HF_HOME to a writable path."))
        return results

    # Symlink creation, not a capability flag: huggingface_hub memoizes support
    # per-directory and a parallel fetch races os.symlink, raising WinError 1314
    # rather than the PermissionError the hub catches. HF_HUB_DISABLE_SYMLINKS
    # does not exist in the installed hub (0.36.2), so priming serially is the
    # remedy -- not an environment variable.
    with tempfile.TemporaryDirectory(dir=cache) as tmp:
        target = pathlib.Path(tmp) / "target.txt"
        target.write_text("x", encoding="utf-8")
        try:
            os.symlink(target, pathlib.Path(tmp) / "link.txt")
            results.append((PASS, "hf.cache.symlinks", "supported"))
        except OSError as exc:
            winerror = getattr(exc, "winerror", None)
            results.append((WARN, "hf.cache.symlinks",
                            f"os.symlink fails (winerror {winerror}). Parallel model "
                            "downloads can crash. Remedy: run this script with --prime "
                            "to fetch serially, or enable Windows Developer Mode."))
    return results


def _check_sandbox_env() -> list[tuple[str, str, str]]:
    # One defect, two symptoms: USERNAME unset breaks getpass.getuser(), and
    # torch._inductor then raises ModuleNotFoundError 'pwd' during an ordinary
    # eager from_pretrained. Setting TORCHINDUCTOR_CACHE_DIR avoids that path.
    username = os.environ.get("USERNAME") or os.environ.get("USER")
    inductor = os.environ.get("TORCHINDUCTOR_CACHE_DIR")
    results = []
    if username:
        results.append((PASS, "env.USERNAME", username))
    elif inductor:
        results.append((PASS, "env.USERNAME",
                        "unset, but TORCHINDUCTOR_CACHE_DIR is set, which avoids the 'pwd' crash"))
    else:
        results.append((FAIL, "env.USERNAME",
                        "unset with no TORCHINDUCTOR_CACHE_DIR; torch._inductor will raise "
                        "ModuleNotFoundError 'pwd' on from_pretrained. Set either variable."))
    results.append((PASS, "env.TORCHINDUCTOR_CACHE_DIR", inductor) if inductor
                   else (WARN, "env.TORCHINDUCTOR_CACHE_DIR",
                         "unset; set it to a writable path if model loading raises 'pwd'"))
    return results


def _check_models(prime: bool) -> list[tuple[str, str, str]]:
    try:
        from huggingface_hub import snapshot_download
    except ImportError:
        return [(WARN, "hf.hub", "huggingface_hub not installed; cannot inspect model cache")]

    results = []
    for model_id in ESMC_MODELS:
        if not prime:
            # local_files_only resolves from cache without touching the network.
            try:
                path = snapshot_download(model_id, local_files_only=True)
                results.append((PASS, f"model.{model_id}", f"cached at {path}"))
            except Exception:
                results.append((WARN, f"model.{model_id}",
                                "not cached; re-run with --prime to fetch serially"))
            continue
        try:
            path = snapshot_download(model_id, max_workers=1)
            results.append((PASS, f"model.{model_id}", f"primed at {path}"))
        except Exception as exc:
            results.append((FAIL, f"model.{model_id}", f"serial download failed: {exc}"))
    return results


def validate_env(prime: bool = False) -> list[tuple[str, str, str]]:
    """Run every check and return (status, name, detail) in report order."""
    results = [_check_python()]
    results += _check_torch()
    results += _check_hf_cache()
    results += _check_sandbox_env()
    results += _check_models(prime)
    return results


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--prime", action="store_true",
                        help="download ESMC-300M and ESMC-600M serially (several GiB)")
    args = parser.parse_args()

    results = validate_env(prime=args.prime)
    width = max(len(name) for _, name, _ in results)
    for status, name, detail in results:
        print(f"[{status}] {name.ljust(width)}  {detail}")

    failed = sum(1 for status, _, _ in results if status == FAIL)
    warned = sum(1 for status, _, _ in results if status == WARN)
    print(f"\n{len(results)} checks: {len(results) - failed - warned} pass, {warned} warn, {failed} fail")
    if failed:
        print("Environment is not ready for local inference; fix the FAIL lines above.")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
