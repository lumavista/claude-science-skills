"""Device selection for torch workloads: GPU when present, CPU as fallback.

This is the canonical copy of a small snippet duplicated verbatim into every
skill in this collection that drives a torch model. Skills are imported
individually, so they cannot share a module; the duplication is deliberate.
Keep copies in sync -- `grep -rl "SCI_DEVICE" skills/` lists them.

Contract: CPU always produces a correct result. The GPU is an accelerator, so
nothing here raises when torch is missing or no GPU is visible -- it returns
"cpu" and the caller proceeds. An explicit SCI_DEVICE is honored verbatim and
is NOT validated against availability, so a forced device fails loudly at the
first tensor op rather than silently running somewhere the caller did not ask
for.
"""

from __future__ import annotations

import os

__all__ = ["select_device", "device_report", "vram_gb", "autocast_dtype", "fits_in_vram"]


def select_device(prefer: str | None = None) -> str:
    """Return the fastest available torch device string; "cpu" if unsure.

    Priority: explicit argument > SCI_DEVICE env var > CUDA > MPS > CPU.
    """
    requested = (prefer or os.environ.get("SCI_DEVICE") or "").strip().lower()
    if requested:
        return requested
    try:
        import torch
    except ImportError:
        return "cpu"
    try:
        if torch.cuda.is_available():
            return "cuda"
        mps = getattr(torch.backends, "mps", None)
        if mps is not None and mps.is_available():
            return "mps"
    except Exception:
        # A broken driver or a CUDA build mismatch must degrade, not crash.
        return "cpu"
    return "cpu"


def vram_gb(device: str | None = None) -> float | None:
    """Total memory of the selected CUDA device in GiB, or None if not CUDA.

    Returns total, not free, memory: it is the stable number to plan a model
    against. Check free memory at allocation time if a run must share the card.
    """
    device = device or select_device()
    if not device.startswith("cuda"):
        return None
    try:
        import torch

        index = int(device.split(":", 1)[1]) if ":" in device else torch.cuda.current_device()
        return torch.cuda.get_device_properties(index).total_memory / (1024**3)
    except Exception:
        return None


def autocast_dtype(device: str | None = None):
    """Preferred reduced precision for inference, or None to stay in fp32.

    bfloat16 on Ampere and newer (compute capability >= 8.0) because it keeps
    fp32's exponent range and so does not need loss scaling; float16 on older
    CUDA cards; fp32 everywhere else, including CPU and MPS, where reduced
    precision is usually slower or unsupported rather than faster.
    """
    device = device or select_device()
    if not device.startswith("cuda"):
        return None
    try:
        import torch

        major = torch.cuda.get_device_capability()[0]
        return torch.bfloat16 if major >= 8 else torch.float16
    except Exception:
        return None


def fits_in_vram(params_billions: float, *, bytes_per_param: int = 2,
                 overhead: float = 1.35, device: str | None = None) -> bool | None:
    """Rough check that a model of this size fits. None when memory is unknown.

    `overhead` covers activations, workspace and allocator fragmentation. This
    is a planning heuristic for picking a checkpoint, not an allocation
    guarantee -- a long sequence can still exceed it through activation growth.
    """
    total = vram_gb(device)
    if total is None:
        return None
    return params_billions * bytes_per_param * overhead <= total


def device_report(device: str | None = None) -> str:
    """One-line human-readable summary for logs and provenance sidecars."""
    device = device or select_device()
    if device.startswith("cuda"):
        try:
            import torch

            index = int(device.split(":", 1)[1]) if ":" in device else torch.cuda.current_device()
            name = torch.cuda.get_device_name(index)
            total = vram_gb(device)
            return f"{device} ({name}, {total:.1f} GiB)" if total else f"{device} ({name})"
        except Exception:
            return device
    if device == "cpu":
        return f"cpu ({os.cpu_count()} logical cores)"
    return device


if __name__ == "__main__":
    print(device_report())
