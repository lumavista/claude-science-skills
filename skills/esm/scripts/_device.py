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

    Returns total, not free, memory: it is the stable number to plan against.
    Free memory is the less useful number than it looks on Windows, where
    `mem_get_info` reports 0 free under pressure while allocations keep
    succeeding -- see `fits_in_vram` for why that matters.
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
                 overhead: float = 1.35, reserve_gb: float = 1.1,
                 device: str | None = None) -> bool | None:
    """Rough check that a model of this size fits. None when memory is unknown.

    `overhead` covers activations, workspace and allocator fragmentation. This
    is a planning heuristic for picking a checkpoint, not an allocation
    guarantee -- a long sequence can still exceed it through activation growth.

    On Windows this heuristic is the ONLY guard available, because the usual
    backstop does not exist. Measured on an 8 GiB consumer GPU (8 GiB) in 2026-10:
    allocating past the card never raised `torch.cuda.OutOfMemoryError` -- an oversized allocation was accepted while `mem_get_info` reported 0 free, because WDDM pages
    GPU allocations into host RAM instead of failing. What breaks is throughput,
    silently: a 4096-cube fp32 matmul held full speed within budget,
    then collapsed once the budget was exceeded, an order-of-magnitude
    cliff with no exception to catch. So returning False here is a reason not to
    proceed, and `try/except OutOfMemoryError` is not a usable fallback on this
    platform.

    `reserve_gb` is subtracted from total memory because a display-attached card
    does not offer all of it: the desktop compositor held ~1.1 GiB of the 8 GiB measured above, leaving ~6 GiB usable. Comparing against total would
    approve a model that then lands in the cliff. Raise it when a browser or a
    second process shares the card; set it to 0 only on a headless device.
    """
    total = vram_gb(device)
    if total is None:
        return None
    budget = total - max(reserve_gb, 0.0)
    if budget <= 0:
        return False
    return params_billions * bytes_per_param * overhead <= budget


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
