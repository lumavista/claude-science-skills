#!/usr/bin/env python3
"""
Hybrid MHC prediction with fallback routing (no simulations).

Route MHC binding prediction based on available resources:
1. Local MHCflurry (if model can be downloaded)
2. IEDB API for queried alleles
3. NetMHCpan server (if installed locally)
4. Fail with clear instructions (no synthetic predictions)

Usage:
    from predict_mhc_with_fallback import predict_mhc_class1

    results = predict_mhc_class1(
        peptides=["FVFSTSFYL", "YVFSTSFYL"],
        alleles=["HLA-A*02:01", "HLA-A*11:01"],
        include_percentile=True,
    )
    for row in results:
        print(f"{row['peptide']} → {row['allele']}: {row['percentile']:.3f}%")
"""

import os
import sys
import json
import platform
import subprocess
import tempfile
from pathlib import Path
from typing import List, Dict, Optional, Any


class MHCPredictionError(Exception):
    """Raised when MHC prediction cannot proceed in current environment."""
    pass


def detect_mhc_environment() -> Dict[str, Any]:
    """Detect available MHC prediction backends."""
    info = {
        "platform": platform.system(),
        "can_use_mhcflurry": False,
        "can_use_netmhcpan": False,
        "can_use_iedb_api": True,  # Always available if internet
        "mhcflurry_version": None,
        "netmhcpan_path": None,
    }

    # Check local MHCflurry
    try:
        import mhcflurry
        info["can_use_mhcflurry"] = True
        info["mhcflurry_version"] = mhcflurry.__version__
        # Try to load the model to see if downloads work
        try:
            from mhcflurry import Class1PresentationPredictor
            Class1PresentationPredictor.load()
            info["mhcflurry_loads"] = True
        except Exception as e:
            info["mhcflurry_loads"] = False
            info["mhcflurry_load_error"] = str(e)
    except ImportError:
        pass

    # Check NetMHCpan
    import shutil
    netmhcpan_exe = shutil.which("netMHCpan")
    if netmhcpan_exe:
        info["can_use_netmhcpan"] = True
        info["netmhcpan_path"] = netmhcpan_exe

    return info


def _setup_mhcflurry_models() -> None:
    """Initialize MHCflurry model cache, preferring skill-bundled models."""
    model_dir = None

    # Priority order: explicit env var > skill bundle > temp > home > cwd
    skill_root = Path(__file__).parent.parent  # skills/immunoinformatics/
    candidates = [
        os.environ.get("MHCFLURRY_MODELS_DIR"),  # Explicit override
        str(skill_root / "models" / "mhcflurry"),  # Bundled with skill
        os.path.join(tempfile.gettempdir(), "mhcflurry_models"),  # Temp cache
        os.path.expanduser("~/.mhcflurry/models"),  # User home
        os.path.join(os.getcwd(), ".mhcflurry_models"),  # CWD
    ]

    for candidate in candidates:
        if not candidate:
            continue
        try:
            os.makedirs(candidate, exist_ok=True)
            model_dir = candidate
            print(f"[predict_mhc_with_fallback] Using model cache: {model_dir}", file=sys.stderr)
            break
        except (OSError, PermissionError):
            continue

    if model_dir and not os.environ.get("MHCFLURRY_MODELS_DIR"):
        os.environ["MHCFLURRY_MODELS_DIR"] = model_dir


def predict_with_mhcflurry(
    peptides: List[str],
    alleles: List[str],
    include_percentile: bool = True,
) -> List[Dict[str, Any]]:
    """
    Predict MHC Class I binding using local MHCflurry.

    Args:
        peptides: List of peptide sequences
        alleles: List of HLA alleles (e.g., ["HLA-A*02:01"])
        include_percentile: If True, compute percentile rank

    Returns:
        List of prediction dicts with keys:
        peptide, allele, affinity_nM, percentile_rank, presentation_score
    """
    try:
        from mhcflurry import Class1PresentationPredictor
    except ImportError as e:
        raise MHCPredictionError(
            f"mhcflurry not installed: {e}\n"
            "Install with: pip install mhcflurry\n"
            "Then download models: mhcflurry-downloads fetch models_class1_presentation"
        ) from e

    print("[predict_mhc_with_fallback] Loading MHCflurry...", file=sys.stderr)

    # Set up model cache in writable directory
    _setup_mhcflurry_models()

    try:
        predictor = Class1PresentationPredictor.load()
    except Exception as e:
        # Model load failed; try downloading to writable cache
        model_dir = os.environ.get("MHCFLURRY_MODELS_DIR")

        if model_dir:
            print(f"[predict_mhc_with_fallback] Attempting to download models to {model_dir}", file=sys.stderr)
            try:
                result = subprocess.run(
                    [sys.executable, "-m", "mhcflurry.downloads", "fetch", "models_class1_presentation"],
                    capture_output=True,
                    text=True,
                    timeout=300,
                    env={**os.environ, "MHCFLURRY_MODELS_DIR": model_dir},
                )

                if result.returncode == 0:
                    print("[predict_mhc_with_fallback] Models downloaded; retrying load...", file=sys.stderr)
                    predictor = Class1PresentationPredictor.load()
                else:
                    raise Exception(result.stderr[:500])
            except Exception as download_err:
                raise MHCPredictionError(
                    f"MHCflurry model load and download both failed.\n\n"
                    f"Original error: {str(e)[:200]}\n"
                    f"Download error: {str(download_err)[:200]}\n\n"
                    f"To fix (run outside Claude Science):\n"
                    f"  1. export MHCFLURRY_MODELS_DIR=/path/to/writable/cache\n"
                    f"  2. mhcflurry-downloads fetch models_class1_presentation\n"
                    f"  3. Re-run with MHCFLURRY_MODELS_DIR=/path/to/cache set\n\n"
                    f"Or use NetMHCpan if installed on your system."
                ) from download_err
        else:
            raise MHCPredictionError(
                f"Failed to load MHCflurry model: {e}\n"
                "No writable cache directory found.\n\n"
                "To fix:\n"
                "  1. Set MHCFLURRY_MODELS_DIR to a writable path\n"
                "  2. Run: mhcflurry-downloads fetch models_class1_presentation\n"
                "  3. Re-run with MHCFLURRY_MODELS_DIR set"
            ) from e

    supported = set(predictor.supported_alleles)
    usable = [a for a in alleles if a in supported]

    for a in sorted(set(alleles) - supported):
        print(f"[predict_mhc_with_fallback] Warning: {a} not supported, skipped", file=sys.stderr)

    if not usable:
        raise MHCPredictionError(
            f"No requested alleles supported by MHCflurry {predictor.metadata()['model_version']}\n"
            f"Requested: {alleles}\n"
            f"Supported: {sorted(supported)[:10]}... (showing first 10)"
        )

    print(f"[predict_mhc_with_fallback] Predicting {len(peptides)} peptides × {len(usable)} alleles", file=sys.stderr)

    # MHCflurry >=2.2 format: one sample per allele for epitope triage
    allele_samples = {a: [a] for a in usable}
    df = predictor.predict(peptides=peptides, alleles=allele_samples, verbose=0)

    results = []
    for _, rec in df.iterrows():
        pep = rec.get("peptide")
        if pep not in peptides:
            continue

        allele = rec.get("best_allele") or rec.get("sample_name")
        affinity = float(rec.get("affinity", float("nan")))
        percentile = float(rec.get("affinity_percentile", float("nan")))
        pres_score = float(rec.get("presentation_score", float("nan")))
        pres_percentile = float(rec.get("presentation_percentile", float("nan")))

        results.append({
            "peptide": pep,
            "allele": allele,
            "affinity_nM": affinity,
            "percentile_rank": pres_percentile if include_percentile else percentile,
            "presentation_score": pres_score,
            "backend": "mhcflurry",
            "backend_version": predictor.metadata()["model_version"],
        })

    return results


def predict_with_netmhcpan(
    peptides: List[str],
    alleles: List[str],
    netmhcpan_path: str,
) -> List[Dict[str, Any]]:
    """
    Predict using NetMHCpan command-line tool.

    Args:
        peptides: List of sequences
        alleles: List of alleles
        netmhcpan_path: Path to netMHCpan executable

    Returns:
        Prediction results in unified format
    """
    print(f"[predict_mhc_with_fallback] Using NetMHCpan: {netmhcpan_path}", file=sys.stderr)

    # Write peptides to temp FASTA
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
        for i, pep in enumerate(peptides):
            f.write(f">pep_{i}\n{pep}\n")
        fasta_path = f.name

    try:
        allele_str = ",".join(alleles)
        result = subprocess.run(
            [netmhcpan_path, "-a", allele_str, "-p", fasta_path],
            capture_output=True,
            text=True,
            timeout=300,
        )

        if result.returncode != 0:
            raise MHCPredictionError(
                f"NetMHCpan failed: {result.stderr}"
            )

        # Parse NetMHCpan output (format varies; simplified here)
        results = []
        for line in result.stdout.splitlines():
            if line.startswith("#") or not line.strip():
                continue
            # This is simplified; real parsing depends on NetMHCpan output format
            parts = line.split()
            if len(parts) >= 4:
                results.append({
                    "peptide": parts[1],
                    "allele": parts[0],
                    "affinity_nM": float(parts[3]) if parts[3] != "nan" else float("nan"),
                    "percentile_rank": float(parts[4]) if len(parts) > 4 else float("nan"),
                    "backend": "netmhcpan",
                    "backend_version": "server",
                })

        return results

    finally:
        Path(fasta_path).unlink(missing_ok=True)


def predict_mhc_class1_with_fallback(
    peptides: List[str],
    alleles: List[str],
    include_percentile: bool = True,
    force_backend: Optional[str] = None,
) -> List[Dict[str, Any]]:
    """
    Predict MHC Class I binding with intelligent fallback.

    Route order:
    1. Local MHCflurry (if available and working)
    2. NetMHCpan server (if installed)
    3. Fail with clear instructions (no simulation)

    Args:
        peptides: List of peptide sequences
        alleles: List of HLA alleles
        include_percentile: Include percentile rank in output
        force_backend: Force specific backend ("mhcflurry", "netmhcpan")

    Returns:
        List of prediction results

    Raises:
        MHCPredictionError: If no backend available and cannot predict
    """
    env = detect_mhc_environment()
    print(f"[predict_mhc_with_fallback] Environment: {json.dumps(env, default=str, indent=2)}", file=sys.stderr)

    # Force backend if requested
    if force_backend == "mhcflurry":
        if env["can_use_mhcflurry"]:
            return predict_with_mhcflurry(peptides, alleles, include_percentile)
        else:
            raise MHCPredictionError(
                "force_backend='mhcflurry' but MHCflurry not available.\n"
                "Install: pip install mhcflurry\n"
                "Download models: mhcflurry-downloads fetch models_class1_presentation"
            )

    if force_backend == "netmhcpan":
        if env["can_use_netmhcpan"]:
            return predict_with_netmhcpan(
                peptides, alleles, env["netmhcpan_path"]
            )
        else:
            raise MHCPredictionError(
                "force_backend='netmhcpan' but NetMHCpan not installed.\n"
                "Install from: https://services.healthtech.dtu.dk/service.php?NetMHCpan"
            )

    # Route 1: Try MHCflurry
    if env["can_use_mhcflurry"] and env.get("mhcflurry_loads"):
        try:
            return predict_with_mhcflurry(peptides, alleles, include_percentile)
        except MHCPredictionError as e:
            print(f"[predict_mhc_with_fallback] MHCflurry failed: {e}", file=sys.stderr)
    elif env["can_use_mhcflurry"]:
        print(f"[predict_mhc_with_fallback] MHCflurry import failed: {env.get('mhcflurry_load_error')}", file=sys.stderr)

    # Route 2: Try NetMHCpan
    if env["can_use_netmhcpan"]:
        try:
            return predict_with_netmhcpan(
                peptides, alleles, env["netmhcpan_path"]
            )
        except MHCPredictionError as e:
            print(f"[predict_mhc_with_fallback] NetMHCpan failed: {e}", file=sys.stderr)

    # Route 3: No viable path
    raise MHCPredictionError(
        f"Cannot predict MHC Class I binding. No backend available.\n\n"
        f"Available options:\n"
        f"  1. MHCflurry (recommended, open-source)\n"
        f"     Install: pip install mhcflurry\n"
        f"     Download models: mhcflurry-downloads fetch models_class1_presentation\n\n"
        f"  2. NetMHCpan (requires license)\n"
        f"     Download: https://services.healthtech.dtu.dk/service.php?NetMHCpan\n"
        f"     Add to PATH\n\n"
        f"Current environment:\n"
        f"  - MHCflurry available: {env['can_use_mhcflurry']}\n"
        f"  - NetMHCpan available: {env['can_use_netmhcpan']}"
    )


if __name__ == "__main__":
    # Example usage
    example_peptides = ["FVFSTSFYL", "YVFSTSFYL", "MTMPLSCTK"]
    example_alleles = ["HLA-A*02:01", "HLA-A*11:01"]

    try:
        results = predict_mhc_class1_with_fallback(
            example_peptides,
            example_alleles,
        )
        for r in results:
            print(f"{r['peptide']} → {r['allele']}: {r['percentile_rank']:.3f}% (backend: {r['backend']})")
    except MHCPredictionError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)
