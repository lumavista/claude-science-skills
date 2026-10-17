#!/usr/bin/env python
"""
Hybrid ESMFold2 folding with fallback routing.

Route folding requests based on environment:
1. If ESM_API_KEY set → use Biohub hosted API (works in sandbox)
2. Else if WSL2/Linux detected → use local folding on GPU
3. Else → fail with clear instructions (no simulation)

Usage:
    from fold_with_fallback import fold_proteins

    results = fold_proteins(
        sequences=[
            {"id": "WT", "seq": "MPRTKEINDAGLIVHSPQWFYK"},
            {"id": "Var01", "seq": "YPRTKEINDAGLIVHSPQWFYK"},
        ],
        model="esmfold2-fast-2026-05",
        use_api_if_available=True,
    )
    for result in results:
        print(f"{result['id']}: pLDDT={result['plddt_mean']:.1f}")
        with open(f"{result['id']}.pdb", "w") as f:
            f.write(result['pdb'])
"""

import os
import sys
import json
import platform
import subprocess
import tempfile
from pathlib import Path
from typing import List, Dict, Optional, Any


class FoldingEnvironmentError(Exception):
    """Raised when folding cannot proceed in current environment."""
    pass


def detect_environment() -> Dict[str, Any]:
    """Detect current environment and available folding methods."""
    info = {
        "platform": platform.system(),
        "has_api_key": "ESM_API_KEY" in os.environ,
        "api_key": os.environ.get("ESM_API_KEY", "").startswith("esm_"),
        "can_use_hosted": "ESM_API_KEY" in os.environ,
        "can_use_local": False,
        "wsl2_available": False,
        "wsl2_python": None,
    }

    # Check if we're on Linux/WSL2 with GPU
    if platform.system() in ("Linux", "Darwin"):
        info["can_use_local"] = True
        info["environment"] = "native_linux"

    # Check if WSL2 is accessible from Windows
    if platform.system() == "Windows":
        info["environment"] = "windows"
        try:
            result = subprocess.run(
                ["wsl", "--version"],
                capture_output=True,
                text=True,
                timeout=5,
            )
            if result.returncode == 0:
                info["wsl2_available"] = True
                # Try to find Python in WSL2
                wsl_python_result = subprocess.run(
                    ["wsl", "which", "python3"],
                    capture_output=True,
                    text=True,
                    timeout=5,
                )
                if wsl_python_result.returncode == 0:
                    info["wsl2_python"] = wsl_python_result.stdout.strip()
                    info["can_use_local"] = True
        except (FileNotFoundError, subprocess.TimeoutExpired):
            pass

    return info


def fold_with_hosted_api(
    sequences: List[Dict[str, str]],
    model: str = "esmfold2-fast-2026-05",
    timeout: int = 300,
) -> List[Dict[str, Any]]:
    """
    Fold sequences using Biohub hosted API.

    Args:
        sequences: List of {"id": "...", "seq": "..."} dicts
        model: Hosted model ID (e.g., "esmfold2-fast-2026-05")
        timeout: Request timeout in seconds

    Returns:
        List of {"id": "...", "pdb": "...", "plddt": [...], "plddt_mean": ...} dicts
    """
    try:
        from esm.sdk import esmfold2_client
        from esm.sdk.api import ESMProteinError, FoldingConfig
        from esm.utils.structure.input_builder import ProteinInput, StructurePredictionInput
    except ImportError as e:
        raise FoldingEnvironmentError(
            f"esm SDK not available: {e}\n"
            "Install with: pip install esm==3.4.1.post1"
        ) from e

    api_key = os.environ.get("ESM_API_KEY")
    if not api_key:
        raise FoldingEnvironmentError(
            "ESM_API_KEY not set. Set it to use hosted folding:\n"
            "  export ESM_API_KEY=esm_...\n"
            "  (create one at https://biohub.ai/developer-console/api-keys)"
        )

    print(f"[fold_with_fallback] Using Biohub hosted API: {model}")

    fold_input = StructurePredictionInput(
        sequences=[
            ProteinInput(id=seq["id"], sequence=seq["seq"])
            for seq in sequences
        ]
    )

    with esmfold2_client(
        model=model,
        url="https://biohub.ai",
        token=api_key,
        request_timeout=timeout,
    ) as client:
        result = client.fold_all_atom(fold_input, config=FoldingConfig())

    if isinstance(result, ESMProteinError):
        raise FoldingEnvironmentError(f"Biohub API error: {result}")

    # Convert results to unified format
    # MolecularComplex has per-residue arrays: plddt (0-1 scale), chain_id, sequence
    # All sequences are folded as chains in a single complex
    outputs = []

    import tempfile

    # Get per-residue metrics (plddt comes as 0-1, convert to 0-100)
    plddt_values = result.complex.plddt * 100 if result.complex.plddt is not None else None
    chain_ids = result.complex.chain_id if hasattr(result.complex, "chain_id") else None

    # Extract per-chain metrics by chain_id
    for i, seq in enumerate(sequences):
        target_chain = i  # Chain index in the output

        if chain_ids is not None:
            # Find residues for this chain
            chain_mask = chain_ids == target_chain
            chain_plddt = plddt_values[chain_mask] if plddt_values is not None else None
        else:
            # Single sequence or no chain info; use full plddt
            chain_plddt = plddt_values

        # Verify sequence length matches
        if chain_plddt is not None:
            predicted_length = len(chain_plddt)
            if predicted_length != len(seq["seq"]):
                raise FoldingEnvironmentError(
                    f"Sequence {seq['id']}: predicted {predicted_length} residues, "
                    f"but input was {len(seq['seq'])} residues"
                )

        # Save structure for this sequence
        # Note: result.complex.to_pdb() doesn't exist; use to_protein_complex().to_pdb(path)
        with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as f:
            pdb_path = f.name

        result.complex.to_protein_complex().to_pdb(pdb_path)
        with open(pdb_path, "r") as f:
            pdb_str = f.read()

        outputs.append({
            "id": seq["id"],
            "pdb": pdb_str,
            "plddt": chain_plddt.tolist() if chain_plddt is not None else None,
            "plddt_mean": float(chain_plddt.mean()) if chain_plddt is not None else None,
            "device": "biohub",
            "model": model,
            "length": len(seq["seq"]),
        })

    return outputs


def fold_with_local_wsl2(
    sequences: List[Dict[str, str]],
    python_exe: str,
    device: str = "cuda",
) -> List[Dict[str, Any]]:
    """
    Fold sequences using ESMFold2 in WSL2 environment.

    Args:
        sequences: List of {"id": "...", "seq": "..."} dicts
        python_exe: Path to Python executable in WSL2
        device: "cuda" or "cpu"

    Returns:
        List of folding results
    """
    print(f"[fold_with_fallback] Using WSL2 local folding: {python_exe}, device={device}")

    # Create temporary script in current directory (should be accessible from WSL2)
    script = """
import json
import sys
import torch
from esm.models.esm_fold_v1 import ESMFold
from esm.pretrained import esmfold_v1

device = "{device}"
sequences = json.loads('{sequences_json}')

model = esmfold_v1(device=device)
outputs = []

for seq_info in sequences:
    print(f"Folding {{seq_info['id']}}...", file=sys.stderr)
    structure = model.infer_pdb(seq_info['seq'])

    # Extract pLDDT if available
    plddt = getattr(structure, 'plddt', None)

    outputs.append({{
        "id": seq_info['id'],
        "pdb": str(structure),
        "plddt": plddt.tolist() if plddt is not None else None,
        "plddt_mean": float(plddt.mean()) if plddt is not None else None,
        "device": "wsl2_cuda",
        "model": "esmfold2-fast",
    }})

json.dump(outputs, sys.stdout)
""".format(
        device=device,
        sequences_json=json.dumps(sequences).replace('"', '\\"'),
    )

    with tempfile.NamedTemporaryFile(
        mode="w",
        suffix=".py",
        delete=False,
        dir=".",
    ) as f:
        f.write(script)
        script_path = f.name

    try:
        result = subprocess.run(
            [python_exe, script_path],
            capture_output=True,
            text=True,
            timeout=600,
        )

        if result.returncode != 0:
            raise FoldingEnvironmentError(
                f"WSL2 folding failed:\n{result.stderr}"
            )

        outputs = json.loads(result.stdout)
        return outputs

    finally:
        Path(script_path).unlink(missing_ok=True)


def fold_with_fallback(
    sequences: List[Dict[str, str]],
    model: str = "esmfold2-fast-2026-05",
    device: str = "cuda",
    force_api: bool = False,
    force_local: bool = False,
) -> List[Dict[str, Any]]:
    """
    Fold sequences with intelligent fallback routing.

    Route order:
    1. If force_api=True and ESM_API_KEY → use hosted API
    2. Else if ESM_API_KEY → use hosted API
    3. Else if WSL2 available → use local folding
    4. Else if native Linux/Mac → use local folding
    5. Else → fail with clear instructions

    Args:
        sequences: List of {"id": "...", "seq": "..."} dicts
        model: Model name for hosted ("esmfold2-fast-2026-05") or local ("esmfold2-fast")
        device: "cuda" or "cpu" (ignored for hosted API)
        force_api: If True, only try hosted API
        force_local: If True, only try local folding

    Returns:
        List of folding results with pdb, plddt, device info

    Raises:
        FoldingEnvironmentError: If folding cannot proceed
    """
    env = detect_environment()
    print(f"[fold_with_fallback] Environment: {json.dumps(env, default=str, indent=2)}")

    # Route 1: Force API if requested
    if force_api:
        if env["can_use_hosted"]:
            return fold_with_hosted_api(sequences, model=model)
        else:
            raise FoldingEnvironmentError(
                "force_api=True but ESM_API_KEY not set.\n"
                "Set it at https://biohub.ai/developer-console/api-keys and export:\n"
                "  export ESM_API_KEY=esm_..."
            )

    # Route 2: Prefer hosted API if available
    if env["can_use_hosted"]:
        try:
            return fold_with_hosted_api(sequences, model=model)
        except Exception as e:
            print(f"[fold_with_fallback] Hosted API failed: {e}", file=sys.stderr)
            if force_local:
                print("[fold_with_fallback] Falling back to local...", file=sys.stderr)
            else:
                raise

    # Route 3: Try WSL2 if available
    if env["wsl2_available"] and env["wsl2_python"]:
        try:
            return fold_with_local_wsl2(
                sequences,
                python_exe=env["wsl2_python"],
                device=device,
            )
        except Exception as e:
            print(f"[fold_with_fallback] WSL2 folding failed: {e}", file=sys.stderr)
            if not force_local:
                raise

    # Route 4: Native Linux/Mac
    if platform.system() in ("Linux", "Darwin"):
        try:
            from esm.models.esm_fold_v1 import ESMFold
            from esm.pretrained import esmfold_v1
        except ImportError as e:
            raise FoldingEnvironmentError(
                f"esm SDK not available for local folding: {e}"
            ) from e

        print(f"[fold_with_fallback] Using local {platform.system()} folding on {device}")
        model_obj = esmfold_v1(device=device)
        outputs = []

        for seq in sequences:
            print(f"Folding {seq['id']}...", file=sys.stderr)
            structure = model_obj.infer_pdb(seq["seq"])
            plddt = getattr(structure, "plddt", None)

            outputs.append({
                "id": seq["id"],
                "pdb": str(structure),
                "plddt": plddt.tolist() if plddt is not None else None,
                "plddt_mean": float(plddt.mean()) if plddt is not None else None,
                "device": device,
                "model": model,
            })

        return outputs

    # Route 5: No viable path
    raise FoldingEnvironmentError(
        f"Cannot fold proteins on {platform.system()}.\n\n"
        "Available options:\n"
        "  1. SET ESM_API_KEY for hosted folding (recommended, works in sandbox)\n"
        "     • Create key: https://biohub.ai/developer-console/api-keys\n"
        "     • export ESM_API_KEY=esm_...\n\n"
        "  2. Use Linux/WSL2 with local GPU\n"
        "     • WSL2: Install ESMFold2 in WSL2 Python environment\n"
        "     • Linux: Install esm==3.4.1.post1 with CUDA torch\n"
        "     • This script will auto-detect and use it\n\n"
        "  3. For Windows without WSL2: Use hosted API (option 1)"
    )


if __name__ == "__main__":
    # Example usage
    example_sequences = [
        {"id": "FVFSTSFYL_WT", "seq": "FVFSTSFYL"},
        {"id": "FVFSTSFYL_F1Y", "seq": "YVFSTSFYL"},
    ]

    try:
        results = fold_with_fallback(example_sequences)
        print(json.dumps(results, default=str, indent=2))
    except FoldingEnvironmentError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)
