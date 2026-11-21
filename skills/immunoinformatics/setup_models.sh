#!/bin/bash
# Download MHCflurry models to skill directory for offline use

set -e

SKILL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODEL_DIR="$SKILL_DIR/models/mhcflurry"

echo "Setting up MHCflurry models in: $MODEL_DIR"

# Create model directory
mkdir -p "$MODEL_DIR"

# Set environment to store models in skill directory
export MHCFLURRY_MODELS_DIR="$MODEL_DIR"

# Download models
echo "Downloading MHCflurry models (this may take a few minutes)..."
python3 -m mhcflurry.downloads fetch models_class1_presentation

echo ""
echo "✓ Models downloaded to: $MODEL_DIR"
echo ""
echo "epitope_triage.py finds this directory on its own, so no further setup is"
echo "needed for it. For anything else that loads MHCflurry directly, export the"
echo "path yourself -- the export above applied only to this script's own shell:"
echo ""
echo "    export MHCFLURRY_MODELS_DIR=\"$MODEL_DIR\""
echo ""
echo "Do NOT commit these weights: they are hundreds of MB and .gitignore"
echo "excludes skills/immunoinformatics/models/ deliberately. Re-run this script"
echo "instead; it is the cache, not a tracked artefact."
