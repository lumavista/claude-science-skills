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
echo "These models are now cached and will be used automatically."
echo "Commit to git: git add $MODEL_DIR && git commit -m 'Cache MHCflurry models'"
