#!/usr/bin/env bash
# setup_wsl.sh - one-time setup of the Linux tools inside WSL (Ubuntu).
#
# Run from the repository root, inside WSL:
#     bash workflow/scripts/setup_wsl.sh
#
# What it does:
#   1. installs Miniforge (conda + mamba) in ~/miniforge3 if no conda is found
#   2. creates the environment 'salmonella-amr-tools' from envs/step4_amr.yml (skipped if it exists)
#   3. downloads the AMRFinderPlus database (amrfinder -u)
#   4. prints the tool versions, which you should keep for the report
#
# Safe to run again: finished parts are skipped.
set -eo pipefail

ENV_NAME="salmonella-amr-tools"
ENV_FILE="envs/step4_amr.yml"
[ -f "$ENV_FILE" ] || { echo "Run this from the repository root (cannot find $ENV_FILE)."; exit 1; }

# 1. conda
if [ -x "$HOME/miniforge3/bin/conda" ]; then
  CONDA_BASE="$HOME/miniforge3"
elif command -v conda >/dev/null 2>&1; then
  CONDA_BASE="$(conda info --base)"
else
  echo ">>> Installing Miniforge into $HOME/miniforge3 ..."
  command -v curl >/dev/null 2>&1 || { echo "curl is missing. Run: sudo apt update && sudo apt install -y curl"; exit 1; }
  curl -L -o /tmp/Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
  bash /tmp/Miniforge3.sh -b -p "$HOME/miniforge3"
  CONDA_BASE="$HOME/miniforge3"
fi
# shellcheck disable=SC1091
source "$CONDA_BASE/etc/profile.d/conda.sh"

# 2. environment
if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo ">>> Environment $ENV_NAME already exists, skipping creation."
else
  echo ">>> Creating environment $ENV_NAME (this can take 5-15 minutes) ..."
  mamba env create -y -f "$ENV_FILE"
fi
conda activate "$ENV_NAME"

# 3. database
echo ">>> Updating the AMRFinderPlus database ..."
amrfinder -u

# 4. versions
echo ">>> Tool versions:"
amrfinder --version || true
amrfinder -l | head -n 5 || true
datasets --version || true
echo
if ! grep -q "conda initialize" "$HOME/.bashrc" 2>/dev/null; then
  "$CONDA_BASE/bin/conda" init bash > /dev/null
  echo "Conda was added to ~/.bashrc. Close this Ubuntu window and open it again before the next step."
fi
echo "Setup finished. In every new WSL window, run:  conda activate $ENV_NAME"
