#!/usr/bin/env bash
set -euo pipefail

# Keep Conda and Poetry env specs in sync:
# - requirements*.txt exported from Poetry (locked)
# - environment.in.yml minimal conda spec (core non-Python deps)
# - conda-lock.yml multi-platform lockfile
# - environment.yml rendered (linux-64 convenience)

here="$(cd "$(dirname "${BASH_SOURCE[0]}")"/.. && pwd)"
cd "$here"

# --- sanity checks
if ! command -v poetry >/dev/null 2>&1; then
  echo "ERROR: poetry not found on PATH" >&2
  exit 1
fi

# ensure export plugin
if ! poetry self show plugins 2>/dev/null | grep -qi "export"; then
  echo "Installing poetry-plugin-export..."
  poetry self add "poetry-plugin-export"
fi

# ensure conda-lock (prefer Poetry-run)
if ! poetry run python -c "import conda_lock" 2>/dev/null; then
  echo "Installing conda-lock into Poetry env (dev group recommended) ..."
  poetry run python -m pip install conda-lock
fi

# --- export Poetry requirements (locked)
echo "Exporting requirements from Poetry..."
poetry export -f requirements.txt -o requirements.txt
poetry export -f requirements.txt --with dev -o requirements-dev.txt

# --- write minimal conda spec (non-Python binaries via conda; Python deps via pip)
echo "Writing environment.in.yml..."
cat > environment.in.yml <<'YAML'
name: stickit
channels:
  - conda-forge
dependencies:
  - python=3.11
  - rdkit
  - openmm
  - openff-toolkit
  - pip
  - pip:
      - -r requirements.txt
YAML

# --- lock for multiple platforms (use micromamba if available)
echo "Locking with conda-lock..."
MAMBA_FLAG=""
if command -v micromamba >/dev/null 2>&1; then
  MAMBA_FLAG="--mamba"
fi

poetry run conda-lock lock $MAMBA_FLAG -f environment.in.yml \
  -p linux-64 -p osx-64 -p osx-arm64 -p win-64

# --- render convenience environment.yml for linux-64
echo "Rendering environment.yml (linux-64)..."
poetry run conda-lock render -p linux-64 -f environment.in.yml > environment.yml

echo "Done. Files updated:"
echo "  requirements.txt"
echo "  requirements-dev.txt"
echo "  environment.in.yml"
echo "  conda-lock.yml"
echo "  environment.yml (linux-64 convenience)"

