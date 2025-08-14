#!/usr/bin/env bash
set -euo pipefail

# Sync environment.yml from Poetry lock without conda-lock

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")"/.. && pwd)"
cd "$repo_root"

# --- ensure poetry export is available
if ! command -v poetry >/dev/null 2>&1; then
  echo "ERROR: poetry not found on PATH" >&2
  exit 1
fi

POETRY_VER_RAW="$(poetry --version | awk '{print $3}')"
POETRY_MAJOR="${POETRY_VER_RAW%%.*}"

# On Poetry 1.x we need the export plugin (<1.9 for Poetry<2)
if [ "$POETRY_MAJOR" -lt 2 ]; then
  if ! poetry self show plugins 2>/dev/null | grep -qi "export"; then
    echo "Installing poetry-plugin-export <1.9.0 for Poetry 1.x ..."
    poetry self add "poetry-plugin-export<1.9.0"
  fi
fi

echo "Exporting requirements from Poetry lock..."
poetry export -f requirements.txt -o requirements.txt
poetry export -f requirements.txt --with dev -o requirements-dev.txt

echo "Writing environment.yml ..."
cat > environment.yml <<'YAML'
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
    - -r requirements-dev.txt
YAML

echo "Done. Updated:"
echo "  requirements.txt"
echo "  requirements-dev.txt"
echo "  environment.yml"
