
# STICkit — AI/CompChem Python Coding Rules (2025)

**Scope:** Computational chemistry + AI/ML repos.  
**You prefer:** Streamlit for prototyping; React (+ FastAPI) for production UIs; seaborn; PostgreSQL + SQLAlchemy (sync); Sphinx + Google‑style docstrings; Mermaid diagrams with auto‑generated class/call graphs; **no Ray** (use `multiprocessing`, optionally PySpark).

## 0) Principles
- Scientist‑first UX; reproducibility; separation of concerns; docs‑as‑code; diagrams‑as‑code.

## 1) Layout
```
STICkit/
  stickit/               # core library
  apps/streamlit/        # Streamlit prototypes
  services/api/          # FastAPI backend for React
  ui/react/              # React app
  scripts/               # utilities (no business logic)
  configs/               # YAML configs (non‑secret) + encrypted or SecretsManager
  docs/                  # Sphinx + Mermaid
  tests/                 # pytest
  data/                  # gitignored
```
Poetry; Python ≥3.11.

## 2) Tooling
Ruff + Black, mypy(strict for new code), Bandit, pre‑commit (also runs diagrams drift check).

## 3) Config & secrets
Defaults in `configs/settings.yaml`. Secrets from **AWS Secrets Manager** (preferred) or `configs/settings.enc.yaml` (SOPS). Avoid env vars for secrets.

## 4) Data/DB
PostgreSQL + SQLAlchemy 2.x (sync). Alembic migrations. Large artifacts in object storage; DB holds metadata/URIs.

## 5) Interfaces
Streamlit for prototypes; React + FastAPI for production. seaborn first; matplotlib as needed.

## 6) Docs
Sphinx + `napoleon` (Google style) + `sphinx_autodoc_typehints` + `myst-parser` + `sphinxcontrib-mermaid`.

## 7) Diagrams
- Class diagrams: `pyreverse -o mmd -p stickit stickit -d docs/_mermaid`
- Call graph: `scripts/generate_callgraph.py` (pyan3 → Mermaid)
- CI **fails** if `docs/_mermaid/*.mmd` drift.

## 8) Testing
`pytest`; `@pytest.mark.gpu` and `--require-gpu`. Hypothesis for invariants; pytest‑benchmark for hot paths.

## 9) ML/DS
MLflow for model storage/load; minimal tracking. Group/scaffold CV; fix seeds; log versions/CUDA.

## 10) CI/CD
Lint → typecheck → tests → docs → diagrams‑check (fail on drift).

## 11) Parallelism (no Ray)
- Node‑local CPU: `multiprocessing` / `ProcessPoolExecutor` (picklable functions).  
- Big ETL/joins: **PySpark** as optional pre‑chemistry stage; chemistry stays pandas/RDKit by default.

## 12) Visualization & trajectories
- 2D: RDKit; stats: seaborn (`context="talk"`, `style="whitegrid"`).  
- 3D/MD:
  - **Notebooks:** **NGLview** (with MDAnalysis/MDTraj).  
  - **Web UIs:** **Mol\*** (React/Streamlit via HTML component).

## 13) Logging
`logging` with optional JSON; redact secrets/proprietary strings.

## 14) Repro & envs
Poetry‑locked; pin native libs (RDKit/OpenMM) by minor release; snapshot runtime env with results.

## 15) Git/release
Conventional commits; changelog; tag releases; diagrams included in diffs.

## 16) Templates
- FastAPI stub; seaborn house style; Mermaid generation commands as above.
