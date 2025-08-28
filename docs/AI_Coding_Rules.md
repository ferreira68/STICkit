
# STICkit — AI/CompChem Python Coding Rules (2025)

**Scope:** Guidelines for STICkit and related repos focused on computational chemistry and AI/ML for drug discovery.  
**Your preferences baked in:** Streamlit for prototyping; React (+ FastAPI) for production UIs; seaborn for plotting; PostgreSQL + SQLAlchemy ORM; Sphinx (Google-style docstrings); Mermaid diagrams with auto‑generated class/call graphs; **no Ray** (prefer `multiprocessing`, optionally PySpark).

---

## 0) Principles
- **Scientist-first UX:** small, inspectable functions; notebooks and CLIs stay thin.
- **Reproducibility:** deterministic seeds, captured environments, data lineage.
- **Separation of concerns:** core algorithms independent from I/O, plotting, and UI.
- **Docs-as-code:** every public function has Google-style docstrings that build into Sphinx.
- **Diagrams-as-code:** Mermaid sources generated in CI from the codebase—fail if stale.

---

## 1) Project layout & packaging
```
repo-root/
  stickit/                 # core library (src or flat layout; prefer src/ for new packages)
  apps/streamlit/          # Streamlit prototypes
  services/api/            # FastAPI service for React frontends
  ui/react/                # React (Vite) app
  scripts/                 # one-off utilities (no business logic)
  configs/                 # YAML configs (non-secret), plus encrypted file or AWS SecretsManager
  docs/                    # Sphinx + Mermaid
  tests/                   # pytest suite
  data/                    # gitignored large assets
```
- **Packaging:** Poetry as the canonical tool. Pin Python to `>=3.11` for modern type system.
- **Imports:** absolute imports inside `stickit/`; avoid circular deps.

---

## 2) Tooling & quality gates
- **Formatting/Linting:** Ruff (lint + import sort) and Black (optional if you prefer Ruff-only formatting).
- **Types:** mypy in `strict` for new/modified modules; legacy modules can be `# type: ignore` gated.
- **Security:** Bandit for code; Pip‑audit/Safety for deps (as CI step).
- **Pre‑commit:** run Ruff/Black/mypy/Bandit and a “diagrams drift” check (see §10).

**pyproject (snippets)**
```toml
[tool.ruff]
line-length = 100
target-version = "py311"
lint.select = ["E","F","I","B","UP","S","SIM"]

[tool.mypy]
python_version = "3.11"
strict = true
```
---

## 3) Configuration & secrets
- **Defaults:** `configs/settings.yaml` (checked into git, non-sensitive).
- **Secrets (preferred):** **AWS Secrets Manager** JSON merged on load.
- **Alternative:** `configs/settings.enc.yaml` encrypted with **sops**.  
- **Avoid env vars for secrets.** Read once at process start and pass down as objects.

**Python loader skeleton**
```python
# stickit/config/settings.py
@dataclass
class DBConfig: url: str
@dataclass
class Settings:
    db: DBConfig
    mlflow_tracking_uri: str = "file:./mlruns"
    secrets_manager: SecretsManagerConfig = SecretsManagerConfig()
    @classmethod
    def load(cls) -> "Settings": ...
```
---

## 4) Data & storage
- **Primary DB:** PostgreSQL + SQLAlchemy 2.x ORM (sync sessions).
- **Migrations:** Alembic with a predictable revision naming scheme.
- **Large artifacts** (SD files, MD trajectories, models): object storage (S3 or on‑prem S3‑compatible). DB stores metadata and URIs only.

**Session pattern (sync)**
```python
engine = create_engine(settings.db.url, future=True)
SessionLocal = sessionmaker(bind=engine, expire_on_commit=False)
with SessionLocal() as s:
    s.add(obj); s.commit()
```

---

## 5) Interfaces (UIs & APIs)
- **Prototyping:** Streamlit (scientist-friendly). Embed images/plots/tables; offload heavy work to library functions.
- **Production UI:** React + FastAPI (OpenAPI schema, Pydantic models). Use REST/JSON; add WebSockets/SSE for streaming progress.
- **Plotting:** **seaborn** first; fallback to matplotlib for low-level control. Share plotting helpers in `stickit/viz/`.

---

## 6) Documentation
- **Sphinx** with `sphinx.ext.autodoc`, `sphinx.ext.napoleon` (Google style), `sphinx_autodoc_typehints`, `myst-parser`, and `sphinxcontrib-mermaid`.
- Every public class/function has a short summary line, Args/Returns/Raises, and Examples.

**Google-style docstring example**
```python
def enumerate_tautomers(mol: Chem.Mol, pH: float = 7.4) -> list[Chem.Mol]:
    \"\"\"Enumerate likely tautomers at a given pH.

    Args:
        mol: Input RDKit molecule.
        pH: Henderson–Hasselbalch target pH.

    Returns:
        List of unique tautomers as RDKit molecules.
    \"\"\"
```

---

## 7) Diagrams & auto-generation (Mermaid)
- **Class diagrams:** `pyreverse -o mermaid -p stickit stickit > docs/_mermaid/classes.mmd`
- **Call graphs:** generate DOT with `pyan3` → convert to Mermaid (`scripts/generate_callgraph.py`).
- **CI must fail** if `docs/_mermaid/*.mmd` drift from source (`make diagrams-check`).

Mermaid snippets are included in docs with:
```rst
.. mermaid::
   graph TD
     A[Enumerate] --> B[Protonate]
```

---

## 8) Testing
- `pytest` with small, deterministic samples. Use `@pytest.mark.gpu` where CUDA is required.
- **GPU-aware CI:** mark GPU tests as skipped unless `--require-gpu` is passed (local/HPC).  
- **Property tests:** Hypothesis for chemistry invariants (charge, valence, stereochemistry).  
- **Benchmarking:** `pytest-benchmark` for hot paths (enumeration, substructure, featurization).

---

## 9) ML/DS conventions
- **MLflow for model storage/retrieval**, minimal tracking by default. Log model flavors (sklearn/pyfunc) and dependencies.
- **Splits:** scaffold or group-aware CV for chemical series; document split logic.
- **Randomness:** fix seeds; log versions and CUDA info on run start.

---

## 10) CI/CD
- **GitHub Actions:** lint → typecheck → tests → docs → `make diagrams-check` (fail on drift).
- **Artifacts:** publish docs HTML; (optional) build/publish Docker images for API/UI.
- **Pre-commit:** run locally; same gates in CI.

---

## 11) Performance & parallelism (no Ray)
- **Node-local parallelism:** `multiprocessing` (or `concurrent.futures.ProcessPoolExecutor`) for CPU‑bound work; keep tasks picklable and small. Chunk iterables; avoid global state.
- **I/O & ETL at scale:** **PySpark** (optional) when the workload is dominated by distributed joins/aggregations/file shuffles and data exceeds a single node’s memory. Keep chemistry steps on pandas/RDKit unless you explicitly adopt RDKit‑on‑Spark.
- **Vectorization:** prefer NumPy/Numba for tight loops before reaching for multiprocessing.

**Multiprocessing template**
```python
from concurrent.futures import ProcessPoolExecutor, as_completed
def _worker(payload):  # keep top-level, pure, picklable
    return heavy_function(payload)

def run_parallel(items, max_workers=None):
    out = []
    with ProcessPoolExecutor(max_workers=max_workers) as ex:
        futures = [ex.submit(_worker, x) for x in items]
        for f in as_completed(futures):
            out.append(f.result())
    return out
```

**PySpark usage (optional)**
```python
from pyspark.sql import SparkSession
spark = SparkSession.builder.appName("stickit").getOrCreate()
df = spark.read.parquet("s3://bucket/ds.parquet")
# Wide ETL/join/aggregate...
pdf = df.toPandas()  # switch to pandas for RDKit chemistry
```

**When to prefer pandas vs Spark?**
- Use **pandas** when:
  - You’re doing chemistry (RDKit) operations or small/medium data (< a few GB in memory).
  - You need quick iteration in notebooks; seamless plotting with seaborn/matplotlib.
- Consider **PySpark** when:
  - Source data live in a data lake (S3/HDFS) and you need scalable joins/group-bys.
  - You’re aggregating >100M rows or multi‑TB inputs **before** chemistry steps.
  - You already have a Spark cluster and need repeatable ETL pipelines.
- Hybrid pattern: ETL in Spark → **narrowed** pandas frame → RDKit/compute → results back to DB/object storage.

---

## 12) Visualization (2D/3D + trajectories)
- **2D:** RDKit drawing utilities; wrap in `stickit/viz/` and return images/figures instead of saving to disk directly.
- **Statistical plots:** **seaborn** with a project “house style” (`context="talk"`, `style="whitegrid"`).
- **3D molecules & MD trajectories:**
  - **Jupyter notebooks:** **NGLview** is recommended. It integrates well with MDTraj/MDAnalysis, supports DCD/XTCTRAJ, and gives interactive widget controls.
  - **Web/React/Streamlit UIs:** **Mol\*** is recommended. It’s a robust WebGL viewer used by RCSB; easy to embed in React and can stream large structures. In Streamlit, embed via `st.components.v1.html(...)`.
  - Keep a small adapter layer that prepares structures/trajectories for the chosen viewer (PDB, mmCIF, topology + frames). Prefer remote object storage URLs for large trajectories.

**NGLview (notebook)**
```python
import nglview as nv, MDAnalysis as mda
u = mda.Universe("topology.psf", "traj.dcd")
nv.show_mdanalysis(u)
```

**Mol\* (React sketch)**
```tsx
// Install molstar; then in a React component:
import { Viewer } from 'molstar/lib/viewer';
useEffect(() => {
  const v = new Viewer(containerRef.current!, { layoutIsExpanded: true });
  v.loadPdb('https://files.rcsb.org/view/1CRN.pdb');
}, []);
```

---

## 13) Logging
- Standard `logging` with optional JSON formatter. No `print` in library code.
- Redact sensitive fields (tokens, secrets, proprietary SMILES) by default.

---

## 14) Reproducibility & environments
- Poetry‑locked dependencies; pin core native libs (RDKit/OpenMM) by minor release.
- Persist a run-time environment snapshot (package versions, CUDA, GPU model) alongside results.

---

## 15) Git & release hygiene
- Conventional commits; one logical change per PR; draft early.
- Changelog maintained; tag releases.
- Docs and diagrams are part of the diff; CI will block if not regenerated.

---

## 16) House templates
**FastAPI stub for React**
```python
from fastapi import FastAPI
from pydantic import BaseModel

class PredictReq(BaseModel): smiles: str
class PredictResp(BaseModel): result: float

app = FastAPI(title="STICkit API")
@app.post("/predict", response_model=PredictResp)
def predict(req: PredictReq) -> PredictResp:
    return PredictResp(result=0.0)
```

**Seaborn style**
```python
import seaborn as sns
def apply_house_style(context="talk", style="whitegrid"): sns.set_theme(context=context, style=style)
```

**Mermaid generation (CI)**
```bash
pyreverse -o mermaid -p stickit stickit > docs/_mermaid/classes.mmd
python scripts/generate_callgraph.py
```

---

## 17) Decisions you can revisit later
- **Async SQLAlchemy:** not needed now (slow‑changing data), but trivial to adopt with FastAPI if concurrency becomes important.
- **Spark for chemistry:** only if you commit to RDKit‑on‑Spark or treat Spark purely as pre‑chemistry ETL.
- **Viewer choice:** NGLview stays the notebook default; Mol\* for any browser‑based UI.
