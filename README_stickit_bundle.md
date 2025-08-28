# STICkit Bundle (Non-breaking Augmentation)

This bundle only **adds** files or uses `*.stickit.*` names so it won't overwrite your repo.
Unzip from the directory that **contains** the `STICkit/` folder:

```bash
unzip stickit_stack_update.zip
```

## What you get
- `Makefile.stickit` with diagram targets; call with `make -f Makefile.stickit ...`
- `.pre-commit-config.yaml` (only if you don't have one; you can rename/merge)
- `scripts/generate_callgraph.py` → Mermaid call graph
- `docs/_mermaid/` placeholders + `docs/diagrams_stickit.md`
- GitHub Actions: `.github/workflows/ci-augment.yml`
- Config/secrets loader: `stickit/config/settings.py` + `configs/settings.yaml`
- SQLAlchemy sync session boilerplate + example model
- MLflow minimal registry wrapper
- Parallel helpers: `stickit/parallel/mp.py`, `stickit/parallel/spark.py` (no Ray)
- Viz: seaborn house style, plus Mol* Streamlit demo and NGLview snippet
- API/UI stubs: FastAPI service, React skeleton

## First-time setup
```bash
poetry install --with dev
pre-commit install
make -f Makefile.stickit diagrams       # generate Mermaid UML & callgraph
git add docs/_mermaid/*.mmd
git commit -m "chore: add Mermaid diagrams"
```

## Optional Spark
Only use PySpark for large ETL/join workloads. Chemistry remains pandas/RDKit-first.
