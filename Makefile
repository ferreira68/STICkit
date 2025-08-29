SHELL := /bin/bash
PY := python

export PYTHONHASHSEED=0

.PHONY: help install test docs env-sync lint fmt typecheck clean uml callgraph diagrams diagrams-check ci-check

help:
	@echo "Targets:"
	@echo "  install         - poetry install with dev deps"
	@echo "  test            - run pytest"
	@echo "  docs            - build Sphinx HTML docs"
	@echo "  env-sync        - sync Poetry -> conda/env files (runs scripts/env-sync.sh if present)"
	@echo "  lint            - run ruff lint checks"
	@echo "  fmt             - auto-fix with ruff (and black if installed)"
	@echo "  typecheck       - run mypy on stickit"
	@echo "  uml             - generate Mermaid class diagrams via pyreverse (.mmd)"
	@echo "  callgraph       - generate Mermaid call graph (PyCG if available, else AST fallback)"
	@echo "  diagrams        - run uml + callgraph"
	@echo "  diagrams-check  - fail if diagrams drift from HEAD"
	@echo "  clean           - remove build and docs artifacts"
	@echo "  ci-check        - lint + typecheck + test + docs + diagrams-check (intended for CI)"

install:
	poetry install --with dev

test:
	pytest -q

docs:
	sphinx-build -b html docs docs/_build/html

env-sync:
	@if [ -x scripts/env-sync.sh ]; then \
	  echo "Running scripts/env-sync.sh..."; \
	  scripts/env-sync.sh; \
	else \
	  echo "No env-sync script found (scripts/env-sync.sh). Skipping."; \
	fi

lint:
	ruff check .

fmt:
	ruff check --fix .
	ruff format .
	@if command -v black >/dev/null 2>&1; then black .; else echo "black not installed; skipped"; fi

typecheck:
	mypy stickit

clean:
	rm -rf dist build docs/_build

uml:
	# Mermaid (.mmd) via pyreverse (requires modern pylint/pyreverse)
	@set -e; \
	tmp=$$(mktemp); \
	if pyreverse -o mmd -p stickit stickit -d docs/_mermaid 2>$$tmp; then \
	  echo "pyreverse: wrote Mermaid .mmd"; \
	else \
	  if grep -q 'Format: "mmd" not recognized' $$tmp; then \
	    echo "pyreverse: falling back to DOT → Mermaid"; \
	    pyreverse -o dot -p stickit stickit -d docs/_mermaid; \
	    $(PY) scripts/dot_to_mermaid.py docs/_mermaid/classes.dot docs/_mermaid/classes.mmd; \
	  else \
	    cat $$tmp; rm -f $$tmp; exit 1; \
	  fi; \
	fi; rm -f $$tmp

callgraph:
	$(PY) scripts/generate_callgraph.py

diagrams: uml callgraph

diagrams-check: diagrams
	@git update-index -q --refresh
	@if git diff --name-only -- docs/_mermaid | grep -q . ; then \
	  echo '\nDiagrams changed. Commit the updates.'; exit 1; \
	else \
	  echo 'Diagrams are up-to-date.'; \
	fi

ci-check: lint typecheck test docs diagrams-check

