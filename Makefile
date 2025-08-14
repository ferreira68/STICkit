.PHONY: help install test docs clean env-sync lint fmt

help:
	@echo "Targets:"
	@echo "  install   - poetry install with dev deps"
	@echo "  test      - run pytest"
	@echo "  docs      - build Sphinx HTML docs"
	@echo "  env-sync  - sync Poetry -> conda-lock/env files"
	@echo "  lint      - placeholder (flake8/ruff if you add)"
	@echo "  fmt       - placeholder (black/ruff format)"
	@echo "  clean     - remove build and docs artifacts"

install:
	poetry install --with dev

test:
	poetry run pytest -q

docs:
	poetry run sphinx-build -b html docs docs/_build/html

env-sync:
	./scripts/sync_env.sh

lint:
	@echo "Add ruff/flake8 if desired"; exit 0

fmt:
	@echo "Add black/ruff if desired"; exit 0

clean:
	rm -rf dist build docs/_build

