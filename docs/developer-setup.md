# Developer setup

This project uses **Poetry** for Python deps and other tools (e.g. RDKit, OpenMM).
Poetry remains the *single source of truth* for Python packages. You can use `poetry build` to create a package
that can be installed with `conda` or `pip`.


## Sphinx setup
Use
```shell
CHROME="$(command -v google-chrome-stable || command -v google-chrome || command -v chromium || command -v chromium-browser)"
```
to set the path to the Chrome executable for documentation builds.