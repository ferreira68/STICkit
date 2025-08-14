import os, sys
sys.path.insert(0, os.path.abspath(".."))
project = "stickit"
extensions = ["sphinx.ext.autodoc", "sphinx_autodoc_typehints", "myst_parser"]
html_theme = "furo"
autodoc_typehints = "description"

