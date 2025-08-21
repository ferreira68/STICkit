import os, sys
sys.path.insert(0, os.path.abspath(".."))
project = "stickit"
extensions = [
        "sphinx.ext.autodoc",
        "sphinx_autodoc_typehints",
        "sphinxcontrib.mermaid",
        "myst_parser"
]
autodoc_typehints = "description"

# Use Mermaid code fences in Markdown files:
myst_fence_as_directive = ["mermaid"]

# --- Mermaid settings ---
mermaid_output_format = "raw"
html_static_path = ["_static"]
mermaid_use_local = "_static/js/mermaid.min.mjs"

# Ensure Mermaid initializes
mermaid_init_js = "mermaid.initialize({ startOnLoad: true });"

# Optional Sphinx theme:
html_theme = "furo"
