import os, sys
sys.path.insert(0, os.path.abspath(".."))
project = "stickit"
extensions = [
        "sphinx.ext.autodoc",
        "sphinx_autodoc_typehints",
        "sphinxcontrib.mermaid",
        "myst_parser"
]
html_theme = "furo"
autodoc_typehints = "description"

# Use Mermaid code fences in Markdown files:
myst_fence_as_directive = ["mermaid"]   # ```mermaid ... ```
# (See MyST + Mermaid note.)  # :contentReference[oaicite:1]{index=1}

# --- Mermaid settings ---
# Option A: simplest (HTML) – render in browser with JS
mermaid_output_format = "raw"           # default; no extra deps

# (Optional) pin versions or use local assets
# mermaid_version = "11.2.0"
mermaid_init_js = "mermaid.initialize({startOnLoad:true});"
# (Config keys reference.)  # :contentReference[oaicite:2]{index=2}

# Optional Sphinx theme:
html_theme = "furo"
