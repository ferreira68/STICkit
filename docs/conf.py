import os, sys
sys.path.insert(0, os.path.abspath(".."))

# Project name
project = "stickit"

# Sphinx extensions
extensions = [
        "sphinx.ext.autodoc",
        "sphinx_autodoc_typehints",
        "sphinxcontrib.mermaid",
        "myst_parser",
]

# Set up theme
html_theme = "furo"

# Autodoc settings
autodoc_typehints = "description"

# Mermaid settings
mermaid_version = "10.3.0"
mermaid_init_js = "mermaid.initialize({startOnLoad:true});"
mermaid_params = ['-p', 'puppeteer-config.json']