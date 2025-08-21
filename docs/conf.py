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
html_theme = "nature"
html_theme_options = {
    "sidebarwidth": "24em",   # widen the sidebar to accommodate long function names
}

# Autodoc settings
autodoc_typehints = "description"

# Mermaid settings
mermaid_version = "11.2.0"
mermaid_d3_zoom = True
mermaid_init_js = "mermaid.initialize({startOnLoad:true});"
