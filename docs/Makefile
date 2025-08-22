SPHINXBUILD ?= sphinx-build
SOURCEDIR   ?= .
BUILDDIR    ?= $(SOURCEDIR)/_build

.PHONY: help clean html

help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)"

clean:
	rm -rf "$(BUILDDIR)"

# Ensure diagrams are generated before building HTML
html: diagrams
	$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)"
