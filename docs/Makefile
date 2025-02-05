# Minimal makefile for Sphinx documentation with Doxygen integration

SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
DOXYGEN       ?= doxygen
SOURCEDIR     = .
BUILDDIR      = _build

.PHONY: help Makefile doxygen html publish

# Default help command
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Run Doxygen to generate XML for Breathe
doxygen:
	@$(DOXYGEN) Doxyfile
	@echo "Doxygen build completed."

# Build HTML docs, ensuring Doxygen runs first
html: doxygen
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Publish target for GitHub Pages
publish: html
	@cp -r $(BUILDDIR)/html/* .
	@rm -rf $(BUILDDIR)
	@git add .
	@git commit -m "Update Sphinx and Doxygen docs"
	@git push origin ghPages
	@echo "Documentation successfully published to ghPages."
