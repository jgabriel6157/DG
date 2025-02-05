# Configuration file for the Sphinx documentation builder.

import os
import sys

# -- Project information -----------------------------------------------------
project = 'Guernica'
copyright = '2024, Jack Gabriel'
author = 'Jack Gabriel'

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',  # Generate documentation from docstrings
    'sphinx.ext.napoleon',  # Support for NumPy and Google docstrings
    'sphinx.ext.viewcode',  # Adds links to source code
    'sphinx.ext.githubpages',  # Enable GitHub Pages integration
    'breathe',  # Integrate Doxygen XML output
]

# Paths for templates and static files
templates_path = ['_templates']
html_static_path = ['_static']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------
html_theme = 'alabaster'

# -- Breathe Configuration (For Doxygen Integration) -------------------------
breathe_projects = {"Guernica": "doxygen_output/xml"}
breathe_default_project = "Guernica"

# -- Doxygen Setup -----------------------------------------------------------
# Automatically run Doxygen when Sphinx is built
def run_doxygen():
    """Run Doxygen to generate XML documentation."""
    print("Running Doxygen...")
    os.system("doxygen Doxyfile")

run_doxygen()
