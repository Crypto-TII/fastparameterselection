# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import sys
import os
project = 'Fast Parameter Selection'
copyright = '2025, Beatrice Biasioli, Elena Kirshanova, Chiara Marcolla, Sergi Rovira'
author = 'Beatrice Biasioli, Elena Kirshanova, Chiara Marcolla, Sergi Rovira'
release = '1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',          # Automatically generate documentation from docstrings
    'sphinx.ext.napoleon',         # Support for NumPy and Google style docstrings
    'sphinx.ext.viewcode',         # Add links to source code
    'sphinx.ext.mathjax',  # Enable MathJax for LaTeX math rendering
]

templates_path = ['_templates']
exclude_patterns = []

# Add the root directory of the project to the Python path
# Adjust the path to your project root
# Adjust the path to point to your src directory
sys.path.insert(0, os.path.abspath('../src'))

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'  # Use the Read the Docs theme
html_static_path = ['_static']

# -- Napoleon settings -------------------------------------------------------
# Configure Napoleon for better docstring parsing
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True

# -- Autodoc settings --------------------------------------------------------
autodoc_default_options = {
    # Include all members (functions, classes, etc.)
    'members': True,
    'undoc-members': True,         # Include members without docstrings
    'show-inheritance': True,      # Show class inheritance
}
# Show type hints in the description instead of the signature
autodoc_typehints = 'description'
