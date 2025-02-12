# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'AlphaFold Analysis'
copyright = '2023, Samuel Murail'
author = 'Samuel Murail'
release = '0.1.4'

# If extensions (or modules to document with autodoc) are in another
# directory, add these directories to sys.path here. If the directory is
# relative to the documentation root, use os.path.abspath to make it
# absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('../../src/'))

# This command doesn't get any sense for me (previous one should be enough)
sys.path.insert(0, os.path.abspath('../../src/af_analysis/'))

print(sys.path)

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.coverage',
    'sphinx.ext.todo',
    'sphinx.ext.viewcode',
    'sphinx.ext.githubpages',
    'sphinx.ext.todo',
    'sphinxarg.ext',
    'sphinx.ext.mathjax',
    'numpydoc',
    'myst_parser',
    'nbsphinx',
    ]

templates_path = ['_templates']
exclude_patterns = []

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = ['.rst', '.md']


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme ='sphinx_rtd_theme'
html_static_path = ['_static']

master_doc = 'index'

man_pages = [
    (master_doc, 'af_analysis',
     'AF-Analysis Documentation',
     [author], 1)
]

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
# The master toctree document.


#autodoc_mock_imports = ["numpy", "scipy", "pytest"]

# Exclude unit pages (tests and data) from the documentation
exclude_patterns = ['af_analysis.test.*']