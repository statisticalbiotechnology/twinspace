#!/usr/bin/env python
#
# msci documentation build configuration file, created by
# sphinx-quickstart on Fri Jun  9 13:47:02 2017.

import os
import sys
sys.path.insert(0, os.path.abspath('..'))
import sphinx_rtd_theme
import MSCI

# -- General configuration ---------------------------------------------

extensions = ['sphinx.ext.autodoc', 'sphinx.ext.viewcode','sphinx_rtd_theme']

templates_path = ['_templates']

source_suffix = '.rst'

master_doc = 'index'

project = 'MSCI'
copyright = "2023, Zahra ELHAMRAOUI"
author = "Zahra ELHAMRAOUI"

version = MSCI.__version__
release = MSCI.__version__

language = 'en'  # Set this to a valid language code

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

pygments_style = 'sphinx'

todo_include_todos = False

# -- Options for HTML output -------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

html_logo = '_static/MSCI_logo-rem.png'

# Ensure the _static directory exists or comment out this line
html_static_path = ['_static']

# -- Options for HTMLHelp output ---------------------------------------

htmlhelp_basename = 'mscidoc'

# -- Options for LaTeX output ------------------------------------------

latex_elements = {}

latex_documents = [
    (master_doc, 'msci.tex',
     'MSCI Documentation',
     'Zahra ELHAMRAOUI', 'manual'),
]

# -- Options for manual page output ------------------------------------

man_pages = [
    (master_doc, 'msci',
     'MSCI Documentation',
     [author], 1)
]

# -- Options for Texinfo output ----------------------------------------

texinfo_documents = [
    (master_doc, 'msci',
     'MSCI Documentation',
     author,
     'msci',
     'One line description of project.',
     'Miscellaneous'),
]
