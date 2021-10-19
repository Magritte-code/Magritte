# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('.'))
import datetime
year = datetime.datetime.now().year

# -- Project information -----------------------------------------------------

project = 'Magritte'
copyright = f'{year}, Frederik De Ceuster'
author = 'Frederik De Ceuster'

# -- Explicitly specify master doc -------------------------------------------
master_doc = 'index'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [ "breathe", 'sphinx_copybutton', 'nbsphinx' ]

breathe_default_project = "Magritte"
breathe_projects = {'Magritte' : 'doxygen/xml'}


# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'
html_logo  = 'images/Magritte_logo_white_plain.svg'
html_theme_options = {'logo_only': True}


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']





# -- Run doxygen on READTHEDOCS server ---------------------------------------

from subprocess import call
# read_the_docs_build = os.environ.get('READTHEDOCS', None) == 'True'
call('doxygen Doxyfile', shell=True)


# Fix to show plotly figures
# https://github.com/readthedocs/sphinx_rtd_theme/issues/788#issuecomment-585785027
nbsphinx_prolog = r"""
.. raw:: html

    <script src='http://cdnjs.cloudflare.com/ajax/libs/require.js/2.1.10/require.min.js'></script>
    <script>require=requirejs;</script>


"""
