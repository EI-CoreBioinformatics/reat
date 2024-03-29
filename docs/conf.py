# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
import subprocess

from importlib import resources as pkg_resources

from annotation import VERSION
import sphinx_rtd_theme

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'REAT'
copyright = '2021, Earlham Institute'
author = 'Luis Yanes'
version = VERSION
release = version
# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx_rtd_theme'
]

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

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']


# Generate a full CLI help for the transcriptome command
transcriptome_help = subprocess.run(['reat', 'transcriptome', '--help'], capture_output=True)
print(transcriptome_help.stdout.decode(), file=open('modules/transcriptome/transcriptome_help.txt', 'w'))

# Generate a full CLI help for the homology command
homology_help = subprocess.run(['reat', 'homology', '--help'], capture_output=True)
print(homology_help.stdout.decode(), file=open('modules/homology/homology_help.txt', 'w'))

# Generate a full CLI help for the prediction command
prediction_help = subprocess.run(['reat', 'prediction', '--help'], capture_output=True)
print(prediction_help.stdout.decode(), file=open('modules/prediction/prediction_help.txt', 'w'))

with pkg_resources.path("annotation.prediction_module", "extrinsic.ei_augustus_generic.cfg") as extrinsic_path:
    with open("modules/prediction/extrinsic.ei_augustus_generic.cfg", 'w') as extrinsic_copy:
        contents = open(extrinsic_path).readlines()
        print(''.join(contents), file=extrinsic_copy)

with pkg_resources.path("annotation.prediction_module", "evm_default_weights.wgt") as evm_path:
    with open("modules/prediction/evm_default_weights.wgt", 'w') as evm_copy:
        contents = open(evm_path).readlines()
        print(''.join(contents), file=evm_copy)
