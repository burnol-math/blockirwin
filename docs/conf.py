# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'block-Irwin'
copyright = '2026, Jean-François Burnol'
author = 'Jean-François Burnol'
release = ''

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# The master toctree document.
master_doc = 'index'


extensions = [
    'sphinx.ext.mathjax',
    'sphinx_math_dollar',
    # 'myst_parser',
]

# source_suffix = {
#     '.rst': 'restructuredtext',
#     '.md': 'markdown',
# }

# myst_enable_extensions = ["dollarmath"]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'

html_sidebars = {
    '**': [
        'about.html',
        'searchfield.html',
        'localtoc.html',
        # 'navigation.html',
        # 'relations.html',
        # 'donate.html',
    ]
}


html_static_path = ['_static']

pygments_style = 'friendly'
