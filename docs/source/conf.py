# Configuration file for the Sphinx documentation builder.
import sys
import os

project = 'WDMwavelet'
copyright = '2025, DrizzleatDusk'
author = 'DrizzleatDusk'
release = '0.0.1'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx.ext.todo'
]

templates_path = ['_templates']
exclude_patterns = []




html_theme = 'sphinx_book_theme'
html_static_path = ['_static']

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(os.path.dirname(current_dir))
target_dir = os.path.join(parent_dir, "wdmwavelet")
sys.path.insert(0, target_dir)

todo_include_todos = True
autosummary_generate = True
autosummary_imported_members = False
autodoc_member_order = 'bysource'
autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'inherited-members': True,
    'show-inheritance': True,
    'exclude-members': '__weakref__',
}