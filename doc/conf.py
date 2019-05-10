# This source code is part of the Gecos package and is distributed
# under the 3-Clause BSD License. Please see 'LICENSE.rst' for further
# information.

__author__ = "Patrick Kunzmann"

from os.path import realpath, dirname, join, isdir, isfile, basename
from os import listdir, makedirs
import sys
import glob
import shutil
import matplotlib
from importlib import import_module
import types
import abc

absolute_path = dirname(realpath(__file__))
package_path = join(dirname(absolute_path), "src")
sys.path.insert(0, package_path)
import gecos


#### General ####

extensions = ["sphinx.ext.autodoc",
              "sphinx.ext.autosummary",
              "sphinx.ext.doctest",
              "sphinx.ext.mathjax",
              "sphinx.ext.viewcode",
              "numpydoc"]

templates_path = ["templates"]
source_suffix = [".rst"]
master_doc = "index"

project = "gecos"
copyright = "Patrick Kunzmann, 2019"
version = gecos.__version__

exclude_patterns = ["build"]

pygments_style = "sphinx"

todo_include_todos = False

# Prevents numpydoc from creating an autosummary which does not work
# due to Gecos' import system
numpydoc_show_class_members = False


#### HTML ####

html_theme = "alabaster"
html_static_path = ["static"]
html_favicon = "static/assets/general/gecos_icon_32p.png"
htmlhelp_basename = "GecosDoc"
html_sidebars = {"**": ["about.html",
                        #"localtoc.html",
                        "navigation.html",
                        "relations.html",
                        "searchbox.html",
                        "donate.html"]}
html_theme_options = {
    "description"   : "Generated color schemes",
    "logo"          : "assets/general/gecos_logo_s.png",
    "logo_name"     : "false",
    "github_user"   : "biotite-dev",
    "github_repo"   : "gecos",
    "github_banner" : "true",
    "page_width"    : "85%",
    "fixed_sidebar" : "true"
}