# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
import inspect

sys.path.insert(0, os.path.abspath(".."))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "spheremap"
copyright = "2023, Fermi Research Alliance, LLC."
author = "Eric Neilsen"
release = "v0.1.0"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.doctest",
    "sphinx.ext.linkcode",
    "sphinx.ext.autosummary",
    "sphinx-prompt",
    "numpydoc",
]

autosummary_generate = True
source_suffix = ".rst"
master_doc = "source/index"

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "classic"
html_static_path = ["_static"]


def linkcode_resolve(domain, info):
    if (domain == "py") and info["module"]:
        filename = info["module"].replace(".", "/")

        # Traverse the nest of objects until we get to the one we
        # are documenting.
        element = sys.modules.get(info["module"])
        for element_name in info["fullname"].split("."):
            parent_element = element
            element = getattr(element, element_name)

        if isinstance(element, property):
            # findsource cannot find the line number for properties,
            # so get the line number for the object of which it is a
            # property instead.
            element = parent_element

        _, line = inspect.findsource(element)
        url = f"https://github.com/lsst/spheremap/blob/tickets/PREOPS-3447/{filename}.py#L{line}"
    else:
        url = None

    return url
