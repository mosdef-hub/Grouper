import os
import sys

# -- Path setup --------------------------------------------------------------
sys.path.insert(0, os.path.abspath(".."))

# -- General configuration ---------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx_rtd_theme",
    # 'breathe',
    # comment out exhale if you only want the bindings.cpp docs, not full C++ tree
    # 'exhale',
    "sphinx.ext.intersphinx",
    "sphinx_autodoc_typehints",
]

# Mock heavy or unavailable imports (so autodoc doesn’t fail)
autodoc_mock_imports = [
    "numpy",
    "pandas",
    "pyarrow",
    "igraph",
    "networkx",
    "cairosvg",
    "matplotlib",
    "rdkit",
    "mbuild",
    "cairo",
    "torch",
    "torch_geometric",
    "sklearn",
    "seaborn",
    "scipy",  # <— prevent scipy/seaborn import errors
    "Grouper.tests",  # <— never import test modules
    # 'Grouper._Grouper'       # <— uncomment if pybind11 module not built
]

templates_path = ["_templates"]

# Exclude build artifacts, OS junk, and all test directories
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "Grouper/tests",
    "Grouper/tests/**",
    "tests",
    "**/tests/**",
    "*/__pycache__/*",
    "modules.rst",
    "Grouper.rst",
]

# -- Options for HTML output -------------------------------------------------
html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

# -- Breathe / Exhale configuration for C++ API (Pybind11 bindings) ----------
# The XML from Doxygen should live at docs/_doxygen/xml
# breathe_projects = {"Grouper": os.path.abspath("_doxygen/xml")}
# breathe_default_project = "Grouper"


# Cross-linking with other docs
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
}

# Optional: if your extension isn’t built, you can mock it dynamically
# try:
#     import Grouper._Grouper
# except ImportError:
#     autodoc_mock_imports += ['Grouper._Grouper']
