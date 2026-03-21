"""Sphinx configuration for the SIAC documentation site."""

from __future__ import annotations

project = "SIAC"
copyright = "2026, SIAC contributors"
author = "SIAC contributors"

extensions = [
    "myst_parser",
    "sphinxcontrib.mermaid",
]

templates_path = ["_templates"]
exclude_patterns = [
    "_build",
    "Thumbs.db",
    ".DS_Store",
    "EARTHACCESS_PLAN.md",
    "PLANS.md",
    "PLANS_CLOUD_MASK.md",
    "PLANS_S2.md",
    "REVIEW_AND_RECOMMENDATIONS.md",
    "TEST_PLAN.md",
    "code_structure.md",
    "config-redesign.md",
    "naming-conventions.md",
]

source_suffix = {
    ".md": "markdown",
}

myst_enable_extensions = [
    "colon_fence",
    "deflist",
]
myst_heading_anchors = 3
myst_fence_as_directive = ["mermaid"]

html_theme = "sphinx_rtd_theme"
html_title = "SIAC Documentation"
html_show_sourcelink = False
html_copy_source = False
