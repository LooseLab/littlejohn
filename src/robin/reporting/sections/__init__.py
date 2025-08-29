"""
Report Sections Package

This package contains the individual sections that make up the ROBIN report.
"""

# IMPORTANT: Avoid importing heavy submodules at package import time to prevent
# side effects (e.g., matplotlib backend initialization) that can interfere with
# signal handling (Ctrl-C) in CLI/GUI contexts. Import submodules explicitly
# where needed, e.g., `from .classification import ClassificationSection`.

__all__ = []
