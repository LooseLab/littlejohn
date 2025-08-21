"""LittleJohn - A Python CLI tool with file watching capabilities."""

__version__ = "0.1.1"
__author__ = "Your Name"
__email__ = "your.email@example.com"

# Import analysis modules via consolidated analysis package
from . import analysis  # noqa: F401

# Export key modules
__all__ = [
	"analysis",
]