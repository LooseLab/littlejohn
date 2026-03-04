"""ROBIN and robin - A Python CLI tool with file watching capabilities."""

# Suppress pkg_resources deprecation warnings from sorted_nearest
import warnings
warnings.filterwarnings(
    "ignore", message="pkg_resources is deprecated", category=UserWarning
)
# Suppress matplotlib tight_layout warnings
warnings.filterwarnings(
    "ignore", message="The figure layout has changed to tight", category=UserWarning
)

__version__ = "0.4.1"
__author__ = "Matt Loose"
__email__ = "matt.loose@nottingham.ac.uk"

# Import analysis modules via consolidated analysis package
from . import analysis  # noqa: F401

# Export key modules
__all__ = [
    "analysis",
]
