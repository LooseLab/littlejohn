"""LittleJohn - A Python CLI tool with file watching capabilities."""

__version__ = "0.1.0"
__author__ = "Your Name"
__email__ = "your.email@example.com"

# Import analysis modules
from . import bam_preprocessor
from . import bed_conversion
from . import nanodx_analysis
from . import fusion_analysis

# Import temp_utilities optionally (may fail if polars is not available)
try:
    from . import temp_utilities
except ImportError:
    # temp_utilities not available, but that's okay
    pass

# Export handler functions
__all__ = [
    'bam_preprocessor',
    'bed_conversion', 
    'nanodx_analysis',
    'fusion_analysis',
    'temp_utilities'
] 