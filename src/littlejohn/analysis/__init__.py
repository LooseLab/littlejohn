"""Analysis package for LittleJohn.

This package consolidates all analytical modules under one namespace.
It re-exports the primary handler entry points for convenience.
"""

from .bam_preprocessor import bam_preprocessing_handler  # noqa: F401
from .cnv_analysis import cnv_handler  # noqa: F401
from .mgmt_analysis import mgmt_handler  # noqa: F401
from .fusion_analysis import fusion_handler  # noqa: F401
from .target_analysis import target_handler  # noqa: F401
from .sturgeon_analysis import sturgeon_handler  # noqa: F401
from .random_forest_analysis import random_forest_handler  # noqa: F401
from .nanodx_analysis import nanodx_handler, pannanodx_handler  # noqa: F401
from .bed_conversion import bed_conversion_handler  # noqa: F401
from .master_csv_manager import MasterCSVManager  # noqa: F401
from .temp_utilities import merge_modkit_files  # noqa: F401
from .methylation_wrapper import (
    locus_figure, save_figure_pickle, load_figure_pickle,  # noqa: F401
)

__all__ = [
	"bam_preprocessing_handler",
	"cnv_handler",
	"mgmt_handler",
	"fusion_handler",
	"target_handler",
	"sturgeon_handler",
	"random_forest_handler",
	"nanodx_handler",
	"pannanodx_handler",
    "bed_conversion_handler",
    "MasterCSVManager",
    "merge_modkit_files",
    "locus_figure",
    "save_figure_pickle",
    "load_figure_pickle",
]

