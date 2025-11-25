# SPDX-FileCopyrightText: 2023-present Matt Loose <matt.loose@nottingham.ac.uk>
#
# SPDX-License-Identifier: MIT

from pathlib import Path

# Directory containing ROBIN model files
# This allows code to access the models directory via: from robin import models; models.DIR
DIR = Path(__file__).parent
