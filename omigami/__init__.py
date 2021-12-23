# -*- coding: utf-8 -*-

"""Top-level package for Omigami."""

__author__ = """Data Revenue GmbH"""
__email__ = "carlos@datarevenue.com"
__version__ = "0.1.0"

from pathlib import Path

from .plotting import MoleculePlotter, SpectraComparisonPlotter
from ._version import get_versions
from .utilities import GnpsHelper
from .utilities import SpectrumDataFrameHelper

__version__ = get_versions()["version"]
del get_versions

ROOT_DIR = Path(__file__).parent.parent
