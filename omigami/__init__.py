# -*- coding: utf-8 -*-

"""Top-level package for Omigami."""

__author__ = """Data Revenue GmbH"""
__email__ = "carlos@datarevenue.com"
__version__ = "0.1.0"

from pathlib import Path

from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions
from .spectra_matching import MS2DeepScore, Spec2Vec
from omigami.plotting import MoleculePlotter

ROOT_DIR = Path(__file__).parent.parent
