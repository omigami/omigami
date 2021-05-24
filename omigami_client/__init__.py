# -*- coding: utf-8 -*-

"""Top-level package for Omigami Client."""

__author__ = """Data Revenue GmbH"""
__email__ = "carlos@datarevenue.com"
__version__ = "0.1.0"

from pathlib import Path

from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions
from omigami_client.spec2vec import Spec2VecClient

ROOT_DIR = Path(__file__).parent.parent
