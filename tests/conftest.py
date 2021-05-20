import pickle
from pathlib import Path

import pytest
from matchms.importing import load_from_mgf

ASSETS_DIR = Path(__file__).parent / "assets"


@pytest.fixture(scope="session")
def mgf_path():
    return str(ASSETS_DIR / "GNPS-COLLECTIONS-MISC.mgf")


@pytest.fixture(scope="session")
def mgf_generator(mgf_path):
    return list(load_from_mgf(mgf_path))[:20]


@pytest.fixture(scope="session")
def sample_response():
    with open(ASSETS_DIR / "sample_response.pickle", "rb") as f:
        response = pickle.load(f)

    return response
