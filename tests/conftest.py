import pickle
from pathlib import Path

import pytest
from matchms.importing import load_from_mgf

from omigami import Spec2Vec
from omigami.config import config

ASSETS_DIR = Path(__file__).parent / "assets"


@pytest.fixture(scope="module")
def spec2vec_client():
    token = config["login"]["dev"]["token"].get()
    client = Spec2Vec(token)
    client._endpoint_url = (
        "https://mlops.datarevenue.com/seldon/seldon/spec2vec/api/v0.1/predictions"
    )
    return client


# please download files below from https://gnps-external.ucsd.edu/gnpslibrary
# and place them inside the assets directory
@pytest.fixture(scope="session")
def small_mgf_path():
    return str(ASSETS_DIR / "gnps_small.mgf")


@pytest.fixture(scope="session")
def mgf_path():
    return str(ASSETS_DIR / "GNPS-COLLECTIONS-MISC.mgf")


@pytest.fixture(scope="session")
def mgf_huge_path():
    return str(ASSETS_DIR / "GNPS-NIST14-MATCHES.mgf")


@pytest.fixture(scope="session")
def mgf_generator(mgf_path):
    return list(load_from_mgf(mgf_path))[:20]


@pytest.fixture(scope="session")
def sample_response():
    with open(ASSETS_DIR / "sample_response.pickle", "rb") as f:
        response = pickle.load(f)

    return response
