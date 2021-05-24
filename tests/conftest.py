import pickle
import yaml
from pathlib import Path

import pytest
from matchms.importing import load_from_mgf

from omigami_client import ROOT_DIR

ASSETS_DIR = Path(__file__).parent / "assets"
DEV_ENV_VARS_PATH = ROOT_DIR / "dev.env"


@pytest.fixture(scope="package")
def token():
    with open(DEV_ENV_VARS_PATH) as yaml_dev_vars_file:
        vars = yaml.safe_load(yaml_dev_vars_file)
        return vars["token"]


# please download files below from https://gnps-external.ucsd.edu/gnpslibrary
# and place them inside the assets directory
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
