import pytest
from matchms.importing import load_from_mgf

from tests.test_omigami_client import ASSETS_DIR


@pytest.fixture(scope="session")
def mgmf_path():
    return str(ASSETS_DIR / "GNPS-COLLECTIONS-MISC.mgf")


@pytest.fixture(scope="session")
def mgmf_generator(mgmf_path):
    return load_from_mgf(mgf_path)