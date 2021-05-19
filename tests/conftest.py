import pytest
from matchms.importing import load_from_mgf

from tests.test_omigami_client import ASSETS_DIR


@pytest.fixture(scope="session")
def mgf_path():
    return str(ASSETS_DIR / "GNPS-COLLECTIONS-MISC.mgf")


@pytest.fixture(scope="session")
def mgmf_generator(mgf_path):
    return list(load_from_mgf(mgf_path))[:20]
