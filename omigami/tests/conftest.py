import pickle
import os
from pathlib import Path

import pytest
from matchms.importing import load_from_mgf

from omigami import Spec2Vec
from omigami.authentication import encrypt_credentials, AUTH
from omigami.omi_settings import config
from omigami.ms2deepscore import MS2DeepScore

ASSETS_DIR = Path(__file__).parent / "assets"


def _set_credentials_and_auth_for_tests():
    """
    This uses the default dev credentials "omigami@dev.org" and setup necessary variable values for testing
    """
    username = os.getenv("OMIGAMI_USERNAME") or config["login"]["dev"]["username"].get()
    pwd = os.getenv("OMIGAMI_PWD") or config["login"]["dev"]["password"].get()
    AUTH.credentials = encrypt_credentials(username, pwd)
    AUTH.self_service_endpoint = (
        "https://mlops.datarevenue.com/.ory/kratos/public/self-service/login/api"
    )


@pytest.fixture(scope="module")
def spec2vec_client():
    _set_credentials_and_auth_for_tests()
    client = Spec2Vec()
    client._PREDICT_ENDPOINT_BASE = "https://mlops.datarevenue.com/seldon/seldon/spec2vec-{ion_mode}/api/v0.1/predictions"
    return client


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


@pytest.fixture(scope="module")
def ms2deepscore_client():
    _set_credentials_and_auth_for_tests()
    client = MS2DeepScore()
    client._PREDICT_ENDPOINT_BASE = "https://mlops.datarevenue.com/seldon/seldon/ms2deepscore-{ion_mode}/api/v0.1/predictions"
    return client


@pytest.fixture()
def ms2deepscore_prediction_endpoints():
    _client = Spec2Vec()
    return {
        "positive": _client._PREDICT_ENDPOINT_BASE.format(ion_mode="positive"),
        "negative": _client._PREDICT_ENDPOINT_BASE.format(ion_mode="negative"),
    }


@pytest.fixture(scope="session")
def mgf_path_of_2_spectra():
    return str(ASSETS_DIR / "gnps_2_spectra.mgf")


@pytest.fixture()
def spec2vec_prediction_endpoints():
    _client = Spec2Vec()
    return {
        "positive": _client._PREDICT_ENDPOINT_BASE.format(ion_mode="positive"),
        "negative": _client._PREDICT_ENDPOINT_BASE.format(ion_mode="negative"),
    }


@pytest.fixture(scope="session")
def spectra_match_data_path():
    return str(ASSETS_DIR / "spectrum_matches.csv")


@pytest.fixture(scope="session")
def spectra_match_data_path_missing_smiles():
    return str(ASSETS_DIR / "spectrum_missing_smiles.csv")


@pytest.fixture(scope="session")
def spectra_match_data_path_web_api_error():
    return str(ASSETS_DIR / "spectrum_matches_error_on_classyfire.csv")
