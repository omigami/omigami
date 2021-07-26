import pickle
from pathlib import Path

import pytest
from matchms.importing import load_from_mgf

from omigami import Spec2Vec
from omigami.config import config
from omigami.ms2deepscore import MS2DeepScore

ASSETS_DIR = Path(__file__).parent / "assets"


@pytest.fixture(scope="module")
def spec2vec_client():
    token = config["login"]["dev"]["token"].get()
    client = Spec2Vec(token)
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
    token = config["login"]["dev"]["token"].get()
    client = MS2DeepScore(token)
    client._endpoint_url = (
        "https://mlops.datarevenue.com/seldon/seldon/ms2deepscore/api/v0.1/predictions"
    )
    return client


@pytest.fixture(scope="session")
def mgf_path_of_pair():
    return str(ASSETS_DIR / "gnps_2_spectra.mgf")


@pytest.fixture()
def spec2vec_prediction_endpoints():
    _client = Spec2Vec("")
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
