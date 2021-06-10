import pandas as pd
import pytest
import requests

from omigami import Spec2Vec
from omigami.spec2vec import InvalidCredentials


def test_build_payload(mgf_generator):
    client = Spec2Vec("token")

    payload = client._build_payload((mgf_generator), {"n_best_spectra": 10})

    assert "data" in payload.keys()
    assert payload["data"]["ndarray"]["parameters"]["n_best_spectra"] == 10
    assert payload["data"]["ndarray"]["data"][0]["Precursor_MZ"] == "240.115"


@pytest.mark.internet_connection
def test_send_request():
    client = Spec2Vec("bad_token")
    small_payload = {
        "data": {
            "ndarray": {
                "parameters": {"n_best_spectra": 10},
                "data": [
                    {
                        "peaks_json": "[[80.060677, 157.0], [337.508301, 230.0]]",
                        "Precursor_MZ": "153.233",
                    }
                ],
            }
        }
    }

    with pytest.raises(InvalidCredentials):
        client._send_request(small_payload)


def test_format_results(sample_response):
    client = Spec2Vec("token")
    requests.Response()

    results = client._format_results(sample_response)

    assert isinstance(results[0], pd.DataFrame)
    assert results[0].index.name == "matches of spectrum-0"
    assert all(results[0]["score"] > 0)


def test_validate_input():
    model_input = {
        "peaks_json": "[[80.060677, 157.0], [337.508301, 230.0]]",
        "Precursor_MZ": "153.233",
    }
    # first validates if the input is correct then we test for errors
    Spec2Vec._validate_input([model_input])

    with pytest.raises(TypeError, match="Spectrum data must be a dictionary."):
        Spec2Vec._validate_input(["not_a_dict"])

    with pytest.raises(KeyError, match="mandatory keys"):
        Spec2Vec._validate_input([{"Precursor_MZ": "1", "peaks_JASON": "[not a list]"}])

    with pytest.raises(
        ValueError, match="peaks_json needs to be a valid python string representation"
    ):
        Spec2Vec._validate_input([{"Precursor_MZ": "1", "peaks_json": "[not a list]"}])

    with pytest.raises(
        ValueError,
        match="peaks_json needs to be a valid python string representation",
    ):
        Spec2Vec._validate_input([{"Precursor_MZ": "1", "peaks_json": 10}])

    with pytest.raises(
        ValueError,
        match="Precursor_MZ needs to be a string representation of a float",
    ):
        Spec2Vec._validate_input([{"Precursor_MZ": "float", "peaks_json": [10]}])


def test_validate_parameters():
    parameters = Spec2Vec._build_parameters(2.0, ["smiles"])
    assert {"n_best_spectra", "include_metadata"} == set(parameters.keys())

    with pytest.raises(ValueError, match="batman"):
        _ = Spec2Vec._build_parameters(2, ["batman"])

    with pytest.raises(ValueError, match="must be an integer"):
        _ = Spec2Vec._build_parameters("robin", ["smiles"])
