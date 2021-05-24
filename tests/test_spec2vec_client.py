from unittest.mock import Mock

import pandas as pd
import pytest
import requests

from omigami_client import Spec2VecClient


def test_match_spectra_from_path_calls(mgf_path):
    client = Spec2VecClient("token")
    client._build_payload = Mock(return_value="payload")
    client._send_request = Mock(return_value="request")
    client._format_results = Mock(return_value="result")

    result = client.match_spectra_from_path(mgf_path, 10)

    assert result == "result"
    client._build_payload.assert_called_once()
    client._send_request.assert_called_once_with("payload")
    client._format_results.assert_called_once_with("request")


def test_build_payload(mgf_generator):
    client = Spec2VecClient("token")

    payload = client._build_payload((mgf_generator), 10)

    assert "data" in payload.keys()
    assert payload["data"]["ndarray"]["parameters"]["n_best_spectra"] == 10
    assert payload["data"]["ndarray"]["data"][0]["Precursor_MZ"] == "240.115"


@pytest.mark.internet_connection
def test_send_request():
    client = Spec2VecClient("bad_token")
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

    response = client._send_request(small_payload)

    assert response.status_code == 401


def test_format_results(sample_response):
    client = Spec2VecClient("token")
    requests.Response()

    results = client._format_results(sample_response)

    assert isinstance(results[0], pd.DataFrame)
    assert results[0].index.name == "matches of spectrum #1"
    assert all(results[0].values > 0)


def test_validate_input():
    model_input = {
        "peaks_json": "[[80.060677, 157.0], [337.508301, 230.0]]",
        "Precursor_MZ": "153.233",
    }
    # first validates if the input is correct then we test for errors
    Spec2VecClient._validate_input([model_input])

    with pytest.raises(TypeError, match="Spectrum data must be a dictionary."):
        Spec2VecClient._validate_input(["not_a_dict"])

    with pytest.raises(KeyError, match="mandatory keys"):
        Spec2VecClient._validate_input(
            [{"Precursor_MZ": "1", "peaks_JASON": "[not a list]"}]
        )

    with pytest.raises(
        ValueError, match="peaks_json needs to be a valid python string representation"
    ):
        Spec2VecClient._validate_input(
            [{"Precursor_MZ": "1", "peaks_json": "[not a list]"}]
        )

    with pytest.raises(
        ValueError,
        match="peaks_json needs to be a valid python string representation",
    ):
        Spec2VecClient._validate_input([{"Precursor_MZ": "1", "peaks_json": 10}])

    with pytest.raises(
        ValueError,
        match="Precursor_MZ needs to be a string representation of a float",
    ):
        Spec2VecClient._validate_input([{"Precursor_MZ": "float", "peaks_json": [10]}])
