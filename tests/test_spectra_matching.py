from io import StringIO
from unittest.mock import Mock

import pandas as pd
import pytest
import requests
from matchms import Spectrum
from matchms.importing import load_from_mgf

from omigami.authentication import AUTH
from omigami.exceptions import InvalidCredentials
from omigami.spectra_matching.spectra_matching import (
    SpectraMatching,
    SPECTRA_LIMIT_PER_REQUEST,
)


def test_build_payload(mgf_generator):
    client = SpectraMatching()

    payload = client._build_payload((mgf_generator), {"n_best_spectra": 10})

    assert "data" in payload.keys()
    assert payload["data"]["ndarray"]["parameters"]["n_best_spectra"] == 10
    assert payload["data"]["ndarray"]["data"][0]["Precursor_MZ"] == "240.115"


@pytest.mark.internet_connection
def test_unauthorized_request(spec2vec_prediction_endpoints):
    AUTH.token = ""
    client = SpectraMatching()
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
        client._send_request(small_payload, spec2vec_prediction_endpoints["positive"])


def test_format_results(sample_response):
    client = SpectraMatching()
    requests.Response()

    results = client._format_results(sample_response)

    assert isinstance(results[0], pd.DataFrame)
    assert results[0].index.name == "matches of spectrum-0"
    assert all(results[0]["score"] > 0)
    assert all(results[0].columns == ["score", "compound_name", "smiles"])


def test_validate_parameters():
    parameters = SpectraMatching._build_parameters(2)
    assert {"n_best_spectra"} == set(parameters.keys())

    with pytest.raises(ValueError, match="must be an integer"):
        _ = SpectraMatching._build_parameters("robin")


def test_validate_input():
    model_input = {
        "peaks_json": "[[80.060677, 157.0], [337.508301, 230.0]]",
        "Precursor_MZ": "153.233",
    }
    client = SpectraMatching()
    # first validates if the input is correct then we test for errors
    client._validate_input([model_input])

    with pytest.raises(TypeError, match="Spectrum data must be a dictionary."):
        client._validate_input(["not_a_dict"])

    with pytest.raises(KeyError, match="mandatory keys"):
        client._validate_input([{"Precursor_MZ": "1", "peaks_JASON": "[not a list]"}])

    with pytest.raises(
        ValueError, match="peaks_json needs to be a valid python string representation"
    ):
        client._validate_input([{"Precursor_MZ": "1", "peaks_json": "[not a list]"}])

    with pytest.raises(
        ValueError,
        match="peaks_json needs to be a valid python string representation",
    ):
        client._validate_input([{"Precursor_MZ": "1", "peaks_json": 10}])

    with pytest.raises(
        ValueError,
        match="Precursor_MZ needs to be a string representation of a float",
    ):
        client._validate_input([{"Precursor_MZ": "float", "peaks_json": [10]}])


def test_create_spectra_generator_file_source(small_mgf_path):
    client = SpectraMatching()

    spectra_generator = client._create_spectra_generator(small_mgf_path)

    assert isinstance(spectra_generator.__next__(), Spectrum)


def test_create_spectra_generator_stream_source(small_mgf_path):
    client = SpectraMatching()
    stream = StringIO(open(small_mgf_path, "r").read())

    spectra_generator = client._create_spectra_generator(stream)

    assert isinstance(spectra_generator.__next__(), Spectrum)


def test_create_spectra_generator_spectrum_source(small_mgf_path):
    client = SpectraMatching()
    spectra = list(load_from_mgf(small_mgf_path))

    spectra_generator = client._create_spectra_generator(spectra)

    assert isinstance(spectra_generator.__next__(), Spectrum)


@pytest.fixture
def spectra_matching_client(mock_send_request):
    client = SpectraMatching(optional_token="token")
    client._ENDPOINT = "endpoint"
    client._send_request = mock_send_request

    return client


def test_make_batch_requests(small_mgf_path, spectra_matching_client):
    spectra_matching_client._cache_results = Mock()
    _46_spectra = list(load_from_mgf(small_mgf_path))

    # Input has 260 spectra -> 3 batch requests of [100, 100, 60]
    input_spectra = _46_spectra * 5
    n_spectra = len(input_spectra)
    n_expected_batches = n_spectra // SPECTRA_LIMIT_PER_REQUEST + 1
    spectra_generator = spectra_matching_client._create_spectra_generator(input_spectra)

    predictions = spectra_matching_client._make_batch_requests(
        spectra_generator, {"n_best_spectra": 2}, "endpoint"
    )

    assert len(predictions) == n_spectra
    assert spectra_matching_client._send_request.call_count == n_expected_batches
    assert spectra_matching_client._cache_results.call_count == n_expected_batches


def test_spectra_matching_caching(spectra_matching_client, small_mgf_path):
    """Tests result caching by running the same function twice"""

    input_spectra = list(load_from_mgf(small_mgf_path))

    spectra = spectra_matching_client.match_spectra(input_spectra, n_best=2)

    assert len(spectra) == 46
    assert len(spectra_matching_client._cached_results) == len(spectra)
    assert spectra_matching_client._send_request.call_count == 1

    _ = spectra_matching_client.match_spectra(input_spectra, n_best=2)

    assert spectra_matching_client._send_request.call_count == 1


def test_reset_cache(spectra_matching_client, small_mgf_path):
    """Resets caches after running the first time"""

    input_spectra = list(load_from_mgf(small_mgf_path))

    spectra = spectra_matching_client.match_spectra(input_spectra, n_best=2)

    assert len(spectra) == 46
    assert len(spectra_matching_client._cached_results) == len(spectra)
    assert spectra_matching_client._send_request.call_count == 1

    spectra_matching_client.reset_cache()
    _ = spectra_matching_client.match_spectra(input_spectra, n_best=2)

    assert spectra_matching_client._send_request.call_count == 2


def test_failed_spectra_caching(spectra_matching_client, small_mgf_path):
    #  TODO
    pass
