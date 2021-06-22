import pytest
import requests

from omigami.ms2deep import MS2Deep, InvalidNumberOfSpectra
from omigami.spec2vec import InvalidCredentials


def test_build_payload(mgf_generator):
    client = MS2Deep("token")

    payload = client._build_payload((mgf_generator[:2]))

    assert "data" in payload.keys()
    assert len(payload["data"]["ndarray"]["data"]) == 2


@pytest.mark.skip(reason="This test should only work when the endpoint is working")
@pytest.mark.internet_connection
def test_send_request():
    client = MS2Deep("bad_token")
    small_payload = {
        "data": {
            "ndarray": {
                "parameters": [],
                "data": [
                    {"intensities": "[80.060677, 337.508301]", "mz": "[157.0, 230.0]"},
                    {"intensities": "[81.060677, 339.508301]", "mz": "[158.0, 240.0]"},
                ],
            }
        }
    }

    with pytest.raises(InvalidCredentials):
        client._send_request(small_payload)


@pytest.mark.skip(
    reason="This test should only work when we can have a sample response of the "
    "endpoint "
)
def test_format_results(sample_response):
    client = MS2Deep("token")
    requests.Response()

    results = client._format_results(sample_response)

    assert isinstance(results, float)


def test_validate_input():
    model_input = [
        {"intensities": "[80.060677, 337.508301]", "mz": "[157.0, 230.0]"},
        {"intensities": "[81.060677, 339.508301]", "mz": "[158.0, 240.0]"},
    ]
    # first validates if the input is correct then we test for errors
    MS2Deep._validate_input(model_input)

    with pytest.raises(
        InvalidNumberOfSpectra,
        match=f"MS2DeepScore can only predict the Tanimoto "
        f"score for a pair of spectra. Please input "
        f"two spectra at a time. Found 1 spectra.",
    ):
        MS2Deep._validate_input(["only_one_spectrum"])

    with pytest.raises(TypeError, match="Spectrum data must be a dictionary."):
        MS2Deep._validate_input(["not_a_dict", "also_not_a_dict"])

    with pytest.raises(KeyError, match="mandatory keys"):
        MS2Deep._validate_input(
            [{"intensi": "[not a list]"}, {"intensities": [], "mz": []}]
        )

    with pytest.raises(
        ValueError, match="intensities needs to be a valid python string representation"
    ):
        MS2Deep._validate_input(
            [
                {"intensities": "[not a list]", "mz": []},
                {"intensities": "[not a list]", "mz": []},
            ]
        )

    with pytest.raises(
        ValueError,
        match="mz needs to be a valid python string representation",
    ):
        MS2Deep._validate_input(
            [{"intensities": [], "mz": 10}, {"intensities": [], "mz": 10}]
        )
