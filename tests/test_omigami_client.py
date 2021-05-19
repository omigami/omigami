from pathlib import Path
from unittest.mock import Mock

import pytest

from omigami_client.omigami_client import OmigamiClient

ASSETS_DIR = Path(__file__).parent / "assets"


def test_match_spectra_from_path(mgf_path):
    client = OmigamiClient("token")
    client._build_payload = Mock(return_value="payload")
    client._send_request = Mock(return_value="request")
    client._format_results = Mock(return_value="result")

    result = client.match_spectra_from_path(mgf_path)

    assert result == "result"
    client._build_payload.assert_called_once()
    client._send_request.assert_called_once_with("payload")
    client._format_results.assert_called_once_with("request")


def test_build_payload(mgmf_generator):
    client = OmigamiClient("token")

    payload = client._build_payload((mgmf_generator), 10)

    assert "data" in payload.keys()
    assert payload["data"]["ndarray"]["parameters"]["n_best_spectra"] == 10
    assert payload["data"]["ndarray"]["data"][0]["Precursor_MZ"] == "240.115"


@pytest.mark.internet_connection
def test_send_request():
    client = OmigamiClient("bad_token")
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


def test_format_results():
    client = OmigamiClient("token")


def test_validate_input():
    pass
