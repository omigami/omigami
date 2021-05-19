from pathlib import Path
from unittest.mock import MagicMock, Mock

from omigami_client.omigami_client import OmigamiClient

ASSETS_DIR = Path(__file__).parent / "assets"


def test_match_spectra_from_path(mgmf_path):
    client = OmigamiClient("token")
    client._build_payload = Mock(return_value="payload")
    client._send_request = Mock(return_value="request")
    client._format_results = Mock(return_value="result")

    result = client.match_spectra_from_path(mgmf_path)

    assert result == "result"
    client._build_payload.assert_called_once()
    client._send_request.assert_called_once_with("payload")
    client._format_results.assert_called_once_with("request")


def test_send_request():
    pass


def test_build_payload():
    pass


def test_format_results():
    pass


def test_validate_input():
    pass
