import pytest

from omigami_client.omigami_client import OmigamiClient


@pytest.mark.internet_connection
@pytest.mark.skip(reason="Need to provide a valid token")
def test_match_spectra_from_path(mgf_path):
    client = OmigamiClient("5iTy7vACUXmlLO9fwGAL8v2WLPbo1SNH")

    result = client.match_spectra_from_path(mgf_path)

    assert result
