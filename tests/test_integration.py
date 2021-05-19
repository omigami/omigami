from omigami_client.omigami_client import OmigamiClient


def test_match_spectra_from_path(mgmf_path):
    client = OmigamiClient("5iTy7vACUXmlLO9fwGAL8v2WLPbo1SNH")

    result = client.match_spectra_from_path(mgmf_path)

    assert result