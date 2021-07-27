import pytest


@pytest.mark.internet_connection
@pytest.mark.skip(reason="Requires a valid token")
def test_match_spectra_from_path_with_2_spectra(mgf_path_of_pair, ms2deepscore_client):
    result = ms2deepscore_client.match_spectra_from_path(
        mgf_path_of_pair, 3, ["compound_name"]
    )
    assert result
    assert set(result[0].columns) == {"score", "compound_name"}


# TODO: add more tests
