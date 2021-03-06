import pytest


@pytest.mark.internet_connection
@pytest.mark.skip(reason="Requires valid credentials")
def test_match_spectra_from_path_with_2_spectra(
    mgf_path_of_2_spectra, ms2deepscore_client
):
    result = ms2deepscore_client.match_spectra_from_path(
        mgf_path_of_2_spectra, 3, ["compound_name"], "positive"
    )
    assert result
    assert set(result[0].columns) == {"score", "compound_name"}
    assert result[0].index[0] == "CCMSLIB00000001547"
    assert result[1].index[0] == "CCMSLIB00000001548"


@pytest.mark.internet_connection
@pytest.mark.skip(reason="Requires valid credentials")
def test_match_spectra_from_path_with_small(small_mgf_path, ms2deepscore_client):
    result = ms2deepscore_client.match_spectra_from_path(
        small_mgf_path, 3, ["compound_name"], "positive"
    )
    assert result
    assert len(result) == 46
