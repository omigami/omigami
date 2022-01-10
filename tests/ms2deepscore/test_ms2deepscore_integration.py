import pytest


@pytest.mark.internet_connection
@pytest.mark.skip(reason="Requires valid credentials")
def test_match_spectra_from_path_with_2_spectra(
    mgf_path_of_2_spectra, ms2deepscore_client
):
    result = ms2deepscore_client.match_spectra(
        source=mgf_path_of_2_spectra, n_best=3, ion_mode="positive"
    )
    assert result
    assert "score" in result[0].columns
    assert result[0].index[0] == "CCMSLIB00000001547"
    assert result[1].index[0] == "CCMSLIB00000001548"


@pytest.mark.internet_connection
@pytest.mark.skip(reason="Requires valid credentials")
def test_match_spectra_from_path_with_small(mgf_46_spectra_path, ms2deepscore_client):
    n_best = 3

    result = ms2deepscore_client.match_spectra(
        source=mgf_46_spectra_path, n_best=n_best, ion_mode="positive"
    )
    assert len(result) == 46
    assert len(result[0]) == n_best
