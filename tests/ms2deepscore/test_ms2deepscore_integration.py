import pytest


@pytest.mark.internet_connection
@pytest.mark.skip(reason="Requires valid credentials")
def test_match_spectra_from_path_with_2_spectra(
    mgf_path_of_2_spectra, ms2deepscore_client
):
    n_best = 3
    result = ms2deepscore_client.match_spectra(
        mgf_path_of_2_spectra, n_best, "positive"
    )

    assert result
    assert len(result) == 2
    assert len(result[0]) == n_best


@pytest.mark.internet_connection
@pytest.mark.skip(reason="Requires valid credentials")
def test_match_spectra_from_path_with_small(mgf_46_spectra_path, ms2deepscore_client):
    n_best = 3

    result = ms2deepscore_client.match_spectra(mgf_46_spectra_path, n_best, "positive")
    assert len(result) > 0
    assert len(result[0]) == n_best
