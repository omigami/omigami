import pytest
from matchms.importing import load_from_mgf


@pytest.mark.internet_connection
@pytest.mark.skip(reason="Requires valid credentials")
def test_match_spectra_from_path_small(mgf_46_spectra_path, spec2vec_client):
    """
    Tests matching spectra against library with a very small ammount of spectra ( < 50 )
    """

    result = spec2vec_client.match_spectra(mgf_46_spectra_path, 10, ion_mode="positive")

    assert result
    assert len(result[0].columns) == 38  # 38 metadata fields from gnps
    assert len(result) == 46


@pytest.mark.internet_connection
@pytest.mark.skip(reason="Requires valid credentials")
def test_match_spectra_from_list(mgf_46_spectra_path, spec2vec_client):
    """
    Tests matching spectra against library with a very small ammount of spectra ( < 50 )
    """

    spectra = list(load_from_mgf(mgf_46_spectra_path))

    result = spec2vec_client.match_spectra(spectra, 10, ion_mode="positive")

    assert result
    assert len(result) == 46


@pytest.mark.internet_connection
@pytest.mark.skip(reason="Requires valid credentials")
def test_match_spectra_from_path(mgf_377_spectra_path, spec2vec_client):
    """
    Tests matching spectra against library with a more substantial ammount of spectra ( > 350 )
    """
    result = spec2vec_client.match_spectra(mgf_377_spectra_path, 10, "positive")

    assert result
    assert len(result) == 377


@pytest.mark.internet_connection
@pytest.mark.skip(reason="Requires valid credentials")
def test_match_spectra_from_path_negative_mode(mgf_377_spectra_path, spec2vec_client):
    result = spec2vec_client.match_spectra(mgf_377_spectra_path, 10, "negative")

    assert result
    assert len(result) == 377


@pytest.mark.internet_connection
@pytest.mark.skip(
    reason="this is a stress test, run only when you know what you are doing"
)
def test_match_spectra_from_path_with_huge_payload(mgf_huge_path, spec2vec_client):
    result = spec2vec_client.match_spectra(mgf_huge_path, 10)

    assert result
    assert len(result) == 5760
