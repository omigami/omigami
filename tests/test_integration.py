import pytest


@pytest.mark.internet_connection
@pytest.mark.skip(reason="Requires a valid token")
def test_match_spectra_from_path_small(small_mgf_path, spec2vec_client):
    result = spec2vec_client.match_spectra_from_path(
        small_mgf_path, 10, ["smiles", "compound_name"]
    )

    assert result
    assert set(result[0].columns) == {"smiles", "score", "compound_name"}
    assert len(result) == 46


@pytest.mark.internet_connection
@pytest.mark.skip(reason="Requires a valid token")
def test_match_spectra_from_path(mgf_path, spec2vec_client):
    result = spec2vec_client.match_spectra_from_path(mgf_path, 10, ["smiles"])

    assert result
    assert len(result) == 377


@pytest.mark.internet_connection
@pytest.mark.skip(
    reason="this is a stress test, run only when you know what you are doing"
)
def test_match_spectra_from_path_with_huge_payload(mgf_huge_path, spec2vec_client):
    result = spec2vec_client.match_spectra_from_path(mgf_huge_path, 10)

    assert result
    assert len(result) == 5760
