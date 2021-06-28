import pytest


@pytest.mark.internet_connection
@pytest.mark.skip(reason="Requires a valid token")
def test_predict_similarity_of_pair(mgf_path_of_pair, ms2deepscore_client):
    result = ms2deepscore_client.predict_similarity_of_pair(
        mgf_path_of_pair,
    )
    assert result


@pytest.mark.internet_connection
@pytest.mark.skip(reason="Requires a valid token")
def test_match_spectra_from_path(mgf_path_of_equal_pair, ms2deepscore_client):
    result = ms2deepscore_client.predict_similarity_of_pair(mgf_path_of_equal_pair)

    assert result
    assert result["Tanimoto Score"] == 1
