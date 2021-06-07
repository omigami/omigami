import pytest

from omigami.spec2vec import Spec2Vec


@pytest.mark.internet_connection
@pytest.mark.skip(reason="Requires a valid token")
def test_match_spectra_from_path(mgf_path, token):
    client = Spec2Vec(token)

    result = client.match_spectra_from_path(mgf_path, 10)

    assert result
    assert len(result) == 377


@pytest.mark.internet_connection
@pytest.mark.skip(
    reason="this is a stress test, run only when you know what you are doing"
)
def test_match_spectra_from_path_with_huge_payload(mgf_huge_path, token):
    client = Spec2Vec(token)

    result = client.match_spectra_from_path(mgf_huge_path, 10)

    assert result
    assert len(result) == 5760
