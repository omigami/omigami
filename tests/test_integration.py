import pytest

from omigami_client.spec2vec import Spec2VecClient


@pytest.mark.internet_connection
@pytest.mark.skip(reason="Need to provide a valid token")
def test_match_spectra_from_path(mgf_path):
    client = Spec2VecClient("token")

    result = client.match_spectra_from_path(mgf_path)

    assert result
