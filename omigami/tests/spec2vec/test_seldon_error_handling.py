import pytest

from spectra_matching import Spec2Vec
from tests.conftest import _set_credentials_and_auth_for_tests


@pytest.mark.skip
def test_seldon_error_handling():
    _set_credentials_and_auth_for_tests()
    mgf_file_path = "" # set here testing mgf file path
    n_best_matches = 10
    ion_mode = "positive"

    client = Spec2Vec()

    client.match_spectra(
        mgf_file_path,
        n_best_matches,
        ion_mode=ion_mode,
    )
