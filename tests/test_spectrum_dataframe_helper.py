from typing import List

import pytest
import pandas as pd
from matchms.importing import load_from_mgf

from omigami.utilities import SpectrumDataFrameHelper


@pytest.fixture
def spectra_dataframe(small_mgf_path) -> List[pd.DataFrame]:
    spectra = list(load_from_mgf(small_mgf_path))
    df = SpectrumDataFrameHelper.from_spectra_list(spectra)
    return df


def test_scale(spectra_dataframe):
    spectrum_df = spectra_dataframe[0]
    new_spectrum_df = SpectrumDataFrameHelper.scale(spectrum_df)

    assert spectrum_df.Intensity.max() > 1
    assert new_spectrum_df.Intensity.max() == 1
    assert new_spectrum_df.Intensity.min() == 0


def test_filter(spectra_dataframe):
    n_peaks = 10
    spectrum_df = spectra_dataframe[0]
    new_spectrum_df = SpectrumDataFrameHelper.filter(spectrum_df, n_peaks)

    assert spectrum_df.shape[0] > n_peaks
    assert new_spectrum_df.shape[0] == n_peaks


def test_crop():
    pass


def test_from_spectra_list():
    pass


def test_from_spectrum():
    pass


def test_from_gnps_id():
    pass


def test_from_gnps_id_list():
    pass
