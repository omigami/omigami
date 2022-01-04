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


def test_crop(spectra_dataframe):
    max_mz = 100
    spectrum_df = spectra_dataframe[0]
    new_spectrum_df = SpectrumDataFrameHelper.crop(spectrum_df, (0, max_mz))

    assert spectrum_df["m/z"].max() > max_mz
    assert new_spectrum_df["m/z"].max() <= max_mz


def test_from_spectra_list(mgf_path):
    spectra = list(load_from_mgf(mgf_path))
    spectra_df = SpectrumDataFrameHelper.from_spectra_list(spectra)

    assert type(spectra_df) == list
    assert isinstance(spectra_df[0], pd.DataFrame)
    assert {"Intensity", "m/z"}.issubset(set(spectra_df[0].columns))


def test_from_spectrum(mgf_path):
    spectra = list(load_from_mgf(mgf_path))
    spectrum_df = SpectrumDataFrameHelper.from_spectrum(spectra[0])

    assert isinstance(spectrum_df, pd.DataFrame)
    assert {"Intensity", "m/z"}.issubset(set(spectrum_df.columns))


@pytest.mark.internet_connection
def test_from_gnps_id():
    spectrum_df = SpectrumDataFrameHelper.from_gnps_id("CCMSLIB00000078899")

    assert isinstance(spectrum_df, pd.DataFrame)
    assert {"Intensity", "m/z"}.issubset(set(spectrum_df.columns))
    assert spectrum_df.shape[0] > 0


@pytest.mark.internet_connection
def test_from_gnps_id_list():
    spectra_df = SpectrumDataFrameHelper.from_gnps_id_list(
        ["CCMSLIB00000078899", "CCMSLIB00000078904"]
    )

    assert type(spectra_df) == list
    assert len(spectra_df) == 2
    assert all(isinstance(spectrum_df, pd.DataFrame) for spectrum_df in spectra_df)
    assert all(
        {"Intensity", "m/z"}.issubset(set(spectrum_df.columns))
        for spectrum_df in spectra_df
    )
    assert all(spectrum_df.shape[0] > 0 for spectrum_df in spectra_df)
