from typing import List

import pytest
import pandas as pd
from matchms.importing import load_from_mgf

from omigami.utilities import SpectrumDataFrameHelper


@pytest.fixture
def spectra_dataframe(mgf_46_spectra_path) -> List[pd.DataFrame]:
    spectra = list(load_from_mgf(mgf_46_spectra_path))
    df = SpectrumDataFrameHelper.from_spectra_list(spectra)
    return df


def test_scale_defaults(spectra_dataframe):
    spectrum_df = spectra_dataframe[0]
    new_spectrum_df = SpectrumDataFrameHelper.scale(spectrum_df)

    assert spectrum_df.Intensity.max() > 1
    assert new_spectrum_df.Intensity.max() == 1
    assert new_spectrum_df.Intensity.min() == 0


def test_scale_custom_values(spectra_dataframe):
    spectrum_df = spectra_dataframe[0]
    new_spectrum_df = SpectrumDataFrameHelper.scale(spectrum_df, 1, 20)

    assert spectrum_df.Intensity.max() > 20
    assert new_spectrum_df.Intensity.max() == 20
    assert new_spectrum_df.Intensity.min() == 1


def test_scale_wrong_values(spectra_dataframe):
    spectrum_df = spectra_dataframe[0]

    with pytest.raises(ValueError):
        new_spectrum_df = SpectrumDataFrameHelper.scale(spectrum_df, max=2, min=5)


def test_filter(spectra_dataframe):
    n_peaks = 10
    spectrum_df = spectra_dataframe[0]
    new_spectrum_df = SpectrumDataFrameHelper.filter(spectrum_df, n_peaks)

    assert spectrum_df.shape[0] > n_peaks
    assert new_spectrum_df.shape[0] == n_peaks


def test_crop(spectra_dataframe):
    max_mz = 100
    spectrum_df = spectra_dataframe[0]
    new_spectrum_df = SpectrumDataFrameHelper.filter_mz(spectrum_df, (0, max_mz))

    assert spectrum_df["m/z"].max() > max_mz
    assert new_spectrum_df["m/z"].max() <= max_mz


def test_from_spectra_list(mgf_46_spectra_path):
    spectra = list(load_from_mgf(mgf_46_spectra_path))
    spectra_dfs = SpectrumDataFrameHelper.from_spectra_list(spectra)

    assert type(spectra_dfs) == list
    assert isinstance(spectra_dfs[0], pd.DataFrame)
    assert {"Intensity", "m/z"} == set(spectra_dfs[0].columns)
    assert len(spectra[0].peaks.intensities) == len(spectra_dfs[0].Intensity)
    assert len(spectra[0].peaks.mz) == len(spectra_dfs[0]["m/z"])
    assert spectra[0].peaks.intensities[0] == spectra_dfs[0].Intensity.iloc[0]
    assert spectra[0].peaks.mz[0] == spectra_dfs[0]["m/z"].iloc[0]


def test_from_spectrum(mgf_46_spectra_path):
    spectra = list(load_from_mgf(mgf_46_spectra_path))
    spectrum_df = SpectrumDataFrameHelper.from_spectrum(spectra[0])

    assert isinstance(spectrum_df, pd.DataFrame)
    assert {"Intensity", "m/z"} == set(spectrum_df.columns)
    assert len(spectra[0].peaks.intensities) == len(spectrum_df.Intensity)
    assert len(spectra[0].peaks.mz) == len(spectrum_df["m/z"])
    assert spectra[0].peaks.intensities[0] == spectrum_df.Intensity.iloc[0]
    assert spectra[0].peaks.mz[0] == spectrum_df["m/z"].iloc[0]


@pytest.mark.internet_connection
def test_from_gnps_id():
    spectrum_df = SpectrumDataFrameHelper.from_gnps_id("CCMSLIB00000078899")

    assert isinstance(spectrum_df, pd.DataFrame)
    assert {"Intensity", "m/z"} == set(spectrum_df.columns)
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
        {"Intensity", "m/z"} == set(spectrum_df.columns) for spectrum_df in spectra_df
    )
    assert all(spectrum_df.shape[0] > 0 for spectrum_df in spectra_df)
