import requests
from typing import List, Optional, Tuple

import pandas as pd
from matchms import Spectrum
from pandas.tests.groupby.test_value_counts import df


class SpectrumDataFrameHelper:
    @staticmethod
    def scale(spectrum_df: pd.DataFrame, min=0, max=1) -> pd.DataFrame:
        """
        Scale the Intensity column according to min and max values. By default, normalizes to 0 and 1.
        """
        if max <= min:
            raise ValueError(
                "'max' argument must have a bigger value than 'min' argument."
            )

        df = spectrum_df.copy()

        min_intensity, max_intensity = df.Intensity.min(), df.Intensity.max()
        df["Intensity"] = (df.Intensity - min_intensity) / (
            max_intensity - min_intensity
        ) * (max - min) + min

        return df

    @staticmethod
    def filter(spectrum_df: pd.DataFrame, num_peaks: int = 5000) -> pd.DataFrame:
        """
        Returns a new dataframe containing the largest n peaks according to num_peaks parameter.
        """
        df = spectrum_df.copy()

        return df.nlargest(num_peaks, "Intensity").sort_values("m/z")

    @staticmethod
    def filter_mz(
        spectrum_df: pd.DataFrame, mz_limits: Tuple[int, int] = (0, 5000)
    ) -> pd.DataFrame:
        """
        Returns a new filtered dataframe according to the provided min and max m/z values.
        """
        constrained_df = spectrum_df.copy()
        mass_low_limit, mass_high_limit = mz_limits
        constrained_df = constrained_df[
            constrained_df["m/z"].between(mass_low_limit, mass_high_limit)
        ]
        return constrained_df

    @staticmethod
    def from_spectra_list(
        spectra_list: List[Spectrum], scale_spectra: bool = False, num_peaks: int = 500
    ) -> List[pd.DataFrame]:
        """
        Creates dataframes from a list of Spectrum objects.
        """
        dfs = []
        for spectrum in spectra_list:
            df = SpectrumDataFrameHelper.from_spectrum(
                spectrum, scale_spectra, num_peaks
            )
            dfs.append(df)
        return dfs

    @staticmethod
    def from_spectrum(
        spectrum: Spectrum, scale_spectrum: bool = False, num_peaks: int = 500
    ) -> pd.DataFrame:
        """
        Creates a dataframe from a Spectrum object.
        """
        merged = list(zip(spectrum.peaks[0].tolist(), spectrum.peaks[1].tolist()))
        df = pd.DataFrame(merged, columns=["m/z", "Intensity"])
        if scale_spectrum:
            df = SpectrumDataFrameHelper.scale(df, num_peaks)
        return df

    @staticmethod
    def from_gnps_id(spectrum_id: str) -> pd.DataFrame:
        """
        Creates a dataframe using a GNPS spectrum ID by fetching GNPS metadata.
        """

        gnps_metadata = requests.get(
            "https://gnps.ucsd.edu/ProteoSAFe/SpectrumCommentServlet?SpectrumID="
            + str(spectrum_id)
        ).json()
        gnps_metadata_df = pd.DataFrame(
            eval(gnps_metadata["spectruminfo"]["peaks_json"]),
            columns=["m/z", "Intensity"],
        )

        return gnps_metadata_df

    @staticmethod
    def from_gnps_id_list(spectrum_id_list: List[str]) -> List[pd.DataFrame]:
        """
        Creates dataframes using a list of GNPS spectrum ID by fetching GNPS metadata.
        """
        dfs = []
        for gnps_id in spectrum_id_list:
            df = SpectrumDataFrameHelper.from_gnps_id(gnps_id)
            dfs.append(df)
        return dfs
