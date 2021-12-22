from typing import List, Optional, Tuple

import pandas as pd
from matchms import Spectrum
from pandas.tests.groupby.test_value_counts import df


class SpectrumDataFrameHelper():

    def __init__(self):
        pass

    @staticmethod
    def scale(self, spectra_df: pd.DataFrame, min=0, max=1, num_peaks=5000) -> pd.DataFrame:
        """
        Scale the Intensity column according to min and max values. By default, normalizes to 0 and 1.
        """
        df = spectra_df.copy()

        min_intensity, max_intensity = df.Intensity.min(), df.Intensity.max()
        df['Intensity'] = (df.Intensity - min_intensity) / (max_intensity - min_intensity) * (max - min) + min

        return df

    @staticmethod
    def filter(self, spectra_df: pd.DataFrame, num_peaks=5000) -> pd.DataFrame:
        """
        Returns a dataframe containing the largest n peaks according to num_peaks parameter
        """
        return df.nlargest(num_peaks, "Intensity").sort_values("m/z")

    @staticmethod
    def crop(
            spectrum_df: pd.DataFrame, display_limits: Optional[Tuple[int, int]] = (0, 5000)
    ) -> pd.DataFrame:
        """
        Constrains the spectrum DF to fit within the display limits.
        """
        constrained_df = spectrum_df.copy()
        mass_low_limit, mass_high_limit = display_limits
        Spectrum_constrained = constrained_df[
            constrained_df["m/z"].between(mass_low_limit, mass_high_limit)
        ]
        return constrained_df

    @staticmethod
    def from_spectra_list(spectra_list: List[Spectrum]) -> pd.DataFrame:
        raise NotImplementedError()

