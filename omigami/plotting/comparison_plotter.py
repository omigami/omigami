from typing import List, Optional, Any, Tuple, Union

import pandas as pd

import matplotlib.pyplot as plt
from omigami.utilities.gnps_helper import GnpsHelper
from omigami.utilities.spectrum_df_helper import SpectrumDataFrameHelper


class SpectraComparisonPlotter:
    def __init__(self):
        pass

    def mirror_plot(
        spectrum_1: Union[pd.DataFrame, str],
        spectrum_2: Union[pd.DataFrame, str],
        labels: Optional[List[str]] = None,
        display_limits: Optional[Tuple[int, int]] = (0, 5000),
    ) -> Any:
        """
        Creates a mirror plot with spectrum 1 on top (blue), and spectrum 2 on the bottom (red).

        Requires:
        spectrum1 : A dataframe version of the spectrum information for top spectrum
        spectrum2 : A dataframe version of the spectrum information for bottom spectrum

        Optional:
        match_id (int): If the dataframes are a
        labels (list): A list containing the names of the spectra.
        Display_Limits (list): left and right boundaries for displaying the spectrum.

        Returns:
        fig (plot): The mirror plot.
        """

        if not labels:
            labels = ["Spectrum 1", "Spectrum 2"]

        gnps_helper = GnpsHelper()
        sdf_helper = SpectrumDataFrameHelper()

        spectra = [spectrum_1, spectrum_2]
        spectra_df_list = []
        for spectrum in spectra:
            if type(spectrum) == str:
                spectrum_df = gnps_helper.fetch_metadata(spectrum)
            elif type(spectrum) == pd.DataFrame:
                if ["m/z", "Intensities"] not in spectrum.columns.values:
                    raise ValueError(
                        "Spectra arguments of type DataFrame must contain 'Intensity' and 'm/z' columns."
                    )
                spectrum_df = spectrum
            else:
                raise ValueError(
                    "Spectra arguments must be either 'str' (gnps ID) or 'pd.DataFrame' type."
                )

            spectrum = sdf_helper.crop(spectrum_df, display_limits)
            spectra_df_list.append(spectrum)

        fig, axes = plt.subplots(figsize=(20, 10), nrows=2, sharex=True)
        axes[0].bar(
            spectra_df_list[0]["m/z"],
            spectra_df_list[0].Intensity,
            color="b",
            label=labels[0],
        )
        axes[1].bar(
            spectra_df_list[1]["m/z"],
            spectra_df_list[1].Intensity,
            color="r",
            label=labels[1],
        )
        axes[1].invert_yaxis()
        axes[1].set_xlabel("m/z")
        fig.legend()
        fig.subplots_adjust(hspace=0)

        return fig
