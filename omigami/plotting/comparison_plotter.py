from typing import List, Optional, Any, Tuple

import pandas as pd

import matplotlib.pyplot as plt
from omigami.utilities.spectrum_df_helper import SpectrumDataFrameHelper


class SpectraComparisonPlotter:
    def mirror_plot(
        self,
        spectrum_1: pd.DataFrame,
        spectrum_2: pd.DataFrame,
        labels: Optional[List[str]] = None,
        display_limits: Optional[Tuple[int, int]] = (0, 5000),
    ) -> Any:
        """
        Creates a mirror plot with spectrum 1 on top (blue), and spectrum 2 on the bottom (red).

        Parameters
        ----------
        spectrum_1 : DataFrame
            A dataframe version of the spectrum information for top spectrum
        spectrum_2 : DataFrame
            A dataframe version of the spectrum information for bottom spectrum
        labels : Optional[List[str]]
            A list containing the names of the spectra to display in the plot.
        display_limits: Optional[Tuple[int, int]]
            Left and right boundaries (in m/z) for displaying the spectrum. Defaults to (0, 5000)

        Returns
        -------
        fig: matplotlib.figure.Figure
            The mirror plot.
        """

        if not labels:
            labels = ["Spectrum 1", "Spectrum 2"]
        elif len(labels) != 2:
            raise ValueError("'labels' arg must be a list of strings with len == 2")

        sdf_helper = SpectrumDataFrameHelper()

        spectra = [spectrum_1, spectrum_2]
        spectra_df_list = []
        for spectrum in spectra:
            if (
                type(spectrum) != pd.DataFrame
                or ["m/z", "Intensity"] not in spectrum.columns.values
            ):
                raise ValueError(
                    "Spectra arguments must be of type DataFrame and must contain 'Intensity' and 'm/z' columns."
                )

            spectrum_df = spectrum
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
