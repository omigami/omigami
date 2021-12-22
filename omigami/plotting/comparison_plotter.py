from typing import List, Optional, Any, Tuple

import pandas as pd

import matplotlib.pyplot as plt
from omigami.utilties.gnps_helper import GnpsHelper
from omigami.utilties.spectrum_df_helper import SpectrumDataFrameHelper


class SpectraComparisonPlotter():

    def __init__(self):
        pass

    def mirror_plot(
        in_spectrum_df: pd.DataFrame,
        spectrum_matches_df: pd.DataFrame,
        match_idx: int,
        labels: Optional[List[str]] = None,
        display_limits: Optional[Tuple[int, int]] = (0, 5000),
    ) -> Any:
        """
        Creates a mirror plot with spectrum 1 on top (blue), and spectrum 2 on the bottom (red).

        Requires:
        spectrum1 : A dataframe version of the spectrum information for top spectrum
        spectrum2 : A dataframe version of the spectrum information for bottom spectrum

        Optional:
        labels (list): A list containing the names of the spectra.
        Display_Limits (list): left and right boundaries for displaying the spectrum.

        Returns:
        fig (plot): The mirror plot.
        """

        if not labels:
            labels = ["Spectrum 1", "Spectrum 2"]

        spectrum_match_id = str(spectrum_matches_df.index.tolist()[match_idx])
        spectrum_match_name = spectrum_matches_df.compound_name.loc[spectrum_match_id]

        sdf_helper = SpectrumDataFrameHelper()
        gnps_helper = GnpsHelper()

        # currently fetching metadata for the match from GNPS, later we should use our own (to get the peaks)
        spectrum_match_df = gnps_helper.fetch_metadata(spectrum_match_id)
        spectrum_match_df = sdf_helper.scale(spectrum_match_df, 100)

        _in_spectrum_df = sdf_helper.crop(in_spectrum_df, display_limits)
        spectrum_match_df = sdf_helper.crop(spectrum_match_df, display_limits)

        fig, axes = plt.subplots(figsize=(20, 10), nrows=2, sharex=True)
        axes[0].bar(
            in_spectrum_df["m/z"], in_spectrum_df.Intensity, color="b", label=labels[0]
        )
        axes[1].bar(
            spectrum_match_df["m/z"],
            spectrum_match_df.Intensity,
            color="r",
            label=labels[1],
        )
        axes[1].invert_yaxis()
        axes[1].set_xlabel("m/z")
        fig.legend()
        fig.subplots_adjust(hspace=0)

        return fig
