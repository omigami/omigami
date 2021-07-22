from typing import List

from matchms import Spectrum
from ms2deepscore import BinnedSpectrum, SpectrumBinner


class MS2DeepScoreSpectrumBinner:
    def __init__(self):
        self.spectrum_binner = SpectrumBinner(number_of_bins=10000)

    def bin_spectra(
        self,
        spectra: List[Spectrum],
    ) -> List[BinnedSpectrum]:

        binned_spectra = self.spectrum_binner.fit_transform(spectra)

        for binned_spectrum, spectrum in zip(binned_spectra, spectra):
            binned_spectrum.set("spectrum_id", spectrum.get("spectrum_id"))
            binned_spectrum.set("inchi", spectrum.get("inchi"))
        return binned_spectra
