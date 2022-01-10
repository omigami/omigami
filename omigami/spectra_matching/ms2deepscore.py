from typing import Optional

from omigami.spectra_matching.spectra_matching import SpectraMatching


class MS2DeepScore(SpectraMatching):
    _algorithm = "ms2deepscore"

    def __init__(self, optional_token: Optional[str] = None):
        super(MS2DeepScore, self).__init__(optional_token)
