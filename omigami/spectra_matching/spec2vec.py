from typing import Optional

from omigami.spectra_matching.spectra_matching import SpectraMatching


class Spec2Vec(SpectraMatching):
    _algorithm = "spec2vec"

    def __init__(self, optional_token: Optional[str] = None):
        super(Spec2Vec, self).__init__(optional_token)
