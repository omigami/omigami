from logging import getLogger
from typing import Dict, Union, List
from matchms.importing import load_from_mgf
from omigami.endpoint import Endpoint

Payload = Dict[str, Dict[str, Dict[str, Union[int, dict]]]]
log = getLogger(__file__)

JSON = Union[List[dict], dict]


class InvalidCredentials(Exception):
    pass


class MS2DeepScore(Endpoint):
    _endpoint_url = "https://omigami.datarevenue.com/seldon/seldon/ms2deepscore/api/v0.1/predictions"

    def __init__(self, token: str):
        super().__init__(token)

    def match_spectra_from_path(
        self,
        mgf_path: str,
        n_best: int,
        include_metadata: List[str] = None,
        **kwargs,
    ) -> float:
        """
        Finds the N best matches for spectra in a local mgf file using MS2DeepScore
        algorithm.

        Parameters
        ----------
        mgf_path: str
            Local path to mgf file with two spectra
        n_best: int
            Number of best matches to select
        include_metadata: List[str]
            Metadata keys to include in the response. Will make response slower. Please
            check the documentation for a list of valid keys.

        Returns
        -------
        A list of pandas dataframes containing the best matches and optionally metadata
        for these matches.

        """
        parameters = self._build_parameters(n_best, include_metadata)

        spectra_generator = load_from_mgf(mgf_path)

        return self._make_batch_requests(
            spectra_generator, parameters, self._endpoint_url
        )
