from logging import getLogger
from typing import Dict, Union, List
from matchms.importing import load_from_mgf

from omigami.authentication import authenticate_client
from omigami.endpoint import Endpoint

Payload = Dict[str, Dict[str, Dict[str, Union[int, dict]]]]
log = getLogger(__file__)

JSON = Union[List[dict], dict]


class InvalidCredentials(Exception):
    pass


class MS2DeepScore(Endpoint):
    _PREDICT_ENDPOINT_BASE = (
        "https://app.omigami.com/seldon/seldon/ms2deepscore-{ion_mode}/api/v0.1/predictions"
    )

    def match_spectra_from_path(
        self,
        mgf_path: str,
        n_best: int,
        include_metadata: List[str] = None,
        ion_mode: str = "positive",
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
        ion_mode: str
            String setting the ion mode. Can either be 'positive' or 'negative'
        include_metadata: List[str]
            Metadata keys to include in the response. Will make response slower. Please
            check the documentation for a list of valid keys.

        Returns
        -------
        A list of pandas dataframes containing the best matches and optionally metadata
        for these matches.

        """
        if ion_mode not in ["positive", "negative"]:
            raise ValueError(
                "Parameter ion_mode should be either set to 'positive' or 'negative'. "
                "Defaults to 'positive'. "
            )

        endpoint = self._PREDICT_ENDPOINT_BASE.format(ion_mode=ion_mode)
        parameters = self._build_parameters(n_best, include_metadata)

        # gets token from user credentials
        authenticate_client()

        # loads spectra
        spectra_generator = load_from_mgf(mgf_path)

        # issue requests respecting the spectra limit per request
        return self._make_batch_requests(spectra_generator, parameters, endpoint)
