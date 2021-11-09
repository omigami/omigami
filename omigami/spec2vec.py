from typing import List

import pandas as pd
from matchms.importing import load_from_mgf

from omigami.authentication import authenticate_client
from omigami.endpoint import Endpoint


class Spec2Vec(Endpoint):
    _PREDICT_ENDPOINT_BASE = (
        "https://app.omigami.com/seldon/seldon/spec2vec-{ion_mode}/api/v0.1/predictions"
    )

    def match_spectra_from_path(
        self,
        mgf_path: str,
        n_best: int,
        include_metadata: List[str] = None,
        ion_mode: str = "positive",
    ) -> List[pd.DataFrame]:
        """
        Finds the N best matches for spectra in a local mgf file using spec2vec algorithm.

        Parameters
        ----------
        mgf_path: str
            Local path to mgf file
        n_best: int
            Number of best matches to select
        include_metadata: List[str]
            Metadata keys to include in the response. Will make response slower. Please
            check the documentation for a list of valid keys.
        ion_mode: str
            Selects which model will be used for the predictions: Either a model trained with
            positive or negative ion mode spectra data. Defaults to positive.

        Returns
        -------
        A list of pandas dataframes containing the best matches and optionally metadata
        for these matches.

        """
        # validates input
        if ion_mode not in ["positive", "negative"]:
            raise ValueError(
                "Parameter ion_mode should be either set to 'positive' or 'negative'. "
                "Defaults to 'positive'. "
            )
        # defines endpoint based on user choice of spectra ion mode
        endpoint = self._PREDICT_ENDPOINT_BASE.format(ion_mode=ion_mode)

        # gets token from user credentials
        authenticate_client()

        parameters = self._build_parameters(n_best, include_metadata)
        # loads spectra
        spectra_generator = load_from_mgf(mgf_path)

        # issue requests respecting the spectra limit per request
        return self._make_batch_requests(spectra_generator, parameters, endpoint)
