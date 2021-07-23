import json
import ast
from abc import abstractmethod
from logging import getLogger
import requests
import pandas as pd
from typing import List, Dict, Any, Union


from omigami.custom_exceptions import InvalidCredentials, NotFoundError
from omigami.visualisation.config import VALID_KEYS, MATCH_LIMIT_PER_REQUEST

JSON = Union[List[dict], dict]
Payload = Dict[str, Dict[str, Dict[str, Union[int, dict]]]]

log = getLogger(__file__)


class Endpoint:

    def __init__(self, token: str):
        self._token = token
        self.mandatory_keys = ["compound_name"]

        # TODO: Find better name. For a key of which you dont need all of but at least one.
        self.semi_optional_keys = []

    @abstractmethod
    def plot(
            self,
            # TODO: Rather from path or as a direct dataframe?
            scores_path: str,
            n_best: int,
            include_metadata: List[str] = VALID_KEYS,
            **kwargs,
    ) -> List[pd.DataFrame]:
        """
        Finds the N best matches for spectra in a local mgf file
        """
        pass

    @staticmethod
    def _build_parameters(n_best: int, include_metadata: List[str]) -> Dict[str, Any]:
        parameters = {}
        try:
            parameters["n_best_spectra"] = int(n_best)
        except ValueError:
            raise ValueError(
                "The number of best features parameter must be an integer."
            )

        if include_metadata:
            for key in include_metadata:
                if key.lower() not in VALID_KEYS:
                    raise ValueError(
                        f"The metadata {key} is not included in the valid keys list. "
                        f"Please check documentation for the list of valid keys."
                    )
            include_metadata = [key.lower() for key in include_metadata]
            parameters["include_metadata"] = include_metadata

        return parameters

    def _send_request(self, payload: Payload, endpoint: str) -> requests.Response:
        api_request = requests.post(
            endpoint,
            json=payload,
            headers={"Authorization": f"Bearer {self._token}"},
            timeout=600,
        )

        if api_request.status_code == 401:
            raise InvalidCredentials(
                "Your credentials are invalid, please revise your API token."
            )
        if api_request.status_code == 404:
            raise NotFoundError("The API endpoint couldn't be reached.")

        return api_request

    def _build_payload(
            self,
            batch: Dict[str, Any],
            parameters: Dict[str, Any],
    ) -> JSON:
        """Extract abundance pairs and Precursor_MZ data, then build the json payload

        Returns:
        --------
        - payload: JSON
            the full request payload with input data
        """

        log.info(f"{len(batch)} spectra found on input file.")
        self._validate_input(batch)

        payload = {
            "data": {
                "ndarray": {
                    "parameters": parameters,
                    "data": batch,
                }
            }
        }
        return payload

    # TODO: Wrtie Validation input
    def _validate_input(self, model_input: List[Dict]):

        for i, spectrum in enumerate(model_input):
            if not isinstance(spectrum, Dict):
                raise TypeError(
                    f"Spectrum data must be a dictionary. Passed type: {type(spectrum)}"
                )

            if any(key not in spectrum.keys() for key in self.mandatory_keys):
                raise KeyError(
                    f"Please include all the mandatory keys in your input data. "
                    f"The mandatory keys are {self.mandatory_keys}",
                )

            if isinstance(spectrum["compound_name"], str):
                try:
                    ast.literal_eval(spectrum["peaks_json"])
                except SyntaxError:
                    raise ValueError(
                        f"compound_name needs to be a valid python string representation "
                        f"of a list or a list. Passed value: {spectrum['compound_name']}",
                        400,
                    )


    def _make_batch_requests(self, spectra_generator, parameters, endpoint):
        batch = []
        requests = []

        for spectrum in spectra_generator:
            batch.append(spectrum)
        if len(batch) == MATCH_LIMIT_PER_REQUEST:
            payload = self._build_payload(batch, parameters)
            requests.append(self._send_request(payload, endpoint))
            batch = []

        if batch:
            payload = self._build_payload(batch, parameters)
            requests.append(self._send_request(payload, endpoint))

        predictions = []
        for r in requests:
            predictions.extend(self._format_results(r))
        return predictions


