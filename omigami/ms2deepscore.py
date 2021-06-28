import ast
import json
import numpy as np
from logging import getLogger
from typing import Dict, Union, List, Generator

import requests
from matchms import Spectrum
from matchms.importing import load_from_mgf


Payload = Dict[str, Dict[str, Dict[str, Union[int, dict]]]]
log = getLogger(__file__)

JSON = Union[List[dict], dict]


class InvalidCredentials(Exception):
    pass


class InvalidNumberOfSpectra(Exception):
    pass


class MS2DeepScore:
    _endpoint_url = "https://omigami.datarevenue.com/seldon/seldon/ms2deepscore/api/v0.1/predictions"

    def __init__(self, token: str):
        self._token = token

    def predict_similarity_of_pair(
        self,
        mgf_path: str,
    ) -> float:
        """
        Predicts molecular structural similarities (Tanimoto scores) from pairs of
        mass spectrometry spectra in local mgf files using MS2DeepScore algorithm

        Parameters
        ----------
        mgf_path: str
            Local path to mgf file with two spectra

        Returns
        -------
        Molecular structural similarities (Tanimoto scores)

        """
        spectra_generator = load_from_mgf(mgf_path)

        payload = self._build_payload(spectra_generator)
        response = self._send_request(payload)

        prediction = self._format_results(response)
        return prediction

    def _build_payload(
        self,
        pair: Generator[Spectrum, None, None],
    ) -> JSON:
        """Build the json payload

        Returns:
        --------
        - payload: JSON
            the full request payload with input data
        """
        spectra = []
        for spectrum in pair:
            spectra.append(
                {
                    "intensities": spectrum.peaks.intensities,
                    "mz": spectrum.peaks.mz,
                }
            )

        self._validate_input(spectra)

        payload = {
            "data": {
                "ndarray": {
                    "parameters": {},  # maybe there will be some later on
                    "data": spectra,
                }
            }
        }
        return payload

    @staticmethod
    def _validate_input(model_input: List[Dict]):
        if len(model_input) != 2:
            raise InvalidNumberOfSpectra(
                f"MS2DeepScore can only predict the Tanimoto score for a pair of "
                f"spectra. Please input two spectra at a time. Found "
                f"{len(model_input)} spectra. "
            )

        for i, spectrum in enumerate(model_input):
            if not isinstance(spectrum, Dict):
                raise TypeError(
                    f"Spectrum data must be a dictionary. Passed type: {type(spectrum)}"
                )

            mandatory_keys = ["intensities", "mz"]
            if any(key not in spectrum.keys() for key in mandatory_keys):
                raise KeyError(
                    f"Please include all the mandatory keys in your input data. "
                    f"The mandatory keys are {mandatory_keys}",
                )
            for key in mandatory_keys:
                if isinstance(spectrum[key], str):
                    try:
                        ast.literal_eval(spectrum[key])
                    except SyntaxError:
                        raise ValueError(
                            f"{key} needs to be a valid python string representation "
                            f"of a list or a list. Passed value: {spectrum[key]}",
                            400,
                        )
                elif not isinstance(spectrum[key], list) and not isinstance(
                    spectrum[key], np.ndarray
                ):
                    raise ValueError(
                        f"{key} needs to be a valid python string representation of a "
                        f"list or a list. Passed value: {spectrum[key]}",
                        400,
                    )

    def _send_request(self, payload: Payload) -> requests.Response:
        api_request = requests.post(
            self._endpoint_url,
            json=payload,
            headers={"Authorization": f"Bearer {self._token}"},
            timeout=600,
        )

        if api_request.status_code == 401:
            raise InvalidCredentials(
                "Your credentials are invalid, please revise your API token."
            )
        return api_request

    @staticmethod
    def _format_results(api_request: requests.Response) -> Dict[str, float]:
        response = json.loads(api_request.text)
        return {"Tanimoto Score": response["jsonData"]}
