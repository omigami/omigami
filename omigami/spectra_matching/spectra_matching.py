from __future__ import annotations

import ast
import json
from io import StringIO
from logging import getLogger
from typing import List, Dict, Any, Union, Generator, Optional

import pandas as pd
import requests
from matchms import Spectrum
from matchms.importing import load_from_mgf

from omigami.authentication import get_session, authenticate_client, set_token
from omigami.exceptions import (
    InvalidCredentials,
    NotFoundError,
    InternalServerError,
    InvalidUsageError,
)
from omigami.omi_settings import HOST_NAME

SPECTRA_LIMIT_PER_REQUEST = 100
VALID_KEYS = {
    "compound_name",
    "inchikey_inchi",
    "inchikey_smiles",
    "instrument",
    "parent_mass",
    "smiles",
    "precursor_mz",
}

log = getLogger(__file__)

SpectraJson = List[Dict[str, str]]
SpectraMatchingParameters = Dict[str, Any]
Payload = Dict[
    "data", Dict["ndarray", Dict[str, Union[SpectraMatchingParameters, SpectraJson]]]
]


class SpectraMatching:
    _PREDICT_ENDPOINT_BASE = (
        f"{HOST_NAME}seldon/seldon/" + "{algorithm}-{ion_mode}/api/v0.1/predictions"
    )
    # "https://dev.omigami.com/seldon/seldon/{algorithm}-{ion_mode}/api/v0.1/predictions"
    _ENDPOINT = None
    _optional_token = None

    mandatory_keys: List[str] = ["peaks_json", "Precursor_MZ"]
    float_keys: List[str] = ["Precursor_MZ"]

    def __init__(self, optional_token: Optional[str] = None):
        self._optional_token = optional_token

    def match_spectra(
        self,
        source: Union[str, list[Spectrum], StringIO],
        n_best: int,
        ion_mode: str = "positive",
    ) -> List[pd.DataFrame]:
        """
        From a spectra source, issues requests to either MS2DeepScore or Spec2Vec endpoints to find the N best library
        matches.

        Parameters
        ----------
        source: str or list[Spectrum] or StringIO
            either a local path to mgf file (str), or a list of preloaded Spectrum objects, or a StringIO obj
            from a loaded mgf file
        n_best: int
            Number of best matches to select
        ion_mode:
            Selects which model will be used for the predictions: Either a model trained with
            positive or negative ion mode spectra data. Defaults to positive.

        Returns
        -------
        A list of pandas dataframes containing the best matches and optionally metadata
        for these matches.

        """

        spectra_generator = self._create_spectra_generator(source)

        if self._ENDPOINT is None:
            raise InvalidUsageError(
                "You should only evoke match_spectra from either "
                "a MS2DeepScore or a Spec2Vec class instance."
            )

        # defines endpoint based on user choice of spectra ion mode
        endpoint = self._PREDICT_ENDPOINT_BASE.format(
            algorithm=self._ENDPOINT, ion_mode=ion_mode
        )

        # validates input
        if ion_mode not in ["positive", "negative"]:
            raise ValueError(
                "Parameter ion_mode should be either set to 'positive' or 'negative'. "
                "Defaults to 'positive'. "
            )

        # gets token from user credentials
        if self._optional_token is not None:
            set_token(self._optional_token)
        else:
            authenticate_client()

        parameters = self._build_parameters(n_best)

        # issue requests respecting the spectra limit per request
        return self._make_batch_requests(spectra_generator, parameters, endpoint)

    @staticmethod
    def _create_spectra_generator(source: Union[str, StringIO, List[Spectrum]]):
        spectra_generator: Generator[Spectrum]
        if type(source) == str or type(source) == StringIO:
            spectra_generator = load_from_mgf(source)
        else:

            def _spectra_generator(spectra_list: List[Spectrum]) -> Generator[Spectrum]:
                for spectrum in spectra_list:
                    yield spectrum

            spectra_generator = _spectra_generator(source)
        return spectra_generator

    def _build_payload(
        self,
        batch: List[Spectrum],
        parameters: SpectraMatchingParameters,
    ) -> Payload:
        """Extract abundance pairs and Precursor_MZ data, then build the json payload

        Returns:
        --------
        - payload: JSON
            the full request payload with input data
        """
        spectra = []
        for spectrum in batch:
            if "pepmass" in spectrum.metadata:
                precursor_mz = str(spectrum.metadata["pepmass"][0])
            elif "precursor_mz" in spectrum.metadata:
                precursor_mz = str(spectrum.metadata["precursor_mz"][0])
            else:
                raise KeyError(
                    "One of ['pepmass' or 'precursor_mz'] must be present in the"
                    "spectrum metadata."
                )

            spectra.append(
                {
                    "peaks_json": str(
                        [
                            [mz, intensity]
                            for mz, intensity in zip(
                                spectrum.peaks.mz, spectrum.peaks.intensities
                            )
                        ]
                    ),
                    "Precursor_MZ": precursor_mz,
                }
            )

        log.info(f"{len(spectra)} spectra found on input file.")
        self._validate_input(spectra)

        payload = {
            "data": {
                "ndarray": {
                    "parameters": parameters,
                    "data": spectra,
                }
            }
        }
        return payload

    def _make_batch_requests(
        self,
        spectra_generator: Generator[Spectrum],
        parameters: SpectraMatchingParameters,
        endpoint: str,
    ) -> List[pd.DataFrame]:
        batch = []
        requests = []

        for spectrum in spectra_generator:
            batch.append(spectrum)
            if len(batch) == SPECTRA_LIMIT_PER_REQUEST:
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

    @staticmethod
    def _send_request(payload: Payload, endpoint: str) -> requests.Response:
        auth = get_session()
        api_request = requests.post(
            endpoint,
            json=payload,
            headers={"Authorization": f"Bearer {auth.session_token}"},
            timeout=600,
        )

        if api_request.status_code == 401:
            raise InvalidCredentials("Your credentials are invalid.")
        if api_request.status_code == 404:
            raise NotFoundError("The API endpoint couldn't be reached.")

        return api_request

    @staticmethod
    def _build_parameters(n_best: int) -> SpectraMatchingParameters:
        parameters = {}
        try:
            parameters["n_best_spectra"] = int(n_best)
        except ValueError:
            raise ValueError(
                "The number of best features parameter must be an integer."
            )

        return parameters

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

            if isinstance(spectrum["peaks_json"], str):
                try:
                    ast.literal_eval(spectrum["peaks_json"])
                except SyntaxError:
                    raise ValueError(
                        f"peaks_json needs to be a valid python string representation "
                        f"of a list or a list. Passed value: {spectrum['peaks_json']}",
                        400,
                    )
            elif not isinstance(spectrum["peaks_json"], list):
                raise ValueError(
                    "peaks_json needs to be a valid python string representation of a "
                    f"list or a list. Passed value: {spectrum['peaks_json']}",
                    400,
                )

            for key in self.float_keys:
                if spectrum.get(key):
                    try:
                        float(spectrum[key])
                    except ValueError:
                        raise ValueError(
                            f"{key} needs to be a string representation of a float. "
                            f"Passed value: {spectrum[key]}",
                            400,
                        )

    @staticmethod
    def _format_results(api_request: requests.Response) -> List[pd.DataFrame]:
        if api_request.status_code == 500:
            raise InternalServerError(
                f"Error in generating the prediction: "
                f"{api_request.json()['status']['info']}.\n"
                f"Please try again later or contact DataRevenue for more information."
            )

        response = json.loads(api_request.text)
        library_spectra_raw = response["jsonData"]

        predicted_spectra = []
        for id_, matches in library_spectra_raw.items():
            scores_df = pd.DataFrame(matches).T.drop(columns="metadata")
            metadata_df = pd.DataFrame(
                pd.DataFrame(matches).loc["metadata"].to_dict()
            ).T
            library_spectra_dataframe = scores_df.join(metadata_df)
            library_spectra_dataframe.index.name = f"matches of {id_}"
            library_spectra_dataframe = library_spectra_dataframe.sort_values(
                by=["score"], ascending=False
            )

            predicted_spectra.append(library_spectra_dataframe)

        return predicted_spectra
