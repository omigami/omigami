import ast
import json
from abc import abstractmethod
from logging import getLogger
from typing import List, Dict, Any, Union
import requests
import pandas as pd
from matchms import Spectrum

from omigami.authentication import get_session
from omigami.exceptions import InvalidCredentials, NotFoundError, InternalServerError

SPECTRA_LIMIT_PER_REQUEST = 100
VALID_KEYS = {
    "compound_name",
    "inchikey_inchi",
    "inchikey_smiles",
    "instrument",
    "parent_mass",
    "smiles",
}

log = getLogger(__file__)

JSON = Union[List[dict], dict]
Payload = Dict[str, Dict[str, Dict[str, Union[int, dict]]]]


def _sort_columns(df: pd.DataFrame):
    sorted_columns = ["score"] + sorted(list(VALID_KEYS))
    return df.reindex(columns=sorted_columns).dropna(axis=1, how="all")


class Endpoint:
    mandatory_keys: List[str] = ["peaks_json", "Precursor_MZ"]
    float_keys: List[str] = ["Precursor_MZ"]

    @abstractmethod
    def match_spectra_from_path(
        self,
        mgf_path: str,
        n_best: int,
        include_metadata: List[str] = None,
        ion_mode: str = "positive",
    ) -> List[pd.DataFrame]:
        """
        Finds the N best matches for spectra in a local mgf file
        """
        pass

    def _build_payload(
        self,
        batch: List[Spectrum],
        parameters: Dict[str, Any],
    ) -> JSON:
        """Extract abundance pairs and Precursor_MZ data, then build the json payload

        Returns:
        --------
        - payload: JSON
            the full request payload with input data
        """
        spectra = []
        for spectrum in batch:
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
                    "Precursor_MZ": str(spectrum.metadata["pepmass"][0]),
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

    def _make_batch_requests(self, spectra_generator, parameters, endpoint):
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

    def _send_request(self, payload: Payload, endpoint: str) -> requests.Response:
        auth = get_session()
        api_request = requests.post(
            endpoint,
            json=payload,
            headers={"Authorization": f"Bearer {auth.session_token}"},
            timeout=600,
        )

        if api_request.status_code == 401:
            raise InvalidCredentials(
                "Your credentials are invalid, please revise your API token."
            )
        if api_request.status_code == 404:
            raise NotFoundError("The API endpoint couldn't be reached.")

        return api_request

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
                "Something went wrong, the requested service is probably unavailable at the moment. "
                "Please try again later or contact DataRevenue for more information."
            )

        response = json.loads(api_request.text)
        library_spectra_raw = response["jsonData"]

        predicted_spectra = []
        for id_, matches in library_spectra_raw.items():
            library_spectra_dataframe = pd.DataFrame(matches).T
            library_spectra_dataframe.index.name = f"matches of {id_}"
            library_spectra_dataframe = _sort_columns(library_spectra_dataframe)
            library_spectra_dataframe = library_spectra_dataframe.sort_values(
                by=["score"], ascending=False
            )

            predicted_spectra.append(library_spectra_dataframe)

        return predicted_spectra
