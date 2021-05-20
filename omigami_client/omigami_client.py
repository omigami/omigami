import ast
import json
from typing import Dict, Union, List, Generator

import pandas as pd
import requests
from matchms import Spectrum
from matchms.importing import load_from_mgf

Payload = Dict[str, Dict[str, Dict[str, Union[int, dict]]]]


class OmigamiClient:
    _endpoint_url = (
        "https://mlops.datarevenue.com/seldon/seldon/spec2vec/api/v0.1/predictions"
    )

    def __init__(self, token: str):
        self._token = token

    def match_spectra_from_path(self, mgf_path: str, n_best: int) -> List[pd.DataFrame]:
        """
        Finds the N best matches for spectra in a local mgf file using spec2vec algorithm.

        Parameters
        ----------
        mgf_path: str
            Local path to mgf file
        n_best:
            Number of best matches to select

        Returns
        -------
        A list of pandas dataframes containing the best matches.

        """
        spectra_generator = load_from_mgf(mgf_path)
        payload = self._build_payload(spectra_generator, n_best)
        api_request = self._send_request(payload)
        prediction = self._format_results(api_request)

        return prediction

    def _build_payload(
        self, spectra_generator: Generator[Spectrum, None, None], n_best_spectra: int
    ) -> Payload:
        """Extract abundance pairs and Precursor_MZ data, then build the json payload"""
        spectra = []
        for spectrum in spectra_generator:
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

        self._validate_input(spectra)

        payload = {
            "data": {
                "ndarray": {
                    "parameters": {
                        "n_best_spectra": n_best_spectra,
                    },
                    "data": spectra,
                }
            }
        }

        return payload

    @staticmethod
    def _validate_input(model_input: List[Dict]):
        for i, spectrum in enumerate(model_input):
            if not isinstance(spectrum, Dict):
                raise TypeError(
                    f"Spectrum data must be a dictionary. Passed type: {type(spectrum)}"
                )

            mandatory_keys = ["peaks_json", "Precursor_MZ"]
            if any(key not in spectrum.keys() for key in mandatory_keys):
                raise KeyError(
                    f"Please include all the mandatory keys in your input data. "
                    f"The mandatory keys are {mandatory_keys}",
                )

            if isinstance(spectrum["peaks_json"], str):
                try:
                    ast.literal_eval(spectrum["peaks_json"])
                except SyntaxError:
                    raise ValueError(
                        "peaks_json needs to be a valid python string representation of "
                        "a list or a list. Passed value: {spectrum['peaks_json']}",
                        400,
                    )
            elif not isinstance(spectrum["peaks_json"], list):
                raise ValueError(
                    "peaks_json needs to be a valid python string representation of a "
                    f"list or a list. Passed value: {spectrum['peaks_json']}",
                    400,
                )

            float_keys = ["Precursor_MZ", "Charge"]
            for key in float_keys:
                if spectrum.get(key):
                    try:
                        float(spectrum[key])
                    except ValueError:
                        raise ValueError(
                            f"{key} needs to be a string representation of a float. "
                            f"Passed value: {spectrum[key]}",
                            400,
                        )

    def _send_request(self, payload: Payload) -> requests.Response:
        api_request = requests.post(
            self._endpoint_url,
            json=payload,
            headers={"Authorization": f"Bearer {self._token}"},
            timeout=600,
        )
        return api_request

    @staticmethod
    def _format_results(api_request: requests.Response) -> List[pd.DataFrame]:
        response = json.loads(api_request.text)
        library_spectra_raw = response["data"]["ndarray"]

        predicted_spectra = []
        for i in range(len(library_spectra_raw)):
            library_spectra_dataframe = pd.DataFrame(
                data=[spectrum["score"] for spectrum in library_spectra_raw[i]],
                index=[
                    spectrum["match_spectrum_id"] for spectrum in library_spectra_raw[i]
                ],
                columns=["score"],
            )
            library_spectra_dataframe.index.name = f"matches of spectrum #{i + 1}"
            predicted_spectra.append(library_spectra_dataframe)

        return predicted_spectra
