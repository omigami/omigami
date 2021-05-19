import json
from collections import Generator
from typing import Dict, Union, List

import pandas as pd
import requests
from matchms import Spectrum
from matchms.importing import load_from_mgf

Payload = Dict[str, Dict[str, Dict[str, Union[int, dict]]]]


class OmigamiClient:
    _endpoint_url = "https://mlops.datarevenue.com/seldon/seldon/spec2vec/api/v0.1/predictions"

    def match_spectra_from_path(self, mgf_path: str) -> List[pd.DataFrame]:
        """TODO"""
        spectra_generator = load_from_mgf(mgf_path)
        payload = self._build_payload(spectra_generator, 10)
        api_request = requests.post(
            self._endpoint_url,
            json=payload,
            headers={"Authorization": "Bearer 5iTy7vACUXmlLO9fwGAL8v2WLPbo1SNH"},
            timeout=600,
        )
        prediction = self._format_results(api_request)

        return prediction

    @staticmethod
    def _build_payload(
            spectra_generator: Generator[Spectrum, None, None], n_best_spectra: int
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

        # build the payload
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
    def _format_results(api_request) -> List[pd.DataFrame]:
        """Formatting of the results"""
        response = json.loads(api_request.text)
        library_spectra_raw = response["data"]["ndarray"]

        predicted_spectra = []
        for i in range(len(library_spectra_raw)):
            library_spectra_dataframe = pd.DataFrame(
                data=[spectrum_id["score"] for spectrum_id in library_spectra_raw[i]],
                index=[
                    spectrum_id["match_spectrum_id"]
                    for spectrum_id in library_spectra_raw[i]
                ],
                columns=["score"],
            )
            library_spectra_dataframe.index.name = f"matches of spectrum #{i + 1}"
            predicted_spectra.append(library_spectra_dataframe)

        return predicted_spectra
