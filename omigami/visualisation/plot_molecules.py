from typing import List, Dict

import pandas as pd

from omigami.visualisation.config import VALID_KEYS
from omigami.visualisation.endpoint import Endpoint
from IPython.core.display import Image
from pandas import read_csv


class PlotMolecules(Endpoint):
    # TODO: URL?

    def __init__(self, token: str):
        super().__init__(token)
        self.endpoint = "https://omigami.datarevenue.com/seldon/seldon/utils/api/v0.1/plotmolecule"
        self.semi_optional_keys = ["smiles", "inchikey_inchi"]

    def _convert_df_to_dict_list(self, df: pd.DataFrame, meta_data: List[str]) -> List[Dict]:
        data = []

        for index, content in df.iterrows():
            data.append(content[meta_data].to_dict())

        return data

    def plot(self,
             scores_path: str,
             n_best: int,
             include_metadata: List[str] = ["smiles", "compound_name"],
             **kwargs,
             ) -> Image:

        # Checks if at least one of the optional keys is found in the included metadata
        if not bool(set(self.semi_optional_keys) & set(include_metadata)):
            raise ValueError(
                "Parameter include_metadata should be either include 'smiles' or 'inchikey_inchi'. "
                f"Defaults is ['smiles', 'compound_name']"
            )

        parameters = self._build_parameters(n_best, include_metadata)

        spectra_generator = read_csv(scores_path)
        spectra_generator = spectra_generator.sort_values(by="score", ascending=False)[:n_best]
        spectra_generator = self._convert_df_to_dict_list(spectra_generator, parameters["include_metadata"])

        self._make_batch_requests(spectra_generator, parameters, self.endpoint)
        return spectra_generator
