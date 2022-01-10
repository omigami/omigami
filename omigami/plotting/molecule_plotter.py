from json.decoder import JSONDecodeError

import itertools
import requests
from PIL.PngImagePlugin import PngImageFile
from matplotlib import pyplot as plt
from matplotlib.container import BarContainer
from rdkit import Chem
from rdkit.Chem import Draw
from typing import Tuple, Dict, List

import pandas as pd

from omigami.omi_settings import CLASSYFIRE_URL, NPCLASSIFIER_URL


class MoleculePlotter:
    def plot_molecule_structure(
        self,
        spectra_matches: pd.DataFrame,
        representation: str = "smiles",
        draw_indices: bool = False,
        img_size: Tuple[int, int] = (200, 200),
        substructure_highlight: str = "",
    ) -> (Dict[str, PngImageFile], List[str]):
        """
        Generate image representations of the hits returned from the spectra matching predictors outputs.
        All structures passed MUST have valid smiles or inchi representations.

        Parameters:
        ----------
        spectra_matches: DataFrame
            DataFrame resulting from spectra matching (e.g. Spec2Vec, MS2DeepScore).
            Needs to feature smiles or inchi, score and compound_name as columns.
        representation: str = 'smiles' or 'inchi'
            The representation of the molecules found in the provided dataframe.
        draw_indices: bool = False
            If true draws the indices of the atoms.
        molecule_image_size: Tuple[int, int] = [200, 200]
            The size of every individual image of a molecule. Need to be provided as a tuple with two ints (x and y).

        Returns:
        -------
            A tuple containing a dict with keys equal to the spectrum GNPS id and value equals to the molecule plot and
            a list of strings containing the legends for each plot.
        """

        img_size_list = list(img_size)
        self._validate_data(spectra_matches, representation)
        spectra_matches = self._clean_matches(spectra_matches, representation)

        try:
            substructure = Chem.MolFromSmarts(substructure_highlight)
            images = {}

            for i, structure in enumerate(spectra_matches[representation]):

                molecule = None
                if representation == "smiles":
                    molecule = Chem.MolFromSmiles(structure)
                elif representation == "inchi":
                    molecule = Chem.MolFromInchi(structure)

                highlight_bonds = self._get_bonds_to_highlight(molecule, substructure)

                if draw_indices:
                    molecule = self._add_index_to_atoms(molecule)

                gnps_id = spectra_matches.index.values[i]
                images[gnps_id] = Draw.MolToImage(
                    molecule,
                    size=img_size_list,
                    highlightBonds=highlight_bonds,
                )

            return images, spectra_matches.compound_name.tolist()

        except NameError:
            raise NameError(
                "You are missing the rdkit module. "
                "Please go to Omigami's README for instructions on how to install it."
            )

    @staticmethod
    def _validate_data(spectra_matches: pd.DataFrame, representation: str = "smiles"):
        if representation not in ["smiles", "inchi"]:
            raise ValueError(
                f"Got unexpected representation string. Needs to be either 'smiles' "
                f"or 'inchi' got {representation} "
            )

        if not isinstance(spectra_matches, pd.DataFrame):
            raise ValueError(
                f"Matches need to be a Pandas DataFrame got {type(spectra_matches)}"
            )

        if "compound_name" not in spectra_matches.columns:
            raise MandatoryColumnMissingError(
                "The provided DataFrame must contain a column named compound_name"
            )

    @staticmethod
    def _clean_matches(
        spectra_matches: pd.DataFrame, representation: str
    ) -> pd.DataFrame:
        """Drops all molecules that are missing a compound_name, are a duplicate
        structure or are missing a structure"""
        spectra_matches = spectra_matches.drop_duplicates("compound_name")
        spectra_matches = spectra_matches.drop_duplicates(representation)
        spectra_matches = spectra_matches[spectra_matches[representation] != ""]
        spectra_matches = spectra_matches.dropna(subset=[representation])

        return spectra_matches

    @staticmethod
    def _get_bonds_to_highlight(molecule, substructure) -> List[int]:
        """Gets the indexes of a molecule substructure as a single list. Where every
        match is represented by two consecutively ints"""
        substructure_matches = molecule.GetSubstructMatches(substructure)
        merged_list = list(itertools.chain(*substructure_matches))
        return merged_list

    # Author: Takayuki Serizawa
    # Original Source: https://iwatobipen.wordpress.com/2017/02/25/draw-molecule-with-atom-index-in-rdkit/
    @staticmethod
    def _add_index_to_atoms(molecule):
        for atom in molecule.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx())
        return molecule

    # Original Source: https://jcheminf.biomedcentral.com/articles/10.1186/s13321-016-0174-y
    @staticmethod
    def plot_classyfire_result(
        spectra_matches: pd.DataFrame, color="g"
    ) -> BarContainer:
        """Uses the ClassyFire API to classify and plot a barchart of the classifications"""
        class_stats = dict()

        smiles_list = spectra_matches["smiles"].to_list()

        for smiles in smiles_list:
            try:
                classyfire_result = requests.get(str(CLASSYFIRE_URL) + smiles).json()
                class_assignment = classyfire_result["class"]["name"]

                if class_assignment in class_stats.keys():
                    class_stats[class_assignment] += 1
                elif class_assignment not in class_stats.keys():
                    class_stats[class_assignment] = 1

            except KeyError:
                pass

        return plt.barh(list(class_stats.keys()), class_stats.values(), color=color)

    # Original Source: https://chemrxiv.org/engage/chemrxiv/article-details/60c74f58702a9ba8dc18bb6b
    @staticmethod
    def plot_NPclassifier_result(
        spectra_matches: pd.DataFrame, color="g"
    ) -> BarContainer:
        """Uses the NP-Classifier API to classify natural products."""
        class_stats = dict()
        class_stats["Cannot_Assign"] = 0

        smiles_list = spectra_matches["smiles"].to_list()

        for smiles in smiles_list:

            try:
                NPclassifier_result = requests.get(
                    str(NPCLASSIFIER_URL) + smiles
                ).json()
                class_assignment = NPclassifier_result["superclass_results"][0]

                if class_assignment in class_stats.keys():
                    class_stats[class_assignment] += 1
                elif class_assignment not in class_stats.keys():
                    class_stats[class_assignment] = 1

            except JSONDecodeError:
                class_stats["Cannot_Assign"] += 1
            except IndexError:
                class_stats["Cannot_Assign"] += 1

        if class_stats["Cannot_Assign"] == 0:
            del class_stats["Cannot_Assign"]

        return plt.barh(list(class_stats.keys()), class_stats.values(), color=color)


class MandatoryColumnMissingError(Exception):
    pass
