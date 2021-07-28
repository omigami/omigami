from typing import List

from PIL.PngImagePlugin import PngImageFile
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.rdchem import Mol
import itertools
import matplotlib.pyplot as plt
import requests

from omigami.config import CLASSYFIRE_URL, NPCLASSIFER_URL


class MandatoryColumnMissingError(Exception):
    pass


# TODO: Dont know about the name though
class MSPlots:

    def plot_molecule_structure_grid(self, spectra_matches: pd.DataFrame, representation: str = 'smiles',
                                     sort_by_score: bool = True, draw_indices: bool = False,
                                     molecule_image_size: List[int] = [200, 200],
                                     substructure_highlight: str = "") -> PngImageFile:
        """
        Generate a grid image representation of the hits returned from Spec2Vec and MS2DeepScore outputs.
        All structures passed MUST have valid smiles or inchi representations.

        Parameters:
        ----------
        spectra_matches: DataFrame
            DataFrame resulting from either Spec2Vec or MS2DeepScore. Need to feature smiles or inchi, score and compound_name as columns.
        representation: str = 'smiles' or 'inchi'
            The representation of the molecules found in the provided dataframe
        sort_by_score: bool = True
            If true sorts the dataframe by the score column
        draw_indices: bool = False
            If true draws the indices of the atoms
        molecule_image_size: List[int, int] = [200, 200]
            The size of every individual image of a molecule. Need to be provided as a list with two ints
        substructure_highlight: str = None
            Needs to be a molecule substructure represented as a smiles.

        Returns:
        -------
            A Plot showing the structure of the passed smiles/inchis
        """

        self._validate_data(spectra_matches, representation)

        if sort_by_score:
            spectra_matches = spectra_matches.sort_values('score', ascending=True)

        spectra_matches = self._clean_matches(spectra_matches, representation)

        substructure = Chem.MolFromSmarts(substructure_highlight)
        mol_render_list = []
        highlight_bonds = []

        for structure in spectra_matches[representation]:

            if representation == 'smiles':
                molecule = Chem.MolFromSmiles(structure)
            elif representation == 'inchi':
                molecule = Chem.MolFromInchi(structure)

            substructure_matches = self._get_bonds_to_highlight(molecule, substructure)
            highlight_bonds.append(substructure_matches)

            if draw_indices:
                molecule = self._mol_with_atom_index(molecule)

            mol_render_list.append(molecule)

        image = Draw.MolsToGridImage(mol_render_list, subImgSize=molecule_image_size,
                                     legends=spectra_matches.compound_name.tolist(),
                                     highlightBondLists=highlight_bonds)
        return image

    @staticmethod
    def plot_classyfire_result(smiles_list: List[str], color="g"):
        """TODO: Ask joe what those classifiers actually classfy"""
        class_stats = dict()
        for smiles in smiles_list:
            try:
                classyfire_result = requests.get(CLASSYFIRE_URL + smiles).json()
                class_assignment = classyfire_result['class']['name']

                if class_assignment in class_stats.keys():
                    class_stats[class_assignment] += 1
                elif class_assignment not in class_stats.keys():
                    class_stats[class_assignment] = 1
            except:
                pass

        return plt.barh(list(class_stats.keys()), class_stats.values(), color=color)

    @staticmethod
    def plot_NPclassifier_result(smiles_list: List[str], color="g"):
        class_stats = dict()
        class_stats['Cannot_Assign'] = 0
        for smiles in smiles_list:

            try:
                NPclassifier_result = requests.get(NPCLASSIFER_URL + smiles).json()
                class_assignment = NPclassifier_result['superclass_results'][0]
                if class_assignment in class_stats.keys():
                    class_stats[class_assignment] += 1
                elif class_assignment not in class_stats.keys():
                    class_stats[class_assignment] = 1

            except:
                class_stats['Cannot_Assign'] += 1
        return plt.barh(list(class_stats.keys()), class_stats.values(), color=color)

    @staticmethod
    def _validate_data(spectra_matches: pd.DataFrame, representation: str = 'smiles'):
        if representation not in ["smiles", "inchi"]:
            raise ValueError(
                f"Got unexpected representation string. Needs to be either 'smiles' or 'inchi' got {representation}"
            )

        if not isinstance(spectra_matches, pd.DataFrame):
            raise ValueError(
                f"Matches need to be a Pandas DataFrame got {type(spectra_matches)}"
            )

        if "compound_name" not in spectra_matches.columns:
            raise MandatoryColumnMissingError("The provided DataFrame must contain a column named compound_name")

    @staticmethod
    def _clean_matches(spectra_matches: pd.DataFrame, representation: str) -> pd.DataFrame:
        spectra_matches = spectra_matches.drop_duplicates('compound_name')
        spectra_matches = spectra_matches.drop_duplicates(representation)
        spectra_matches = spectra_matches[spectra_matches[representation] != ""]
        spectra_matches = spectra_matches.dropna()

        return spectra_matches

    @staticmethod
    def _get_bonds_to_highlight(molecule: Mol, substructure: Mol) -> List[int]:
        substructure_matches = molecule.GetSubstructMatches(substructure)
        merged_list = list(itertools.chain(*substructure_matches))
        return merged_list

    # Author: Takayuki Serizawa
    # Original Source: https://iwatobipen.wordpress.com/2017/02/25/draw-molecule-with-atom-index-in-rdkit/
    @staticmethod
    def _mol_with_atom_index(molecule: Mol) -> Mol:
        for atom in molecule.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx())
        return molecule
