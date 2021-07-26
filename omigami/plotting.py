from PIL.PngImagePlugin import PngImageFile
from pandas import DataFrame
from rdkit import Chem
from rdkit.Chem import Draw


def _raise_errors(spectra_matches: DataFrame, representation: str = 'smiles'):
    if representation not in ["smiles", "inchi"]:
        raise ValueError(
            f"Got unexpected representation string. Needs to be either 'smiles' or 'inchi' got {representation}"
        )

    if not isinstance(spectra_matches, type(DataFrame())):
        raise ValueError(
            f"Matches need to be a Pandas Dataframe got {type(spectra_matches)}"
        )


def _clean_matches(spectra_matches: DataFrame, representation: str) -> DataFrame:
    spectra_matches = spectra_matches.drop_duplicates('compound_name')
    spectra_matches = spectra_matches.drop_duplicates(representation)
    spectra_matches = spectra_matches[spectra_matches[representation] != ""]
    spectra_matches = spectra_matches[spectra_matches[representation].notna()]

    return spectra_matches


def plot_molecule_structure_grid(spectra_matches: DataFrame, representation: str = 'smiles',
                                 sort_values: bool = True) -> PngImageFile:
    """
    Generate a grid image representation of the hits returned from Spec2Vec and MS2DeepScore outputs.
    All structures passed MUST have valid smiles or inchi representations.

    Parameters:
    ----------
    spectra_matches: DataFrame
    representation: str = 'smiles' or 'inchi'
    sort_values: bool = True

    Returns:
        A Plot showing the structure of the passed smiles/inchis
    """

    _raise_errors(spectra_matches, representation)

    if sort_values:
        spectra_matches = spectra_matches.sort_values('score', ascending=True)

    spectra_matches = _clean_matches(spectra_matches, representation)

    structure_list = [spectra_matches.loc[x][representation] for x in spectra_matches.index.tolist()]

    mol_render_list = None
    if representation == 'smiles':
        mol_render_list = [Chem.MolFromSmiles(structure) for structure in structure_list]
    elif representation == 'inchi':
        mol_render_list = [Chem.MolFromInchi(structure) for structure in structure_list]

    image = Draw.MolsToGridImage(mol_render_list, legends=spectra_matches.compound_name.tolist())
    return image
