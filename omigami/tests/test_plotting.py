import matplotlib
import pandas as pd
import pytest
from PIL.PngImagePlugin import PngImageFile
from rdkit import Chem
from omigami import MSPlotter


def test_plot_molecule_structure_grid(spectra_match_data_path):
    dataset = pd.read_csv(spectra_match_data_path)
    plotter = MSPlotter()
    image = plotter.plot_molecule_structure_grid(dataset)

    assert image.size == (600, 400)
    assert isinstance(image, PngImageFile)


def test_plot_molecule_highlighting(spectra_match_data_path):
    dataset = pd.read_csv(spectra_match_data_path)
    molecule = Chem.MolFromSmiles(dataset["smiles"][2])
    substructure = Chem.MolFromSmarts("C(=O)")

    substructure_match = MSPlotter._get_bonds_to_highlight(molecule, substructure)

    assert substructure_match == [2, 3]


def test_plot_clean_matches(spectra_match_data_path_missing_smiles):
    dataset = pd.read_csv(spectra_match_data_path_missing_smiles)

    results = MSPlotter._clean_matches(dataset, "smiles")

    assert len(results) == 5
    assert isinstance(results, pd.DataFrame)

@pytest.mark.internet_connection
def test_plot_classyfire_result(spectra_match_data_path):

    spectra_match_data_path = pd.read_csv(spectra_match_data_path)
    smiles = spectra_match_data_path["smiles"]

    plot = MSPlotter.plot_classyfire_result(smiles)

    assert isinstance(plot, matplotlib.container.BarContainer)


@pytest.mark.internet_connection
def test_plot_NPclassifier_result(spectra_match_data_path):

    spectra_match_data_path = pd.read_csv(spectra_match_data_path)
    smiles = spectra_match_data_path["smiles"]

    plot = MSPlotter.plot_NPclassifier_result(smiles)

    assert isinstance(plot, matplotlib.container.BarContainer)
