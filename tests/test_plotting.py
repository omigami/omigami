import matplotlib
import pandas as pd
import pytest
from PIL.PngImagePlugin import PngImageFile
from matchms.importing import load_from_mgf
from rdkit import Chem
from omigami.plotting import MoleculePlotter, SpectraComparisonPlotter
from omigami.utilities import SpectrumDataFrameHelper


def test_plot_molecule_structure(spectra_match_data_path):
    dataset = pd.read_csv(spectra_match_data_path)
    plotter = MoleculePlotter()
    images, legends = plotter.plot_molecule_structure(dataset, img_size=(200, 200))

    assert images[0].size == (200, 200)
    assert isinstance(images[0], PngImageFile)
    assert isinstance(legends[0], str)


def test_plot_molecule_highlighting(spectra_match_data_path):
    dataset = pd.read_csv(spectra_match_data_path)
    molecule = Chem.MolFromSmiles(dataset["smiles"][2])
    substructure = Chem.MolFromSmarts("C(=O)")

    substructure_match = MoleculePlotter._get_bonds_to_highlight(molecule, substructure)

    assert substructure_match == [2, 3]


def test_plot_clean_matches(spectra_match_data_path_missing_smiles):
    dataset = pd.read_csv(spectra_match_data_path_missing_smiles)

    results = MoleculePlotter._clean_matches(dataset, "smiles")

    assert len(results) == 5
    assert isinstance(results, pd.DataFrame)


@pytest.mark.internet_connection
def test_plot_classyfire_result(spectra_match_data_path_web_api_error):
    spectra_match_data_path = pd.read_csv(spectra_match_data_path_web_api_error)

    plot = MoleculePlotter.plot_classyfire_result(spectra_match_data_path)

    assert isinstance(plot, matplotlib.container.BarContainer)


@pytest.mark.internet_connection
def test_plot_NPclassifier_result(spectra_match_data_path_web_api_error):
    spectra_match_data_path = pd.read_csv(spectra_match_data_path_web_api_error)

    plot = MoleculePlotter.plot_NPclassifier_result(spectra_match_data_path)

    assert isinstance(plot, matplotlib.container.BarContainer)


def test_mirror_plot(mgf_46_spectra_path):
    spectra = list(load_from_mgf(mgf_46_spectra_path))
    spectra_df = SpectrumDataFrameHelper.from_spectra_list(spectra)

    plotter = SpectraComparisonPlotter()
    fig = plotter.mirror_plot(spectra_df[0], spectra_df[1])

    assert isinstance(fig, matplotlib.pyplot.Figure)
