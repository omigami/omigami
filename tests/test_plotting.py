import pytest
import pandas as pd
# TODO: In the Omigami-core repo the test folder is found inside of the omigami folder and is called test instead of
#  testS. Should we change that
from PIL.PngImagePlugin import PngImageFile

from tests.conftest import spectra_match_data_path, spectra_match_data_path_missing_smiles

from omigami import plotting


def test_plot_molecule_structure_grid(spectra_match_data_path):

    dataset = pd.read_csv(spectra_match_data_path)

    image = plotting.plot_molecule_structure_grid(dataset)

    assert image.size == (600, 400)
    assert isinstance(image, PngImageFile)
    # TODO: Data a pixel comparison?


def test_plot_clean_matches(spectra_match_data_path_missing_smiles):

    dataset = pd.read_csv(spectra_match_data_path_missing_smiles)

    results = plotting._clean_matches(dataset, "smiles")

    assert len(results) == 5
    assert isinstance(results, pd.DataFrame)