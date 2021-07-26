import pytest

from omigami.visualisation.plot_molecules import PlotMolecules

from tests.conftest import prediction_dataset_missing_smiles


@pytest.mark.internet_connection
# @pytest.mark.skip(reason="Requires a valid token")
def test_plot(prediction_dataset_missing_smiles):
    plot_request = PlotMolecules("key")

    response = plot_request.plot(dataset_path=prediction_dataset_missing_smiles, n_best=5)

    assert response
