import pytest

from omigami.visualisation.plot_molecules import PlotMolecules

from tests.conftest import prediction_set_path


@pytest.mark.internet_connection
# @pytest.mark.skip(reason="Requires a valid token")
def test_plot(prediction_set_path):
    plot_request = PlotMolecules("key")

    response = plot_request.plot(scores_path=prediction_set_path, n_best=5)

    assert response
