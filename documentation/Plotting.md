# Plotting
Plotting graphs works the same way for both Spec2Vec and MS2DeepScore.


The following example will plot the structures of the molecules
```python
from omigami.plotting import MoleculePlotter
plotter = MoleculePlotter()
result = plotter.plot_molecule_structure_grid(
    spectra_matches=best_matches,
    representation="smiles",
    draw_indices=True,
    img_size=(600, 600),
    substructure_highlight="C(=O)"
)
# result is a tuple of 
# a plot identified by spectrum_id and 
# a list of compound names of the best matches
```

The following code allows us to plot the results of the [Classyfire](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-016-0174-y) model API.
```python
from omigami.plotting import MoleculePlotter
plotter = MoleculePlotter()
plotter.plot_classyfire_result(best_matches, color="green")
```
<img src="docs/images/classyfire_plot.png" width="500">

Furthermore, Omigami provides the possibility to use the [NPClassifier](https://www.researchgate.net/publication/344008670_NPClassifier_A_Deep_Neural_Network-Based_Structural_Classification_Tool_for_Natural_Products) API.

```python
from omigami.plotting import MoleculePlotter
plotter = MoleculePlotter()
plotter.plot_NPclassifier_result(best_matches, color="orange")
```
<img src="docs/images/NP_classifier_plot.png" width="500">
