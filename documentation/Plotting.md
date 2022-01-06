# Plotting
Plotting graphs works the same way for both Spec2Vec and MS2DeepScore.

## Molecule Plot
The following example will plot the structures of the molecules.
```python
from omigami.plotting import MoleculePlotter
from omigami.spectra_matching import Spec2Vec

client = Spec2Vec()
mgf_file_path = "path_to_file.mgf"
best_matches = client.match_spectra(mgf_file_path, n_best=10)

plotter = MoleculePlotter()
plots, legends = plotter.plot_molecule_structure(
    spectra_matches=best_matches[1],
    representation="smiles",
    draw_indices=True,
    img_size=(600, 600),
    substructure_highlight="C(=O)"
)
first_match = list(plots.values())[0]
first_match
print(legends[0])
```
![](https://raw.githubusercontent.com/omigami/omigami/release/0.3.0/documentation/images/molecule_plot.png)

## ClassyFire Plot
The following code allows us to plot the results of the [Classyfire](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-016-0174-y) model API.
```python
from omigami.plotting import MoleculePlotter
plotter = MoleculePlotter()
plotter.plot_classyfire_result(best_matches[1], color="green")
```
![](https://raw.githubusercontent.com/omigami/omigami/release/0.3.0/documentation/images/classyfire_plot.png)

## NPClassifier Plot
Furthermore, Omigami provides the possibility to use the [NPClassifier](https://www.researchgate.net/publication/344008670_NPClassifier_A_Deep_Neural_Network-Based_Structural_Classification_Tool_for_Natural_Products) API.

```python
from omigami.plotting import MoleculePlotter
plotter = MoleculePlotter()
plotter.plot_NPclassifier_result(best_matches[1], color="orange")
```
![](https://raw.githubusercontent.com/omigami/omigami/release/0.3.0/documentation/images/NP_classifier_plot.png)


## Spectra Comparison Mirror Plot
The following code snippet will plot the Spectra Comparison Mirror Plot. It 
takes two spectra returned and  creates a mirror plot with spectrum 1 on top 
(blue), and spectrum 2 on the bottom (red).
```python
from omigami.plotting import SpectraComparisonPlotter
from omigami.utilities import SpectrumDataFrameHelper
from matchms.importing import load_from_mgf

plotter = SpectraComparisonPlotter()

# lets load in memory our input spectra so we can compare input and output (matches)
input_spectra = list(load_from_mgf(mgf_file_path))

input_df = SpectrumDataFrameHelper.from_spectrum(input_spectra[0])

spectrum_match_id = best_matches[0].index[0]
output_df = SpectrumDataFrameHelper.from_gnps_id(spectrum_match_id)

plot = plotter.mirror_plot(
    spectrum_1=input_df,
    spectrum_2=output_df,
    labels=["Input Spectrum", "First Spectrum Match"],
    display_limits=(0,500)
)
```
![](https://raw.githubusercontent.com/omigami/omigami/release/0.3.0/documentation/images/mirror_plot.png)
