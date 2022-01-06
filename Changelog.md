# Changelog

## Releases
***
v0.3.0
***
Release Date: 2022-01-10

### New Features
- **Mirror plot**: Add `SpectraComparisonPlotter().mirror_plot()` to `omigami.plotting` 
module. It takes two spectra returned from `match_spectra()` and  creates a mirror 
plot with spectrum 1 on top (blue), and spectrum 2 on the bottom (red).
- **Caching functionality**: Both `Spec2Vec` and `MS2DeepScore`objects cache 
successful predictions, also there are 2 new methods:
  - if predictions for spectra fails, they can be called by `faied_spectra()`.
  - cache contains successful predictions of spectra, implemented `reset_cache()`
    to reset it.
- **Implement SpectrumDataFrameHelper:** Accessible from `omigami.utilities`, 
this class serves as a helper to create and manipulate Pandas DataFrame versions 
of Spectrum objects, useful for some plotting functions (such as `mirror_plot()`).

### Bugfixes
- **More detailed error messages**: Omigami was unable to catch errors from 
backend during `match_spectra()` call. Omigami properly handles errors and 
displays them.

Previous error messages lacking the cause of the error. E.g. 
> InternalServerError: Something went wrong, the requested service is 
> probably unavailable at the moment. Please try again later or contact DataRevenue 
> for more information.

Current error message displays the cause. E.g.
> InternalServerError: Error in generating the prediction: SpectraMatchingError
> 500: 'No data found from filtering with precursor MZ for spectra at indices 
> [5, 6, 7, 39, 40, 42]. Try increasing mz_range filtering.'. Please try again
> later or contact DataRevenue for more information.


### Breaking Changes
- **Rename`match_spectra_from_path()`**: available as `match_spectra()`
- **Update `match_spectra()` arguments**: Support passing a list of 
`matchms.Spectrum` objects, in addition to the path to the mgf file. Current 
arguments are:
  - `source`: path to mgf file or a list of `matchms.Spectrum` objects
  - `n_best`: number of best matches
  - `ion_mode`: `positive` or `negative`
- **Move plotting objects to `plotting` module**: `MoleculePlotter` is available
in `omigami.plotting` instead of `omigami`. To import it please do 
`from omigami.plotting import MoleculePlotter`
- **Rename `plot_molecule_structure_grid()`**: available as `plot_molecule_structure()`
- **Update `plot_molecule_structure()` arguments**: Rename `molecule_image_size` 
argument to `img_size`, changed `img_size` type from list to tuple.