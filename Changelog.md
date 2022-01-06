# Changelog

## Omigami 0.3.0
### Summary of changes
- renamed `match_spectra_from_path` to `match_spectra`
- updated function arguments of `match_spectra`
- added caching functionality to predictions
- moved plotting functionality to `omigami.plotting`
- added mirror plot

  
### Details

- renamed `match_spectra_from_path` to `match_spectra`, it takes following as 
arguments:
  - `source`: path to mgf file or a list of `matchms.Spectrum` objects
  - `n_best`: number of best matches
  - `ion_mode`: `positive` or `negative`
- caching of successful spectrum predictions:
  - if predictions for spectra fails, they can be called by `faied_spectra()` 
  - cache contains successful predictions of spectra, implemented `reset_cache()` 
to reset it
- moved plotting functionality to `omigami.plotting`
  - `MoleculePlotter` is available in `omigami.plotting` instead of `omigami`
- added `SpectraComparisonPlotter().mirror_plot` available in `omigami.plotting`:
  - It takes two spectra returned from `match_spectra` and creates a mirror plot 
with spectrum 1 on top (blue), and spectrum 2 on the bottom (red)