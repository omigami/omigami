# MS2DeepScore
## Quickstart

```python
from omigami.spectra_matching import MS2DeepScore

client = MS2DeepScore()

mgf_file_path = "path_to_file.mgf" # or a list of matchms.Spectrum objects
n_best_matches = 10
ion_mode = "positive"  # either positive or negative

result = client.match_spectra(
    mgf_file_path, n_best_matches, ion_mode,
)

# Spectra is sent in batches to the predictor. 
# If any of the batches fail, the following command will return the list of spectra in the failed batch. 
# The problem may have been caused by just one or more spectrum inside the failed batch.
failed_spectra = client.failed_spectra

# successful spectrum predictions will be saved to cache, to reset the cache use:
client.reset_cache()
```


## Notebook
You can find a [tutorial](https://github.com/omigami/omigami/blob/master/notebooks/ms2deepscore/tutorial.ipynb) notebook in the `/notebooks/` folder.
