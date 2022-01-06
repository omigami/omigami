# Spec2Vec
## Quickstart

```python
from omigami.spectra_matching import Spec2Vec

client = Spec2Vec()

mgf_file_path = "path_to_file.mgf" # or a list of matchms.Spectrum objects
n_best_matches = 10
ion_mode = "positive"  # either positive or negative

result = client.match_spectra(
    mgf_file_path, n_best_matches, ion_mode,
)

# Spectra is sent in batches to the predictor. 
# If any of the batches fail, the following command will return the list of spectra in the failed batch. 
# The problem may have been caused by just one or more spectrum inside the failed batch.
failed_spectra = client.failed_spectra()

# successful spectrum predictions will be saved to cache, to reset the cache use:
client.reset_cache()
```

The supported metadata keys for omigami are (case insensitive):
- "smiles",
- "compound_name",
- "instrument",
- "parent_mass",
- "inchikey_smiles",
- "inchikey_inchi"
- "precursor_mz" or "pepmass"


## Notebook
You can find a spec2vec tutorial [here](https://github.com/omigami/omigami/blob/master/notebooks/spec2vec/tutorial.ipynb) and ms2deepscore tutorial [here](https://github.com/omigami/omigami/blob/master/notebooks/ms2deepscore/tutorial.ipynb).
