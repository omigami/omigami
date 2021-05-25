<img src="./omigami-gradient.png" width="50%">

# Omigami Client

[![PyPI version shields.io](https://img.shields.io/pypi/v/omigami_client.svg)](https://pypi.python.org/pypi/omigami_client) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

<!-- image:: https://img.shields.io/travis/datarevenue-berlin/omigami.svg :target: https://travis-ci.org/datarevenue-berlin/omigami -->

<!-- image:: https://readthedocs.org/projects/omigami/badge/?version=latest :target: https://omigami.readthedocs.io/en/latest/?badge=latest :alt: Documentation Status -->


A client package to access Omigami services.

## Installation

```sh
pip install omigami_client
```

## Acknowledgement

### Spec2Vec 
Huber F, Ridder L, Verhoeven S, Spaaks JH, Diblen F, Rogers S, et al. (2021) Spec2Vec: Improved mass spectral similarity scoring through learning of structural relationships. PLoS Comput Biol 17(2): e1008724. https://doi.org/10.1371/journal.pcbi.1008724

## Motivation

TODO

## Features

- [x] Spec2Vec spectra matching
- [ ] MS2Deep score

## Usage

### Spec2Vec
#### Quickstart

```python
from omigami_client import Spec2VecClient

client = Spec2VecClient(token="your_token")
mgf_file_path = "path_to_file.mgf"
n_best_matches = 10

result = client.match_spectra_from_path(mgf_file_path, n_best_matches)
```

#### Notebooks
You can find a [tutorial](https://github.com/omigami/omigami_client/blob/master/notebooks/tutorial.ipynb) notebook and a [minimal example](https://github.com/omigami/omigami_client/blob/master/notebooks/minimal_example.ipynb) notebook in the `/notebooks/` folder.

## How it works

### Spec2Vec
1. Save your spectra data in a MGF file locally
2. Create an Spec2VecClient with your user token
3. Call `match_spectra_from_path` with the location of your mgf file.
4. The MGF spectra data will be processed and sent to the spec2vec model that will convert it into embeddings. 
5. These embeddings will be compared against the reference embeddings around the Precursor MZ.
6. The N best matches per spectrum are returned on the response as pandas dataframes.  

## Contribute to Omigami

1. Fork it (https://github.com/omigami/omigami_client/fork)
2. Create your feature branch (git checkout -b feature/fooBar)
3. Commit your changes (git commit -am 'Add some fooBar')
4. Push to the branch (git push origin feature/fooBar)
5. Create a new Pull Request

## License
MIT license - free software.