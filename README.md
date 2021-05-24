<img src="./omigami-gradient.png" width="50%">

# Omigami Client

[![PyPI version shields.io](https://img.shields.io/pypi/v/omigami_client.svg)](https://pypi.python.org/pypi/omigami_client) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

<!-- image:: https://img.shields.io/travis/datarevenue-berlin/omigami.svg :target: https://travis-ci.org/datarevenue-berlin/omigami -->

<!-- image:: https://readthedocs.org/projects/omigami/badge/?version=latest :target: https://omigami.readthedocs.io/en/latest/?badge=latest :alt: Documentation Status -->


TODO

## Installation

```sh
pip install omigami_client
```

## Acknowledgement

TODO
This package is based on an algorithm first introduced by Carl Brunius in [Variable selection and validation in multivariate modelling (2019)](https://academic.oup.com/bioinformatics/article/35/6/972/5085367).

**Citation**: Data Revenue, based on *Variable selection and validation in multivariate modelling (2019) [DOI:10.1093/bioinformatics/bty710](https://doi.org/10.1093/bioinformatics/bty710)*

## Motivation

- TODO

## Features

- [x] Spec2vec spectra matching
- [ ] MS2Deep score

## Usage

### A minimal example

TODO

```python
from omigami_client import OmigamiClient

client = OmigamiClient(token="your_token")
mgf_file_path = "path_to_mgf"
n_best_matches = 10

result = client.match_spectra_from_path(mgf_file_path, n_best_matches)
```

## How it works

1. TODO
2. TODO

## Contribute to Omigami

1. Fork it (https://github.com/omigami/omigami_client/fork)
2. Create your feature branch (git checkout -b feature/fooBar)
3. Commit your changes (git commit -am 'Add some fooBar')
4. Push to the branch (git push origin feature/fooBar)
5. Create a new Pull Request

## License
MIT license - free software.