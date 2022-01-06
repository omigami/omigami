![image](https://raw.githubusercontent.com/omigami/omigami/master/omigami-gradient.png)

# Omigami

[![PyPI version shields.io](https://img.shields.io/pypi/v/omigami-client.svg)](https://pypi.python.org/pypi/omigami-client) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

<!-- image:: https://img.shields.io/travis/datarevenue-berlin/omigami.svg :target: https://travis-ci.org/datarevenue-berlin/omigami -->

<!-- image:: https://readthedocs.org/projects/omigami/badge/?version=latest :target: https://omigami.readthedocs.io/en/latest/?badge=latest :alt: Documentation Status -->


A package to access Omigami services.

## Installation
Omigami officially supports python 3.7 up to 3.9. Support for newer versions remains to be tested. To install it, simply:

```sh
pip install omigami-client
```
To be able to visualise the molecule plots, you need to install `rdkit`.
For Mac and Linux users, please do:

```sh
pip install "omigami-client[plots]"
```
For Windows users, you need to do the following:

```sh
conda install -c conda-forge rdkit
```

## Acknowledgement

### Spec2Vec
Huber F, Ridder L, Verhoeven S, Spaaks JH, Diblen F, Rogers S, et al. (2021) Spec2Vec: Improved mass spectral similarity scoring through learning of structural relationships. PLoS Comput Biol 17(2): e1008724. https://doi.org/10.1371/journal.pcbi.1008724

### MS2DeepScore
Florian Huber, Sven van der Burg, Justin J.J. van der Hooft, Lars Ridder. (2021) MS2DeepScore - a novel deep learning similarity measure for mass fragmentation spectrum comparisons. bioRxiv 2021, doi: https://doi.org/10.1101/2021.04.18.440324

## Motivation
Motivation
We aim to support metabolomics research by providing the following :
- Easy-to-use tools
- Access and scalability to the newest algorithms
- Maintenance, support and documentation of software, models and data
- A community of metabolomics researchers via our [Slack](https://join.slack.com/t/ml4metabolomics/shared_invite/zt-r39udtdg-G36YE6GQt1YdVIwTdeT8aw)

## Features

- [x] Spec2Vec spectra matching
- [x] MS2Deep score

## Usage

### Configuration
Before start using Omigami, you will need to signup at omigami.com and verify your email to confirm your account.
Then, you should configure your credentials in your working machine.
You can do so by using the command line utility `omigami credentials-helper`
and providing the parameters `--username <yourusername> --password <yourpassword>`

Your encrypted credentials will be safely stored in your machine, if you want to revoke the configuration at any point,
you may use `omigami credentials-helper --unset`.

## How it works

### Spec2Vec
1. Save your spectra data in a MGF file locally
2. Create an Spec2Vec with your user token
3. Call `match_spectra` with the location of your mgf file.
4. The MGF spectra data will be processed and sent to the spec2vec model that will convert it into embeddings.
5. These embeddings will be compared against the reference embeddings around the Precursor MZ.
6. The N best matches per spectrum are returned on the response as pandas dataframes.

### MS2DeepScore
1. Save your pair of spectra data in a MGF file locally
2. Create an MS2DeepScore object with your user token
3. Call `match_spectra` with the location of your mgf file.
4. The MGF spectra data will be processed and sent to the trained neural network that will predict the molecular structural similarity.
5. The prediction is returned on the response as a pandas dataframe.  
