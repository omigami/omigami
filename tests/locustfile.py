from random import randrange
from typing import List

from locust import HttpUser, task
from matchms import Spectrum
from matchms.importing import load_from_mgf

from config import ENV
from omigami_client import Spec2VecClient
from omigami_client.spec2vec import JSON
from tests.conftest import ASSETS_DIR

# initial randomized payload setup

# settings
MAX_DIFFERENT_BATCHES = 5
BATCH_SIZE = 1000
TOKEN = ENV["token"].get()
MGF_PATH = str(ASSETS_DIR / "GNPS-NIST14-MATCHES.mgf")
PREDICT_ENDPOINT = (
    "https://mlops.datarevenue.com/seldon/seldon/spec2vec/api/v0.1/predictions"
)

# spec2vec client
spec2vec = Spec2VecClient("no-need-for-token-here")
spectra_generator = load_from_mgf(MGF_PATH)


def create_payloads_from_file() -> List[JSON]:
    # construct initial batches nested list
    batches: List[List[Spectrum]] = []
    for i in range(0, MAX_DIFFERENT_BATCHES):
        batches.append([])

    # assigns the next spectrum from generator to a random batch
    # until all batches are full
    full_batches = [-1]
    while True:
        # gets the list index of the batch to put spectrum in at random
        batch_ind = -1
        while batch_ind in full_batches:
            batch_ind = randrange(MAX_DIFFERENT_BATCHES)

        batches[batch_ind].append(next(spectra_generator))
        if len(batches[batch_ind]) == BATCH_SIZE:
            full_batches.append(batch_ind)
            if len(full_batches) - 1 == MAX_DIFFERENT_BATCHES:
                break

    payloads = []
    for batch in batches:
        payloads.append(spec2vec._build_payload(batch, 10))

    del batches
    return payloads


PAYLOADS = create_payloads_from_file()


class QuickstartUser(HttpUser):
    payload: JSON

    @task
    def predict_377(self):
        self.client.post(
            PREDICT_ENDPOINT,
            json=self.payload,
            headers={"Authorization": f"Bearer {TOKEN}"},
            timeout=600,
        )

    def on_start(self):
        self.payload = PAYLOADS[randrange(MAX_DIFFERENT_BATCHES)]
