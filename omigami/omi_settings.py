import os

import confuse

config = confuse.Configuration("omigami", __name__)
client_config = config["client"]

CLASSYFIRE_URL = client_config["plotting"]["urls"]["classyfire"]
NPCLASSIFIER_URL = client_config["plotting"]["urls"]["NPclassifier"]


def get_credentials_path():
    return os.path.expanduser("~") + "/.omigami/credentials"


def get_credentials_folder_path():
    return os.path.expanduser("~") + "/.omigami"


class ConfigurationError(Exception):
    pass


OMIGAMI_ENV = os.getenv("OMIGAMI_ENV", "prod")
HOST_NAME = client_config["host"][OMIGAMI_ENV].get(str)
