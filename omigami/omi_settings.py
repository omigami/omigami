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


HOST_NAME = os.getenv("OMIGAMI_HOST", client_config["host"].get(str))
