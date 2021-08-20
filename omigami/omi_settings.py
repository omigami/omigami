import confuse
import os

config = confuse.Configuration("omigami", __name__)

CLASSYFIRE_URL = config["plotting"]["urls"]["classyfire"]
NPCLASSIFIER_URL = config["plotting"]["urls"]["NPclassifier"]


def get_credentials_path():
    return os.path.expanduser("~") + "/.omigami/credentials"


def get_credentials_folder_path():
    return os.path.expanduser("~") + "/.omigami"


class ConfigurationError(Exception):
    pass
