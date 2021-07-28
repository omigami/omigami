import confuse

config = confuse.Configuration("omigami", __name__)

CLASSYFIRE_URL = config["plotting"]["urls"]["classyfire"]
NPCLASSIFER_URL = config["plotting"]["url"]["NPclassifier"]