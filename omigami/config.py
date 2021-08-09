import confuse

config = confuse.Configuration("omigami", __name__)

CLASSYFIRE_URL = config["plotting"]["urls"]["classyfire"]
NPCLASSIFIER_URL = config["plotting"]["urls"]["NPclassifier"]