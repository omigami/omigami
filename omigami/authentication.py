import pickle
from typing import Dict

import requests
from cryptography.fernet import Fernet

from omigami.exceptions import InvalidCredentials, NotFoundError, ServerAuthError
from omigami.omi_settings import get_credentials_path, ConfigurationError
from datetime import datetime

"""
 > AUTH CLASS
"""


class Auth:
    """
    Class necessary to hold mutable session-related values
    Only one instance of this class is expected to exist
    """

    credentials: Dict[str, bytes]
    session_token: str
    session_expiration: datetime
    self_service_endpoint: str

    def __init__(self):
        self.credentials = {}
        self.session_token = ""
        self.session_expiration = datetime.now()
        self.self_service_endpoint = (
            "https://omigami.datarevenue.com/.ory/kratos/public/self-service/login/api"
        )


"""
 > AUTH: 
 > global singleton auth information, to be used inside this module by all endpoints 
"""
AUTH = Auth()


"""
 > PUBLIC FUNCTIONS
"""


def authenticate_client():
    """
    Will find the user's credentials as stored in his machine (or raise exception if not set) and
    use them to get a session token. If token already present under expiration date, do nothing.
    If credentials are already fetched, use them when a new token is necessary.
    """
    if not AUTH.credentials:
        AUTH.credentials = _get_configured_credentials()

    if not AUTH.session_token or AUTH.session_expiration <= datetime.now():
        creds = _decrypt_credentials()
        _get_session_token_using_credentials(creds)
        del creds


def encrypt_credentials(username: str, password: str) -> Dict[str, bytes]:
    """
    Utilizes Fernet library to generate a random key to encrypt credentials and save both creds and key to disk.
    """
    key = Fernet.generate_key()
    f = Fernet(key)
    creds = {
        "u": f.encrypt(username.encode("ascii")),
        "p": f.encrypt(password.encode("ascii")),
        "k": key,
    }
    return creds


def get_session() -> Auth:
    """
    Returns the singleton AUTH object
    """
    return AUTH


"""
 > PRIVATE FUNCTIONS
"""


def _get_configured_credentials() -> Dict[str, bytes]:
    """
    Get the encryupted credentials stored in disk
    """
    path = get_credentials_path()

    credentials: Dict[str, bytes]
    with open(path, "rb") as file_handle:
        credentials = pickle.load(file_handle)
        if len(credentials) == 0:
            raise ConfigurationError(
                "You have not setup your credentials yet. "
                "Please do so by using 'omigami credentials-helper' CLI functionality and try again."
            )
        if not all(key in ["k", "u", "p"] for key in credentials.keys()):
            raise ConfigurationError(
                "Something seems wrong with your credentials. "
                "Please, run 'omigami credentials-helper --unset' to remove them and then set them again."
            )

    return credentials


def _decrypt_credentials() -> Dict[str, str]:
    """
    Uses Fernet + key to decrypt encrypted credentials
    """
    f = Fernet(AUTH.credentials["k"])
    decripted = {
        "u": f.decrypt(AUTH.credentials["u"]).decode("ascii"),
        "p": f.decrypt(AUTH.credentials["p"]).decode("ascii"),
    }
    return decripted


def _get_session_token_using_credentials(credentials: Dict[str, str]) -> None:
    """
    Uses the AUTH object and decrypted credentials to get the self service endpoint, authenticate,
    and update the AUTH object with session information
    """
    flow = requests.get(AUTH.self_service_endpoint)
    action_url = flow.json()["ui"]["action"]
    api_request = requests.post(
        action_url,
        headers={"Content-Type": "application/json", "Accept": "application / json"},
        json={
            "password_identifier": credentials["u"],
            "password": credentials["p"],
            "method": "password",
        },
    )

    if api_request.status_code == 401:
        raise InvalidCredentials(
            "Your credentials are invalid, please revise your API token."
        )
    if api_request.status_code == 404:
        raise NotFoundError("The API endpoint couldn't be reached.")

    if "session_token" not in api_request.json().keys():
        raise ServerAuthError(
            "The server did not responded accordingly. Please try again or contact DataRevenue for assistance."
        )

    AUTH.session_token = api_request.json()["session_token"]
    expiration = api_request.json()["session"]["expires_at"]
    # removes microseconds from the string
    expiration = expiration.split(".")[0]
    AUTH.expiration_date = datetime.strptime(expiration, "%Y-%m-%dT%H:%M:%S")
