from typing import List, Dict

import requests
from cryptography.fernet import Fernet

from omigami.exceptions import InvalidCredentials, NotFoundError, ServerAuthError
from omigami.omi_settings import get_credentials_path, ConfigurationError
from datetime import datetime


class Auth:
    credentials: Dict[str]
    session_token: str
    session_expiration: datetime

    def __init__(self):
        self.credentials = {}
        self.session_token = ""
        self.session_expiration = datetime.now()


# global auth information, to be used by all endpoints
AUTH = Auth()


# authentication related methods
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
        AUTH.token = _get_session_token_using_credentials(creds)
        del creds


def _get_configured_credentials() -> Dict[str]:
    path = get_credentials_path()

    credentials: List[str]
    with open(path, "r") as cred_file:
        credentials = cred_file.read().splitlines()
        if len(credentials) == 0:
            raise ConfigurationError(
                "You have not setup your credentials yet. "
                "Please do so by using 'omigami credentials-helper' CLI functionality and try again."
            )

    return {"u": credentials[0], "p": credentials[1], "k": credentials[2]}


def _decrypt_credentials() -> Dict[str]:
    f = Fernet(AUTH.credentials["k"])
    decripted = {
        "u": f.decrypt(AUTH.credentials["u"]),
        "p": f.decrypt(AUTH.credentials["p"]),
    }
    return decripted


def _get_session_token_using_credentials(credentials):
    flow = requests.get(
        "https://omigami.datarevenue.com/.ory/kratos/public/self-service/login/api"
    )
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

    if "session_token" not in api_request.json().keys:
        raise ServerAuthError(
            "The server did not responded accordingly. Please try again or contact DataRevenue for assistance."
        )

    AUTH.token = api_request.json()["session_token"]
    expiration = api_request.json()["session_token"]
    # removes microseconds from the string
    expiration = expiration.split(".", 0)[0]
    AUTH.expiration_date = datetime.strptime(expiration, "%Y-%m-%dT%H:%M:%S")
