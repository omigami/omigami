import pickle
from pathlib import Path

import click
from click import ClickException

from omigami.authentication import encrypt_credentials
from omigami.omi_settings import get_credentials_path, get_credentials_folder_path


@click.command()
@click.option("--username", help="Your Omigami.com username")
@click.option("--password", help="Your Omigami.com password")
@click.option(
    "--unset",
    is_flag=True,
    help="Remove previously setup credentials for the current user",
)
def credentials_helper(username, password, unset):
    """
    CLI Helper for configuring the user machine to access the Omigami endpoints
    """
    path = get_credentials_path()
    Path(get_credentials_folder_path()).mkdir(parents=True, exist_ok=True)

    if unset:
        open(path, "w").close()
        print("Crendetials successfully unset.")
        return

    if not username or not password:
        raise ClickException(
            "Please provide username and password using --username and --password arguments, "
            "placing values between single quotes is recommended"
        )

    creds = encrypt_credentials(username, password)

    with open(path, "wb") as file_handle:
        pickle.dump(creds, file_handle, protocol=pickle.HIGHEST_PROTOCOL)

    print("Crendetials successfully saved.")


@click.group()
def omigami():
    pass


omigami.add_command(credentials_helper)
