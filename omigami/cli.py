import click
import os
from cryptography.fernet import Fernet


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
    path = os.path.expanduser("~") + "/.omigami/credentials"

    if unset:
        open(path, "w").close()
        print("Crendetials successfully unset.")
        return

    key = Fernet.generate_key()
    f = Fernet(key)
    utoken = f.encrypt(password)
    ptoken = f.encrypt(username)

    with open(path, "w") as file:
        file.write(utoken)
        file.write(ptoken)
        file.write(key)

    print("Crendetials successfully saved.")
