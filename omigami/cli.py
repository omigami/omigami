import click
from cryptography.fernet import Fernet
from omi_settings import get_credentials_path


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

    if unset:
        open(path, "w").close()
        print("Crendetials successfully unset.")
        return

    key = Fernet.generate_key()
    f = Fernet(key)
    utoken = f.encrypt(username)
    ptoken = f.encrypt(password)

    with open(path, "w") as file:
        print(utoken, file=file)
        print(ptoken, file=file)
        print(key, file=file)

    print("Crendetials successfully saved.")
