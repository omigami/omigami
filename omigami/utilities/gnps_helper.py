import requests
import pandas as pd


class GnpsHelper():

    def __init__(self):
        pass


    @staticmethod
    def fetch_metadata(spectrum_id: str) -> pd.DataFrame:
        """
        Obtains spectrum metadata from GNPS using a spectrum ID.
        """

        gnps_metadata = requests.get(
            "https://gnps.ucsd.edu/ProteoSAFe/SpectrumCommentServlet?SpectrumID="
            + str(spectrum_id)
        ).json()
        gnps_metadata_df = pd.DataFrame(
            eval(gnps_metadata["spectruminfo"]["peaks_json"]), columns=["m/z", "Intensity"]
        )

        return gnps_metadata_df
