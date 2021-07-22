from rdkit import Chem


def visualize_hits(
    spectrum_number: int, representation: str = "smiles", sort_values: bool = True
):
    """
    spectrum_number = int
    representation = 'smiles' or 'inchi'
    Generate a grid image representation of the hits returned from Spec2Vec and MS2DeepScore outputs.
    for RDkit to work, all structures passed MUST have valid smiles or inchi representations.

    """

    sorted_results = spectra_matches[spectrum_number].sort_values(
        "score", ascending=sort_values
    )  # sort for top matches
    sorted_results = sorted_results.drop_duplicates(
        "compound_name"
    )  # remove duplicate structures whos 'compound_name' is the same.
    sorted_results = sorted_results.drop_duplicates(
        "smiles"
    )  # remove duplicate structures whose 'compound_name' varies
    structure_list = [
        sorted_results.loc[x].smiles for x in sorted_results.index.tolist()
    ]
    structure_list.remove(
        ""
    )  # remove cases where spectra are not annotated with a SMILES string
    if representation == "smiles":
        mol_render_list = [
            Chem.MolFromSmiles(structure) for structure in structure_list
        ]
    elif representation == "inchi":
        mol_render_list = [Chem.MolFromInchi(structure) for structure in structure_list]
    for m in mol_render_list:
        tmp = Chem.AllChem.Compute2DCoords(m)
    image = Chem.Draw.MolsToGridImage(
        mol_render_list, legends=sorted_results.compound_name.tolist()
    )
    return image
