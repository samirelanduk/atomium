"""Contains the function for creating Models from PdbDataFiles."""

from ..structures.models import Model
from ..structures.atoms import Atom

def pdb_data_file_to_model(data_file):
    """Converts a :py:class:`.PdbDataFile` to a :py:class:`.Model`

    :param PdbDataFile data_file: The PdbDataFile.
    :rtype: ``Model``"""

    model = Model()
    return model


def atom_dict_to_atom(d):
    """Turns an atom ``dict`` into an actual :py:class:`.Atom`.

    :param dict d: The ``dict`` to convert.
    :rtype: ``Atom``"""
    
    atom = Atom(
     d["element"] if d["element"] else "X",
     d["x"], d["y"], d["z"],
     atom_id=d["atom_id"], name=d["atom_name"],
     charge=d["charge"] if d["charge"] else 0
    )
    atom.temp_chain = d["chain_id"]
    atom.temp_residue_id = d["chain_id"] + str(d["residue_id"])
    if d["insert_code"]: atom.temp_residue_id += d["insert_code"]
    atom.temp_residue_name = d["residue_name"]
    return atom
