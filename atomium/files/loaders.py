from requests import get
from .pdbstring2pdbdict import pdb_string_to_pdb_dict
from .pdbdict2pdb import pdb_dict_to_pdb

def string_from_file(path):
    """Opens a file from the given path and returns the contents as a string.

    :param str path: The path to the file.
    :rtype: ``str``"""

    with open(path) as f:
        return f.read()


def fetch_string(code, pdbe=False):
    """Gets a the filestring of a PDB from the RCSB web services.

    :param str code: The PDB code to fetch.
    :param bool pdbe: If ``True``, the PDB will instead be fetched from PDBe.
    :raises TypeError: if the code is not a string.
    :raises ValueError: if the code is not four caracters long.
    :rtype: ``str``"""

    if not isinstance(code, str):
        raise TypeError("PDB code {} is not string".format(code))
    if len(code) != 4:
        raise ValueError("PDB code {} is not of length 4".format(code))
    url = "https://files.rcsb.org/view/{}.pdb"
    if pdbe:
        url = "https://www.ebi.ac.uk/pdbe/entry-files/pdb{}.ent"
    response = get(url.format(code.lower()))
    if response.status_code == 200:
        return response.text


def pdb_data_from_file(path):
    filestring = string_from_file(path)
    return pdb_string_to_pdb_dict(filestring)


def fetch_data(code, **kwargs):
    filestring = fetch_string(code, **kwargs)
    return pdb_string_to_pdb_dict(filestring)


def pdb_from_file(path):
    pdb_dict = pdb_data_from_file(path)
    return pdb_dict_to_pdb(pdb_dict)


def fetch(code, **kwargs):
    pdb_dict = fetch_data(code, **kwargs)
    return pdb_dict_to_pdb(pdb_dict)
