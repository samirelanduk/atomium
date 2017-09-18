"""This module contains tools for processing strings into more manageable
forms, and back into strings."""

def string2lines(s):
    """Takes a string and breaks it into lines on line breaks. It will first
    convert any windows line breaks to unix line breaks.

    :param str s: The string to break up.
    :rtype: ``str``."""

    lines = s.replace("\r\n", "\n").split("\n")
    while not lines[-1].strip():
        lines.pop()
    return lines


def string_from_file(path):
    """Opens a file from the given path and returns the contents as a string.

    :param str path: The path to the file.
    :rtype: ``str``"""

    with open(path) as f:
        return f.read()



def string_to_file(string, path):
    """Saves a string to a given path as a file.

    :param str string: The string to save.
    :param str path: The file to save it in."""

    with open(path, "w") as f:
        f.write(string)


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
    from requests import get
    url = "https://files.rcsb.org/view/{}.pdb"
    if pdbe:
        url = "https://www.ebi.ac.uk/pdbe/entry-files/pdb{}.ent"
    response = get(url.format(code.lower()))
    if response.status_code == 200:
        return response.text
