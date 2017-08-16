"""Contains the function for creating strings from PdbFiles."""

def pdb_file_to_pdb_string(pdb_file):
    """Takes an :py:class:`.PdbFile` and turns it into a string in
    .pdb format.

    :param PdbFile pdb_file: the PdbFile to convert.
    :rtype: ``str``"""

    return "\n".join([record.text().ljust(80) for record in pdb_file.records()])
