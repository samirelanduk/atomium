"""Contains the functions for creating PdbDataFiles from Pdbs."""

from .structure2pdbdatafile import structure_to_pdb_data_file

def pdb_to_pdb_data_file(pdb):
    """Converts a :py:class:`.Pdb` to a :py:class:`.PdbDataFile`

    :param Pdb pdb: The Pdb.
    :rtype: ``PdbDataFile``"""

    data_files = [structure_to_pdb_data_file(m, model=i) for i, m in enumerate(pdb.models(), start=1)]
    for data_file in data_files[1:]:
        data_files[0].atoms += data_file.atoms
        data_files[0].heteroatoms += data_file.heteroatoms
    return data_files[0]
