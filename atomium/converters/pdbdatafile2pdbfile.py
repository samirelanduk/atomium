"""Contains the function for creating PdbFiles from PdbDataFiles."""

from math import ceil
from ..parse.pdbfile import PdbRecord

def pack_structure(data_file, pdb_file):
    """Adds the :py:class:`.PdbRecord` objetcs to a :py:class:`.PdbFile`
    pertaining to structure, from a :py:class:`.PdbDataFile`.

    :param PdbDataFile data_file: The source PdbDataFile.
    :param PdbFile pdb_file: The PdbFile to update."""

    for atom in data_file.atoms:
        pdb_file._records.append(atom_dict_to_record(atom))
    for atom in data_file.heteroatoms:
        pdb_file._records.append(atom_dict_to_record(atom, hetero=True))
    pdb_file._records += conections_list_to_records(data_file.connections)


def atom_dict_to_record(d, hetero=False):
    """Converts an atom ``dict`` to an ATOM or HETATM :py:class:`.PdbRecord`.

    :param dict d: The ``dict`` to pack.
    :param bool hetero: If ``True``, a HETATM record will be made. If ``False``\
    (the default), an ATOM record will be made.
    :rtype: ``PdbRecord``"""

    line = "{:6}{:5} {:4}{:1}{:3} {:1}{:4}{:1}   "
    line += "{:8}{:8}{:8} {:.3f}{:6}          {:>2}{:2d}"
    atom_name = d.get("atom_name", "") if d.get("atom_name") else ""
    atom_name = " " + atom_name if len(atom_name) < 4 else atom_name
    line = line.format(
     "HETATM" if hetero else "ATOM",
     d.get("atom_id", "") if d.get("atom_id") else "",
     atom_name,
     d.get("alt_loc", "") if d.get("alt_loc") else "",
     d.get("residue_name", "") if d.get("residue_name") else "",
     d.get("chain_id", "") if d.get("chain_id") else "",
     d.get("residue_id", "") if d.get("residue_id") else "",
     d.get("insert_code", "") if d.get("insert_code") else "",
     d.get("x", "") if d.get("x") else "",
     d.get("y", "") if d.get("y") else "",
     d.get("z", "") if d.get("z") else "",
     d.get("occupancy", 1) if d.get("occupancy") else 1,
     d.get("temperature_factor", "") if d.get("temperature_factor") else "",
     d.get("element", "") if d.get("element") else "",
     d.get("charge", 0) if d.get("charge") else 0
    )
    return PdbRecord(line)


def conections_list_to_records(l):
    """Converts a list of connections to CONECT :py:class:`.PdbRecord` objects.

    :param list l: The list of connections.
    :returns: ``list`` of ``PdbRecord``"""

    records = []
    for connection in l:
        for line_num in range(ceil(len(connection["bond_to"]) / 4)):
            bonded_ids = connection["bond_to"][line_num * 4: (line_num + 1) * 4]
            line = "CONECT" + "{:5}" * (len(bonded_ids) + 1)
            records.append(PdbRecord(line.format(
             connection["atom"], *bonded_ids
            )))
    return records
