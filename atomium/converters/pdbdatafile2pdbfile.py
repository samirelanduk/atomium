"""Contains the function for creating PdbFiles from PdbDataFiles."""

from ..parse.pdbfile import PdbRecord

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
