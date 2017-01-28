"""This module handles the logic of converting a :py:class:`.PdbDataFile` to a
:py:class:`.PdbFile`"""

import math
from ..pdb.pdbfile import PdbFile, PdbRecord
from ..pdb.pdbdatafile import PdbDataFile

def pdb_file_from_pdb_data_file(data_file):
    """Takes a :py:class:`.PdbDataFile`, converts it to a :py:class:`.PdbFile`,
    and returns it.

    :param PdbDataFile data_file: The :py:class:`.PdbDataFile` to convert.
    :rtype: :py:class:`.PdbFile`"""

    if not isinstance(data_file, PdbDataFile):
        raise TypeError(
         "pdb_file_from_pdb_data_file can only convert PdbDataFiles"
        )
    pdb_file = PdbFile()
    pdb_file._source = data_file

    create_compnd_records(pdb_file, data_file)
    create_atom_records(pdb_file, data_file)
    create_atom_records(pdb_file, data_file, hetero=True)
    create_conect_records(pdb_file, data_file)

    return pdb_file


def create_compnd_records(pdb_file, data_file):
    """Takes a :py:class:`.PdbFile` and creates COMPND records in it based on
    the data in the provided :py:class:`.PdbDataFile`

    :param PdbFile pdb_file: the PDB File to update.
    :param PdbDataFile data_file: The source Pdb Data File"""

    lines = []
    for compound in data_file.compounds():
        segments = []
        if "MOL_ID" in compound:
            segments.append("MOL_ID: %s;" % str(compound["MOL_ID"]))
        if "MOLECULE" in compound:
            segments.append("MOLECULE: %s;" % str(compound["MOLECULE"]))
        if "CHAIN" in compound:
            segments.append("CHAIN: %s;" % ", ".join(compound["CHAIN"]))
        if "SYNONYM" in compound:
            segments.append("SYNONYM: %s;" % ", ".join(compound["SYNONYM"]))
        if "EC" in compound:
            segments.append("EC: %s;" % str(compound["EC"]))
        if "ENGINEERED" in compound:
            segments.append(
             "ENGINEERED: %s;" % "YES" if compound["ENGINEERED"] else "NO"
            )
        for segment in segments:
            if len(segment) <= 69:
                lines.append(segment)
            else:
                chunks = segment.split(" ")[::-1]
                growing_line = ""
                while chunks:
                    chunk = chunks.pop()
                    if len(growing_line + chunk) > 69:
                        lines.append(growing_line)
                        growing_line = chunk + " "
                    else:
                        growing_line += chunk + " "
                    if not chunks:
                        lines.append(growing_line)

    for index, line in enumerate(lines):
        pdb_file.add_record(PdbRecord("COMPND %s%s" % (
         str(index + 1).rjust(3) + " " if index != 0 else "   ",
         line
        )))


def create_atom_records(pdb_file, data_file, hetero=False):
    """Takes a :py:class:`.PdbFile` and creates ATOM and HETATM records in it
    based on the data in the provided :py:class:`.PdbDataFile`

    :param PdbFile pdb_file: the PDB File to update.
    :param PdbDataFile data_file: The source Pdb Data File
    :param bool hetero: if True, the function will create HETATM records, and\
    if False, ATOM records will be created. Default is False."""

    atoms = data_file.heteroatoms() if hetero else data_file.atoms()
    for atom in atoms:
        record_fragments = []
        record_fragments.append("HETATM" if hetero else "ATOM  ")
        record_fragments.append("%5i" % atom["atom_id"] + " ")
        record_fragments.append("%-4s" % atom["atom_name"][:4])
        record_fragments.append(atom["alt_loc"][0] if atom["alt_loc"] else " ")
        record_fragments.append(
         atom["residue_name"][0:3] + " " if atom["residue_name"] else "    "
        )
        record_fragments.append(
         atom["chain_id"][0] if atom["chain_id"] else " "
        )
        record_fragments.append(
         ("%4i" % atom["residue_id"]) if atom["residue_id"] else "    "
        )
        record_fragments.append(
         atom["insert_code"][0] + "   " if atom["insert_code"] else "    "
        )
        record_fragments.append(_number_to_8_char_string(atom["x"]))
        record_fragments.append(_number_to_8_char_string(atom["y"]))
        record_fragments.append(_number_to_8_char_string(atom["z"]))
        record_fragments.append(_number_to_6_char_string(atom["occupancy"]))
        record_fragments.append(
         _number_to_6_char_string(atom["temperature_factor"])
        )
        record_fragments.append(" " * 10)
        record_fragments.append("%-2s" % atom["element"])
        record_fragments.append(
         ("%-2i" % atom["charge"]) if atom["charge"] else "  "
        )
        pdb_file.add_record(PdbRecord("".join(record_fragments)))


def create_conect_records(pdb_file, data_file):
    """Takes a :py:class:`.PdbFile` and creates CONECT records in it based on
    the data in the provided :py:class:`.PdbDataFile`

    :param PdbFile pdb_file: the PDB File to update.
    :param PdbDataFile data_file: The source Pdb Data File"""

    for connection in data_file.connections():
        record_count = math.ceil(len(connection["bonded_atoms"]) / 4)
        for n in range(record_count):
            pdb_file.add_record(PdbRecord("CONECT%5i%5s%5s%5s%5s" % (
             connection["atom_id"],
             str(connection["bonded_atoms"][(n * 4)])\
              if (n * 4) < len(connection["bonded_atoms"]) else "",
             str(connection["bonded_atoms"][(n * 4) + 1])
              if (n * 4) + 1 < len(connection["bonded_atoms"]) else "",
             str(connection["bonded_atoms"][(n * 4) + 2])
              if (n * 4) + 2 < len(connection["bonded_atoms"]) else "",
             str(connection["bonded_atoms"][(n * 4) + 3])
              if (n * 4) + 3 < len(connection["bonded_atoms"]) else ""
            )))


def _number_to_8_char_string(number):
    return _number_to_n_char_string(number, 8)


def _number_to_6_char_string(number):
    return _number_to_n_char_string(number, 6)


def _number_to_n_char_string(number, n):
    if number is None:
        return " " * n
    else:
        number = round(number, n)
        int_component = str(int(number))
        if number < 0 and int_component[0] != "-":
            int_component = "-" + int_component
        float_component = number - int(number)
        if len(int_component) >= 6 or float_component == 0 or "e" in str(float_component):
            float_component = ".0"
        else:
            float_component = str(round(float_component, (n - 1) - len(int_component)))
            if float_component[0] == "1":
                float_component = ".0"
                int_component = str(int(int_component) + 1)
            else:
                float_component = "." + float_component.split(".")[-1]
        return (int_component + float_component).ljust(n)
