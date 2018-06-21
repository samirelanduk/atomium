"""This module handles the conversion of PDB data dictionaries to .pdb
filestrings."""

from datetime import datetime
from math import ceil

def pdb_dict_to_pdb_string(pdb_dict):
    """Converts a data ``dict`` to a .pdb filestring.

    :param dict pdb_dict: The data dictionary to pack.
    :rtype: ``str``"""

    from .utilities import lines_to_string
    lines = []
    pack_annotation(lines, pdb_dict)
    pack_structure(lines, pdb_dict)
    return lines_to_string(lines)


def pack_annotation(lines, pdb_dict):
    """Adds non-structural records to a list of lines.

    :param list lines: The record lines to add to.
    :param dict pdb_dict: The data dictionary to pack."""

    pack_header(lines, pdb_dict)
    pack_title(lines, pdb_dict)
    pack_source(lines, pdb_dict)
    pack_keywords(lines, pdb_dict)
    pack_technique(lines, pdb_dict)
    pack_resolution(lines, pdb_dict)
    pack_rfactor(lines, pdb_dict)
    pack_biomolecules(lines, pdb_dict)
    pack_sequences(lines, pdb_dict)


def pack_header(lines, pdb_dict):
    """Adds a HEADER record to a list of lines.

    :param list lines: The record lines to add to.
    :param dict pdb_dict: The data dictionary to pack."""

    if (pdb_dict["deposition_date"] or pdb_dict["code"]
     or pdb_dict["classification"]):
        lines.append("HEADER    {}{}   {}".format(
         pdb_dict["classification"].ljust(40) if
          pdb_dict["classification"] else " " * 40,
         pdb_dict["deposition_date"].strftime("%d-%b-%y").upper() if
          pdb_dict["deposition_date"] else " " * 9,
         pdb_dict["code"] if pdb_dict["code"] else "    "
        ).ljust(80))


def pack_title(lines, pdb_dict):
    """Adds TITLE records to a list of lines.

    :param list lines: The record lines to add to.
    :param dict pdb_dict: The data dictionary to pack."""

    if pdb_dict["title"]:
        lines += split_string(pdb_dict["title"], "TITLE", 11)


def pack_resolution(lines, pdb_dict):
    """Adds REMARK records to a list of lines for resolution.

    :param list lines: The record lines to add to.
    :param dict pdb_dict: The data dictionary to pack."""

    if pdb_dict["resolution"] is not None:
        lines.append("REMARK   2".ljust(80))
        lines.append("REMARK   2 RESOLUTION.    {:.2f} ANGSTROMS.".format(
         pdb_dict["resolution"]
        ).ljust(80))


def pack_rfactor(lines, pdb_dict):
    """Adds REMARK records to a list of lines for rfactor.

    :param list lines: The record lines to add to.
    :param dict pdb_dict: The data dictionary to pack."""

    if any([pdb_dict["rfactor"], pdb_dict["rfree"], pdb_dict["rcount"]]):
        lines.append("REMARK   3".ljust(80))
    if pdb_dict["rfactor"] is not None:
        lines.append("REMARK   3   R VALUE            (WORKING SET) : {}".format(
         pdb_dict["rfactor"]
        ).ljust(80))
    if pdb_dict["rfree"] is not None:
        lines.append("REMARK   3   FREE R VALUE                     : {}".format(
         pdb_dict["rfree"]
        ).ljust(80))
    if pdb_dict["rcount"] is not None:
        lines.append("REMARK   3   FREE R VALUE TEST SET COUNT      : {}".format(
         int(pdb_dict["rcount"])
        ).ljust(80))


def pack_source(lines, pdb_dict):
    """Adds SOURCE records to a list of lines for source organism and
    expression system.

    :param list lines: The record lines to add to.
    :param dict pdb_dict: The data dictionary to pack."""

    source = ""
    if pdb_dict["organism"]:
        source += "ORGANISM_SCIENTIFIC: {};".format(
         pdb_dict["organism"]
        ).ljust(70)
    if pdb_dict["expression_system"]:
        source += "EXPRESSION_SYSTEM: {};".format(
         pdb_dict["expression_system"]
        ).ljust(70)
    if source: lines += split_string(source, "SOURCE", 11)


def pack_technique(lines, pdb_dict):
    """Adds a EXPDTA record to a list of lines for experimental technique.

    :param list lines: The record lines to add to.
    :param dict pdb_dict: The data dictionary to pack."""

    if pdb_dict["technique"] is not None:
        lines.append("EXPDTA    {}".format(pdb_dict["technique"]).ljust(80))


def pack_keywords(lines, pdb_dict):
    """Adds a KEYWDS record to a list of lines for PDB keyword tags.

    :param list lines: The record lines to add to.
    :param dict pdb_dict: The data dictionary to pack."""

    if pdb_dict["keywords"]:
        lines += split_string(", ".join(pdb_dict["keywords"]), "KEYWDS", 11)


def pack_sequences(lines, pdb_dict):
    """Adds SEQRES records to a list of lines for sequence annotation.

    :param list lines: The record lines to add to.
    :param dict pdb_dict: The data dictionary to pack."""

    if pdb_dict["sequences"]:
        for chain in sorted(pdb_dict["sequences"].keys()):
            residues = pdb_dict["sequences"][chain]
            length = len(residues)
            line_count = ceil(length / 13)
            for line_num in range(line_count):
                lines += ["SEQRES {:>3} {} {:>4}  {}".format(
                 line_num + 1, chain, length,
                 " ".join(residues[line_num * 13: (line_num + 1) * 13])
                ).ljust(80)]


def pack_biomolecules(lines, pdb_dict):
    """Adds REMARK 350 records to a list of lines for assembly annotation.

    :param list lines: The record lines to add to.
    :param dict pdb_dict: The data dictionary to pack."""

    for bm in pdb_dict["biomolecules"]:
        lines.append("REMARK 350")
        lines.append("REMARK 350 BIOMOLECULE: {}".format(bm["id"]))
        if bm["software"]:
            lines.append("REMARK 350 SOFTWARE USED: {}".format(bm["software"]))
        if bm["buried_surface_area"]:
            lines.append("REMARK 350 TOTAL BURIED SURFACE AREA:" +
             " {} ANGSTROM**2".format(int(bm["buried_surface_area"])))
        if bm["surface_area"]:
            lines.append("REMARK 350 SURFACE AREA OF THE COMPLEX:" +
             " {} ANGSTROM**2".format(int(bm["surface_area"])))
        if bm["delta_energy"]:
            lines.append("REMARK 350 CHANGE IN SOLVENT FREE ENERGY:" +
             " {} KCAL/MOL".format(bm["delta_energy"]))
        current_chains = ""
        for index, transformation in enumerate(bm["transformations"]):
            if transformation["chains"] != current_chains:
                lines.append("REMARK 350 APPLY THE FOLLOWING TO CHAINS: "
                 + ", ".join(transformation["chains"]))
                current_chains = transformation["chains"]
            for line in (1, 2, 3):
                matrix_line = transformation["matrix"][line - 1]
                lines.append("REMARK 350   BIOMT" +
                 "{}   {} {}{:.6f} {}{:.6f} {}{:.6f}       {}{:.5f}".format(
                  line, index + 1,
                  " " if matrix_line[0] >= 0 else "", matrix_line[0],
                  " " if matrix_line[1] >= 0 else "", matrix_line[1],
                  " " if matrix_line[2] >= 0 else "", matrix_line[2],
                  " " if transformation["vector"][line - 1] >= 0 else "",
                  transformation["vector"][line - 1],
                 )
                )


def pack_structure(lines, pdb_dict):
    """Adds structure records to a list of lines.

    :param list lines: The record lines to add to.
    :param dict pdb_dict: The data dictionary to pack."""

    if len(pdb_dict["models"]) == 1:
        pack_model(lines, pdb_dict["models"][0], multi=0)
    else:
        for index, model in enumerate(pdb_dict["models"], start=1):
            pack_model(lines, model, multi=index)
    pack_connections(lines, pdb_dict)


def pack_model(lines, model_list, multi=0):
    """Adds structure records for a single model to a list of lines.

    :param list lines: The record lines to add to.
    :param dict model_list: The model list to pack.
    :param int multi: The model number, if applicable"""

    if multi > 0:
        lines.append("MODEL        {}".format(multi).ljust(80))
    atom_dicts, heteroatom_dicts = [], []
    for chain in model_list:
        for residue in chain["residues"]:
            atom_dicts += residue["atoms"]
        for ligand in chain["ligands"]:
            heteroatom_dicts += ligand["atoms"]
    for atom in sorted(atom_dicts, key=lambda a: a["atom_id"]):
        lines.append(atom_dict_to_atom_line(atom, hetero=False))
        if atom["anisotropy"] != [0] * 6:
            lines.append(atom_dict_to_anisou_line(atom))
    for atom in sorted(heteroatom_dicts, key=lambda a: a["atom_id"]):
        lines.append(atom_dict_to_atom_line(atom, hetero=True))
        if atom["anisotropy"] != [0] * 6:
            lines.append(atom_dict_to_anisou_line(atom))
    if multi > 0:
        lines.append("ENDMDL".ljust(80))


def atom_dict_to_atom_line(d, hetero=False):
    """Converts an atom ``dict`` to an ATOM or HETATM record.

    :param dict d: The atom dictionary to pack.
    :param bool hetero: if ``True`` a HETATM record will be made, not ATOM.
    :rtype: ``str``"""

    line = "{:6}{:5} {:4}{:1}{:3} {:1}{:4}{:1}   "
    line += "{:>8}{:>8}{:>8}{:6}{:6}          {:>2}{:2}"
    atom_name = d["atom_name"] if d["atom_name"] else ""
    atom_name = " " + atom_name if len(atom_name) < 4 else atom_name
    occupancy = "  1.00" if (
     d["occupancy"] == 1
    ) else "{:.2f}".format(d["occupancy"]).rjust(6)
    line = line.format(
     "HETATM" if hetero else "ATOM",
     d["atom_id"],
     atom_name,
     d["alt_loc"] if d["alt_loc"] else "",
     d["residue_name"] if d["residue_name"] else "",
     d["chain_id"],
     d["residue_id"] if d["residue_id"] else "",
     d["insert_code"],
     "{:.3f}".format(d["x"]) if d["x"] is not None else "",
     "{:.3f}".format(d["y"]) if d["y"] is not None else "",
     "{:.3f}".format(d["z"]) if d["z"] is not None else "",
     occupancy,
     d["temp_factor"] if d["temp_factor"] is not None else "",
     d["element"] if d["element"] else "",
     str(d["charge"])[::-1] if d["charge"] else "",
    )
    return line


def atom_dict_to_anisou_line(d):
    """Converts an atom ``dict`` to an ANISOU record.

    :param dict d: The atom dictionary to pack."""

    line = "{:6}{:5} {:4}{:1}{:3} {:1}{:4}{:1} "
    line += "{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}      {:>2}{:2}"
    atom_name = d["atom_name"] if d["atom_name"] else ""
    atom_name = " " + atom_name if len(atom_name) < 4 else atom_name
    anisotropy = [round(x * 10000 )for x in d["anisotropy"]]
    line = line.format(
     "ANISOU",
     d["atom_id"],
     atom_name,
     d["alt_loc"] if d["alt_loc"] else "",
     d["residue_name"] if d["residue_name"] else "",
     d["chain_id"],
     d["residue_id"] if d["residue_id"] else "",
     d["insert_code"],
     anisotropy[0] if anisotropy[0] is not 0 else "",
     anisotropy[1] if anisotropy[1] is not 0 else "",
     anisotropy[2] if anisotropy[2] is not 0 else "",
     anisotropy[3] if anisotropy[3] is not 0 else "",
     anisotropy[4] if anisotropy[4] is not 0 else "",
     anisotropy[5] if anisotropy[5] is not 0 else "",
     d["element"] if d["element"] else "",
     str(d["charge"])[::-1] if d["charge"] else "",
    )
    return line


def pack_connections(lines, pdb_dict):
    """Adds CONECT records to a list of lines.

    :param list lines: The record lines to add to.
    :param dict pdb_dict: The data dictionary to pack."""

    for connection in pdb_dict["connections"]:
        for line_num in range(ceil(len(connection["bond_to"]) / 4)):
            bonded_ids = connection["bond_to"][line_num * 4: (line_num + 1) * 4]
            line = "CONECT" + "{:5}" * (len(bonded_ids) + 1)
            lines.append(line.format(connection["atom"], *bonded_ids).ljust(80))



def split_string(string, record, start):
    """Takes a string and splits it into multple PDB records, with number
    continuation lines.

    :param str string: The string to split.
    :param str record: The record name.
    :param int start: The position to start the string at on each line."""

    space = 81 - start
    if len(string) <= space:
        return ["{}{}".format(record.ljust(start - 1), string).ljust(80)]
    else:
        words = string.split(" ")
        lines, line = [], record.ljust(start - 1)
        while words:
            if len(line) + len(words[0]) > 79:
                lines.append(line[:80].ljust(80))
                line = "{}{:<2}{} ".format(
                 record.ljust(start - 2), len(lines) + 1, words.pop(0).lstrip()
                )
            else:
                line += words.pop(0).lstrip() + " "
    if len(line.rstrip()) > 10: lines.append(line[:80].ljust(80))
    return lines
