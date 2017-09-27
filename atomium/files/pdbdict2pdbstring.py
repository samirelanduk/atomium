"""This module handles the conversion of PDB data dictionaries to .pdb
filestrings."""

from datetime import datetime
from math import ceil

def pdb_dict_to_pdb_string(pdb_dict):
    lines = []
    pack_header(lines, pdb_dict)
    pack_structure(lines, pdb_dict)
    return lines_to_string(lines)


def pack_header(lines, pdb_dict):
    if pdb_dict["deposition_date"] or pdb_dict["code"]:
        lines.append("HEADER{}{}   {}".format(
         " " * 44,
         pdb_dict["deposition_date"].strftime("%d-%b-%y").upper() if
          pdb_dict["deposition_date"] else " " * 9,
         pdb_dict["code"] if pdb_dict["code"] else "    "
        ).ljust(80))
    if pdb_dict["title"]:
        chunks_needed = (len(pdb_dict["title"]) - 1) // 70 + 1
        title_chunks = [
         pdb_dict["title"][i * 70:i * 70 + 70] for i in range(chunks_needed)
        ]
        title_records = ["TITLE    {}{}{}".format(
         number if number > 1 else " ",
         " " if number > 1 else "",
         chunk
        ).ljust(80) for number, chunk in enumerate(title_chunks, start=1)]
        lines += title_records


def pack_structure(lines, pdb_dict):
    if len(pdb_dict["models"]) == 1:
        pack_model(lines, pdb_dict["models"][0], sole=True)
    else:
        for model in pdb_dict["models"]:
            pack_model(lines, model)
    pack_connections(lines, pdb_dict)


def pack_model(lines, model_dict, sole=False):
    if not sole:
        lines.append("MODEL".ljust(80))
    for chain in model_dict["chains"]:
        for residue in chain["residues"]:
            for atom in residue["atoms"]:
                lines.append(atom_dict_to_atom_line(atom, hetero=False))
    for molecule in model_dict["molecules"]:
        for atom in molecule["atoms"]:
            lines.append(atom_dict_to_atom_line(atom, hetero=True))
    if not sole:
        lines.append("ENDMDL".ljust(80))


def atom_dict_to_atom_line(d, hetero=False):
    line = "{:6}{:5} {:4}{:1}{:3} {:1}{:4}{:1}   "
    line += "{:8}{:8}{:8}  {:5}{:6}         {:>2}{:2}"
    atom_name = d["atom_name"] if d["atom_name"] else ""
    atom_name = " " + atom_name if len(atom_name) < 4 else atom_name
    temp_factor = "{:.2f}".format(d["temp_factor"]) if d["temp_factor"] else ""
    line = line.format(
     "HETATM" if hetero else "ATOM",
     d["atom_id"] if d["atom_id"] else "",
     atom_name,
     d["alt_loc"] if d["alt_loc"] else "",
     d["residue_name"] if d["residue_name"] else "",
     d["chain_id"],
     d["residue_id"] if d["residue_id"] else "",
     d["insert_code"],
     d["x"] if d["x"] else "",
     d["y"] if d["y"] else "",
     d["z"] if d["z"] else "",
     "{:.2f}".format(d["occupancy"]) if d["occupancy"] and d["occupancy"] != 1 else "",
     "{:.2f}".format(d["temp_factor"]) if d["temp_factor"] else "",
     d["element"] if d["element"] else "",
     str(d["charge"])[::-1] if d["charge"] else "",
    )
    return line


def pack_connections(lines, pdb_dict):
    for connection in pdb_dict["connections"]:
        for line_num in range(ceil(len(connection["bond_to"]) / 4)):
            bonded_ids = connection["bond_to"][line_num * 4: (line_num + 1) * 4]
            line = "CONECT" + "{:5}" * (len(bonded_ids) + 1)
            lines.append(line.format(connection["atom"], *bonded_ids).ljust(80))


def lines_to_string(lines):
    return "\n".join(lines)
