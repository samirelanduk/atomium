"""This module handles the conversion of PDB data dictionaries to .pdb
filestrings."""

from datetime import datetime

def pdb_dict_to_pdb_string(pdb_dict):
    lines = []
    pack_header(lines, pdb_dict)
    pack_structure(lines, pdb_dict)
    return pdb_lines_to_string(lines)


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
    pass


def pdb_lines_to_string(lines):
    pass
