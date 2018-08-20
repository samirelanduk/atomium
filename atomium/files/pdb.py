"""Contains functions for dealing with the .pdb file format."""

from copy import deepcopy
from itertools import groupby, chain
from math import ceil
from datetime import datetime
import re
from .data import generate_higher_structures
from .data import DATA_DICT, MODEL_DICT, CHAIN_DICT, RESIDUE_DICT, ATOM_DICT

def pdb_string_to_pdb_dict(filestring):
    """Takes the filecontents of a .pdb file and produces an pdb file
    dictionary from them.

    :param str filestring: The contents of a .pdb file.
    :rtype: ``dict``"""

    pdb_dict = {}
    lines = list(filter(lambda l: bool(l.strip()), filestring.split("\n")))
    model_records = ("ATOM", "HETATM", "ANISOU", "MODEL", "TER", "ENDMDL")
    model = {}
    for line in lines:
        record, value = line[:6].rstrip(), line[6:].rstrip()
        if record == "REMARK":
            if "REMARK" not in pdb_dict: pdb_dict["REMARK"] = {}
            number = value.lstrip().split()[0]
            try:
                pdb_dict["REMARK"][number].append(value[4:])
            except: pdb_dict["REMARK"][number] = [value[4:]]
        elif record not in model_records:
            try:
                pdb_dict[record].append(value)
            except: pdb_dict[record] = [value]
        else:
            if "MODEL" not in pdb_dict: pdb_dict["MODEL"] = []
            if record == "MODEL":
                model = {}
            elif record == "ENDMDL":
                pdb_dict["MODEL"].append(model)
            else:
                try:
                    model[record].append(value)
                except: model[record] = [value]
    if not pdb_dict["MODEL"]: pdb_dict["MODEL"].append(model)
    return pdb_dict


def pdb_dict_to_data_dict(pdb_dict):
    """Takes a basic .pdb dict and turns it into a standard atomium data
    dictionary.

    :param dict pdb_dict: The .pdb dictionary.
    :rtype: ``dict``"""

    d = deepcopy(DATA_DICT)
    update_description_dict(pdb_dict, d)
    update_experiment_dict(pdb_dict, d)
    update_quality_dict(pdb_dict, d)
    update_geometry_dict(pdb_dict, d)
    update_models_list(pdb_dict, d)
    return d


def update_description_dict(pdb_dict, data_dict):
    """Creates the description component of a standard atomium data dictionary
    from a .pdb dictionary.

    :param dict pdb_dict: The .pdb dictionary to read.
    :param dict data_dict: The data dictionary to update."""

    extract_header(pdb_dict, data_dict["description"])
    extract_title(pdb_dict, data_dict["description"])
    extract_keywords(pdb_dict, data_dict["description"])
    extract_authors(pdb_dict, data_dict["description"])


def update_experiment_dict(pdb_dict, data_dict):
    """Creates the experiment component of a standard atomium data dictionary
    from a .pdb dictionary.

    :param dict pdb_dict: The .pdb dictionary to update.
    :param dict data_dict: The data dictionary to update."""

    extract_technique(pdb_dict, data_dict["experiment"])
    extract_source(pdb_dict, data_dict["experiment"])


def update_quality_dict(pdb_dict, data_dict):
    """Creates the quality component of a standard atomium data dictionary
    from a .pdb dictionary.

    :param dict pdb_dict: The .pdb dictionary to update.
    :param dict data_dict: The data dictionary to update."""

    extract_resolution_remark(pdb_dict, data_dict["quality"])
    extract_rvalue_remark(pdb_dict, data_dict["quality"])


def update_geometry_dict(pdb_dict, data_dict):
    """Creates the geometry component of a standard atomium data
    dictionary from a .pdb dictionary.

    :param dict pdb_dict: The .pdb dictionary to update.
    :param dict data_dict: The data dictionary to update."""

    extract_assembly_remark(pdb_dict, data_dict["geometry"])


def update_models_list(pdb_dict, data_dict):
    """Creates the models component of a standard atomium data dictionary
    from a .pdb dictionary.

    :param dict pdb_dict: The .pdb dictionary to update.
    :param dict data_dict: The data dictionary to update."""

    for model_dict in pdb_dict["MODEL"]:
        model = deepcopy(MODEL_DICT)
        for line in model_dict.get("ATOM", []):
            model["atoms"].append(atom_line_to_atom_dict(line))
        for line in model_dict.get("HETATM", []):
            model["atoms"].append(atom_line_to_atom_dict(line, polymer=False))
        assign_anisou(model_dict, model)
        data_dict["models"].append(model)
    generate_higher_structures(data_dict["models"])
    extract_sequence(pdb_dict, data_dict["models"])
    extract_connections(pdb_dict, data_dict["models"])


def extract_header(pdb_dict, description_dict):
    """Takes a ``dict`` and adds header information to it by parsing the HEADER
    line.

    :param dict pdb_dict: the ``dict`` to read.
    :param dict description_dict: the ``dict`` to update."""

    if pdb_dict.get("HEADER"):
        line = pdb_dict["HEADER"][0]
        if line[44:53].strip():
            description_dict["deposition_date"] = datetime.strptime(
             line[44:53], "%d-%b-%y"
            ).date()
        if line[56:60].strip(): description_dict["code"] = line[56:60]
        if line[4:44].strip():
            description_dict["classification"] = line[4:44].strip()


def extract_title(pdb_dict, description_dict):
    """Takes a ``dict`` and adds title information to it by parsing TITLE lines.

    :param dict pdb_dict: the ``dict`` to read.
    :param dict description_dict: the ``dict`` to update."""

    if pdb_dict.get("TITLE"):
        description_dict["title"] = merge_lines(pdb_dict["TITLE"], 4)


def extract_keywords(pdb_dict, description_dict):
    """Takes a ``dict`` and adds keyword information to it by parsing KEYWD
    lines.

    :param dict pdb_dict: the ``dict`` to read.
    :param dict description_dict: the ``dict`` to update."""

    if pdb_dict.get("KEYWDS"):
        text = merge_lines(pdb_dict["KEYWDS"], 4)
        description_dict["keywords"] = [w.strip() for w in text.split(",")]


def extract_authors(pdb_dict, description_dict):
    """Takes a ``dict`` and adds author information to it by parsing AUTHOR
    lines.

    :param dict pdb_dict: the ``dict`` to read.
    :param dict description_dict: the ``dict`` to update."""

    if pdb_dict.get("AUTHOR"):
        text = merge_lines(pdb_dict["AUTHOR"], 4)
        description_dict["authors"] = [w.strip() for w in text.split(",")]


def extract_technique(pdb_dict, experiment_dict):
    """Takes a ``dict`` and adds technique information to it by parsing EXPDTA
    lines.

    :param dict pdb_dict: the ``dict`` to read.
    :param dict experiment_dict: the ``dict`` to update."""

    if pdb_dict.get("EXPDTA"):
        if pdb_dict["EXPDTA"][0].strip():
            experiment_dict["technique"] = pdb_dict["EXPDTA"][0].strip()


def extract_source(pdb_dict, experiment_dict):
    """Takes a ``dict`` and adds source information to it by parsing SOURCE
    lines.

    :param dict pdb_dict: the ``dict`` to read.
    :param dict experiment_dict: the ``dict`` to update."""

    if pdb_dict.get("SOURCE"):
        data = merge_lines(pdb_dict["SOURCE"], 4)
        patterns = {
         "source_organism": r"ORGANISM_SCIENTIFIC\: (.+?);",
         "expression_system": r"EXPRESSION_SYSTEM\: (.+?);"
        }
        for attribute, pattern in patterns.items():
            matches = re.findall(pattern, data)
            if matches:
                experiment_dict[attribute] = matches[0]


def extract_resolution_remark(pdb_dict, quality_dict):
    """Takes a ``dict`` and adds resolution information to it by parsing REMARK
    2 lines.

    :param dict pdb_dict: the ``dict`` to read.
    :param dict quality_dict: the ``dict`` to update."""

    if pdb_dict.get("REMARK") and pdb_dict["REMARK"].get("2"):
        for remark in pdb_dict["REMARK"]["2"]:
            try:
                quality_dict["resolution"] = float(remark.strip().split()[1])
                break
            except: pass


def extract_rvalue_remark(pdb_dict, quality_dict):
    """Takes a ``dict`` and adds resolution information to it by parsing REMARK
    3 lines.

    :param dict pdb_dict: the ``dict`` to read.
    :param dict quality_dict: the ``dict`` to update."""

    if pdb_dict.get("REMARK") and pdb_dict["REMARK"].get("3"):
        patterns = {
         "rvalue": r"R VALUE[ ]{2,}\(WORKING SET\) : (.+)",
         "rfree": r"FREE R VALUE[ ]{2,}: (.+)",
        }
        for attribute, pattern in patterns.items():
            for remark in pdb_dict["REMARK"]["3"]:
                matches = re.findall(pattern, remark.strip())
                if matches:
                    try:
                        quality_dict[attribute] = float(matches[0].strip())
                    except: pass
                    break


def extract_assembly_remark(pdb_dict, geometry_dict):
    """Takes a ``dict`` and adds assembly information to it by parsing REMARK
    350 lines.

    :param dict pdb_dict: the ``dict`` to read.
    :param dict geometry_dict: the ``dict`` to update."""

    if pdb_dict.get("REMARK") and pdb_dict["REMARK"].get("350"):
        groups = [list(g) for k, g in groupby(
         pdb_dict["REMARK"]["350"], lambda x: "ECULE:" in x
        )][1:]
        assemblies = [list(chain(*a)) for a in zip(groups[::2], groups[1::2])]
        for a in assemblies:
            geometry_dict["assemblies"].append(
             assembly_lines_to_assembly_dict(a)
            )


def assembly_lines_to_assembly_dict(lines):
    """Takes the lines representing a single biological assembly and turns
    them into an assembly dictionary.

    :param list lines: The REMARK lines to read.
    :rtype: ``dict``"""

    assembly = {
     "transformations": [], "software": None, "buried_surface_area": None,
     "surface_area": None, "delta_energy": None, "id": 0
    }
    patterns = [[r"(.+)SOFTWARE USED: (.+)", "software", lambda x: x],
     [r"(.+)BIOMOLECULE: (.+)", "id", int],
     [r"(.+)SURFACE AREA: (.+) [A-Z]", "buried_surface_area", float],
     [r"(.+)AREA OF THE COMPLEX: (.+) [A-Z]", "surface_area", float],
     [r"(.+)FREE ENERGY: (.+) [A-Z]", "delta_energy", float]]
    t = None
    for line in lines:
        for p in patterns:
            matches = re.findall(p[0], line)
            if matches: assembly[p[1]] = p[2](matches[0][1].strip())
        if "APPLY THE FOLLOWING" in line:
            if t: assembly["transformations"].append(t)
            t = {"chains": [], "matrix": [], "vector": []}
        if "CHAINS:" in line:
            t["chains"] += [c.strip() for c in
             line.split(":")[-1].strip().split(",") if c.strip()]
        if "BIOMT" in line:
            values = [float(x) for x in line.split()[2:]]
            if len(t["matrix"]) == 3:
                assembly["transformations"].append(t)
                t = {"chains": t["chains"], "matrix": [], "vector": []}
            t["matrix"].append(values[:3])
            t["vector"].append(values[-1])
    if t: assembly["transformations"].append(t)
    return assembly


def atom_line_to_atom_dict(line, polymer=True):
    """Takes an ATOM or HETATM line and converts it to an atom ``dict``.

    :param str line: the atom record to parse.
    :param bool polymer: is this atom in a chain or not?
    :rtype: ``dict``"""

    a = deepcopy(ATOM_DICT)
    if line[:5].strip(): a["id"] = int(line[:5].strip())
    if line[6:10].strip(): a["name"] = line[6:10].strip()
    if line[10].strip(): a["alt_loc"] = line[10]
    if line[11:14].strip(): a["residue_name"] = line[11:14].strip()
    if line[15].strip(): a["chain_id"] = line[15]
    if line[16:20].strip(): a["residue_id"] = int(line[16:20].strip())
    a["residue_insert"] = line[20].strip()
    if line[24:32].strip(): a["x"] = float(line[24:32].strip())
    if line[32:40].strip(): a["y"] = float(line[32:40].strip())
    if line[40:48].strip(): a["z"] = float(line[40:48].strip())
    if line[48:54].strip(): a["occupancy"] = float(line[48:54].strip())
    if line[54:60].strip(): a["bfactor"] = float(line[54:60].strip())
    if line[70:72].strip(): a["element"] = line[70:72].strip()
    if line[72:76].strip():
        try:
            a["charge"] = int(line[72:76].strip())
        except: a["charge"] = int(line[72:76][::-1].strip())
    a["full_res_id"] = "{}{}{}".format(
     a["chain_id"] or "", str(a["residue_id"] or "") or "", a["residue_insert"]
    ) or None
    a["polymer"] = polymer
    return a


def assign_anisou(model_dict, model):
    """Aassigns anisotropy information the atoms in a model dictionary, by
    pasrsing ANISOU lines.

    :param dict pdb_dict: The ``dict`` to read.
    :param list model_list: The list of model dictionaries to update."""

    anisotropy = {int(line[:5].strip()): [
     int(line[n * 7 + 22:n * 7 + 29]) / 10000 for n in range(6)
    ] for line in model_dict.get("ANISOU", [])}
    for atom in model["atoms"]:
        atom["anisotropy"] = anisotropy.get(atom["id"], [0, 0, 0, 0, 0, 0])


def extract_sequence(pdb_dict, model_list):
    """Adds sequence information to the chains of each model in a list of model
    dictionaries, by parsing SEQRES lines.

    :param dict pdb_dict: The ``dict`` to read from.
    :param list model_list: The list of model dictionaries to update."""

    sequences = {}
    if pdb_dict.get("SEQRES"):
        text = merge_lines(pdb_dict["SEQRES"], 4)
        blocks, code = text.split(), None
        blocks = [b for b in blocks if not b.isdigit()]
        for block in blocks:
            if len(block) == 1 and block not in sequences:
                sequences[block], code = [], block
            elif code and not len(block) == 1:
                sequences[code].append(block)
    for model in model_list:
        for chain in model["chains"]:
            chain["full_sequence"] = sequences.get(chain["id"], [])


def extract_connections(pdb_dict, model_list):
    """Adds connectivity information to each model in a list of model
    dictionaries, by parsing CONECT lines.

    :param dict pdb_dict: The ``dict`` to read from.
    :param list model_list: The list of model dictionaries to update."""

    lines = pdb_dict.get("CONECT", [])
    atom_ids = sorted(list(set([int(r[:5].strip()) for r in lines])))
    connections = [{"atom": id_, "bond_to": []} for id_ in atom_ids]
    for connection in connections:
        relevant_lines = [line[5:] for line in lines
          if int(line[:5].strip()) == connection["atom"]]
        line_numbers = [[
         line[n * 5:n * 5 + 5].strip() for n in range(4)
        ] for line in relevant_lines]
        for line in line_numbers:
            for number in line:
                if number:
                    connection["bond_to"].append(int(number))
    for model in model_list: model["connections"] = connections.copy()


def merge_lines(lines, start, join=" "):
    """Gets a single continuous string from a sequence of lines.

    :param list lines: The lines to merge.
    :param int start: The start point in each record.
    :param str join: The string to join on.
    :rtype: ``str``"""

    string = join.join([line[start:].strip() for line in lines])
    return string
