"""Contains functions for dealing with the .mmtf file format."""

from datetime import datetime
import re
from itertools import groupby, chain
from .data import CODES

def pdb_string_to_pdb_dict(filestring):
    """Takes a .pdb filestring and turns into a ``dict`` which represents its
    record structure. Only lines which aren't empty are used.

    The resultant dictionary has line types as the keys, which point to the
    lines as its value. So ``{"TITLE": ["TITLE line 1", "TITLE line 2"]}`` etc.

    The exceptions are the REMARK records, where there is a sub-dictionary with
    REMARK numbers as keys, and the structure records themselves which are just
    arranged into lists - one for each model.

    :param str filestring: the .pdb filestring to process.
    :rtype: ``dict``"""

    pdb_dict = {}
    lines = list(filter(lambda l: bool(l.strip()), filestring.split("\n")))
    lines = [[line[:6].rstrip(), line.rstrip()] for line in lines]
    model_recs = ("ATOM", "HETATM", "ANISOU", "MODEL", "TER", "ENDMDL")
    model = []
    in_model = False
    for head, line in lines:
        if head == "REMARK":
            if "REMARK" not in pdb_dict: pdb_dict["REMARK"] = {}
            number = line.lstrip().split()[1]
            update_dict(pdb_dict["REMARK"], number, line)
        elif head in model_recs:
            if "MODEL" not in pdb_dict: pdb_dict["MODEL"] = [[]]
            if head == "ENDMDL":
                pdb_dict["MODEL"].append([])
            elif head != "MODEL":
                pdb_dict["MODEL"][-1].append(line)
        else:
            update_dict(pdb_dict, head, line)
    if "MODEL" in pdb_dict and not pdb_dict["MODEL"][-1]: pdb_dict["MODEL"].pop()
    return pdb_dict


def update_dict(d, key, value):
    """Takes a dictionary where the values are lists, and adds a value to one of
    the lists at the specific key. If the list doesn't exist, it creates it
    first.

    The dictionary is changed in place.

    :param dict d: the dictionary to update.
    :param str key: the location of the list.
    :param str value: the value to add to the list."""

    try:
        d[key].append(value)
    except: d[key] = [value]


def pdb_dict_to_data_dict(pdb_dict):
    """Converts an .pdb dictionary into an atomium data dictionary, with the
    same standard layout that the other file formats get converted into.

    :param dict pdb_dict: the .pdb dictionary.
    :rtype: ``dict``"""

    data_dict = {
     "description": {
      "code": None, "title": None, "deposition_date": None,
      "classification": None, "keywords": [], "authors": []
     }, "experiment": {
      "technique": None, "source_organism": None, "expression_system": None
     }, "quality": {"resolution": None, "rvalue": None, "rfree": None},
      "geometry": {"assemblies": []}, "models": []
    }
    update_description_dict(pdb_dict, data_dict)
    update_experiment_dict(pdb_dict, data_dict)
    update_quality_dict(pdb_dict, data_dict)
    update_geometry_dict(pdb_dict, data_dict)
    update_models_list(pdb_dict, data_dict)
    return data_dict


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
    """Creates the geometry component of a standard atomium data dictionary
    from a .pdb dictionary.

    :param dict pdb_dict: The .pdb dictionary to update.
    :param dict data_dict: The data dictionary to update."""

    extract_assembly_remark(pdb_dict, data_dict["geometry"])


def update_models_list(pdb_dict, data_dict):
    sequences = make_sequences(pdb_dict)
    for model_lines in pdb_dict["MODEL"]:
        aniso = make_aniso(model_lines)
        last_ter = get_last_ter_line(model_lines)
        model = {"polymer": {}, "non-polymer": {}, "water": {}}
        for index, line in enumerate(model_lines):
            if line[:6] in ["ATOM  ", "HETATM"]:
                chain_id = line[21] if index < last_ter else id_from_line(line)
                res_id = id_from_line(line)
                if index < last_ter:
                    add_atom_to_polymer(line, model, chain_id, res_id, aniso)
                else:
                    add_atom_to_non_polymer(line, model, res_id, aniso)
            for chain_id, chain in model["polymer"].items():
                chain["sequence"] = sequences.get(chain_id, "")
        data_dict["models"].append(model)


def extract_header(pdb_dict, description_dict):
    """Takes a ``dict`` and adds header information to it by parsing the HEADER
    line.

    :param dict pdb_dict: the ``dict`` to read.
    :param dict description_dict: the ``dict`` to update."""

    if pdb_dict.get("HEADER"):
        line = pdb_dict["HEADER"][0]
        if line[50:59].strip():
            description_dict["deposition_date"] = datetime.strptime(
             line[50:59], "%d-%b-%y"
            ).date()
        if line[62:66].strip(): description_dict["code"] = line[62:66]
        if line[10:50].strip():
            description_dict["classification"] = line[10:50].strip()


def extract_title(pdb_dict, description_dict):
    """Takes a ``dict`` and adds header information to it by parsing the TITLE
    lines.

    :param dict pdb_dict: the ``dict`` to read.
    :param dict description_dict: the ``dict`` to update."""

    if pdb_dict.get("TITLE"):
        description_dict["title"] = merge_lines(pdb_dict["TITLE"], 10)


def extract_keywords(pdb_dict, description_dict):
    """Takes a ``dict`` and adds header information to it by parsing the KEYWDS
    line.

    :param dict pdb_dict: the ``dict`` to read.
    :param dict description_dict: the ``dict`` to update."""

    if pdb_dict.get("KEYWDS"):
        text = merge_lines(pdb_dict["KEYWDS"], 10)
        description_dict["keywords"] = [w.strip() for w in text.split(",")]


def extract_authors(pdb_dict, description_dict):
    """Takes a ``dict`` and adds header information to it by parsing the AUTHOR
    line.

    :param dict pdb_dict: the ``dict`` to read.
    :param dict description_dict: the ``dict`` to update."""

    if pdb_dict.get("AUTHOR"):
        text = merge_lines(pdb_dict["AUTHOR"], 10)
        description_dict["authors"] = [w.strip() for w in text.split(",")]


def extract_technique(pdb_dict, experiment_dict):
    """Takes a ``dict`` and adds technique information to it by parsing EXPDTA
    lines.

    :param dict pdb_dict: the ``dict`` to read.
    :param dict experiment_dict: the ``dict`` to update."""

    if pdb_dict.get("EXPDTA"):
        if pdb_dict["EXPDTA"][0].strip():
            experiment_dict["technique"] = pdb_dict["EXPDTA"][0][6:].strip()


def extract_source(pdb_dict, experiment_dict):
    """Takes a ``dict`` and adds source information to it by parsing SOURCE
    lines.

    :param dict pdb_dict: the ``dict`` to read.
    :param dict experiment_dict: the ``dict`` to update."""

    if pdb_dict.get("SOURCE"):
        data = merge_lines(pdb_dict["SOURCE"], 10)
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
                quality_dict["resolution"] = float(remark[10:].strip().split()[1])
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
            values = [float(x) for x in line.split()[4:]]
            if len(t["matrix"]) == 3:
                assembly["transformations"].append(t)
                t = {"chains": t["chains"], "matrix": [], "vector": []}
            t["matrix"].append(values[:3])
            t["vector"].append(values[-1])
    if t: assembly["transformations"].append(t)
    return assembly


def make_sequences(pdb_dict):
    """Creates a mapping of chain IDs to sequences, by parsing SEQRES records.

    :param dict pdb_dict: the .pdb dictionary to read.
    :rtype: ``dict``"""

    seq = {}
    if pdb_dict.get("SEQRES"):
        for line in pdb_dict["SEQRES"]:
            chain, residues = line[11], line[19:].strip().split()
            if chain not in seq:
                seq[chain] = []
            seq[chain] += residues
    return {k: "".join([CODES.get(r, "X") for r in v]) for k, v in seq.items()}


def make_aniso(model_lines):
    """Creates a mapping of chain IDs to anisotropy, by parsing ANISOU records.

    :param dict pdb_dict: the .pdb dictionary to read.
    :rtype: ``dict``"""

    return {int(line[6:11].strip()): [
     int(line[n * 7 + 28:n * 7 + 35]) / 10000 for n in range(6)
    ] for line in model_lines if line[:6] == "ANISOU"}


def get_last_ter_line(model_lines):
    """Gets the index of the last TER record in a list of records. 0 will be
    returned if there are none.

    :param list model_lines: the lines to search.
    :rtype: ``int``"""

    last_ter = 0
    for index, line in enumerate(model_lines[::-1]):
        if line[:3] == "TER":
            last_ter = len(model_lines) - index - 1
            break
    return last_ter


def id_from_line(line):
    """Creates a residue ID from an atom line.

    :param str line: the ATOM or HETATM line record.
    :rtype: ``str``"""

    return "{}.{}{}".format(line[21], line[22:26].strip(), line[26].strip())


def add_atom_to_polymer(line, model, chain_id, res_id, aniso_dict):
    """Takes an .pdb ATOM or HETATM record, converts it, and adds it to a
    polymer dictionary.

    :param dict line: the line to read.
    :param dict model: the model to update.
    :param str chain_id: the chain ID to add to.
    :param str res_id: the molecule ID to add to.
    :param dict aniso_dict: lookup dictionary for anisotropy information."""

    try:
        model["polymer"][chain_id]["residues"][res_id]["atoms"][
         int(line[6:11])
        ] = atom_line_to_dict(line, aniso_dict)
    except:
        try:
            model["polymer"][chain_id]["residues"][res_id] = {
             "name": line[17:20].strip(),
             "atoms": {int(line[6:11]): atom_line_to_dict(line, aniso_dict)},
             "number": len(model["polymer"][chain_id]["residues"]) + 1
            }
        except:
            model["polymer"][chain_id] = {
             "internal_id": chain_id,
             "residues": {res_id: {
              "name": line[17:20].strip(),
              "atoms": {int(line[6:11]): atom_line_to_dict(line, aniso_dict)},
              "number": 1
             }}
            }


def add_atom_to_non_polymer(line, model, res_id, aniso_dict):
    """Takes an .pdb ATOM or HETATM record, converts it, and adds it to a
    non-polymer dictionary.

    :param dict line: the line to read.
    :param dict model: the model to update.
    :param str res_id: the molecule ID to add to.
    :param dict aniso_dict: lookup dictionary for anisotropy information."""

    key = "water" if line[17:20] == "HOH" else "non-polymer"
    try:
        model[key][res_id]["atoms"][
         int(line[6:11])
        ] = atom_line_to_dict(line, aniso_dict)
    except:
        model[key][res_id] = {
         "name": line[17:20].strip(),
         "internal_id": line[21], "polymer": line[21],
         "atoms": {int(line[6:11]): atom_line_to_dict(line, aniso_dict)}
        }


def atom_line_to_dict(line, aniso_dict):
    """Converts an ATOM or HETATM record to an atom dictionary.

    :param str line: the record to convert.
    :param dict aniso_dict: the anisotropy dictionary to use.
    :rtype: ``dict``"""

    a = {
     "occupancy": 1, "bvalue": None, "charge": 0,
     "anisotropy": aniso_dict.get(int(line[6:11].strip()), [0, 0, 0, 0, 0, 0])
    }
    a["name"] = line[12:16].strip() or None
    a["alt_loc"] = line[16].strip() or None
    a["x"] = float(line[30:38].strip())
    a["y"] = float(line[38:46].strip())
    a["z"] = float(line[46:54].strip())
    if line[54:60].strip(): a["occupancy"] = float(line[54:60].strip())
    if line[60:66].strip(): a["bvalue"] = float(line[60:66].strip())
    a["element"] = line[76:78].strip() or None
    if line[78:80].strip():
        try:
            a["charge"] = int(line[78:80].strip())
        except: a["charge"] = int(line[78:80][::-1].strip())
    return a


def merge_lines(lines, start, join=" "):
    """Gets a single continuous string from a sequence of lines.

    :param list lines: The lines to merge.
    :param int start: The start point in each record.
    :param str join: The string to join on.
    :rtype: ``str``"""

    string = join.join([line[start:].strip() for line in lines])
    return string
