"""Contains functions for dealing with the .pdb file format."""

from copy import deepcopy
from itertools import groupby, chain
from math import ceil
from datetime import datetime
import re
from .data import generate_higher_structures
from .data import DATA_DICT, MODEL_DICT, CHAIN_DICT, RESIDUE_DICT, ATOM_DICT
from ..models.data import CODES

def pdb_string_to_pdb_dict(filestring):
    """Takes the filecontents of a .pdb file and produces an pdb file
    dictionary from them.

    :param str filestring: The contents of a .pdb file.
    :rtype: ``dict``"""

    pdb_dict = {}
    lines = list(filter(lambda l: bool(l.strip()), filestring.split("\n")))
    lines = [[line[:6].rstrip(), line[6:].rstrip()] for line in lines]
    pre_m, model, post_m = [], [], []
    model_recs = ("ATOM", "HETATM", "ANISOU", "MODEL", "TER", "ENDMDL")
    for l in lines:
        (model if l[0] in model_recs else post_m if model else pre_m).append(l)
    for line in pre_m:
        if line[0] == "REMARK":
            if "REMARK" not in pdb_dict: pdb_dict["REMARK"] = {}
            number = line[1].lstrip().split()[0]
            try:
                pdb_dict["REMARK"][number].append(line[1][4:])
            except: pdb_dict["REMARK"][number] = [line[1][4:]]
        else:
            try:
                pdb_dict[line[0]].append(line[1])
            except: pdb_dict[line[0]] = [line[1]]
    pdb_dict["MODEL"], m, het = [], {"ATOM": []}, True
    for line in model[::-1]:
        if line[0] == "ENDMDL":
            m, het = {"ATOM": []}, True
        elif line[0] == "MODEL":
            pdb_dict["MODEL"].append(m)
        else:
            if line[0] == "TER": het = False
            if het and line[0] == "ATOM": line[0] = "HETATM"
            if not het and line[0] == "HETATM": line[0] = "ATOM"
            try:
                m[line[0]].append(line[1])
            except: m[line[0]] = [line[1]]
    if not pdb_dict["MODEL"]: pdb_dict["MODEL"].append(m)
    pdb_dict["MODEL"].reverse()
    for m in pdb_dict["MODEL"]:
        if "ATOM" in m: m["ATOM"].reverse()
        if "HETATM" in m: m["HETATM"].reverse()
    for line in post_m:
        try:
            pdb_dict[line[0]].append(line[1])
        except: pdb_dict[line[0]] = [line[1]]
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


def file_to_pdb_string(file):
    """Takes a :py:class:`.File` and turns it into a .pdb filestring that
    represents it.

    :param File file: the File to convert.
    :rtype: ``str``"""

    lines = []
    pack_annotation(file, lines)
    pack_structure(file, lines)
    return "\n".join(lines)


def pack_annotation(file, lines):
    """Adds non-structural records to a list of lines.

    :param File file: The :py:class:`.File` to read.
    :param list lines: The record lines to add to."""

    pack_header(file, lines)
    pack_title(file, lines)
    pack_source(file, lines)
    pack_keywords(file, lines)
    pack_technique(file, lines)
    pack_resolution(file, lines)
    pack_rvalue(file, lines)
    pack_assemblies(file, lines)
    pack_sequences(file, lines)


def pack_structure(file, lines):
    """Adds structural records to a list of lines.

    :param File file: The :py:class:`.File` to read.
    :param list lines: The record lines to add to."""

    for i, model in enumerate(file.models, start=1):
        if len(file.models) > 1: lines.append("MODEL        {}".format(i))
        for atom in sorted(model.atoms(), key=lambda a: a.id):
            atom_to_atom_line(atom, lines)
        if len(file.models) > 1: lines.append("ENDMDL")
    pack_connections(file.model, lines)


def pack_header(file, lines):
    """Adds a HEADER record to a list of lines.

    :param File file: The :py:class:`.File` to read.
    :param list lines: The record lines to add to."""

    if (file.deposition_date or file.code
     or file.classification):
        lines.append("HEADER    {}{}   {}".format(
         file.classification.ljust(40) if
          file.classification else " " * 40,
         file.deposition_date.strftime("%d-%b-%y").upper() if
          file.deposition_date else " " * 9,
         file.code if file.code else "    "
        ).ljust(80))


def pack_title(file, lines):
    """Adds TITLE records to a list of lines.

    :param File file: The :py:class:`.File` to read.
    :param list lines: The record lines to add to."""

    if file.title:
        lines += split_string(file.title, "TITLE", 11)


def pack_resolution(file, lines):
    """Adds REMARK records to a list of lines for resolution.

    :param File file: The :py:class:`.File` to read.
    :param list lines: The record lines to add to."""

    if file.resolution is not None:
        lines.append("REMARK   2".ljust(80))
        lines.append("REMARK   2 RESOLUTION.    {:.2f} ANGSTROMS.".format(
         file.resolution
        ).ljust(80))


def pack_rvalue(file, lines):
    """Adds REMARK records to a list of lines for rvalue.

    :param File file: The :py:class:`.File` to read.
    :param list lines: The record lines to add to."""

    if any([file.rvalue, file.rfree]):
        lines.append("REMARK   3".ljust(80))
    if file.rvalue is not None:
        lines.append("REMARK   3   R VALUE            (WORKING SET) : {}".format(
         file.rvalue
        ).ljust(80))
    if file.rfree is not None:
        lines.append("REMARK   3   FREE R VALUE                     : {}".format(
         file.rfree
        ).ljust(80))


def pack_source(file, lines):
    """Adds SOURCE records to a list of lines for source organism and
    expression system.

    :param File file: The :py:class:`.File` to read.
    :param list lines: The record lines to add to."""

    source = ""
    if file.source_organism:
        source += "ORGANISM_SCIENTIFIC: {};".format(
         file.source_organism
        ).ljust(70)
    if file.expression_system:
        source += "EXPRESSION_SYSTEM: {};".format(
         file.expression_system
        ).ljust(70)
    if source: lines += split_string(source, "SOURCE", 11)


def pack_technique(file, lines):
    """Adds a EXPDTA record to a list of lines for experimental technique.
    :param list lines: The record lines to add to.

    :param File file: The :py:class:`.File` to read.
    :param list lines: The record lines to add to."""

    if file.technique is not None:
        lines.append("EXPDTA    {}".format(file.technique).ljust(80))


def pack_keywords(file, lines):
    """Adds a KEYWDS record to a list of lines for PDB keyword tags.

    :param File file: The :py:class:`.File` to read.
    :param list lines: The record lines to add to."""

    if file.keywords:
        lines += split_string(", ".join(file.keywords), "KEYWDS", 11)


def pack_sequences(file, lines):
    """Adds SEQRES records to a list of lines for sequence annotation.

    :param File file: The :py:class:`.File` to read.
    :param list lines: The record lines to add to."""

    lookup = {v: k for k, v in CODES.items() if k not in ["HIP", "HIE"]}
    try:
        for chain in sorted(file.model.chains(), key=lambda c: c.id):
            if chain.rep_sequence:
                residues = [lookup.get(c, "???") for c in chain.rep_sequence]
                length = len(residues)
                line_count = ceil(length / 13)
                for line_num in range(line_count):
                    lines += ["SEQRES {:>3} {} {:>4}  {}".format(
                     line_num + 1, chain.id, length,
                     " ".join(residues[line_num * 13: (line_num + 1) * 13])
                    ).ljust(80)]
    except AttributeError: pass


def pack_assemblies(file, lines):
    """Adds REMARK 350 records to a list of lines for assembly annotation.

    :param File file: The :py:class:`.File` to read.
    :param list lines: The record lines to add to."""

    for bm in file.assemblies:
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


def atom_to_atom_line(a, lines):
    """Converts an :py:class:`.Atom` to an ATOM or HETATM record and adds it to
    a list of lines. ANISOU records will also be made.

    :param Atom a: The Atom to pack.
    :param list l: The lines to add to."""

    line = "{:6}{:5} {:4} {:3} {:1}{:4}{:1}   "
    line += "{:>8}{:>8}{:>8}  1.00{:6}          {:>2}{:2}"
    id_, residue_name, chain_id, residue_id, insert_code = "", "", "", "", ""
    if a.residue:
        id_, residue_name = a.residue.id, a.residue.name
        chain_id = a.chain.id if a.chain is not None else ""
        residue_id = int("".join([c for c in id_ if c.isdigit()]))
        insert_code = id_[-1] if id_ and id_[-1].isalpha() else ""
    elif a.ligand:
        id_, residue_name = a.ligand.id, a.ligand.name
        chain_id = a.chain.id if a.chain is not None else ""
        residue_id = int("".join([c for c in id_ if c.isdigit()]))
    atom_name = a.name if a.name else ""
    atom_name = " " + atom_name if len(atom_name) < 4 else atom_name
    occupancy = "  1.00"
    line = line.format(
     "HETATM" if a.residue is None else "ATOM",
     a.id, atom_name, residue_name, chain_id, residue_id, insert_code,
     "{:.3f}".format(a.x) if a.x is not None else "",
     "{:.3f}".format(a.y) if a.y is not None else "",
     "{:.3f}".format(a.z) if a.z is not None else "",
     a.bfactor if a.bfactor else "",
     a.element if a.element else "",
     str(int(a.charge))[::-1] if a.charge else "",
    )
    lines.append(line)
    if a.anisotropy != [0, 0, 0, 0, 0, 0]:
        lines.append(atom_to_anisou_line(a, atom_name,
         residue_name, chain_id, residue_id, insert_code))


def atom_to_anisou_line(a, name, res_name, chain_id, res_id, insert):
    """Converts an :py:class:`.Atom` to an ANISOU record.

    :param Atom a: The Atom to pack.
    :param str name: The atom name to use.
    :param str res_name: The residue name to use.
    :param str chain_id: The chain ID to use.
    :param str res_id: The residue ID to use.
    :param str insert: The residue insert code to use.
    :rtype: ``str``"""

    line = "{:6}{:5} {:4} {:3} {:1}{:4}{:1} "
    line += "{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}      {:>2}{:2}"
    anisotropy = [round(x * 10000 )for x in a.anisotropy]
    line = line.format(
     "ANISOU",
     a.id, name, res_name, chain_id, res_id, insert,
     anisotropy[0] if anisotropy[0] is not 0 else "",
     anisotropy[1] if anisotropy[1] is not 0 else "",
     anisotropy[2] if anisotropy[2] is not 0 else "",
     anisotropy[3] if anisotropy[3] is not 0 else "",
     anisotropy[4] if anisotropy[4] is not 0 else "",
     anisotropy[5] if anisotropy[5] is not 0 else "",
     a.element if a.element else "",
     str(int(a.charge))[::-1] if a.charge else "",
    )
    return line


def pack_connections(structure, lines):
    """Adds CONECT records to a list of lines.

    :param AtomStructure structure: The :py:class:`.AtomStructure` to read.
    :param list lines: The record lines to update."""

    for atom in sorted(structure.atoms(), key=lambda a: a.id):
        if not atom.residue:
            bond_ids = sorted([a.id for a in atom.bonded_atoms])
            for line_num in range(ceil(len(bond_ids) / 4)):
                bonded_ids = bond_ids[line_num * 4: (line_num + 1) * 4]
                line = "CONECT" + "{:5}" * (len(bonded_ids) + 1)
                lines.append(line.format(atom.id, *bonded_ids).ljust(80))


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
