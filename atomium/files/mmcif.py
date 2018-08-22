"""Contains functions for dealing with the .mmcif file format."""

import shlex
import re
from copy import deepcopy
from datetime import datetime
from .data import generate_higher_structures, DATA_DICT, MODEL_DICT, ATOM_DICT

def mmcif_string_to_mmcif_dict(filestring):
    """Takes the filecontents of a .cif file and produces an atomium data
    dictionary from them.

    :param str filestring: The contents of a .cif file.
    :rtype: ``dict``"""

    lines = list(filter(lambda l: l and l[0] != "#", filestring.split("\n")))
    lines = consolidate_strings(lines)
    blocks = mmcif_lines_to_mmcif_blocks(lines)
    mmcif_dict = {}
    for block in blocks:
        if block["lines"][0] == "loop_":
            mmcif_dict[block["category"]] = loop_block_to_list(block)
        else:
            mmcif_dict[block["category"]] = category_block_to_dict(block)
    strip_quotes(mmcif_dict)
    return mmcif_dict


def consolidate_strings(lines):
    """In the .mmcif format, long strings are placed on a line of their own and
    marked with semi-colons. This function goes through the lines of such a file
    and creates new lines where these strings are on the preceding line where
    they semantically belong.

    :param list lines: The lines of a .mmcif file.
    :rtype: ``list``"""

    new_lines = []
    while lines:
        line = lines.pop(0)
        if line.startswith(";"):
            string = [line[1:].strip()]
            while not lines[0].startswith(";"):
                string.append(lines.pop(0))
            lines.pop(0)
            new_lines[-1] += " \"{}\"".format(" ".join(string))
        else:
            new_lines.append(line)
    return new_lines


def mmcif_lines_to_mmcif_blocks(lines):
    """Takes a list if .mmcif lines and groups them into blocks based on
    category name.

    :param list lines: The lines of a .mmcif file.
    :rtype: ``list``"""

    category = None
    block, blocks = [], []
    while lines:
        line = lines.pop(0)
        if line.startswith("data_"): continue
        if line.startswith("_"):
            line_category = line.split(".")[0]
            if line_category != category:
                if category:
                    blocks.append({"category": category[1:], "lines": block})
                category = line_category
                block = []
        if line.startswith("loop_"):
            if category:
                blocks.append({"category": category[1:], "lines": block})
            category = lines[0].split(".")[0]
            block = []
        block.append(line)
    if block: blocks.append({"category": category[1:], "lines": block})
    return blocks


def category_block_to_dict(block):
    """Takes a category block from a .mmcif file and turns it into a dictionary.

    :param dict block: The block to convert.
    :rtype: ``dict``"""

    d = {}
    category = block["lines"][0].split(".")[0]
    for index in range(len(block["lines"]) - 1):
        if block["lines"][index + 1][0] != "_":
            block["lines"][index] += " " + block["lines"][index + 1]
    block["lines"] = [l for l in block["lines"] if l[0] == "_"]
    for line in block["lines"]:
        name = line.split(".")[1].split()[0]
        value = line
        if line.startswith("_"):
            value = " ".join(line.split()[1:])
        if value[0] == "'" and value[-1] == "'": value = value[1:-1]
        d[name] = value
    return [d]


def loop_block_to_list(block):
    """Takes a loop block from a .mmcif file and turns it into a table list.

    :param dict block: The block to convert.
    :rtype: ``list``"""

    names = [l for l in block["lines"] if l.startswith("_" + block["category"])]
    lines = [l for l in block["lines"][1:] if l not in names]
    names = [l.split(".")[1].split()[0] for l in names]
    lines = [split_values(l) for l in lines]
    l = []
    for n in range(len(lines) - 1):
        while n < len(lines) - 1 and\
         len(lines[n]) + len(lines[n + 1]) <= len(names):
            lines[n] += lines[n + 1]
            lines.pop(n + 1)
    for line in lines:
        l.append({
         name: value for name, value in zip(names, line)
        })
    return l


def split_values(line):
    """Takes a line and splits it on strings, ignoring those in strings. Like
    ``shlex.split()`` only faster and it handles quotes inside strings.

    :param str line: The line to split.
    :rtype: ``list``"""
    if not re.search("[\'\"]", line): return line.split()
    chars = list(line.strip())
    values, value, in_string = [], [], False
    while chars:
        char = chars.pop(0)
        if char == " " and not in_string:
            values.append(value)
            value = []
        elif char in "'\"":
            if in_string and chars and chars[0] != " ":
                value.append(char)
            else:
                in_string = not in_string
        else:
            value.append(char)
    values.append(value)
    return ["".join(v) for v in values if v]


def strip_quotes(mmcif_dict):
    """Takes a finished mmcif dictionary, and trims any strings that have now
    unneeded quote marks around them.

    :param dict mmcif_dict: The ``dict`` to clean up."""

    for key, value in mmcif_dict.items():
        dicts = [value] if isinstance(value, dict) else value
        for d in dicts:
            for key2, value2 in d.items():
                if value2[0] == '"' and value2[-1] == '"':
                    d[key2] = value2[1:-1]


def mmcif_dict_to_data_dict(mmcif_dict):
    """Takes a basic .mmcif dict and turns it into a standard atomium data
    dictionary.

    :param dict mmcif_dict: The .mmcif dictionary.
    :rtype: ``dict``"""

    d = deepcopy(DATA_DICT)
    update_description_dict(mmcif_dict, d)
    update_experiment_dict(mmcif_dict, d)
    update_quality_dict(mmcif_dict, d)
    update_geometry_dict(mmcif_dict, d)
    update_models_list(mmcif_dict, d)
    return d


def update_description_dict(mmcif_dict, data_dict):
    """Creates the description component of a standard atomium data dictionary
    from a .mmcif dictionary.

    :param dict mmcif_dict: The .mmcif dictionary.
    :param dict data_dict: The data dictionary to update."""

    try:
        data_dict["description"]["code"] = mmcif_dict["entry"][0]["id"]
    except: pass
    try:
        data_dict["description"]["title"] = mmcif_dict["struct"][0]["title"]
    except: pass
    try:
        data_dict["description"]["deposition_date"] = datetime.strptime(
         mmcif_dict["pdbx_database_status"][0]["recvd_initial_deposition_date"],
         "%Y-%m-%d"
        ).date()
    except: pass
    try:
        data_dict["description"]["classification"] =\
         mmcif_dict["struct_keywords"][0]["pdbx_keywords"]
    except: pass
    try:
        data_dict["description"]["keywords"] = [
         w.strip() for w in mmcif_dict["struct_keywords"][0]["text"].split(",")
        ]
    except: pass
    try:
        data_dict["description"]["authors"] = [
         a["name"] for a in mmcif_dict["audit_author"]
        ]
    except: pass


def update_experiment_dict(mmcif_dict, data_dict):
    """Creates the experiment component of a standard atomium data dictionary
    from a .mmcif dictionary.

    :param dict mmcif_dict: The .mmcif dictionary.
    :param dict data_dict: The data dictionary to update."""

    try:
        data_dict["experiment"]["technique"] = mmcif_dict["exptl"][0]["method"]
    except: pass
    try:
        data_dict["experiment"]["source_organism"] =\
         mmcif_dict["entity_src_gen"][0]["pdbx_gene_src_scientific_name"]
    except: pass
    try:
        data_dict["experiment"]["expression_system"] =\
         mmcif_dict["entity_src_gen"][0]["pdbx_host_org_scientific_name"]
    except: pass


def update_quality_dict(mmcif_dict, data_dict):
    """Creates the quality component of a standard atomium data dictionary
    from a .mmcif dictionary.

    :param dict mmcif_dict: The .mmcif dictionary.
    :param dict data_dict: The data dictionary to update."""

    try:
        data_dict["quality"]["resolution"] = float(
         mmcif_dict["reflns"][0]["d_resolution_high"]
        )
    except: pass
    try:
        data_dict["quality"]["rvalue"] = float(
         mmcif_dict["refine"][0]["ls_R_factor_R_work"]
        )
    except: pass
    try:
        data_dict["quality"]["rfree"] = float(
         mmcif_dict["refine"][0]["ls_R_factor_R_free"]
        )
    except: pass


def update_geometry_dict(mmcif_dict, data_dict):
    """Creates the geometry component of a standard atomium data dictionary
    from a .mmcif dictionary.

    :param dict mmcif_dict: The .mmcif dictionary.
    :param dict data_dict: The data dictionary to update."""

    assembly_ids = [a["id"] for a in mmcif_dict.get("pdbx_struct_assembly", [])]
    for a_id in assembly_ids:
        assembly = {
         "id": int(a_id), "surface_area": None, "buried_surface_area": None,
         "software": None, "delta_energy": None, "transformations": []
        }
        assign_software_to_assembly(assembly, mmcif_dict)
        assign_metrics_to_assembly(assembly, mmcif_dict)
        assign_transformations_to_assembly(assembly, mmcif_dict)
        data_dict["geometry"]["assemblies"].append(assembly)


def update_models_list(mmcif_dict, data_dict):
    """Creates the models component of a standard atomium data dictionary
    from a .mmcif dictionary.

    :param dict mmcif_dict: The .mmcif dictionary.
    :param dict data_dict: The data dictionary to update."""

    model = deepcopy(MODEL_DICT)
    model_num = int(mmcif_dict["atom_site"][0]["pdbx_PDB_model_num"])
    for line in mmcif_dict["atom_site"]:
        if model["atoms"] and int(line["pdbx_PDB_model_num"]) > model_num:
            model_num = int(line["pdbx_PDB_model_num"])
            if model["atoms"]:
                data_dict["models"].append(model)
                model = deepcopy(MODEL_DICT)
        model["atoms"].append(atom_line_to_atom_dict(line))
    if model["atoms"]: data_dict["models"].append(model)
    generate_higher_structures(data_dict["models"])
    for model in data_dict["models"]:
        for res in mmcif_dict.get("pdbx_poly_seq_scheme", []):
            for chain in model["chains"]:
                if res["asym_id"] == chain["id"]:
                    chain["full_sequence"].append(res["mon_id"])
        for atom in model["atoms"]:
            for anisou in mmcif_dict.get("atom_site_anisotrop", []):
                if float(anisou["id"]) == atom["id"]:
                    for x in ("11", "12", "13", "22", "23", "33"):
                        atom["anisotropy"].append(
                         float(anisou["U[{}][{}]".format(*list(x))]) / 10000
                        )
                    break
            else: atom["anisotropy"] = [0, 0, 0, 0, 0, 0]


def assign_software_to_assembly(assembly, mmcif_dict):
    """Takes an assembly dict, and goes through an mmcif dictionary looking for
    relevant software information to update it with.

    :param dict assembly: The assembly to update.
    :param dict mmcif_dict: The dictionary to read."""

    for a in mmcif_dict.get("pdbx_struct_assembly", []):
        if a["id"] == str(assembly["id"]):
            assembly["software"] = None if a["method_details"] == "?" else a["method_details"]
            break


def assign_metrics_to_assembly(assembly, mmcif_dict):
    """Takes an assembly dict, and goes through an mmcif dictionary looking for
    relevant energy etc. information to update it with.

    :param dict assembly: The assembly to update.
    :param dict mmcif_dict: The dictionary to read."""

    for a in mmcif_dict.get("pdbx_struct_assembly_prop", []):
        if a["biol_id"] == str(assembly["id"]):
            if a["type"] == "MORE":
                assembly["delta_energy"] = float(a["value"])
            elif a["type"] == "SSA (A^2)":
                assembly["surface_area"] = float(a["value"])
            elif a["type"] == "ABSA (A^2)":
                assembly["buried_surface_area"] = float(a["value"])


def assign_transformations_to_assembly(assembly, mmcif_dict):
    """Takes an assembly dict, and goes through an mmcif dictionary looking for
    relevant transformation information to update it with.

    :param dict assembly: The assembly to update.
    :param dict mmcif_dict: The dictionary to read."""

    for a in mmcif_dict.get("pdbx_struct_assembly_gen", []):
        if a["assembly_id"] == str(assembly["id"]):
            t_ids = get_transformation_ids(a["oper_expression"])
            for t_id in t_ids:
                for oper in mmcif_dict["pdbx_struct_oper_list"]:
                    if oper["id"] == t_id:
                        assembly["transformations"].append({
                         "chains": a["asym_id_list"].split(","),
                         "matrix": [[float(oper["matrix[{}][{}]".format(x, y)])
                          for y in (1, 2, 3)] for x in (1, 2, 3)],
                         "vector": [float(oper["vector[{}]".format(y)])
                          for y in (1, 2, 3)]
                        })


def get_transformation_ids(expression):
    """Takes an operator expression from an .mmcif transformation dict, and
    works out what transformation IDs it is referring to.

    :param str expression: The expression to parse.
    :rtype: ``list``"""

    matches = re.findall(r"\((.+?)\)", expression)
    if matches:
        ids = []
        for match in matches:
            if "," in match:
                ids += match.split(",")
            elif "-" in match:
                bounds = [int(x) for x in match.split("-")]
                ids += list(range(bounds[0], bounds[1] + 1))
            else: ids.append(match)
        return [str(x) for x in ids]
    else:
        return expression.split(",")


def atom_line_to_atom_dict(d):
    """Takes an atom dictionary from an mmcif dictionary and turns it into an
    atomium atom dictionary.

    :param dict d: The mmcif atom dict.
    :rtype: ``dict``"""

    a = deepcopy(ATOM_DICT)
    table = [
     ["id", int, "id"], ["element", lambda v: v, "type_symbol"],
     ["name", lambda v: v, "label_atom_id"],
     ["alt_loc", lambda v: v, "label_alt_id"],
     ["x", float, "Cartn_x"], ["y", float, "Cartn_y"], ["z", float, "Cartn_z"],
     ["residue_id", int, "auth_seq_id"],
     ["residue_name", lambda v: v, "label_comp_id"],
     ["residue_insert", lambda v: v, "pdbx_PDB_ins_code"],
     ["chain_id", lambda v: v, "auth_asym_id"],
     ["bfactor", float, "B_iso_or_equiv"], ["occupancy", float, "occupancy"],
     ["charge", float, "pdbx_formal_charge"],
    ]
    for row in table:
        if d[row[2]] not in "?.":
            a[row[0]] = row[1](d[row[2]])
    a["full_res_id"] = "{}{}{}".format(
     a["chain_id"] or "", str(a["residue_id"] or "") or "", a["residue_insert"]
    ) or None
    a["polymer"] = d["label_seq_id"] != "."
    return a
