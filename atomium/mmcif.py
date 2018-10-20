"""Contains functions for dealing with the .cif file format."""

from collections import deque
import re
from datetime import datetime
import numpy as np

def mmcif_string_to_mmcif_dict(filestring):
    """Takes a .cif filestring and turns into a ``dict`` which represents its
    table structure. Only lines which aren't empty and which don't begin with
    ``#`` are used.

    Multi-line strings are consolidated onto one line, and the whole thing is
    then split into the blocks that will become table lists. At the end, quote
    marks are removed from any string which retains them.

    :param str filestring: the .cif filestring to process.
    :rtype: ``dict``"""

    lines = deque(filter(lambda l: l and l[0] != "#", filestring.split("\n")))
    lines = consolidate_strings(lines)
    blocks = mmcif_lines_to_mmcif_blocks(lines)
    mmcif_dict = {}
    for block in blocks:
        if block["lines"][0] == "loop_":
            mmcif_dict[block["category"]] = loop_block_to_list(block)
        else:
            mmcif_dict[block["category"]] = non_loop_block_to_list(block)
    strip_quotes(mmcif_dict)
    return mmcif_dict


def consolidate_strings(lines):
    """Generally, .cif files have a one file line to one table row
    correspondence. Sometimes however, a string cell is given a line of its own,
    breaking the row over several lines. This function takes the lines of a .cif
    file and puts all table rows on a single line.

    :param deque lines: the .cif file lines.
    :rtype: ``deque``"""

    new_lines = deque()
    while lines:
        line = lines.popleft()
        if line.startswith(";"):
            string = [line[1:].strip()]
            while not lines[0].startswith(";"):
                string.append(lines.popleft())
            lines.popleft()
            new_lines[-1] += " \"{}\"".format(" ".join(string))
        else:
            new_lines.append(line)
    return new_lines


def mmcif_lines_to_mmcif_blocks(lines):
    """A .cif file is ultimately a list of tables. This function takes a list of
    .cif file lines and splits them into these table blocks. Each block will be
    a ``dict`` containing a category name and a list of lines.

    :param deque lines: the .cif file lines.
    :rtype: ``list``"""

    category = None
    block, blocks = [], []
    while lines:
        line = lines.popleft()
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


def non_loop_block_to_list(block):
    """Takes a simple block ``dict`` with no loop and turns it into a table
    ``list``.

    :param dict block: the .cif block to process.
    :rtype: ``list``"""

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
        d[name] = value
    return [d]


def loop_block_to_list(block):
    """Takes a loop block ``dict`` where the initial lines are table headers and
    turns it into a table ``list``. Sometimes a row is broken over several lines
    so this function deals with that too.

    :param dict block: the .cif block to process.
    :rtype: ``list``"""

    names, lines, header = [], [], True
    body_start = 0
    for index, line in enumerate(block["lines"][1:], start=1):
        if not line.startswith("_" + block["category"]):
            body_start = index
            break
    names = [l.split(".")[1].rstrip() for l in block["lines"][1:body_start]]
    lines = [split_values(l) for l in block["lines"][body_start:]]
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
    """The body of a .cif table is a series of lines, with each cell divided by
    whitespace. This function takes a string line and breaks it into cells.

    There are a few peculiarities to handle. Sometimes a cell is a string
    enclosed in quote marks, and spaces within this string obviously shouldn't
    be used to break the line. This function handles all of that.

    :param str line: the .cif line to split.
    :rtype: ``list``"""

    if not re.search("[\'\"]", line): return line.split()
    chars = deque(line.strip())
    values, value, in_string = [], [], False
    while chars:
        char = chars.popleft()
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
    """Goes through each table in the mmcif ``dict`` and removes any unneeded
    quote marks from the cells.

    :param dict mmcif_dict: the almost finished .mmcif dictionary to clean."""

    for name, table in mmcif_dict.items():
        for row in table:
            for key, value in row.items():
                for char in "'\"":
                    if value[0] == char and value[-1] == char:
                        row[key] = value[1:-1]


def mmcif_dict_to_data_dict(mmcif_dict):
    """Converts an .mmcif dictionary into an atomium data dictionary, with the
    same standard layout that the other file formats get converted into.

    :param dict mmcif_dict: the .mmcif dictionary.
    :rtype: ``dict``"""

    data_dict = {
     "description": {
      "code": None, "title": None, "deposition_date": None,
      "classification": None, "keywords": [], "authors": []
     }, "experiment": {
      "technique": None, "source_organism": None, "expression_system": None
     }, "quality": {"resolution": None, "rvalue": None, "rfree": None},
     "geometry": {"assemblies": []}
    }
    update_description_dict(mmcif_dict, data_dict)
    update_experiment_dict(mmcif_dict, data_dict)
    update_quality_dict(mmcif_dict, data_dict)
    update_geometry_dict(mmcif_dict, data_dict)
    return data_dict


def update_description_dict(mmcif_dict, data_dict):
    """Takes a data dictionary and updates its description sub-dictionary with
    information from a .mmcif dictionary.

    :param dict mmcif_dict: the .mmcif dictionary to read.
    :param dict data_dict: the data dictionary to update."""

    mmcif_to_data_transfer(mmcif_dict, data_dict,
     "description", "code", "entry", "id")
    mmcif_to_data_transfer(mmcif_dict, data_dict,
     "description", "title", "struct", "title")
    mmcif_to_data_transfer(
     mmcif_dict, data_dict, "description", "deposition_date",
     "pdbx_database_status", "recvd_initial_deposition_date", date=True
    )
    mmcif_to_data_transfer(mmcif_dict, data_dict,
     "description", "classification", "struct_keywords", "pdbx_keywords")
    mmcif_to_data_transfer(mmcif_dict, data_dict,
     "description", "keywords", "struct_keywords", "text", split=True)
    mmcif_to_data_transfer(mmcif_dict, data_dict,
     "description", "authors", "audit_author", "name", multi=True)


def update_experiment_dict(mmcif_dict, data_dict):
    """Takes a data dictionary and updates its experiment sub-dictionary with
    information from a .mmcif dictionary.

    :param dict mmcif_dict: the .mmcif dictionary to read.
    :param dict data_dict: the data dictionary to update."""

    mmcif_to_data_transfer(mmcif_dict, data_dict,
     "experiment", "technique", "exptl", "method")
    mmcif_to_data_transfer(mmcif_dict, data_dict, "experiment",
     "source_organism", "entity_src_gen", "pdbx_gene_src_scientific_name")
    mmcif_to_data_transfer(mmcif_dict, data_dict, "experiment",
     "expression_system", "entity_src_gen", "pdbx_host_org_scientific_name")


def update_quality_dict(mmcif_dict, data_dict):
    """Takes a data dictionary and updates its quality sub-dictionary with
    information from a .mmcif dictionary.

    :param dict mmcif_dict: the .mmcif dictionary to read.
    :param dict data_dict: the data dictionary to update."""

    mmcif_to_data_transfer(mmcif_dict, data_dict,
     "quality", "resolution", "reflns", "d_resolution_high", func=float)
    mmcif_to_data_transfer(mmcif_dict, data_dict,
     "quality", "rvalue", "refine", "ls_R_factor_R_work", func=float)
    mmcif_to_data_transfer(mmcif_dict, data_dict,
     "quality", "rfree", "refine", "ls_R_factor_R_free", func=float)


def update_geometry_dict(mmcif_dict, data_dict):
    """Takes a data dictionary and updates its geometry sub-dictionary with
    information from a .mmcif dictionary.

    :param dict mmcif_dict: the .mmcif dictionary to read.
    :param dict data_dict: the data dictionary to update."""

    data_dict["geometry"]["assemblies"] = [{
     "id": int(a["id"]), "software": a.get("method_details", None),
     "delta_energy": None, "buried_surface_area": None, "surface_area": None,
     "transformations": []
    } for a in mmcif_dict.get("pdbx_struct_assembly", [])]
    operations = {o["id"]: [
     [float(o["matrix[{}][{}]".format(r, c)]) for c in [1, 2, 3]
    ] + [float(o["vector[{}]".format(r)])] for r in [1, 2, 3]] + [[0, 0, 0, 1]]
     for o in mmcif_dict.get("pdbx_struct_oper_list", [])}
    for assembly in data_dict["geometry"]["assemblies"]:
        if assembly["software"] == "?": assembly["software"] = None
        assign_metrics_to_assembly(mmcif_dict, assembly)
        assign_transformations_to_assembly(mmcif_dict, operations, assembly)


def assign_metrics_to_assembly(mmcif_dict, assembly):
    """Takes an assembly dict, and goes through an mmcif dictionary looking for
    relevant energy etc. information to update it with.

    :param dict mmcif_dict: The dictionary to read.
    :param dict assembly: The assembly to update."""

    for a in mmcif_dict.get("pdbx_struct_assembly_prop", []):
        if a["biol_id"] == str(assembly["id"]):
            if a["type"] == "MORE":
                assembly["delta_energy"] = float(a["value"])
            elif a["type"] == "SSA (A^2)":
                assembly["surface_area"] = float(a["value"])
            elif a["type"] == "ABSA (A^2)":
                assembly["buried_surface_area"] = float(a["value"])


def assign_transformations_to_assembly(mmcif_dict, operations, assembly):
    """Takes an assembly dict, and goes through an mmcif dictionary looking for
    relevant transformation information to update it with.

    :param dict mmcif_dict: the .mmcif dictionary to read.
    :param dict operations: the processed operations matrices.
    :param dict assembly: the assembly to update."""

    for gen in mmcif_dict.get("pdbx_struct_assembly_gen", []):
        if gen["assembly_id"] == str(assembly["id"]):
            op_ids_groups = get_operation_id_groups(gen["oper_expression"])
            operations = operation_id_groups_to_operations(operations, op_ids_groups)
            for operation in operations:
                assembly["transformations"].append({
                 "chains": gen["asym_id_list"].split(","),
                 "matrix": [row[:3] for row in operation[:3]],
                 "vector": [row[-1] for row in operation[:3]]
                })


def get_operation_id_groups(expression):
    """Takes an operator expression from an .mmcif transformation dict, and
    works out what transformation IDs it is referring to. For example, (1,2,3)
    becomes [[1, 2, 3]], (1-3)(8-11,17) becomes [[1, 2, 3], [8, 9, 10, 11, 17]],
    and so on.

    :param str expression: The expression to parse.
    :rtype: ``list``"""

    if expression[0] != "(": expression = "({})".format(expression)
    groups = re.findall(r"\((.+?)\)", expression)
    group_ids = []
    for group in groups:
        ids = []
        elements = group.split(",")
        for element in elements:
            if "-" in element:
                bounds = [int(x) for x in element.split("-")]
                ids += [str(n) for n in list(range(bounds[0], bounds[1] + 1))]
            else:
                ids.append(element)
        group_ids.append(ids)
    return group_ids


def operation_id_groups_to_operations(operations, operation_id_groups):
    """Creates a list of operation matrices for an assembly, from a list of
    operation IDs - cross multiplying as required.

    :param dict operations: the parsed .mmcif operations.
    :param list operation_id_groups: the operation IDs."""

    operation_groups = [[operations[i] for i in ids] for ids in operation_id_groups]
    while len(operation_groups) and len(operation_groups) != 1:
        operations = []
        for op1 in operation_groups[0]:
            for op2 in operation_groups[1]:
                operations.append(
                 np.matmul(op1, op2),
                )
        operation_groups[0] = operations
        operation_groups.pop(1)
    return operation_groups[0]


def mmcif_to_data_transfer(mmcif_dict, data_dict, d_cat, d_key, m_table, m_key,
                           date=False, split=False, multi=False, func=None):
    """A function for transfering a bit of data from a .mmcif dictionary to a
    data dictionary, or doing nothing if the data doesn't exist.

    :param dict mmcif_dict: the .mmcif dictionary to read.
    :param dict data_dict: the data dictionary to update.
    :param str d_cat: the top-level key in the data dictionary.
    :param str d_key: the data dictionary field to update.
    :param str m_table: the name of the .mmcif table to look in.
    :param str m_key: the .mmcif field to read.
    :param bool date: if True, the value will be converted to a date.
    :param bool split: if True, the value will be split on commas.
    :param bool multi: if True, every row in the table will be read.
    :param function func: if given, this will be applied to the value."""

    try:
        if multi:
            value = [row[m_key] for row in mmcif_dict[m_table]]
        else:
            value = mmcif_dict[m_table][0][m_key]
        if date: value = datetime.strptime(value, "%Y-%m-%d").date()
        if split: value = value.replace(", ", ",").split(",")
        if func: value = func(value)
        data_dict[d_cat][d_key] = value
    except: pass
