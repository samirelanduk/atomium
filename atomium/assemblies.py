import re
import numpy as np

def extract_assemblies(mmcif_dict):
    assemblies = [{
        "id": int(a["id"]), "software": a.get("method_details", None),
        "delta_energy": None, "buried_surface_area": None, "surface_area": None,
        "transformations": []
    } for a in mmcif_dict.get("pdbx_struct_assembly", [])]
    operations = {o["id"]: [
        [float(o["matrix[{}][{}]".format(r, c)]) for c in [1, 2, 3]
    ] + [float(o["vector[{}]".format(r)])] for r in [1, 2, 3]] + [[0, 0, 0, 1]]
        for o in mmcif_dict.get("pdbx_struct_oper_list", [])}
    for assembly in assemblies:
        if assembly["software"] == "?": assembly["software"] = None
        assign_metrics_to_assembly(mmcif_dict, assembly)
        assign_transformations_to_assembly(mmcif_dict, operations, assembly)
    return assemblies



def assign_metrics_to_assembly(mmcif_dict, assembly):
    """Takes an assembly dict, and goes through an mmcif dictionary looking for
    relevant energy etc. information to update it with.

    :param dict mmcif_dict: The dictionary to read.
    :param dict assembly: The assembly to update."""

    for a in mmcif_dict.get("pdbx_struct_assembly_prop", []):
        if a["biol_id"] == str(assembly["id"]):
            if a["type"] == "MORE":
                assembly["delta_energy"] = float(a["value"].split("/")[0])
            elif a["type"] == "SSA (A^2)":
                assembly["surface_area"] = float(a["value"].split("/")[0])
            elif a["type"] == "ABSA (A^2)":
                assembly["buried_surface_area"] = float(a["value"].split("/")[0])


def assign_transformations_to_assembly(mmcif_dict, operations, assembly):
    """Takes an assembly dict, and goes through an mmcif dictionary looking for
    relevant transformation information to update it with.

    :param dict mmcif_dict: the .mmcif dictionary to read.
    :param dict operations: the processed operations matrices.
    :param dict assembly: the assembly to update."""

    for gen in mmcif_dict.get("pdbx_struct_assembly_gen", []):
        if gen["assembly_id"] == str(assembly["id"]):
            op_ids_groups = get_operation_id_groups(gen["oper_expression"])
            ops = operation_id_groups_to_operations(operations, op_ids_groups)
            for operation in ops:
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

    operation_groups = [[
     operations[i] for i in ids
    ] for ids in operation_id_groups]
    while len(operation_groups) and len(operation_groups) != 1:
        operations = []
        for op1 in operation_groups[0]:
            for op2 in operation_groups[1]:
                operations.append(np.matmul(op1, op2))
        operation_groups[0] = operations
        operation_groups.pop(1)
    return operation_groups[0]