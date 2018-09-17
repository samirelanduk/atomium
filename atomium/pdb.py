"""Contains functions for dealing with the .mmtf file format."""

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
