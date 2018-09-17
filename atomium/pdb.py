"""Contains functions for dealing with the .mmtf file format."""

from datetime import datetime
import re

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
     }, "quality": {"resolution": None, "rvalue": None, "rfree": None}
    }
    update_description_dict(pdb_dict, data_dict)
    update_experiment_dict(pdb_dict, data_dict)
    update_quality_dict(pdb_dict, data_dict)
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


def merge_lines(lines, start, join=" "):
    """Gets a single continuous string from a sequence of lines.

    :param list lines: The lines to merge.
    :param int start: The start point in each record.
    :param str join: The string to join on.
    :rtype: ``str``"""
    
    string = join.join([line[start:].strip() for line in lines])
    return string
