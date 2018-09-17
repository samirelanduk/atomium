"""Contains various file handling helper functions."""

import builtins
from requests import get
from .mmcif import mmcif_string_to_mmcif_dict, mmcif_dict_to_data_dict
from .mmtf import mmtf_bytes_to_mmtf_dict
from .pdb import pdb_string_to_pdb_dict

def open(path, *args, **kwargs):
    try:
        with builtins.open(path) as f: filestring = f.read()
    except:
        with builtins.open(path, "rb") as f: filestring = f.read()
    return parse_string(filestring, path, *args, **kwargs)


def fetch(code, *args, **kwargs):
    if code.startswith("http"):
        url = code
    elif code.endswith(".mmtf"):
        url = "https://mmtf.rcsb.org/v1.0/full/{}".format(code[:-5])
    else:
        if "." not in code: code += ".cif"
        url = "https://files.rcsb.org/view/" + code.lower()
    response = get(url, stream=True)
    if response.status_code == 200:
        text = response.content if code.endswith(".mmtf") else response.text
        return parse_string(text, code, *args, **kwargs)
    raise ValueError("Could not find anything at {}".format(url))


def parse_string(filestring, path, file_dict=False, data_dict=False):
    file_func, data_func = get_parse_functions(filestring, path)
    parsed = file_func(filestring)
    if not file_dict:
        parsed = data_func(parsed)
    return parsed


def get_parse_functions(filestring, path):
    if "." in path:
        ending = path.split(".")[-1]
        if ending in ("mmtf", "cif", "pdb"):
            return {
             "cif": (mmcif_string_to_mmcif_dict, mmcif_dict_to_data_dict),
             "mmtf": (mmtf_bytes_to_mmtf_dict, None),
             "pdb": (pdb_string_to_pdb_dict, None)
            }[ending]
