"""Contains various file handling helper functions."""

import builtins
from .mmcif import mmcif_string_to_mmcif_dict

def open(path, *args, **kwargs):
    with builtins.open(path) as f:
        filestring = f.read()
    return parse_string(filestring, path, *args, **kwargs)


def parse_string(filestring, path, file_dict=True, data_dict=True):
    file_func, data_func = get_parse_functions(filestring, path)
    parsed = file_func(filestring)
    return parsed


def get_parse_functions(filestring, path):
    return (mmcif_string_to_mmcif_dict, None)
