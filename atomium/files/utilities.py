"""This module contains various utility functions for dealing with files."""

import builtins
from requests import get
import paramiko
from .pdb import pdb_string_to_pdb_dict, pdb_dict_to_data_dict
from .mmcif import mmcif_string_to_mmcif_dict, mmcif_dict_to_data_dict
from .xyz import xyz_string_to_xyz_dict, xyz_dict_to_data_dict
from .data import data_dict_to_file

def determine_file_type(path, filestring):
    """Takes a file path and contents, and uses them to work out which of the
    known filetypes it should be interpreted as. File extensions will always be
    believed, and if that can't be used, contents will be examined and a guess
    made.

    :param str path: The structure file path.
    :param str filestring: The structure file contents.
    :rtype: ``str``"""

    if "." in path:
        ending = path.split(".")[-1]
        if ending in ("pdb", "cif", "xyz"):
            return ending
    if "\nATOM" in filestring:
        if "loop_\n" in filestring:
            return "cif"
        else:
            return "pdb"
    return "xyz"


def open(path, *args, **kwargs):
    """Opens a structure file at the given path on disk. Supported filetypes are
    .pdb, .cif, and .xyz - if another file extension (or no extension) is given,
    atomium will use the filecontents to try and guess the format.

    :param str path: The location of the file on disk.
    :rtype: ``Container``"""

    with builtins.open(path) as f:
        filestring = f.read()
    return parse_string(filestring, path, *args, **kwargs)



def fetch(identifier, *args, **kwargs):
    if identifier.startswith("http"):
        url = identifier
    else:
        if "." not in identifier: identifier += ".pdb"
        url = "https://files.rcsb.org/view/" + identifier.lower()
    response = get(url)
    if response.status_code == 200:
        return parse_string(response.text, identifier, *args, **kwargs)
    raise ValueError("Could not find anything at {}".format(url))


def fetch_over_ssh(hostname, username, path, *args, password=None, **kwargs):
    client = paramiko.SSHClient()
    try:
        client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        if not password:
            client.load_system_host_keys()
            client.connect(hostname=hostname, username=username)
        else:
            client.connect(hostname=hostname, username=username, password=password)
        stdin, stdout, stderr = client.exec_command("less " + path)
        filestring = stdout.read().decode()
    finally:
        client.close()
    return parse_string(filestring, path, *args, **kwargs)


def parse_string(filestring, path, file_dict=False, data_dict=False):
    filetype = determine_file_type(path, filestring)
    parsed = {
     "pdb": pdb_string_to_pdb_dict,
     "cif": mmcif_string_to_mmcif_dict,
     "xyz": xyz_string_to_xyz_dict
    }[filetype](filestring)
    if not file_dict:
        parsed = {
         "pdb": pdb_dict_to_data_dict,
         "cif": mmcif_dict_to_data_dict,
         "xyz": xyz_dict_to_data_dict
        }[filetype](parsed)
        if not data_dict:
            parsed = data_dict_to_file(parsed)
            parsed._filetype = filetype
    return parsed
