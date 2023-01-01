"""Functions for getting and saving files."""

import gzip
import builtins
import requests
from atomium.mmcif import mmcif_string_to_mmcif_dict
from atomium.mmcif import save_mmcif_dict as save_mmcif_dict_to_mmcif
from atomium.bcif import bcif_string_to_mmcif_dict
from atomium.bcif import save_mmcif_dict as save_mmcif_dict_to_bcif
from atomium.pdb import pdb_string_to_mmcif_dict
from atomium.pdb import save_mmcif_dict as save_mmcif_dict_to_pdb
from atomium.mmtf import mmtf_string_to_mmcif_dict
from atomium.mmtf import save_mmcif_dict as save_mmcif_dict_to_mmtf
from atomium.file import File

def open(path, dictionary=False):
    """Opens a local structure file. The file can be gzipped or uncompressed,
    and you can opt to only get a dictionary representation if you wish. The
    file extension is used to determine how to parse it.
    
    :param str path: the location of the file:
    :param bool dictionary: whether to stop at dictionary parsing.
    :returns: ``File`` or ``dict``"""

    path = str(path)
    if path[-3:] == ".gz":
        try:
            with gzip.open(path) as f: filestring = f.read().decode()
        except:
            with gzip.open(path, "rb") as f: filestring = f.read()
        return parse_filestring(filestring, path[:-3], dictionary=dictionary)
    try:
        with builtins.open(path) as f: filestring = f.read()
    except:
        with builtins.open(path, "rb") as f: filestring = f.read()
    return parse_filestring(filestring, path, dictionary=dictionary)


def fetch(code, dictionary=False):
    """Fetches a remote structure file. You can provide a PDB code or a generic
    URL, and you can opt to only get a dictionary representation if you wish.
    The file extension is used to determine what format to get.
    
    :param str code: the PDB code (or URL).
    :param bool dictionary: whether to stop at dictionary parsing.
    :returns: ``File`` or ``dict``"""

    if code.startswith("http"):
        url = code
    elif code.endswith(".bcif"):
        url = "https://models.rcsb.org/{}.bcif".format(code[:-5].lower())
    elif code.endswith(".mmtf"):
        url = "https://mmtf.rcsb.org/v1.0/full/{}".format(code[:-5].lower())
    else:
        if "." not in code: code += ".cif"
        url = "https://files.rcsb.org/view/" + code.lower()
    response = requests.get(url, stream=True)
    if response.status_code == 200:
        if code.endswith(".bcif") or code.endswith(".mmtf"):
            text = response.content
        else:
            text = response.text
        return parse_filestring(text, code, dictionary=dictionary)
    raise ValueError("Could not find anything at {}".format(url))


def parse_filestring(filestring, filename, dictionary=False):
    """Converts a filestring to a Python dictionary or :py:class:`File`, with
    the specific parsing done based on file extension and file contents.
    
    :param filestring: the contents of a structure file.
    :param str filename: the name of the structure (with extension).
    :param bool dictionary: whether to stop at dictionary parsing.
    :returns: ``File`` or ``dict``"""

    filetype = determine_filetype(filestring, filename)
    mmcif_dict = {
        "mmcif": mmcif_string_to_mmcif_dict,
        "bcif": bcif_string_to_mmcif_dict,
        "mmtf": mmtf_string_to_mmcif_dict,
        "pdb": pdb_string_to_mmcif_dict,
    }[filetype](filestring)
    if dictionary: return mmcif_dict
    return File(mmcif_dict)


def determine_filetype(filestring, filename):
    """Works out the filetype based on filestring and extension.
    
    :param filestring: the contents of a structure file.
    :param str filename: the name of the structure (with extension).
    :rtype: ``str``"""

    if isinstance(filestring, bytes):
        return "mmtf" if filename.endswith("mmtf") else "bcif"
    return "pdb" if filename.endswith("pdb") else "mmcif"


def save_dictionary(mmcif, path):
    """Saves an mmCIF dictionary to disk, based on the extension given.
    
    :param dict mmcif: the mmIF dictionary representation.
    :param path: the location to save to."""

    ext = str(path).split(".")[-1]
    {
        "cif": save_mmcif_dict_to_mmcif,
        "bcif": save_mmcif_dict_to_bcif,
        "pdb": save_mmcif_dict_to_pdb,
        "mmtf": save_mmcif_dict_to_mmtf,
    }[ext](mmcif, path)