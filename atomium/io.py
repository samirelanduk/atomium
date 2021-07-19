import builtins
import requests
from .mmcif import mmcif_string_to_mmcif_dict
from .file import File

def open(path, dictionary=False):
    try:
        with builtins.open(path) as f:
            filestring = f.read()
    except:
        with builtins.open(path, "rb") as f:
            filestring = f.read()
    return parse_filestring(filestring, dictionary=dictionary)


def fetch(code, dictionary=False):
    if code.startswith("http"):
        url = code
    elif code.endswith(".mmtf"):
        url = "https://mmtf.rcsb.org/v1.0/full/{}".format(code[:-5].lower())
    else:
        if "." not in code: code += ".cif"
        url = "https://files.rcsb.org/view/" + code.lower()
    response = requests.get(url, stream=True)
    if response.status_code == 200:
        text = response.content if code.endswith(".mmtf") else response.text
        return parse_filestring(text, dictionary=dictionary)
    raise ValueError("Could not find anything at {}".format(url))


def parse_filestring(filestring, dictionary=False):
    filetype = determine_filetype(filestring)
    mmcif_dict = {
        "mmcif": mmcif_string_to_mmcif_dict
    }[filetype](filestring)
    if dictionary: return mmcif_dict
    return File(mmcif_dict)


def determine_filetype(filestring):
    return "mmcif"