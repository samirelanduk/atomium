import builtins
import requests
import gzip
import paramiko
from .mmcif import mmcif_string_to_mmcif_dict
from .bcif import bcif_string_to_mmcif_dict
from .mmtf import mmtf_string_to_mmcif_dict
from .pdb import pdb_string_to_mmcif_dict
from .mmcif import save_mmcif_dict as save_mmcif_dict_to_mmcif
from .bcif import save_mmcif_dict as save_mmcif_dict_to_bcif
from .pdb import save_mmcif_dict as save_mmcif_dict_to_pdb
from .file import File

def open(path, dictionary=False):
    if str(path)[-3:] == ".gz":
        try:
            with gzip.open(path) as f: filestring = f.read().decode()
        except:
            with gzip.open(path, "rb") as f: filestring = f.read()
        return parse_filestring(filestring, path[:-3], dictionary=dictionary)
    try:
        with builtins.open(path) as f:
            filestring = f.read()
    except:
        with builtins.open(path, "rb") as f:
            filestring = f.read()
    return parse_filestring(filestring, path, dictionary=dictionary)


def fetch(code, dictionary=False):
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
        text = response.content if code.endswith(".bcif") or code.endswith(".mmtf") else response.text
        return parse_filestring(text, code, dictionary=dictionary)
    raise ValueError("Could not find anything at {}".format(url))


def fetch_over_ssh(hostname, username, path, password=None, dictionary=False):
    client = paramiko.SSHClient()
    try:
        client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        if not password:
            client.load_system_host_keys()
            client.connect(hostname=hostname, username=username)
        else:
            client.connect(
                hostname=hostname, username=username, password=password 
            )
        stdout = client.exec_command("less " + path)[1]
        filestring = stdout.read().decode()
    finally:
        client.close()
    return parse_filestring(filestring, path, dictionary=dictionary)


def parse_filestring(filestring, filename, dictionary=False):
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
    if isinstance(filestring, bytes):
        return "mmtf" if filename.endswith("mmtf") else "bcif"
    return "pdb" if filename.endswith("pdb") else "mmcif"


def save_dictionary(d, path):
    ext = str(path).split(".")[-1]
    {
        "cif": save_mmcif_dict_to_mmcif,
        "bcif": save_mmcif_dict_to_bcif,
        "pdb": save_mmcif_dict_to_pdb,
    }[ext](d, path)