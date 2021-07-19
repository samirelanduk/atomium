import builtins
import requests
import gzip
import paramiko
from .mmcif import mmcif_string_to_mmcif_dict
from .file import File

def open(path, dictionary=False):
    if str(path)[-3:] == ".gz":
        try:
            with gzip.open(path) as f: filestring = f.read().decode()
        except:
            with gzip.open(path, "rt") as f: filestring = f.read()
        return parse_filestring(filestring, dictionary=dictionary)
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
    return parse_filestring(filestring, dictionary=dictionary)


def parse_filestring(filestring, dictionary=False):
    filetype = determine_filetype(filestring)
    mmcif_dict = {
        "mmcif": mmcif_string_to_mmcif_dict
    }[filetype](filestring)
    if dictionary: return mmcif_dict
    return File(mmcif_dict)


def determine_filetype(filestring):
    return "mmcif"