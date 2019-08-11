import requests
import random
import subprocess
import sys
import os
import json
from datetime import datetime
sys.path.append(os.path.join("..", "atomium"))
import atomium
import io
from Bio.PDB import *
import matplotlib.pyplot as plt

def get_string(code):
    if code.endswith(".mmtf"):
        url = "https://mmtf.rcsb.org/v1.0/full/{}".format(code[:-5].lower())
    else:
        url = "https://files.rcsb.org/view/" + code.lower()
    response = requests.get(url, stream=True)
    if response.status_code == 200:
        text = response.content if code.endswith(".mmtf") else response.text
    return text


if len(sys.argv) > 1 and sys.argv[1] == "--rebuild":
    query = "<orgPdbQuery>"\
    "<queryType>org.pdb.query.simple.ChemCompFormulaQuery</queryType>"\
    "<formula>ZN</formula></orgPdbQuery>"
    url = "https://www.rcsb.org//pdb/rest/search/"
    response = requests.post(url, data=query.encode(), headers={
     "Content-Type": "application/x-www-form-urlencoded"
    })
    codes = response.text.split()

    data = []
    while len (data) != 1000:
        code = random.choice(codes)
        d = {"code": code}
        print(len(data) + 1, code)
        for ext in [".cif", ".mmtf", ".pdb"]:
            try:
                string = get_string(code + ext)
                with open("temp" + ext, "wb" if ext == ".mmtf" else "w") as f:
                    f.write(string)
                start = datetime.now()
                pdb = atomium.open("temp" + ext)
                end = datetime.now()
                delta = end - start
                d[ext] = delta.total_seconds()
                if "atoms" not in d:
                    d["atoms"] = len(pdb.model.atoms())
                if "models" not in d:
                    d["models"] = len(pdb.models)
            except:
                d[ext] = None

        if d[".pdb"]:
            try:
                parser = PDBParser(QUIET=True, get_header=True)
                start = datetime.now()
                parser.get_structure("", "temp.pdb")
                end = datetime.now()
                delta = end - start
                d["biopython"] = delta.total_seconds()
            except:
                d["biopython"] = None
             
            try:
                with open("temp.pdb", "w") as f: f.write(string)
                start = datetime.now()
                subprocess.check_output("scripts/readapdb temp.pdb", shell=True)
                end = datetime.now()
                delta = end - start
                d["bioplib"] = delta.total_seconds()
            except:
                d["bioplib"] = None

        else:
            d["biopython"] = None
            d["bioplib"] = None

        data.append(d)
        codes.remove(code)

        for f in ["temp.pdb", "temp.cif", "temp.mmtf"]:
            if os.path.exists(f):
                os.remove(f)
        with open("scripts/speed.json", "w") as f:
            json.dump(data, f)


with open("scripts/speed.json") as f:
    data = json.load(f)

print("There are {} data points".format(len(data)))

pdbs_x, pdbs_y = zip(*[[d["atoms"], d[".pdb"]] for d in data if d[".pdb"] and d["models"] == 1 and d[".pdb"] < 20])
cifs_x, cifs_y = zip(*[[d["atoms"], d[".cif"]] for d in data if d[".cif"] and d["models"] == 1 and d[".cif"] < 20])
mmtfs_x, mmtfs_y = zip(*[[d["atoms"], d[".mmtf"]] for d in data if d[".mmtf"] and d["models"] == 1 and d[".mmtf"] < 20])
bios_x, bios_y = zip(*[[d["atoms"], d["biopython"]] for d in data if d["biopython"] and d["models"] == 1 and d["biopython"] < 20])
bioplib_x, bioplib_y = zip(*[[d["atoms"], d["bioplib"]] for d in data if d["bioplib"] and d["models"] == 1 and d["bioplib"] < 20])

def best_fit(X, Y, label):
    xbar = sum(X)/len(X)
    ybar = sum(Y)/len(Y)
    n = len(X) # or len(Y)
    numer = sum([xi*yi for xi,yi in zip(X, Y)]) - n * xbar * ybar
    denum = sum([xi**2 for xi in X]) - n * xbar**2
    b = numer / denum
    a = ybar - b * xbar
    print("{} best fit line:\ny = {:.6f} + {:.6f}x".format(label, a, b))
    return a, b

plt.xscale("log")
plt.yscale("log")
best_fit(cifs_x, cifs_y, "cif")
best_fit(pdbs_x, pdbs_y, "pdb")
best_fit(mmtfs_x, mmtfs_y, "mmtf")
best_fit(bios_x, bios_y, "biopython")
best_fit(bioplib_x, bioplib_y, "bioplib")
plt.scatter(cifs_x, cifs_y, s=6, c="#FD7272", label=".cif", alpha=0.3, linewidths=0)
plt.scatter(pdbs_x, pdbs_y, s=6, c="#58B19F", label=".pdb", alpha=0.3, linewidths=0)
plt.scatter(mmtfs_x, mmtfs_y, s=6, c="#182C61", label=".mmtf", alpha=0.3, linewidths=0)
plt.xlabel("Atom Count")
plt.ylabel("Parse time (s)")
plt.xlim([100, 1000000])
plt.ylim([0.001, 100])
plt.legend(loc=2)
plt.savefig("scripts/format-speed.png", dpi=1000)
plt.clf()


plt.xscale("log")
plt.yscale("log")
plt.scatter(pdbs_x, pdbs_y, s=6, c="#58B19F", label=".pdb (atomium)", alpha=0.5, linewidths=0)
plt.scatter(bios_x, bios_y, s=6, c="#D6A2E8", label=".pdb (biopython)", alpha=0.5, linewidths=0)
plt.scatter(bioplib_x, bioplib_y, s=6, c="#333333", label=".pdb (bioplib)", alpha=0.5, linewidths=0)
plt.xlabel("Atom Count")
plt.ylabel("Parse time (s)")
plt.xlim([100, 100000])
plt.ylim([0.001, 10])
plt.legend(loc=2)
plt.savefig("scripts/library-speed.png", dpi=1000)
plt.clf()

