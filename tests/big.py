import sys
sys.path.insert(0, ".")
import atomium
import requests
from random import shuffle
from tqdm import tqdm
import xml.etree.ElementTree as ET
print()

SUBSET = 5000

# Get all codes
response = requests.get("https://www.rcsb.org/pdb/rest/getCurrent")
codes = [child.attrib["structureId"] for child in ET.fromstring(response.text)]
print(f"There are {len(codes)} codes")

# Go through them
print(f"Processing a random {SUBSET} of them...")
shuffle(codes)
sub_codes = codes[:SUBSET]
results = {}
for code in tqdm(sub_codes):
    results[code] = {}
    for ext in ("cif", "mmtf", "pdb"):
        try:
            pdb = atomium.fetch(f"{code}.{ext}")
            results[code][ext] = str(pdb.model)
        except Exception as e:
            results[code][ext] = str(e)
    models = [results[code][ext] for ext in ("cif", "mmtf", "pdb")
     if results[code][ext][0] == "<"]
    results[code]["match"] = len(set(models)) == 1

no_pdb = [c for c in sub_codes if "could not find" in results[c]["pdb"].lower()]
print(f"{len(no_pdb)} codes had no .pdb representation:")
print((" ".join(no_pdb) + "\n") if no_pdb else "")


mismatch = [code for code in sub_codes if not results[code]["match"]]
print(f"{len(mismatch)} codes had different models for different extensions:")
print((" ".join(mismatch) + "\n") if mismatch else "")

print("Other errors:")
for code in sub_codes:
    for ext in ("cif", "mmtf", "pdb"):
        if results[code][ext][0] != "<":
            if ext != "pdb" or code not in no_pdb:
                print(" ", code, ext, results[code][ext])
print()
