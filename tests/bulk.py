import os
import sys
sys.path.append(".")
import random
import atomium
import argparse
import pdbsearch
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument("count", help="Number of structures", type=int)
args = parser.parse_args()

print("Getting PDB codes in random order...")
codes = pdbsearch.search(limit=None)
random.shuffle(codes)
codes = codes[:args.count]

print("Checking...")
for code in tqdm(codes):
    # mmCIF parse/save
    try:
        mmcif = atomium.fetch(code, dictionary=True)
        atomium.save_dictionary(mmcif, "saved.cif")
        saved = atomium.open("saved.cif", dictionary=True)
        assert mmcif == saved
    except Exception as e:
        print(code)
        raise e
    finally:
        if os.path.exists("saved.cif"): os.remove("saved.cif")