import os
import sys
sys.path.append(".")
import random
import atomium
import argparse
import pdbsearch
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument("--count", help="Number of structures", type=int, required=False)
parser.add_argument("--code", help="Specific code to check", type=str, required=False)
args = parser.parse_args()

if args.code:
    codes = [args.code]
else:
    print("Getting PDB codes in random order...")
    codes = pdbsearch.search(limit=None)
    random.shuffle(codes)
    codes = codes[:args.count or len(codes)]

print("Checking...")
for code in tqdm(codes):
    # mmCIF parse/save
    try:
        mmcif = atomium.fetch(code, dictionary=True)
        atomium.save_dictionary(mmcif, "saved.cif")
        saved = atomium.open("saved.cif", dictionary=True)
        assert mmcif == saved
    except Exception as e:
        print("mmCIF", code)
        raise e
    finally:
        if os.path.exists("saved.cif"): os.remove("saved.cif")
    
    # BinaryCIF parse/save
    try:
        mmcif = atomium.fetch(f"{code}.bcif", dictionary=True)
        atomium.save_dictionary(mmcif, "saved.bcif")
        saved = atomium.open("saved.bcif", dictionary=True)
        assert mmcif == saved
    except Exception as e:
        print("BinaryCIF", code)
        raise e
    finally:
        if os.path.exists("saved.bcif"): os.remove("saved.bcif")
    
    # PDB parse/save
    try:
        mmcif = atomium.fetch(f"{code}.pdb", dictionary=True)
        atomium.save_dictionary(mmcif, "saved.pdb")
    except Exception as e:
        print("PDB", code)
        raise e
    finally:
        if os.path.exists("saved.pdb"): os.remove("saved.pdb")