
import sys
sys.path.insert(0, ".")
import atomium
from random import shuffle
from big import get_all_codes

print("Getting PDB codes...")
codes = get_all_codes()
print("There are {} codes.".format(len(codes)))

print("Parsing...")
shuffle(codes)
for code in codes:
    print("\tParsing {}.pdb...".format(code))
    pdb = atomium.fetch(code + ".pdb")
    print("\tParsing {}.cif...".format(code))
    pdb = atomium.fetch(code + ".cif")
    print()
