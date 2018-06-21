
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
    print("\tParsing {}...".format(code))
    pdb = atomium.fetch(code)
    print("\tSuccess ({} model{}).".format(
     len(pdb.models), "" if len(pdb.models) == 1 else "s"
    ))
    print("\t", *[b["software"] for b in pdb.biomolecules])
    print()
