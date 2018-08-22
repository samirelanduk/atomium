
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
    try:
        pdb = atomium.fetch(code + ".pdb")
    except ValueError:
        print("    Doesn't exist.")
        pdb = None
    print("\tParsing {}.cif...".format(code))
    cif = atomium.fetch(code + ".cif")

    if pdb:
        assert len(pdb.model.chains()) == len(cif.model.chains())
        assert len(pdb.model.residues()) == len(cif.model.residues())
        assert len(pdb.model.ligands()) == len(cif.model.ligands())
        assert len(pdb.model.atoms()) == len(cif.model.atoms())
        assert len(pdb.assemblies) == len(cif.assemblies)
    print()
