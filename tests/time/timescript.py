import sys
sys.path.insert(0, ".")
import atomium

pdb = atomium.pdb_from_file("tests/time/{}.pdb")
pdb.save("tests/time/temp.pdb")
