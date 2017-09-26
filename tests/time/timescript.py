import sys
sys.path.insert(0, ".")
import atomium

pdb = atomium.fetch("{}", pdbe=True)
pdb.save("tests/time/temp.pdb")
