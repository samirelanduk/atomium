import sys
sys.path.insert(0, ".")
import atomium

pdb = atomium.open("tests/time/{}.{}")
#pdb.save("tests/time/temp.pdb")
