from atomium.files.pdbstring2pdbdict import pdb_string_to_pdb_dict
from pprint import pprint

with open("tests/integration/files/1lol.pdb") as f:
    filestring = f.read()

pdb_dict = pdb_string_to_pdb_dict(filestring)
pdb_dict["connections"] = pdb_dict["connections"][:10]
pdb_dict["models"][0]["molecules"] = pdb_dict["models"][0]["molecules"][:5]
pprint(pdb_dict["models"][0]["chains"][1]["residues"][:2])
