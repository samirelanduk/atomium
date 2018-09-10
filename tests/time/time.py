import os
import subprocess
import datetime
import sys
sys.path.insert(0, ".")
import atomium
import requests

here = dir_path = os.path.dirname(os.path.realpath(__file__))
version = atomium.__version__.replace(".", "-")
today = datetime.datetime.now().date().strftime("%Y-%m-%d")

pdbs = ["1lol", "5xme", "4f5x"]
for pdb in pdbs:
    for ext in ["cif", "pdb"]:
        try:
            id = pdb + ".cif" if ext == "cif" else "pdb{}.ent".format(pdb)
            filestring = requests.get(
             "http://www.ebi.ac.uk/pdbe/entry-files/{}".format(id)
            ).text
            with open("{}/{}.{}".format(here, pdb, ext), "w") as f:
                f.write(filestring)
            with open("{}/timescript.py".format(here)) as f:
                script = f.read().format(pdb, ext)
            with open("{}/{}.py".format(here, pdb), "w") as f:
                f.write(script)
            try:
                subprocess.call(
                 "python -m cProfile -o {}/profiles/{}.{}-{}-{}.prof {}/{}.py".format(
                  here, pdb, ext, version, today, here, pdb
                 ), shell=True
                )
            finally:
                try:
                    os.remove("{}/{}.py".format(here, pdb))
                    os.remove("{}/{}.{}".format(here, pdb, ext))
                    os.remove("{}/temp.{}".format(here, ext))
                except: pass
        except: pass
