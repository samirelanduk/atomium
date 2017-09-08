import os
import subprocess
import atomium

here = dir_path = os.path.dirname(os.path.realpath(__file__))
version = atomium.__version__

pdbs = ["1LOL"]
for pdb in pdbs:
    with open("{}/timescript.py".format(here)) as f:
        script = f.read().format(pdb)
    with open("{}/{}.py".format(here, pdb), "w") as f:
        f.write(script)
    try:
        subprocess.call(
         "python -m cProfile -o {}/profiles/{}-{}.prof {}/{}.py".format(
          here, pdb, version, here, pdb
         ), shell=True
        )
    finally:
        os.remove("{}/{}.py".format(here, pdb))
