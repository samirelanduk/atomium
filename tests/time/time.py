from datetime import datetime
import subprocess
import os

# What files do we want to time?
files = [
 "1lol.cif", "5xme.cif", "4v6x.cif",
 "1lol.mmtf", "4v6x.mmtf"
]

# Make a directory for these profiles
now = datetime.now().strftime("%Y-%m-%d_%H%M%S")
subprocess.call("mkdir tests/time/profiles/{}".format(now), shell=True)

# Go through each file to time
for filename in files:
    # Make a script to time the run
    with open("tests/time/savescript.py") as f:
        script = f.read().format(filename)
    with open("tests/time/{}.py".format(filename), "w") as f:
        f.write(script)

    # Run script and time it
    subprocess.call(
     "python -m cProfile -o tests/time/profiles/{}/{}.prof tests/time/{}.py".format(
      now, filename.replace(".", "-"), filename
     ), shell=True
    )

    # Tidy up
    try:
        os.remove("tests/time/{}.py".format(filename))
    except: pass
