#! /usr/bin/env python3

"""This script converts a protein structure file to JSON. The file will be saved
in the same location as the original file."""

import sys
import os
import atomium
import json

if len(sys.argv) < 2:
    print("Please provide a file to convert")
    sys.exit()

path = sys.argv[1]
data = atomium.open(path, data_dict=True)
location = os.path.sep.join(path.split(os.path.sep)[:-1])
filename = path.split(os.path.sep)[-1].split(".")[0]

with open(f"{location}/{filename}.json", "w") as f:
    json.dump(data, f, default=str)
