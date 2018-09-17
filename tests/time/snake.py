import os
import sys
from fnmatch import fnmatch
import subprocess
import signal
from time import sleep

# Get profiles
profile_dirs = [d for d in os.listdir("tests/time/profiles") if "." not in d]
profile_dirs = sorted(profile_dirs)

# Which set of profiles?
index = int(sys.argv[1]) if len(sys.argv) >= 2 else -1
profile_dir = profile_dirs[index]

# Which profiles should be used within this set?
profiles = [f for f in os.listdir("tests/time/profiles/{}".format(profile_dir))
 if f.endswith(".prof")]
if len(sys.argv) >= 3:
    profiles = [f for f in profiles if fnmatch(f, sys.argv[2] + ".prof")]

# Run each profile
for profile in profiles:
    pro = subprocess.Popen(
     "snakeviz tests/time/profiles/{}/{}".format(profile_dir, profile),
     stdout=subprocess.PIPE,
     shell=True,
     preexec_fn=os.setsid
    )
    sleep(1.5)
    os.killpg(os.getpgid(pro.pid), signal.SIGTERM)
