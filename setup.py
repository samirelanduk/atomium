from setuptools import setup

with open("README.rst") as f:
    long_description = f.read()

setup(
 name="atomium",
 version="1.0.0",
 description="A molecular modeller and file parser.",
 long_description=long_description,
 long_description_content_type="text/x-rst",
 url="https://atomium.samireland.com",
 author="Sam Ireland",
 author_email="mail@samireland.com",
 license="MIT",
 classifiers=[
  "Development Status :: 4 - Beta",
  "Intended Audience :: Science/Research",
  "License :: OSI Approved :: MIT License",
  "Topic :: Scientific/Engineering :: Chemistry",
  "Programming Language :: Python :: 3",
  "Programming Language :: Python :: 3.6",
  "Programming Language :: Python :: 3.7",
 ],
 keywords="chemistry bioinformatics proteins biochemistry molecules PDB MMCIF CIF MMTF",
 packages=["atomium"],
 python_requires="!=2.*, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, !=3.4.*, !=3.5.*",
 install_requires=["numpy", "requests", "rmsd", "paramiko", "msgpack", "valerius"]
)
