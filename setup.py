from setuptools import setup

setup(
 name="molecupy",
 version="2.0.0",
 description="A Python molecular modeller with PDB parsing.",
 url="https://molecupy.readthedocs.io",
 author="Sam Ireland",
 author_email="mail@samireland.com",
 license="MIT",
 classifiers=[
  "Development Status :: 4 - Beta",
  "Intended Audience :: Science/Research",
  "License :: OSI Approved :: MIT License",
  "Topic :: Scientific/Engineering :: Chemistry",
  "Programming Language :: Python :: 3",
  "Programming Language :: Python :: 3.0",
  "Programming Language :: Python :: 3.1",
  "Programming Language :: Python :: 3.2",
  "Programming Language :: Python :: 3.3",
  "Programming Language :: Python :: 3.4",
  "Programming Language :: Python :: 3.5",
 ],
 keywords="chemistry bioinformatics proteins biochemistry molecules",
 packages=["molecupy", "molecupy.structures", "molecupy.pdb"],
 install_requires=["requests", "omnicanvas"]
)
