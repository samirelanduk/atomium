from setuptools import setup

setup(
 name="molecupy",
 version="1.0.3",
 description="A Python PDB parser.",
 url="https://molecupy.readthedocs.org",
 author="Sam Ireland",
 author_email="Sam.Ireland@ed.ac.uk",
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
 packages=["molecupy", "molecupy.structures"],
 install_requires=["requests", "omnicanvas"]
)
