molecuPy Documentation
======================

molecuPy is a Python parser for Protein Data Bank (PDB) files. It provides
utilities for reading and analysing the structural data contained therein.

Example
-------

  >>> import molecupy
  >>> pdb = molecupy.get_pdb_remotely("1LOL")
  >>> pdb.title
  'CRYSTAL STRUCTURE OF OROTIDINE MONOPHOSPHATE DECARBOXYLASE COMPLEX WITH XMP'
  >>> pdb.model
  '<Model (3431 atoms)>'
  >>> pdb.model.get_chain_by_id("A").get_mass()
  20630.8656


Table of Contents
-----------------

.. toctree::
   :maxdepth: 2
