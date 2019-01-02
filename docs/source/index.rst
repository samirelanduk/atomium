atomium
========

atomium is a molecular modeller and file parser, capable of reading from and
writing to .pdb, .cif and .mmtf files.

Example
-------

    >>> import atomium
    >>> pdb = atomium.fetch("5HVD")
    >>> pdb.model
    <Model (1 chain, 6 ligands)>
    >>> pdb.model.chain("A")
    <Chain A (255 residues)>

Table of Contents
-----------------

.. toctree ::
  installing
  overview
  api
  contributing
  changelog
