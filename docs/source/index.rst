atomium
=======

atomium is a molecular modeller and file parser.

Example
-------

  >>> import atomium
  >>> pdb = atomium.fetch("5HVD")
  >>> pdb.model()
  <Model (2156 atoms)>
  >>> pdb.model().chain("A")
  <Chain (204 residues)>



Table of Contents
-----------------

.. toctree ::

    installing
    overview
    api
    changelog
