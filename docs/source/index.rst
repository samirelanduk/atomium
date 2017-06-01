atomium
=======

atomium is a molecular modeller and file parser.

Example
-------

  >>> import atomium
  >>> glucose = atomium.xyz_from_file("glucose.xyz")
  >>> glucose.comment()
  'glucose from 2gbp'
  >>> glucose.model()
  <Model (12 atoms)>
  >>> glucose.model().mass()
  168.0606



Table of Contents
-----------------

.. toctree ::

    installing
    api
    changelog
