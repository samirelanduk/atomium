Overview
--------

atomium allows you to open .xyz files, manipulate the model within, and save
them as new .xyz files.

From .xyz
~~~~~~~~~

You can load a .xyz file as follows:

  >>> import atomium
  >>> glucose = atomium.xyz_from_file("glucose.xyz")
  >>> glucose.comment()
  'glucose from 2gbp'
  >>> glucose.model()
  <Model (12 atoms)>

The :py:class:`.Xyz` object you get has a :py:meth:`~.Xyz.comment` property,
which describes the file, and a :py:meth:`~.Xyz.model` property, which returns
the :py:class:`.Model` the file describes.

The Model
~~~~~~~~~

A :py:class:`.Model` is a representation of some molecular system and every
:py:class:`.Atom` within it, as described by a file.

As an :py:class:`.AtomicStructure` you can query its atoms, transform it in
space, get its mass or formula, and get its centre of mass and radius of
gyration:

  >>> model = glucose.model()
  >>> model.atoms()
  {<C Atom at (38.553, 30.4, 50.259)>, <C Atom at (35.884, 30.895, 49.12)>, <C A
  tom at (36.177, 29.853, 50.124)>, <C Atom at (37.296, 30.296, 51.074)>, <O Ato
  m at (39.261, 32.018, 46.92)>, <C Atom at (38.357, 31.29, 49.044)>, <C Atom at
   (39.559, 31.209, 48.082)>, <O Atom at (37.441, 29.265, 52.113)>, <O Atom at (
  34.923, 29.775, 50.91)>, <O Atom at (34.968, 30.34, 48.234)>, <O Atom at (37.1
  55, 30.858, 48.364)>, <O Atom at (39.572, 30.954, 51.086)>}
  >>> model.atoms(element="O")
  {<O Atom at (37.441, 29.265, 52.113)>, <O Atom at (39.261, 32.018, 46.92)>, <O
   Atom at (37.155, 30.858, 48.364)>, <O Atom at (34.968, 30.34, 48.234)>, <O At
  om at (34.923, 29.775, 50.91)>, <O Atom at (39.572, 30.954, 51.086)>}
  >>> model.atom(element="O")
  <O Atom at (37.441, 29.265, 52.113)>
  >>> model.mass()
  168.0606
  >>> model.formula()
  Counter({'C': 6, 'O': 6})
  >>> model.translate(34, -12, 3.5)
  >>> model.rotate("x", 45)
  >>> model.atom(element="O")
  <O Atom at (71.441, -27.11613084494172, 51.53252799931321)>
  >>> model.center_of_mass()
  (71.39909500620611, -24.411126748628675, 50.69765860848817)
  >>> model.radius_of_gyration()
  2.3076405766875925

:py:meth:`~.AtomicStructure.atoms` returns all matching elements as a ``set``
while :py:meth:`~.AtomicStructure.atom` returns the first matching atom.

The atoms themselves have properties for their coordinates and elements, and
also for finding the distance between them:

  >>> atom = model.atom(element="C")
  >>> atom.x(), atom.y(), atom.z()
  (72.553, -25.00258867597513, 51.02411822364008)
  >>> atom.element()
  'C'
  >>> atom.distance_to(model.atom(element="O"))
  2.4417381104450953

Instead of an atom, you can also provide a coordinate and get the atom's
distance to that:

  >>> atom.distance_to(model.center_of_mass())
  1.3371237139950765


Saving
~~~~~~

A model can be saved to file using:

  >>> model.save("new.xyz", description="Modifed glucose")

The ``Xyz`` object itself can also be saved:

  >>> glucose.comment("Modified glucose")
  >>> glucose.save("new.xyz")
