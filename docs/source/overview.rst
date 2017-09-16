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

From .pdb
~~~~~~~~~

A .pdb can also be loaded from a file, but they can also be fetched directly
from the RCSB over the internet using the PDB code:

  >>> pdb = atomium.pdb_from_file("1LOL.pdb")
  >>> pdb2 = atomium.fetch("5HVD")
  >>> pdb2.model()
  <Model (2156 atoms)>

If the PDB has multiple models, these can be accessed using the
:py:meth:`~.Pdb.models` method.

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

Atoms can be bonded to one another using the :py:meth:`~.Atom.bond` method:

  >>> other_atom = model.atom(element="O")
  >>> atom.bond(other_atom)
  >>> atom.bonds()
  {"<C-O Bond>"}
  >>> atom.bonded_atoms()
  {<O Atom at (37.441, 29.265, 52.113)>}
  >>> atom.bond_with(other_atom)
  <C-O Bond>
  >>> atom.unbond(other_atom)
  >>> atom.bonds()
  {}
  >>> atom.bonded_atoms()
  {}


Sub-Structures
~~~~~~~~~~~~~~

Molecules
#########

PDB files contain descriptions of the various molecular units within the model.
The simplest way to access these is to get the :py:class:`.Molecule` objects in
the model:

  >>> pdb.model().molecules(water=False)
  {<Molecule A2001 (XMP, 24 atoms)>, <Molecule B5002 (BU2, 6 atoms)>, <Molecule A5
  001 (BU2, 6 atoms)>, <Chain (204 residues)>, <Molecule B2002 (XMP, 24 atoms)>, <
  Chain (214 residues)>}
  >>> pdb.model().molecules(water=False, generic=True)
  {<Molecule B2002 (XMP, 24 atoms)>, <Molecule B5002 (BU2, 6 atoms)>, <Molecule A2
  001 (XMP, 24 atoms)>, <Molecule A5001 (BU2, 6 atoms)>}

In the first case all molecules (excluding water molecules) are returned - these
include generic :py:class:`.Molecule` objects, used to represent the small
molecules in the PDB, and also :py:class:`.Chain` objects, which are the main
macromolecular unit of the PDB.

Other criteria can be used:

  >>> pdb.model().molecules(name="XMP")
  {<Molecule B2002 (XMP, 24 atoms)>, <Molecule A2001 (XMP, 24 atoms)>}
  >>> pdb.model().molecule(name="XMP")
  <Molecule B2002 (XMP, 24 atoms)>
  >>> pdb.model().molecule("B5002")
  <Molecule B5002 (BU2, 6 atoms)>

Here, all XMP molecules are returned, then the first matching XMP molecule, then
the molecule with ID 'B5002'.

Chains
######

You can specifically get chains in much the same way:

  >>> pdb.model().chains()
  {<Chain (214 residues)>, <Chain (204 residues)>}
  >>> pdb.model().chain("A")
  <Chain (204 residues)>
  >>> pdb.model().chain("B")
  <Chain (214 residues)>

A :py:class:`.Chain` is a useful object in its own right:

  >>> pdb.model().chain("A").length()
  204

Residues
########

Both models and chains are :py:class:`.ResidueStructure` objects, which allows
you to access their :py:class:`.Residue` objects:

  >>> pdb.model().residues(name="SER")
  {<Residue B1221 (SER, 6 atoms)>, <Residue B1204 (SER, 6 atoms)>, <Residue B112
  7 (SER, 6 atoms)>, <Residue A221 (SER, 6 atoms)>, <Residue A204 (SER, 6 atoms)
  >, <Residue A179 (SER, 6 atoms)>, <Residue B1165 (SER, 6 atoms)>, <Residue B11
  75 (SER, 6 atoms)>, <Residue A127 (SER, 6 atoms)>, <Residue B1050 (SER, 6 atom
  s)>, <Residue B1158 (SER, 6 atoms)>, <Residue A158 (SER, 6 atoms)>, <Residue B
  1105 (SER, 6 atoms)>, <Residue A165 (SER, 6 atoms)>, <Residue A175 (SER, 6 ato
  ms)>, <Residue A50 (SER, 6 atoms)>, <Residue B1179 (SER, 6 atoms)>, <Residue A
  105 (SER, 6 atoms)>}
  >>> pdb.model().residue("A23")
  <Residue A23 (ASN, 8 atoms)>

Residues are also a kind of Molecule, and have other useful properties:

  >>> pdb.model().residue("A23").name()
  'ASN'
  >>> pdb.model().residue("A23").chain()
  <Chain (204 residues)>
  >>> pdb.model().residue("A23").next()
  <Residue A24 (ARG, 11 atoms)>
  >>> pdb.model().residue("A23").previous()
  <Residue A22 (MET, 8 atoms)>


Saving
~~~~~~

A model can be saved to file using:

  >>> model.save("new.xyz", description="Modifed glucose")
  >>> model.save("new.pdb")

Any structure can be saved in this way, so you can save chains or molecules to
their own seperate files if you so wish.

  >>> model.chain("A").save("chainA.pdb")
  >>> model.chain("B").save("chainB.pdb")
  >>> model.molecule(name="XMP").save("ligand.xyz")

The ``Xyz`` or ``Pdb`` object itself can also be saved:

  >>> glucose.comment("Modified glucose")
  >>> glucose.save("new.xyz")
  >>> pdb.save("new.pdb")
