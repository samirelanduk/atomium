|travis| |coveralls| |pypi|

.. |travis| image:: https://api.travis-ci.org/samirelanduk/atomium.svg?branch=0.9.1
  :target: https://travis-ci.org/samirelanduk/atomium/

.. |coveralls| image:: https://coveralls.io/repos/github/samirelanduk/atomium/badge.svg?branch=0.9.1
  :target: https://coveralls.io/github/samirelanduk/atomium/

.. |pypi| image:: https://img.shields.io/pypi/pyversions/atomium.svg
  :target: https://pypi.org/project/atomium/

atomium
=======

atomium is a molecular modeller and file parser.

Example
-------

  >>> import atomium
  >>> pdb = atomium.fetch("5HVD")
  >>> pdb.model
  <Model (2156 atoms)>
  >>> pdb.model.chain("A")
  <Chain (204 residues)>





Installing
----------

pip
~~~

atomium can be installed using pip:

``$ pip3 install atomium``

atomium is written for Python 3, and does not support Python 2. It currently
requires Python 3.5 and above.

If you get permission errors, try using ``sudo``:

``$ sudo pip3 install atomium``


Development
~~~~~~~~~~~

The repository for atomium, containing the most recent iteration, can be
found `here <http://github.com/samirelanduk/atomium/>`_. To clone the
atomium repository directly from there, use:

``$ git clone git://github.com/samirelanduk/atomium.git``


Requirements
~~~~~~~~~~~~

atomium requires the Python library
`NumPy <http://www.numpy.org/>`_ (for now) and
`requests <https://docs.python-requests.org/>`_ - pip will install these
automatically when it installs atomium.


Overview
--------

atomium allows you to open .pdb and .xyz files, manipulate the model within,
and save them as new files.

From .xyz
~~~~~~~~~

You can load a .xyz file as follows:

  >>> import atomium
  >>> glucose = atomium.xyz_from_file("glucose.xyz")
  >>> glucose.title
  'glucose from 2gbp'
  >>> glucose.model
  <Model (12 atoms)>

The ``Xyz`` object you get has a ``Xyz.title`` property,
which describes the file, and a ``Xyz.model`` property, which returns
the ``Model`` the file describes.


From .pdb
~~~~~~~~~

A .pdb can also be loaded from a file, but they can also be fetched directly
from the RCSB over the internet using the PDB code:

  >>> pdb = atomium.pdb_from_file("1LOL.pdb")
  >>> pdb2 = atomium.fetch("5HVD")
  >>> pdb.deposition_date
  datetime.date(2002, 5, 6)
  >>> pdb.resolution
  1.9
  >>> pdb.rfactor
  0.193
  >>> pdb.classification
  'LYASE'
  >>> pdb.technique
  'X-RAY DIFFRACTION'
  >>> pdb2.model
  <Model (2156 atoms)>

If the PDB has multiple models, these can be accessed using the
``Pdb.models`` property. They also have ``Pdb.title``,
``Pdb.code`` and ``Pdb.deposition_date`` properties, as well
as other parsed properties - see the full documentation for details.


The Model
~~~~~~~~~

A ``Model`` is a representation of some molecular system and every
``Atom`` within it, as described by a file.

As an ``AtomicStructure`` you can query its atoms, transform it in
space, get its mass or formula, and get its centre of mass and radius of
gyration:

  >>> model = glucose.model
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
  >>> model.mass
  168.0606
  >>> model.formula
  Counter({'C': 6, 'O': 6})
  >>> model.translate(34, -12, 3.5)
  >>> model.rotate("x", 45)
  >>> model.atom(element="O")
  <O Atom at (71.441, -27.11613084494172, 51.53252799931321)>
  >>> model.center_of_mass
  (71.39909500620611, -24.411126748628675, 50.69765860848817)
  >>> model.radius_of_gyration
  2.3076405766875925

``AtomicStructure.atoms`` returns all matching elements as a ``set``
while ``AtomicStructure.atom`` returns the first matching atom.

For pairwise comparisons, structures also have the
``AtomicStructure.pairwise_atoms`` generator which will yield all
unique atom pairs in the structure. These can obviously get very big indeed - a
5000 atom PDB file would have about 12 million unique pairs.

The atoms themselves have properties for their coordinates and elements, and
also for finding the distance between them:

  >>> atom = model.atom(element="C")
  >>> atom.x, atom.y, atom.z
  (72.553, -25.00258867597513, 51.02411822364008)
  >>> atom.location
  (72.553, -25.00258867597513, 51.02411822364008)
  >>> atom.element()
  'C'
  >>> atom.distance_to(model.atom(element="O"))
  2.4417381104450953

Instead of an atom, you can also provide a coordinate and get the atom's
distance to that:

  >>> atom.distance_to(model.center_of_mass)
  1.3371237139950765

Atoms can be bonded to one another using the ``Atom.bond_to`` method:

  >>> other_atom = model.atom(element="O")
  >>> atom.bond_to(other_atom)
  >>> atom.bonds()
  {"<C-O Bond>"}
  >>> atom.bonded_atoms()
  {<O Atom at (37.441, 29.265, 52.113)>}
  >>> atom.bond_with(other_atom)
  <C-O Bond>
  >>> atom.unbond_from(other_atom)
  >>> atom.bonds
  {}
  >>> atom.bonded_atoms()
  {}


Sub-Structures
~~~~~~~~~~~~~~

Molecules
#########

PDB files contain descriptions of the various molecular units within the model.
The simplest way to access these is to get the ``Molecule`` objects in
the model:

  >>> pdb.model.molecules(water=False)
  {<Molecule A2001 (XMP, 24 atoms)>, <Molecule B5002 (BU2, 6 atoms)>, <Molecule A5
  001 (BU2, 6 atoms)>, <Chain (204 residues)>, <Molecule B2002 (XMP, 24 atoms)>, <
  Chain (214 residues)>}
  >>> pdb.model.molecules(water=False, generic=True)
  {<Molecule B2002 (XMP, 24 atoms)>, <Molecule B5002 (BU2, 6 atoms)>, <Molecule A2
  001 (XMP, 24 atoms)>, <Molecule A5001 (BU2, 6 atoms)>}

In the first case all molecules (excluding water molecules) are returned - these
include generic ``Molecule`` objects, used to represent the small
molecules in the PDB, and also ``Chain`` objects, which are the main
macromolecular unit of the PDB.

Other criteria can be used:

  >>> pdb.model.molecules(name="XMP")
  {<Molecule B2002 (XMP, 24 atoms)>, <Molecule A2001 (XMP, 24 atoms)>}
  >>> pdb.model.molecule(name="XMP")
  <Molecule B2002 (XMP, 24 atoms)>
  >>> pdb.model.molecule("B5002")
  <Molecule B5002 (BU2, 6 atoms)>

Here, all XMP molecules are returned, then the first matching XMP molecule, then
the molecule with ID 'B5002'.

Any molecule can try and determine its binding site with the
``Molecule.site`` method:

  >>> pdb.model.molecule("B5002").site()
  <'B5002' Site (8 residues)>
  >>> pdb.model.molecule("B5002").site().residues()
  {<Residue B1096 (ILE, 8 atoms)>, <Residue B1157 (PRO, 7 atoms)>, <Residue B1
  123 (LEU, 8 atoms)>, <Residue B1070 (ASP, 8 atoms)>, <Residue B1042 (LYS, 9
  atoms)>, <Residue B1072 (LYS, 9 atoms)>, <Residue B1156 (GLY, 4 atoms)>, <Re
  sidue B1155 (VAL, 7 atoms)>}

These are all the residues with a non-hydrogen atom within 4 Angstroms of a
non-hydrogen atom in the molecule. The crirteria - cutoff distance, whether to
include main chain atoms etc., can be modified using arguments. See the full
docs for details.

You can also get RMSD with, and superimpose onto, other molecules. Again the
details are in the full docs.

Chains
######

You can specifically get chains in much the same way:

  >>> pdb.model.chains()
  {<Chain (214 residues)>, <Chain (204 residues)>}
  >>> pdb.model.chain("A")
  <Chain (204 residues)>
  >>> pdb.model.chain("B")
  <Chain (214 residues)>

A ``Chain`` is a useful object in its own right:

  >>> pdb.model.chain("A").length()
  204

Residues
########

Both models and chains are made of residues objects, which allows
you to access their ``Residue`` objects:

  >>> pdb.model.residues(name="SER")
  {<Residue B1221 (SER, 6 atoms)>, <Residue B1204 (SER, 6 atoms)>, <Residue B112
  7 (SER, 6 atoms)>, <Residue A221 (SER, 6 atoms)>, <Residue A204 (SER, 6 atoms)
  >, <Residue A179 (SER, 6 atoms)>, <Residue B1165 (SER, 6 atoms)>, <Residue B11
  75 (SER, 6 atoms)>, <Residue A127 (SER, 6 atoms)>, <Residue B1050 (SER, 6 atom
  s)>, <Residue B1158 (SER, 6 atoms)>, <Residue A158 (SER, 6 atoms)>, <Residue B
  1105 (SER, 6 atoms)>, <Residue A165 (SER, 6 atoms)>, <Residue A175 (SER, 6 ato
  ms)>, <Residue A50 (SER, 6 atoms)>, <Residue B1179 (SER, 6 atoms)>, <Residue A
  105 (SER, 6 atoms)>}
  >>> pdb.model.residue("A23")
  <Residue A23 (ASN, 8 atoms)>

Residues are also a kind of Molecule, and have other useful properties:

  >>> pdb.model.residue("A23").name()
  'ASN'
  >>> pdb.model.residue("A23").chain()
  <Chain (204 residues)>
  >>> pdb.model.residue("A23").next
  <Residue A24 (ARG, 11 atoms)>
  >>> pdb.model.residue("A23").previous
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

  >>> glucose.title("Modified glucose")
  >>> glucose.save("new.xyz")
  >>> pdb.title("Modified PDB")
  >>> pdb.save("new.pdb")


Changelog
---------

Release 0.9.1
~~~~~~~~~~~~~

`17 May 2018`

* Added Residue one-letter codes.
* Fixed stray print statement.


Release 0.9.0
~~~~~~~~~~~~~

`10 April 2018`

* Turned many methods into properties.
* Added full residue name generation.
* Made bind site detection more picky.
* Added coordinate rounding to deal with floating point rounding errors.
* Atomic structures now 'copy'able.
* Refactored atom querying.
* Added grid generation.
* Implemented Kabsch superposition/rotation.
* Implemented RMSD comparison.
* Created Complex class (for later).


Release 0.8.0
~~~~~~~~~~~~~

`2 December 2017`

* Added option to get water residues in binding sites.
* Added extra PDB meta information parsing, such as:

	* Classification
	* Experimental Technique
	* Source Organism
	* Expression Organism
	* R-factor


Release 0.7.0
~~~~~~~~~~~~~

`2 November 2017`

* PDBs with multiple occupancy can now be parsed correctly.
* Added pairwise atom generator.
* PDB parser now extracts resolution.
* Further speed increased to PDB parser.
* Miscellaneous bug fixes.
* Implemented Continuous Integration.


Release 0.6.0
~~~~~~~~~~~~~

`3 October 2017`

* Now allows for fetching and opening of PDB data dictionaries.
* Added parsing/saving of HEADER and TITLE records in PDB files.
* Added ability to exclude elements from atom search.
* Added ability to get nearby atoms in a model.
* Added bind site identification.
* Fixed chain length bottleneck in PDB model saving.
* Overhauled PDB parsing by replacing classes with built in Python types.
* Fixed bug where numerical residue names were interpreted as integers.
* Changed atoms so that they can allow negative B factors.
* Added loading of .xyz data dictionaries.
* Miscellaneous speed increases.

Release 0.5.0
~~~~~~~~~~~~~

`16 September 2017`

* Added atom temperature factors.
* Added bond vector production.
* Added parse time tests and reduced parse time by over a half.
* Changed way atoms are stored in structures to make ID lookup orders of \
  magnitude faster.
* Made IDs immutable.
* Added multiple model parsing and saving.
* Added option to fetch PDBs from PDBe rather than RCSB.


Release 0.4.0
~~~~~~~~~~~~~

`26 August 2017`

* Added PDB parsing.
* Added PDB saving.
* Gave atoms ability to get specific bond with other atom.
* Added bond angle calculation.
* Added ability to filter out water molecules.

Release 0.3.0
~~~~~~~~~~~~~

`11 August 2017`

* Added classes for Molecules, Chains, Residues, and their interfaces.
* Added charges to atoms and structures.
* Add ability to create AtomicStructures from AtomicStructures.


Release 0.2.0
~~~~~~~~~~~~~

`14 June 2017`

* Made all Atomic Structures savable.
* Added Atom IDs and uniqueness constraints.
* Added Atom Bonds.


Release 0.1.1
~~~~~~~~~~~~~

`1 June 2017`

* Fixed setup.py
* Minor typos


Release 0.1.0
~~~~~~~~~~~~~

`1 June 2017`

* Added basic Model and Atom classes.
* Added .xyz parsing.
