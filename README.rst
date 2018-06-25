|travis| |coveralls| |pypi|

.. |travis| image:: https://api.travis-ci.org/samirelanduk/atomium.svg?branch=0.10.1
  :target: https://travis-ci.org/samirelanduk/atomium/

.. |coveralls| image:: https://coveralls.io/repos/github/samirelanduk/atomium/badge.svg?branch=0.10.1
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

atomium is a Python library for opening and saving .pdb and .xyz files, and
presenting and manipulating the information contained within.

Loading Data
~~~~~~~~~~~~

Whist you can use atomium to create models from scratch to build an entirely
*de novo* structure, in practice you would generally use it to load molecular
data from an existing file...

	>>> import atomium
	>>> pdb1 = atomium.pdb_from_file('../1LOL.pdb')
	>>> xyz1 = atomium.xyz_from_file('/structures/glucose.xyz')
	>>> pdb2 = atomium.fetch('5HVD')

In that latter case, you don't need the file to be saved locally - it will just
go and grab the PDB with that code from the RCSB.

The rest of this guide will focus on .pdb files. .xyz files are very simple
structures, and the only annotation they really contain is a ``.title``.

Using Data
~~~~~~~~~~

Once you've got your ``Pdb`` object, what can you do with it?

Annotation
##########

There is various meta information contained within the ``Pdb`` object.

    >>> pdb1.title
    'CRYSTAL STRUCTURE OF OROTIDINE MONOPHOSPHATE DECARBOXYLASE COMPLEX WITH XMP'
    >>> pdb1.deposition_date
    datetime.date(2002, 5, 6)
    >>> pdb1.keywords
    ['TIM BARREL', 'LYASE']
    >>> pdb1.classification
    'LYASE'
    >>> pdb1.organism
    'METHANOTHERMOBACTER THERMAUTOTROPHICUS STR. DELTA H'
    >>> pdb1.resolution
    1.9
    >>> pdb1.rfactor
    0.193
    >>> pdb1.rfree
    0.229

atomium doesn't currently parse *every* bit of information from .pdb files, but
there is more than those shown above. See `the full API docs <api/pdb.html>`_
for more details.

Models and Assembly
###################

All .pdb files contain one or more models - little universes containing a
molecular scene.

    >>> pdb1.model
    <Model (3431 atoms)>
    >>> pdb1.models
    (<Model (3431 atoms)>,)

Most just contain one - it's generally those that come from NMR experiments
which contain multiple models.

This model contains the 'asymmetric unit' - this is one or more protein
(usually) chains arranged in space, which may not be how the molecule arranges
itself in real life. It might just be how they arranged themselves in the
experiment. To create the 'real thing' from the asymmetric unit, you use
**biological assemblies.**

Most .pdb files contain one or more biological assemblies - instructions for how
to create a more realistic structure from the chains present, which in atomium
are accessed using ``Pdb.biomolecules``.

In practice, what you need to know is that you can create a new model - not the
one already there containing the asymmetric unit - as follows...

    >>> pdb3 = atomium.fetch('1XDA')
    >>> pdb3.model
    <Model (1842 atoms)>
    >>> pdb3.generate_assembly(1)
    <Model (924 atoms)>
    >>> pdb3.generate_assembly(10)
    <Model (2730 atoms)>
    >>> pdb3.generate_best_assembly()
    <Model (5550 atoms)>

Here you load a .pdb with multiple possible assemblies, have a quick look at
the asymmetric unit with 1,842 atoms, generate two of its possible biological
assemblies by passing in their IDs, and then generate the 'best' of the
assemblies, which is the one with the lowest (that is, most negative) delta
free energy change as described in the .pdb file. In this case it is a hexameric
formation.

Model Contents
##############

The basic structures within a model are chains, residues, ligands, and atoms.

    >>> pdb1.model.chains()
    {<Chain (B, 1748 atoms)>, <Chain (A, 1683 atoms)>}
    >>> pdb1.model.chain('B')
    <Chain (B, 1748 atoms)>
    >>> pdb1.model.residues(name='TYR')
    {<Residue TYR (A206, 12 atoms)>, <Residue TYR (A45, 12 atoms)>, <Residue TYR
     (A37, 12 atoms)>, <Residue TYR (B1154, 12 atoms)>, <Residue TYR (B1206, 12
    atoms)>, <Residue TYR (A154, 12 atoms)>, <Residue TYR (B1045, 12 atoms)>, <R
    esidue TYR (B1037, 12 atoms)>}
    >>> pdb1.model.residues(name_regex='TYR|PRO')
    {<Residue PRO (B1046, 7 atoms)>, <Residue TYR (A37, 12 atoms)>, <Residue PRO
     (A157, 7 atoms)>, <Residue TYR (B1206, 12 atoms)>, <Residue PRO (B1228, 7 a
    toms)>, <Residue PRO (A211, 7 atoms)>, <Residue PRO (B1077, 7 atoms)>, <Resi
    due PRO (B1129, 7 atoms)>, <Residue TYR (A45, 12 atoms)>, <Residue TYR (A154
    , 12 atoms)>, <Residue PRO (A180, 7 atoms)>, <Residue PRO (B1157, 7 atoms)>,
    <Residue TYR (B1037, 12 atoms)>, <Residue TYR (A206, 12 atoms)>, <Residue PR
    O (B1189, 7 atoms)>, <Residue PRO (A161, 7 atoms)>, <Residue PRO (A101, 7 at
    oms)>, <Residue PRO (A46, 7 atoms)>, <Residue TYR (B1045, 12 atoms)>, <Resid
    ue PRO (A77, 7 atoms)>, <Residue PRO (A129, 7 atoms)>, <Residue PRO (B1211,
    7 atoms)>, <Residue TYR (B1154, 12 atoms)>, <Residue PRO (B1180, 7 atoms)>,
    <Residue PRO (B1101, 7 atoms)>, <Residue PRO (B1161, 7 atoms)>}
    >>> pdb1.model.chain('B').residue('B1206')
    <Residue TYR (B1206, 12 atoms)>
    >>> pdb1.model.ligands(water=False)
    {<Ligand XMP (B2002, 24 atoms)>, <Ligand BU2 (A5001, 6 atoms)>, <Ligand XMP
    (A2001, 24 atoms)>, <Ligand BU2 (B5002, 6 atoms)>}
    >>> pdb1.model.ligand(name='BU2').atoms()
    {<C (C1) Atom 3194 at (2.646, 45.112, 48.995)>, <C (C4) Atom 3199 at (-0.456
    , 44.629, 51.162)>, <C (C3) Atom 3197 at (0.706, 44.197, 50.309)>, <O (O3) A
    tom 3198 at (1.101, 42.889, 50.701)>, <O (O1) Atom 3195 at (1.781, 45.484, 4
    7.929)>, <C (C2) Atom 3196 at (1.922, 45.088, 50.288)>}
    >>> pdb1.model.ligand(name='BU2').atoms(mass__gt=12)
    {<C (C4) Atom 3199 at (-0.456, 44.629, 51.162)>, <O (O3) Atom 3198 at (1.101
    , 42.889, 50.701)>, <C (C2) Atom 3196 at (1.922, 45.088, 50.288)>, <C (C1) A
    tom 3194 at (2.646, 45.112, 48.995)>, <C (C3) Atom 3197 at (0.706, 44.197, 5
    0.309)>, <O (O1) Atom 3195 at (1.781, 45.484, 47.929)>}
    >>> pdb1.model.ligand(name='BU2').atoms(mass__gt=14)
    {<O (O3) Atom 3198 at (1.101, 42.889, 50.701)>, <O (O1) Atom 3195 at (1.781,
     45.484, 47.929)>}

The examples above demonstrate atomium's selection language. In the case of the
molecules - ``Model``, ``Chain``, ``Residue`` and
``Ligand`` - you can pass in an ``id`` or ``name``, or search by regex
pattern with ``id_regex`` or ``name_regex``.

Atoms have an even more powerful syntax. You can pass in *any* property of atoms
such as ``charge=1``, any comparitor of a property such as ``mass__lt=100``, or
any regex of a property such as ``name_regex='[^C]'``.

For pairwise comparisons, structures also have the
``AtomStructure.pairwise_atoms`` generator which will yield all
unique atom pairs in the structure. These can obviously get very big indeed - a
5000 atom PDB file would have about 12 million unique pairs.

Structures can be moved around and otherwise compared with each other...

    >>> pdb1.model.ligand(id='B2002').mass
    351.1022
    >>> pdb1.model.ligand(id='B2002').formula
    Counter({'C': 10, 'O': 9, 'N': 4, 'P': 1})
    >>> pdb1.model.ligand(id='B2002').nearby_atoms(2.8)
    {<O (O) Atom 3377 at (-24.077, 59.423, 53.919)>, <O (O) Atom 3418 at (-14.53
    5, 62.938, 57.757)>, <O (OD1) Atom 1636 at (-22.92, 57.72, 52.315)>}
    >>> pdb1.model.ligand(id='B2002').nearby_atoms(2.8, name='OD1')
    {<O (OD1) Atom 1636 at (-22.92, 57.72, 52.315)>}
    >>> pdb1.model.ligand(id='B2002').nearby_residues(2.8)
    {<Residue ASP (B1020, 8 atoms)>}
    >>> pdb1.model.ligand(id='B2002').nearby_residues(2.8, ligands=True)
    {<Ligand HOH (B3155, 1 atom)>, <Ligand HOH (B3059, 1 atom)>, <Residue ASP (B
    1020, 8 atoms)>}
    >>> import math
    >>> pdb1.model.ligand(id='B2002').rotate(math.pi / 2, 'x')
    >>> pdb1.model.ligand(id='B2002').translate(10, 10, 15)
    >>> pdb1.model.ligand(id='B2002').center_of_mass
    (-9.886734282781484, -42.558415679537184, 77.33400578435568)
    >>> pdb1.model.ligand(id='B2002').radius_of_gyration
    3.6633506511540825
    >>> pdb1.model.ligand(id='B2002').rmsd_with(pdb1.model.ligand(id='A2001'))
    90.55588214099254
    >>> other_ligand = pdb1.model.ligand(id='A2001')
    >>> pdb1.model.ligand(id='B2002').rmsd_with(other_ligand)
    90.55588214099254
    >>> pdb1.model.ligand(id='B2002').rmsd_with(other_ligand, superimpose=True)
    0.13325557235580035

Here we look at one of the ligands, identify its mass and molecular formula,
look at what atoms are within 2.8 Angstroms of it, and what residues are within
that same distance, rotate it and translate it through space, see where its new
center of mass is, and then finally get its RMSD with the other similar ligand
in the model - first using their locations 'as is', and then by seeing what the
RMSD would be if they were superimposed in such a way as to minimise RMSD.

The ``Atom`` objects themselves have their own useful properties.

    >>> pdb1.model.atom(97)
    <C (CA) Atom 97 at (-12.739, 31.201, 43.016)>
    >>> pdb1.model.atom(97).mass
    12.0107
    >>> pdb1.model.atom(97).anisotropy
    [0, 0, 0, 0, 0, 0]
    >>> pdb1.model.atom(97).bfactor
    24.87
    >>> pdb1.model.atom(97).location
    (-12.739, 31.201, 43.016)
    >>> pdb1.model.atom(97).distance_to(pdb1.model.atom(1))
    26.18289982030257
    >>> pdb1.model.atom(97).bonded_atoms
    {<N (N) Atom 96 at (-11.649, 32.148, 42.889)>, <C (C) Atom 98 at (-12.515, 3
    0.319, 44.247)>, <C (CB) Atom 100 at (-12.897, 30.387, 41.732)>}
    >>> pdb1.model.atom(97).nearby_atoms(2)
    {<N (N) Atom 96 at (-11.649, 32.148, 42.889)>, <C (C) Atom 98 at (-12.515, 3
    0.319, 44.247)>, <C (CB) Atom 100 at (-12.897, 30.387, 41.732)>}
    >>> pdb1.model.atom(97).is_metal
    False
    >>> pdb1.model.atom(97).residue
    <Residue ASN (A23, 8 atoms)>
    >>> pdb1.model.atom(97).chain
    <Chain (A, 1683 atoms)>

Chains are a bit different from other structures in that they are iterable,
indexable, and return their residues as a tuple, not a set...

    >>> pdb1.model.atom(97).chain
    <Chain (A, 1683 atoms)>
    >>> pdb1.model.chain('A')
    <Chain (A, 1683 atoms)>
    >>> len(pdb1.model.chain('A'))
    204
    >>> pdb1.model.chain('A')[10]
    <Residue LEU (A21, 8 atoms)>
    >>> pdb1.model.chain('A').residues()[:5]
    (<Residue VAL (A11, 7 atoms)>, <Residue MET (A12, 8 atoms)>, <Residue ASN (A
    13, 8 atoms)>, <Residue ARG (A14, 11 atoms)>, <Residue LEU (A15, 8 atoms)>)
    >>> pdb1.model.chain('A').sequence
    'VMNRLILAMDLMNRDDALRVTGEVREYIDTVKIGYPLVLSEGMDIIAEFRKRFGCRIIADFKVADIPETNEKICR
    ATFKAGADAIIVHGFPGADSVRACLNVAEEMGREVFLLTEMSHPGAEMFIQGAADEIARMGVDLGVKNYVGPSTRP
    ERLSRLREIIGQDSFLISPGGETLRFADAIIVGRSIYLADNPAAAAAGIIESI'
    >>> pdb1.model.chain('A').rep_sequence
    'LRSRRVDVMDVMNRLILAMDLMNRDDALRVTGEVREYIDTVKIGYPLVLSEGMDIIAEFRKRFGCRIIADFKVAD
    IPETNEKICRATFKAGADAIIVHGFPGADSVRACLNVAEEMGREVFLLTEMSHPGAEMFIQGAADEIARMGVDLGV
    KNYVGPSTRPERLSRLREIIGQDSFLISPGVGAQGGDPGETLRFADAIIVGRSIYLADNPAAAAAGIIESIKDLLI
    PE'

In those latter two cases, two different sequences are returned. The first just
returns the sequence of residues actually present in the model, whereas the
second is the 'real' sequence that exists in nature. Some of them will be
missing from the model for practical reasons.

Residues can generate name information based on their three letter code, and are
aware of their immediate neighbors.

    >>> pdb1.model.residue('A100')
    <Residue PHE (A100, 11 atoms)>
    >>> pdb1.model.residue('A100').name
    'PHE'
    >>> pdb1.model.residue('A100').code
    'F'
    >>> pdb1.model.residue('A100').full_name
    'phenylalanine'
    >>> pdb1.model.residue('A100').next
    <Residue PRO (A101, 7 atoms)>
    >>> pdb1.model.residue('A100').previous
    <Residue GLY (A99, 4 atoms)>

Saving Data
~~~~~~~~~~~

A model can be saved to file using:

  >>> model.save("new.xyz", description="Modifed glucose")
  >>> model.save("new.pdb")

Any structure can be saved in this way, so you can save chains or molecules to
their own seperate files if you so wish.

  >>> model.chain("A").save("chainA.pdb")
  >>> model.chain("B").save("chainB.pdb")
  >>> model.ligand(name="XMP").save("ligand.xyz")

The ``Pdb`` or ``Xyz`` object itself can also be saved:

  >>> pdb.title = "Modified PDB"
  >>> pdb.save("new.pdb")

Note that if the model you are saving is one from a biological assembly, it will
likely have many duplicated IDs, so saving to file may create unexpected
results.


Changelog
---------

Release 0.10.1
~~~~~~~~~~~~~~

`25 June 2018`

* Added function for returning best biological assembly.
* Fixed bug with sorting None energy assemblies.
* Fixed bug pertaining to excessive atom duplication when creating assembly.


Release 0.10.0
~~~~~~~~~~~~~~

`22 June 2018`

* Parsing of .pdb keywords.
* Parsing of atom anisotropy.
* Parsing of .pdb sequence information.
* More R-factor information.
* Biological assembly parsing and generation.
* More powerful transformations rather than just simple rotation.
* Backend simplifications.
* Powerful new atom querying syntax.


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
