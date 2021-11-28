atomium
=======

|travis| |coveralls| |pypi| |version| |commit| |downloads|

.. |travis| image:: https://api.travis-ci.org/samirelanduk/atomium.svg?branch=master
  :target: https://travis-ci.org/samirelanduk/atomium/

.. |coveralls| image:: https://coveralls.io/repos/github/samirelanduk/atomium/badge.svg?branch=master
  :target: https://coveralls.io/github/samirelanduk/atomium/

.. |pypi| image:: https://img.shields.io/pypi/pyversions/atomium.svg
  :target: https://pypi.org/project/atomium/

.. |version| image:: https://img.shields.io/pypi/v/atomium.svg
  :target: https://pypi.org/project/atomium/

.. |commit| image:: https://img.shields.io/github/last-commit/samirelanduk/atomium/master.svg
  :target: https://github.com/samirelanduk/atomium/tree/master/

.. |downloads| image:: https://img.shields.io/pypi/dm/atomium.svg
  :target: https://pypi.org/project/atomium/


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



Installing
----------

pip
~~~

atomium can be installed using pip:

``$ pip3 install atomium``

atomium is written for Python 3, and does not support Python 2.

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

atomium requires `requests <http://docs.python-requests.org/>`_ for fetching
structures from the RCSB, `paramiko <http://www.paramiko.org//>`_ for
fetching structures over SSH,
`msgpack <https://github.com/msgpack/msgpack-python>`_ for parsing .mmtf files,
and `valerius <https://valerius.samireland.com>`_ for dealing with sequences.


Testing
~~~~~~~

To test a local version of atomium, cd to the atomium directory and run:

``$ python -m unittest discover tests``

You can opt to only run unit tests or integration tests:

``$ python -m unittest discover tests.unit``
``$ python -m unittest discover tests.integration``

You can run the 'big test' to get a random 1000 structures, parse them all, and
report any problems:

``$ python tests/big.py``

Finally, to perform speed profiles you can run:

``$ python tests/time/time.py``

...which creates various profiles that SnakeViz can visualise.



Overview
--------

atomium is a Python library for opening and saving .pdb, .cif and .mmtf files,
and presenting and manipulating the information contained within.


Loading Data
~~~~~~~~~~~~

While you can use atomium to create models from scratch to build an entirely
*de novo* structure, in practice you would generally use it to load molecular
data from an existing file...

	>>> import atomium
	>>> pdb1 = atomium.open('../1LOL.pdb')
	>>> mmtf1 = atomium.open('/structures/glucose.mmtf')
	>>> cif1 = atomium.open('/structures/1XDA.cif')
	>>> pdb3 = atomium.open('./5CPA.pdb.gz')
	>>> pdb2 = atomium.fetch('5XME.pdb')
	>>> cif2 = atomium.fetch('5XME')

In that latter case, you don't need the file to be saved locally - it will just
go and grab the PDB with that code from the RCSB.

atomium will use the file extension you provide to decide how to parse it. If
there isn't one, or it doesn't recognise the extension, it will peek at the
file contents and try and guess whether it should be interpreted as .pdb, .cif
or .mmtf.


Using Data
~~~~~~~~~~

Once you've got your ``File`` object, what can you do with it?

Annotation
##########

There is meta information contained within the ``File`` object:

    >>> pdb1.title
    'CRYSTAL STRUCTURE OF OROTIDINE MONOPHOSPHATE DECARBOXYLASE COMPLEX WITH XMP'
    >>> pdb1.deposition_date
    datetime.date(2002, 5, 6)
    >>> pdb1.keywords
    ['TIM BARREL', 'LYASE']
    >>> pdb1.classification
    'LYASE'
    >>> pdb1.source_organism
    'METHANOTHERMOBACTER THERMAUTOTROPHICUS STR. DELTA H'
    >>> pdb1.resolution
    1.9
    >>> pdb1.rvalue
    0.193
    >>> pdb1.rfree
    0.229

atomium doesn't currently parse *every* bit of information from these
files, but there is more than those shown above. See
`the full API docs <api/pdb.html>`_ for more details. In particular, you can
access the processed intermediate MMCIF dictionary to get *any* attribute of
these structures.

Models and Assembly
###################

All .pdb files contain one or more models - little universes containing a
molecular scene.

    >>> pdb1.model
    <Model (2 chains, 4 ligands)>
    >>> pdb1.models
    (<Model (2 chains, 4 ligands)>,)

Most just contain one - it's generally those that come from NMR experiments
which contain multiple models. You can easily iterate through these to get their
individual metrics:

    >>> for model in pdb2.models:
            print(model.center_of_mass)

This model contains the 'asymmetric unit' - this is one or more protein
(usually) chains arranged in space, which may not be how the molecule arranges
itself in real life. It might just be how they arranged themselves in the
experiment. To create the 'real thing' from the asymmetric unit, you use
**biological assemblies.**

Most .pdb files contain one or more biological assemblies - instructions for how
to create a more realistic structure from the chains present, which in atomium
are accessed using ``File.assemblies``.

In practice, what you need to know is that you can create a new model - not the
one already there containing the asymmetric unit - as follows...

    >>> pdb3 = atomium.fetch('1XDA')
    >>> pdb3.model
    <Model (8 chains, 16 ligands)>
    >>> pdb3.generate_assembly(1)
    <Model (2 chains, 4 ligands)>
    >>> pdb3.generate_assembly(10)
    <Model (6 chains, 12 ligands)>
    >>> [pdb.generate_assembly(n + 1) for n in range(len(pdb.assemblies))]
    [<Model (2 chains, 4 ligands)>, <Model (2 chains, 4 ligands)>, <Model (2 cha
    ins, 4 ligands)>, <Model (2 chains, 4 ligands)>, <Model (12 chains, 24 ligan
    ds)>, <Model (12 chains, 24 ligands)>, <Model (6 chains, 12 ligands)>, <Mode
    l (6 chains, 12 ligands)>, <Model (6 chains, 12 ligands)>, <Model (6 chains,
     12 ligands)>, <Model (4 chains, 8 ligands)>, <Model (4 chains, 8 ligands)>]

Here you load a .pdb with multiple possible assemblies, have a quick look at
the asymmetric unit with 1,842 atoms, and then generate first , and then all,
of its possible biological assemblies by passing in their IDs.


Model Contents
##############

The basic structures within a model are chains, residues, ligands, and atoms.

    >>> pdb1.model.chains()
    {<Chain A (204 residues)>, <Chain B (214 residues)>}
    >>> pdb1.model.chain('B')
    <Chain B (214 residues)>
    >>> pdb1.model.residues(name='TYR')
    {<Residue TYR (A.37)>, <Residue TYR (B.1037)>, <Residue TYR (A.45)>, <Residu
    e TYR (A.154)>, <Residue TYR (B.1206)>, <Residue TYR (B.1154)>, <Residue TYR
     (B.1045)>, <Residue TYR (A.206)>}
    >>> pdb1.model.residues(name__regex='TYR|PRO')
    {<Residue PRO (A.101)>, <Residue PRO (A.46)>, <Residue PRO (A.161)>, <Residu
    e TYR (A.45)>, <Residue PRO (B.1046)>, <Residue TYR (A.154)>, <Residue TYR (
    B.1206)>, <Residue TYR (B.1045)>, <Residue PRO (B.1189)>, <Residue TYR (A.37
    )>, <Residue PRO (B.1129)>, <Residue PRO (B.1077)>, <Residue PRO (A.211)>, <
    Residue PRO (B.1180)>, <Residue PRO (B.1157)>, <Residue PRO (B.1211)>, <Resi
    due PRO (B.1228)>, <Residue PRO (B.1101)>, <Residue TYR (B.1154)>, <Residue
    PRO (A.157)>, <Residue PRO (A.77)>, <Residue PRO (A.180)>, <Residue TYR (B.1
    037)>, <Residue PRO (A.129)>, <Residue PRO (B.1161)>, <Residue TYR (A.206)>}
    >>> pdb1.model.chain('B').residue('B.1206')
    <Residue TYR (B.1206)>
    >>> pdb1.model.chain('B').residue('B.1206').helix
    True
    >>> pdb1.model.ligands()
    {<Ligand BU2 (A.5001)>, <Ligand XMP (A.2001)>, <Ligand BU2 (B.5002)>, <Ligan
    d XMP (B.2002)>}
    >>> pdb1.model.ligand(name='BU2').atoms()
    {<Atom 3196 (O3)>, <Atom 3192 (C1)>, <Atom 3193 (O1)>, <Atom 3197 (C4)>, <At
    om 3194 (C2)>, <Atom 3195 (C3)>}
    >>> pdb1.model.ligand(name='BU2').atoms(mass__gt=12)
    {<Atom 3196 (O3)>, <Atom 3192 (C1)>, <Atom 3193 (O1)>, <Atom 3197 (C4)>, <At
    om 3194 (C2)>, <Atom 3195 (C3)>}
    >>> pdb1.model.ligand(name='BU2').atoms(mass__gt=14)
    {<Atom 3196 (O3)>, <Atom 3193 (O1)>}

The examples above demonstrate atomium's selection language. In the case of the
molecules - ``Model``, ``Chain``, ``Residue`` and
``Ligand`` - you can pass in an ``id`` or ``name``, or search by regex
pattern with ``id__regex`` or ``name__regex``.

These structures have an even more powerful syntax too - you can pass in *any*
property such as ``charge=1``, any comparitor of a property such as
``mass__lt=100``, or any regex of a property such as ``name__regex='[^C]'``.

For pairwise comparisons, structures also have the
``AtomStructure.pairwise_atoms`` generator which will yield all
unique atom pairs in the structure. These can obviously get very big indeed - a
5000 atom PDB file would have about 12 million unique pairs.

Structures can be moved around and otherwise compared with each other...

    >>> pdb1.model.ligand(id='B:2002').mass
    351.1022
    >>> pdb1.model.ligand(id='B.2002').formula
    Counter({'C': 10, 'O': 9, 'N': 4, 'P': 1})
    >>> pdb1.model.ligand(id='B:2002').nearby_atoms(2.8)
    {<Atom 3416 (O)>, <Atom 3375 (O)>, <Atom 1635 (OD1)>}
    >>> pdb1.model.ligand(id='B.2002').nearby_atoms(2.8, name='OD1')
    {<Atom 1635 (OD1)>}
    >>> pdb1.model.ligand(id='B.2002').nearby_residues(2.8)
    {<Residue ASP (B.1020)>}
    >>> pdb1.model.ligand(id='B.2002').nearby_structures(2.8, waters=True)
    {<Residue ASP (B.1020)>, <Water HOH (B.3155)>, <Water HOH (B.3059)>}
    >>> import math
    >>> pdb1.model.ligand(id='B.2002').rotate(math.pi / 2, 'x')
    >>> pdb1.model.ligand(id='B.2002').translate(10, 10, 15)
    >>> pdb1.model.ligand(id='B.2002').center_of_mass
    (-9.886734282781484, -42.558415679537184, 77.33400578435568)
    >>> pdb1.model.ligand(id='B.2002').radius_of_gyration
    3.6633506511540825
    >>> pdb1.model.ligand(id='B.2002').rmsd_with(pdb1.model.ligand(id='A.2001'))
    0.133255572356

Here we look at one of the ligands, identify its mass and molecular formula,
look at what atoms are within 2.8 Angstroms of it, and what residues are within
that same distance, rotate it and translate it through space, see where its new
center of mass is, and then finally get its RMSD with the other similar ligand
in the model.

Any operation which involves identifying nearby structures or atoms can be sped
up - dramatically in the case of very large structures - by calling
``Model.optimise_distances`` on the ``Model`` first. This
prevents atomium from having to compare every atom with every other atom every
time a proximity check is made.

The ``Atom`` objects themselves have their own useful properties.

    >>> pdb1.model.atom(97)
    <Atom 97 (CA)>
    >>> pdb1.model.atom(97).mass
    12.0107
    >>> pdb1.model.atom(97).anisotropy
    [0, 0, 0, 0, 0, 0]
    >>> pdb1.model.atom(97).bvalue
    24.87
    >>> pdb1.model.atom(97).location
    (-12.739, 31.201, 43.016)
    >>> pdb1.model.atom(97).distance_to(pdb1.model.atom(1))
    26.18289982030257
    >>> pdb1.model.atom(97).nearby_atoms(2)
    {<Atom 96 (N)>, <Atom 98 (C)>, <Atom 100 (CB)>}
    >>> pdb1.model.atom(97).is_metal
    False
    >>> pdb1.model.atom(97).structure
    <Residue ASN (A.23)>
    >>> pdb1.model.atom(97).chain
    <Chain A (204 residues)>

Chains are a bit different from other structures in that they are iterable,
indexable, and return their residues as a tuple, not a set...

    >>> pdb1.model.atom(97).chain
    <Chain A (204 residues)>
    >>> pdb1.model.chain('A')
    <Chain A (204 residues)>
    >>> len(pdb1.model.chain('A'))
    204
    >>> pdb1.model.chain('A')[10]
    <Residue LEU (A.21)>
    >>> pdb1.model.chain('A').residues()[:5]
    (<Residue VAL (A.11)>, <Residue MET (A.12)>, <Residue ASN (A.13)>, <Residue
    ARG (A.14)>, <Residue LEU (A.15)>)
    >>> pdb1.model.chain('A').sequence
    'LRSRRVDVMDVMNRLILAMDLMNRDDALRVTGEVREYIDTVKIGYPLVLSEGMDIIAEFRKRFGCRIIADFKVAD
    IPETNEKICRATFKAGADAIIVHGFPGADSVRACLNVAEEMGREVFLLTEMSHPGAEMFIQGAADEIARMGVDLGV
    KNYVGPSTRPERLSRLREIIGQDSFLISPGVGAQGGDPGETLRFADAIIVGRSIYLADNPAAAAAGIIESIKDLLI
    PE'

The sequence is
the 'real' sequence that exists in nature. Some of them will be
missing from the model for practical reasons.

Residues can generate name information based on their three letter code, and are
aware of their immediate neighbors.

    >>> pdb1.model.residue('A.100')
    <Residue PHE (A.100)>
    >>> pdb1.model.residue('A.100').name
    'PHE'
    >>> pdb1.model.residue('A.100').code
    'F'
    >>> pdb1.model.residue('A.100').full_name
    'phenylalanine'
    >>> pdb1.model.residue('A.100').next
    <Residue PRO (A.101)>
    >>> pdb1.model.residue('A.100').previous
    <Residue GLY (A.99)>

Saving Data
~~~~~~~~~~~

A model can be saved to file using:

  >>> model.save("new.cif")
  >>> model.save("new.pdb")

Any structure can be saved in this way, so you can save chains or molecules to
their own seperate files if you so wish.


  >>> model.chain("A").save("chainA.pdb")
  >>> model.chain("B").save("chainB.cif")
  >>> model.ligand(name="XMP").save("ligand.mmtf")

Note that if the model you are saving is one from a biological assembly, it will
likely have many duplicated IDs, so saving to file may create unexpected
results.


Changelog
---------

Release 1.0.11
~~~~~~~~~~~~~~

`27 November 2021`

* Optimised distance lookup for finding atoms within sphere.


Release 1.0.10
~~~~~~~~~~~~~~

`29 May 2021`

* Fixed secondary structure parsing for multi character asym IDs in mmCIF.


Release 1.0.9
~~~~~~~~~~~~~

`4 February 2021`

* Fixed temperature factor zero-padding in PDB saving.
* Fixed MMTF decode bug in Ubuntu.


Release 1.0.8
~~~~~~~~~~~~~

`9 December 2020`

* HETATM identity now preserved when parsing PDB files


Release 1.0.7
~~~~~~~~~~~~~

`5 November 2020`

* Fixed blank ANISOU values in PDB saving.
* Fixed negative residue IDs in PDB saving.
* Fixed SyntaxWarning messages on PDB saving.

Release 1.0.6
~~~~~~~~~~~~~

`8 September 2020`

* Added handling of new branched entities in MMCIF/MMTF.

Release 1.0.5
~~~~~~~~~~~~~

`21 July 2020`

* Added ability to open compressed .gz files.


Release 1.0.4
~~~~~~~~~~~~~

`1 May 2020`

* Made TER records more compliant in saved PDB files.
* Specified required msgpack version to fix MMTF parsing issue.


Release 1.0.3
~~~~~~~~~~~~~

`5 December 2019`

* Made quality information detection more broad.
* Improved documentation. 


Release 1.0.2
~~~~~~~~~~~~~

`1 October 2019`

* Added distance optimiser for proximity checks.
* Improved test coverage.


Release 1.0.1
~~~~~~~~~~~~~

`26 September 2019`

* Added a pdb2json script for converting local structure files to JSON.
* Improved speed comparison checks.


Release 1.0.0
~~~~~~~~~~~~~

`23 June 2019`

* Saving now issues warning if the stucture has duplicate IDs.
* Missing residues parsed for all three file types.
* Crystallographic information now parsed.
* Refactor of atomic structures.
* Refactor of .mmtf parsing.
* Structure copying now retains all properties.
* Fixed bug in parsing .cif expression systems.
* Full names of ligands and modified residues now parsed.
* Secondary structure information parsed and available now.
* Atoms now have covalent radius property for calculating bond cutoffs.
* .pdb parsing can now handle heavy water (DOD).
* General speed improvements.


Release 0.12.2
~~~~~~~~~~~~~~

`4 February 2019`

* Angle between superimposed atoms now possible.
* Fixed source speices lookup in .cif files.
* Fixed bug relating to embedded quotes in .cif files.


Release 0.12.1
~~~~~~~~~~~~~~

`13 January 2019`

* Fixed assembly parsing bug in small number of .cif files.


Release 0.12.0
~~~~~~~~~~~~~~

`2 January 2019`

* Refactored parse utilities to improve speed.
* Added support for .mmtf files.
* Added file writing for all three file types (.pdb, .cif, .mmtf).
* Made .cif the default file type.
* General library restructuring.


Release 0.11.1
~~~~~~~~~~~~~~

`13 September 2018`

* Fixed bug pertaining to residues with ID 0.
* Fixed bug pertaining to SEQRES parsing when chain ID is numeric.
* Changed format of residue IDs to include colon.
* Considerable speed improvements in .mmcif parsing.


Release 0.11.0
~~~~~~~~~~~~~~

`22 August 2018`

* Added .mmcif parsing.
* Changed how parsing in general is done under the hood.
* Added atom angle calculation.
* Fixed bug where modified residues were treated as ligands if authors used HETATM records.


Release 0.10.2
~~~~~~~~~~~~~~

`29 July 2018`

* Added function for getting PDBs over SSH.
* Fixed biological assembly parsing bug.
* Fixed chain copying of sequence bug.


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
