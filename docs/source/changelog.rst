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
