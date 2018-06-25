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
