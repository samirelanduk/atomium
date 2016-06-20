Changelog
---------

Release 0.4.0
~~~~~~~~~~~~~

`20 June 2016`

* Secondary Structure

    * Added Alpha Helix class.

    * Added Beta Strand class.

* Residue distance matrices

    * Chains can now generate SVG distance matrices showing the distances between residues.

* Missing residues

    * Chains can now produce a combined list of all residue IDs, missing and present.


Release 0.3.0
~~~~~~~~~~~~~

`1 June 2016`

* Atom connectivity

    * Covalent bonds are now added, and atoms now know about their neighbours.

* Residue connectivity

    * Residues are now aware of which residue they are covalently bound to in their chain.

* Atomic contacts

    * Added methods for calculating the internal and external atomic contacts of any atomic structure.

* Bug fixes

    * Fixed bug where PDB files could not have site mapping parsed where there was no space between the chain ID and residue ID.


Release 0.2.0
~~~~~~~~~~~~~

`19 May 2016`

* Protein Sequences

    * Residuic Sequences can now return their amino acid sequence as a string

* Binding Sites

    * Added a class for binding sites
    * Mapped sites to ligands
    * Added methods for getting sites for ligands

* Insert codes

    * Incorporated insert codes into residue IDs


Release 0.1.0
~~~~~~~~~~~~~

`16 May 2016`

* Basic PDB parsing

  * Models
  * Chains
  * Residues
  * Atoms
  * Small Molecules
