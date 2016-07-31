Overview
--------

Creating Pdb objects
~~~~~~~~~~~~~~~~~~~~

There are two main ways to create a Pdb object from a PDB file. The first is
from a local PDB file:

    >>> import molecupy
    >>> pdb = molecupy.get_pdb_from_file("path/to/file.pdb")

This is the quickest way, though it is not always convenient to store PDB files
locally. The second way is to fetch the PDB file over the internet:

    >>> import molecupy
    >>> pdb = molecup.get_pdb_remotely("1LOL")

This takes longer, but it means you can access any published PDB file without
needing to manually download them first.

However the text of the PDB file is obtained, the process of parsing it is
always the same:

    1. First a :py:class:`.PdbFile` object is created, which is a
    representation of the file itself. This is essentially a list of records,
    with methods for getting records of a certain name.

    2. This is used to make a :py:class:`.PdbDataFile` object. This is the
    object which extracts the data from the file, and is essentially an
    unstructured list of values.

    3. This is used to make a :py:class:`.Pdb` object, by using the values in
    the data file to create a user-friendly handle to the information. This is
    the object returned by the above two methods.


Accessing Pdb properties
~~~~~~~~~~~~~~~~~~~~~~~~

Aside from structural information, PDB files also contain many other pieces of
information about the file, such as its title, experimental techniques used to
create it, publication information etc.

    >>> pdb.pdb_code()
    '1LOL'
    >>> pdb.deposition_date()
    datetime.date(2002, 5, 6)
    >>> pdb.authors()
    ['N.WU', 'E.F.PAI']

molecuPy is a reasonably forgiving parser. If records are missing from the PDB
file - even records which the PDB specification insists *must* be present - the
file will still parse, and any missing properties will just be set to ``None``
or an empty list, whichever is appropriate.


Pdb Models
~~~~~~~~~~

The heart of a Pdb is its model. A :py:class:`.PdbModel` represents the
structure contained in that PDB file, and is the environment in which all other
molecules and structures are based.

All Pdb objects have a list of models, which in most cases will contain a single
model. Structures created from NMR will often have multiple models - each
containing the same structures but with slightly different coordinates. For ease
of use, all Pdb objects also have a ``model`` method, which points to the
first model in the list.

    >>> pdb.models()
    [<Model (3431 atoms)>]
    >>> pdb.model()
    <Model (3431 atoms)>

The PdbModel class is an atomic structure (i.e. it inherits from
:py:class:`.AtomicStructure`) which means you can get certain atomic properties
directly from the model, such as mass, formula, and the atoms
themselves:

    >>> pdb.model().mass()
    20630.8656
    >>> pdb.model().formula()
    Counter({'C': 2039, 'O': 803, 'N': 565, 'S': 22, 'P': 2})
    >>> len(pdb.model().atoms)
    3431
    >>> pdb.model().get_atoms_by_element("P")
    {<PdbAtom 3200 (P)>, <PdbAtom 3230 (P)>}
    >>> pdb.get_atom_by_id(23)
    <PdbAtom 23 (N)>


The chains and small molecules of the model exist as sets, and can be queried
by ID or name:

    >>> pdb.model().chains()
    {<Chain B (214 residues)>, <Chain A (204 residues)>}
    >>> len(pdb.model().small_molecules()) # Includes solvent molecules
    184
    >>> pdb.model().get_chain_by_id("B")
    <Chain B (214 residues)>
    >>> pdb.model().get_small_molecules_by_name("XMP")
    {<SmallMolecule (XMP)>, <SmallMolecule (XMP)>}


.. note::

   PDB files are not always perfect representations of the real molecular
   structures they are created from. Sometimes there are missing atoms, and
   sometimes there are missing residues. Future versions of molecuPy will flag
   these and maybe even fill them in, but for now simply bear in mind that there
   may be missing atoms and disconnected chains.


Chains
~~~~~~

A :py:class:`.Chain` object is an ordered sequence of Residue objects, and
they are the macromolecular structures which constitute the bulk of the model.

    >>> pdb.model().get_chain_by_id("A")
    <Chain A (204 residues)>
    >>> pdb.model().get_chain_by_id("A").chain_id
    'A'
    >>> pdb.model().get_chain_by_id("A").residues[0]
    <Residue (VAL)>

Chains inherit from :py:class:`.ResiduicStructure` and
:py:class:`.ResiduicSequence` and so have methods for retrieving residues:

    >>> pdb.model().get_chain_by_id("A").get_residue_by_id("A23")
    <Residue (ASN)>
    >>> pdb.model().get_chain_by_id("A").get_residue_by_name("ASP")
    <Residue (ASP)>
    >>> pdb.model().get_chain_by_id("A").get_residues_by_name("ASN")
    {<Residue A5 (ASN)>, <Residue A23 (ASN)>, <Residue A23A (ASN)>, <Residue A10
    1(ASN)>, <Residue A141 (ASN)>, <Residue A199 (ASN)>}
    >>> pdb.model().get_chain_by_id("A").sequence_string()
    'VMNRLILAMDLMNRDDALRVTGEVREYIDTVKIGYPLVLSEGMDIIAEFRKRFGCRIIADFKVADIPETNEKICR
    ATFKAGADAIIVHGFPGADSVRACLNVAEEMGREVFLLTEMSHPGAEMFIQGAADEIARMGVDLGVKNYVGPSTRP
    ERLSRLREIIGQDSFLISPGGETLRFADAIIVGRSIYLADNPAAAAAGIIESI'

Like pretty much everything else in molecuPy, chains are ultimately atomic
structures, and have the usual atomic structure methods for getting mass,
retrieving atoms etc.

The :py:class:`.PdbResidue` objects themselves are also atomic structures, and
behave very similar to small molecules. They have ``downstream_residue`` and
``upstream_residue`` methods for getting the next and previous residue in
their chain respectively.


Small Molecules
~~~~~~~~~~~~~~~

Many PDB files also contain non-macromolecular objects, such as ligands, and
solvent molecules. In molecuPy, these are represented as
:py:class:`.SmallMolecule` objects.

There's not a great deal to be said about small molecules. They are atomic
structures, so you can get their mass, get atoms by name/ID etc.

    >>> pdb.model().get_small_molecule_by_name("BU2")
    <SmallMolecule A500 (BU2)>
    >>> pdb.model().get_small_molecule_by_name("XMP").atoms
    {<PdbAtom 3240 (C)>, <PdbAtom 3241 (N)>, <PdbAtom 3242 (N)>, <PdbAtom 3243 (
    C)>, <PdbAtom 3244 (O)>, <PdbAtom 3245 (C)>, <PdbAtom 3246 (O)>, <PdbAtom 32
    47 (C)>, <PdbAtom 3248 (N)>, <PdbAtom 3249 (C)>, <PdbAtom 3250 (C)>, <PdbAto
    m 3251 (O)>, <PdbAtom 3252 (C)>, <Atom 3253 (O)>, <PdbAtom 3230 (P)>, <PdbAt
    om 3231 (O)>, <PdbAtom 3232 (O)>, <PdbAtom 3233 (O)>, <PdbAtom 3234 (O)>, <P
    dbAtom 3235 (C)>, <PdbAtom 3236 (C)>, <PdbAtom 3237 (O)>, <Atom 3238 (C)>, <
    PdbAtom 3239 (N)>}
    >>> pdb.model().get_small_molecule_by_name("XMP").get_atom_by_id(3252)
    <PdbAtom 3252 (C)>

The :py:class:`.BindSite` binding site of the molecule, if there is one, can be
determined in one of two ways. If the PDB file already defines the site, it can
be found with:

    >>> pdb.model().get_small_molecule_by_name("XMP").bind_site()
    <BindSite AC3 (11 residues)>

If there isn't one defined, you can try to predict it using atomic distances:

    >>> pdb.model().get_small_molecule_by_name("XMP").predict_bind_site()
    <BindSite CALC (5 residues)>

All atomic structures can do this, but it is perhaps most useful with small
molecules.


Atoms
~~~~~

Pdb structures - like everything else in the universe really - are ultimately
collections of Atom - :py:class:`.Atom` - objects. They possess a few key
properties from which much of everything else is created:

    >>> pdb.model().get_atom_by_id(28)
    <PdbAtom 28 (C)>
    >>> pdb.model().get_atom_by_id(28).atom_id()
    28
    >>> pdb.model().get_atom_by_id(28).atom_name()
    'CB'
    >>> pdb.model().get_atom_by_id(28).element()
    'C'
    >>> pdb.model().get_atom_by_id(28).mass()
    12.0107

  molecuPy draws a distinction between generic atom objects, and
  :py:class:`.PdbAtom`, which have coordinates. These are the atoms listed in
  the PDB file as being observed in the experiment that produced it.

  Why the distinction? PDB files also list *missing* atoms - atoms known to be
  present in the structure depicted but which were not observed in the data. For
  those the generic :py:class:`.Atom` class is used.

  There are also missing residues, which are represented here as ordinary
  residues composed entirely of missing atoms. All residues have a
  ``is_missing`` method to make this clear.

The distance between any two PDB atoms can be calculated easily:

    >>> atom1 = pdb.model().get_atom_by_id(23)
    >>> atom2 = pdb.model().get_atom_by_id(28)
    >>> atom1.distance_to(atom2)
    7.931296047935668

Bonds will be assigned where possible - the bonds between atoms in
standard residues are inferred from atom names, and PDB files contain
annotations for other covalent bonds. These are assigned to the atoms as
:py:class:`.tBond` objects.

    >>> pdb.model().get_atom_by_id(27).bonds()
    {<Bond between Atom 27 and Atom 101>, <Bond between Atom 100 and Atom 27>}

The atoms directly bonded to any atom can be obtained with
``bonded_atoms``, and the set of all atoms that are `accessible` is accessed
with ``accessible_atoms``.

    >>> pdb.model().get_atom_by_id(3201)
    <PdbAtom 3200 (P)>
    >>> pdb.model().get_atom_by_id(3201).bonded_atoms()
    {<PdbAtom 3200 (P)>}
    >>> pdb.model().get_atom_by_id(3200).bonded_atoms()
    {<PdbAtom 3203 (O)>, <PdbAtom 3201 (O)>, <PdbAtom 3204 (O)>, <PdbAtom 3202 (
    O)>}
    >>> pdb.model().get_atom_by_id(3200).accessible_atoms()
    {<PdbAtom 3214 (O)>, <PdbAtom 3215 (C)>, <PdbAtom 3216 (O)>, <PdbAtom 3217 (
    C)>, <PdbAtom 3218 (N)>, <PdbAtom 3219 (C)>, <PdbAtom 3201 (O)>, <PdbAtom 32
    20 (C)>, <PdbAtom 3202 (O)>, <PdbAtom 3221 (O)>, <PdbAtom 3203 (O)>, <PdbAto
    m 3222 (C)>, <PdbAtom 3204 (O)>, <PdbAtom 3223 (O)>, <PdbAtom 3205 (C)>, <Pd
    bAtom 3206 (C)>, <PdbAtom 3207 (O)>, <PdbAtom 3208 (C)>, <PdbAtom 3209 (N)>,
     <PdbAtom 3210 (C)>, <PdbAtom 3211 (N)>, <PdbAtom 3212 (N)>, <PdbAtom 3213 (
    C)>}


Similarly, all atoms have a ``model`` method which refers back to their Model,
and as long as this is the case, they can use their ``local_atoms`` method to
return a set of all atoms within a given distance.

    >>> pdb.model().get_atom_by_id(3201).local_atoms(5) # Atoms within 5A
    {<PdbAtom 3214 (O)>, <PdbAtom 3215 (C)>, <PdbAtom 3216 (O)>, <PdbAtom 3217 (
    C)>, <PdbAtom 3218 (N)>}


Binding Sites
~~~~~~~~~~~~~

:py:class:`.BindSite` objects represent binding sites. They are residuic
structures, with the usual residuic structure methods, as well as a ``ligand``
property.

    >>> pdb.model().sites()
    {<BindSite AC2 (5 residues)>, <BindSite AC1 (4 residues)>, <BindSite AC4 (11
     residues)>, <BindSite AC3 (11 residues)>}
    >>> pdb.model().get_site_by_id("AC1").residues()
    {<Residue A10 (ASP)>, <Residue A11 (LEU)>, <Residue A34 (LYS)>}
    >>> pdb.model().get_site_by_id("AC1").ligand()
    <SmallMolecule A1000 (BU2)>



Secondary Structure
~~~~~~~~~~~~~~~~~~~

:py:class:`.Chain` objects have a ``alpha_helices`` property and a
``beta_strands`` property, which are sets of :py:class:`.AlphaHelix` objects
and :py:class:`.BetaStrand` objects respectively.
