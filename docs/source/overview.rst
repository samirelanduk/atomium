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

Once you've got your :py:class:`.File` object, what can you do with it?

Annotation
##########

There is meta information contained within the :py:class:`.File` object:

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
are accessed using :py:attr:`~.File.assemblies`.

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
molecules - :py:class:`.Model`, :py:class:`.Chain`, :py:class:`.Residue` and
:py:class:`.Ligand` - you can pass in an ``id`` or ``name``, or search by regex
pattern with ``id__regex`` or ``name__regex``.

These structures have an even more powerful syntax too - you can pass in *any*
property such as ``charge=1``, any comparitor of a property such as
``mass__lt=100``, or any regex of a property such as ``name__regex='[^C]'``.

For pairwise comparisons, structures also have the
:py:meth:`~.AtomStructure.pairwise_atoms` generator which will yield all
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
:py:meth:`~.Model.optimise_distances` on the :py:class:`.Model` first. This
prevents atomium from having to compare every atom with every other atom every
time a proximity check is made.

The :py:class:`.Atom` objects themselves have their own useful properties.

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
