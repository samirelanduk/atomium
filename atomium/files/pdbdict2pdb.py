"""This module handles the conversion of PDB data dictionaries to Pdb
objects."""

from .pdb import Pdb
from ..structures import Model, Chain, Residue, Molecule, Atom
from ..structures.reference import bonds

def pdb_dict_to_pdb(pdb_dict):
    """Converts a data ``dict`` to a :py:class:`.Pdb`

    :param dict pdb_dict: The data dictionary to load.
    :rtype: :py:class:`.Pdb`"""

    pdb = Pdb()
    pdb._deposition_date = pdb_dict["deposition_date"]
    pdb._code = pdb_dict["code"]
    pdb._title = pdb_dict["title"]
    pdb._resolution = pdb_dict["resolution"]
    pdb._organism = pdb_dict["organism"]
    pdb._expression_system = pdb_dict["expression_system"]
    pdb._technique = pdb_dict["technique"]
    pdb._classification = pdb_dict["classification"]
    pdb._rfactor = pdb_dict["rfactor"]
    pdb._models = [model_dict_to_model(
     d, pdb_dict["connections"]
    ) for d in pdb_dict["models"]]
    return pdb


def model_dict_to_model(model_dict, connections):
    """Converts a model ``dict`` to a :py:class:`.Model`

    :param dict model_dict: The model dictionary to load.
    :rtype: :py:class:`.Model`"""

    model = Model()
    for chain in model_dict["chains"]:
        model.add_chain(chain_dict_to_chain(chain))
    for molecule in model_dict["molecules"]:
        model.add_molecule(residue_dict_to_residue(molecule, molecule=True))
    bond_atoms(model, connections)
    return model


def chain_dict_to_chain(chain_dict):
    """Converts a chain ``dict`` to a :py:class:`.Chain`

    :param dict chain_dict: The chain dictionary to load.
    :rtype: :py:class:`.Chain`"""

    residues = [residue_dict_to_residue(res) for res in chain_dict["residues"]]
    for index, residue in enumerate(residues[:-1]):
        residue.next(residues[index + 1])
    return Chain(*residues, chain_id=chain_dict["chain_id"])


def residue_dict_to_residue(residue_dict, molecule=False):
    """Converts a residue ``dict`` to a :py:class:`.Residue` or
    :py:class:`.Molecule`.

    :param dict residue_dict: The residue dictionary to load.
    :param bool molecule: if ``True``, a :py:class:`.Molecule` will be returned\
    instead of a :py:class:`.Residue`.
    :rtype: :py:class:`.Residue` or :py:class:`.Molecule`"""

    alt_loc = None
    if any([atom["occupancy"] < 1 for atom in residue_dict["atoms"]]):
        if any([atom["alt_loc"] for atom in residue_dict["atoms"]]):
            alt_loc = sorted([
             atom["alt_loc"] for atom in residue_dict["atoms"] if atom["alt_loc"]
            ])[0]
    atoms = [
     atom_dict_to_atom(atom) for atom in residue_dict["atoms"]
      if atom["occupancy"] ==1 or atom["alt_loc"] is None
       or atom["alt_loc"] == alt_loc
    ]
    if not atoms: print(residue_dict["id"])
    MolClass = Molecule if molecule else Residue
    residue = MolClass(*atoms, name=residue_dict["name"])
    residue._id = residue_dict["id"]
    return residue


def atom_dict_to_atom(atom_dict):
    """Converts an atom ``dict`` to a :py:class:`.Atom`

    :param dict atom_dict: The atom dictionary to load.
    :rtype: :py:class:`.Atom`"""

    return Atom(
     atom_dict["element"], atom_dict["x"], atom_dict["y"], atom_dict["z"],
     atom_id=atom_dict["atom_id"], name=atom_dict["atom_name"],
     charge=atom_dict["charge"],
     bfactor=atom_dict["temp_factor"] if atom_dict["temp_factor"] else 0
    )


def bond_atoms(model, connections):
    """Bonds the atoms of a :py:class:`.Model` in a sensible way.

    :param Model model: The ``Model`` to be connected up."""

    make_intra_residue_bonds(model.residues(), bonds)
    make_inter_residue_bonds(model.residues())
    make_connections_bonds(model, connections)


def make_intra_residue_bonds(residues, bonds):
    """Takes some :py:class:`.Residue` objects and bonds together its atoms
    internally, using a ``dict`` and the residue names as a reference.

    :param residues: A collection of Residues.
    :param dict d: The reference ``dict``"""

    for residue in residues:
        res_ref = bonds.get(residue.name())
        if res_ref:
            for atom in residue.atoms():
                atom_ref = res_ref.get(atom.name())
                if atom_ref:
                    for atom2 in residue.atoms():
                        if atom2.name() in atom_ref:
                            atom.bond(atom2)


def make_inter_residue_bonds(residues):
    """Takes some :py:class:`.Residue` objects and bonds them together with
    peptide bonds. If the relevant atoms are more than 5 Angstroms apart, no
    bond will be made.

    :param residues: A collection of Residues.`"""

    for residue in residues:
        if residue.next():
            c = residue.atom(name="C")
            n = residue.next().atom(name="N")
            if c and n and c.distance_to(n) < 5:
                c.bond(n)


def make_connections_bonds(model, connections):
    """Takes a :py:class:`.Model` and a connections ``list`` and connects the
    atoms according to the specifications in the ``list``.

    :param Model model: A Model to connect up.
    :param list connections: The connections list from a data dictionary"""

    for connection in connections:
        try:
            connection["bond_to"].remove(connection["atom"])
        except: pass
        atom = model.atom(atom_id=connection["atom"])
        if atom:
            for other in connection["bond_to"]:
                other_atom = model.atom(atom_id=other)
                if other_atom:
                    atom.bond(other_atom)
