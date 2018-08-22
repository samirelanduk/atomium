"""Contains the code for dealing with and defining atomium data dictionaries."""

from copy import deepcopy
from itertools import groupby
from .file import File
from ..models import *
from ..models.data import CODES, BONDS

def generate_higher_structures(models):
    """Takes a list of model dictionaries, and for each it goes through the
    atoms, and creates dictionaries for chains, residues and ligands based on
    the IDs it finds in the atoms. The model dictionaries are updated in place.

    :param list models: The list of model dictionaries."""

    for model in models:
        chain_ids, residue_id, ligand_id = [], None, None
        for atom in model["atoms"]:
            if atom["chain_id"] not in chain_ids:
                chain_ids.append(atom["chain_id"])
                chain = deepcopy(CHAIN_DICT)
                chain["id"] = atom["chain_id"]
                model["chains"].append(chain)
            if (atom["polymer"] and residue_id != atom["full_res_id"]) or\
             (not atom["polymer"] and ligand_id != atom["full_res_id"]):
                if atom["polymer"]:
                    residue_id = atom["full_res_id"]
                else:
                    ligand_id = atom["full_res_id"]
                het = deepcopy(RESIDUE_DICT)
                het["id"] = atom["full_res_id"]
                het["name"] = atom["residue_name"]
                het["chain_id"] = atom["chain_id"]
                model["residues" if atom["polymer"] else "ligands"].append(het)


def data_dict_to_file(d):
    """Takes an atomium data dictionary and turns it into a :py:class:`.File`
    object.

    :param dict d: The atomium data dictionary.
    :rtype: ``File``"""

    file_ = File()
    for key, value in d.items():
        if key != "models":
            for key2, value2 in value.items():
                setattr(file_, "_" + key2, value2)
    file_._models = [model_dict_to_model(m) for m in d["models"]]
    return file_


def model_dict_to_model(m):
    """Takes an atomium model dictionary, and turns it into a :py:class:`.Model`
    object.

    :param dict m: The dictionary to read.
    :rtype: ``Model``"""

    atoms = [atom_dict_to_atom(a) for a in m["atoms"]]
    model = Model(*atoms)
    for chain in m["chains"]:
        residues = [r for r in m["residues"] if r["chain_id"] == chain["id"]]
        ligands = [l for l in m["ligands"] if l["chain_id"] == chain["id"]]
        residue_objects, ligand_objects = [], []
        for residue in residues:
            residue_objects.append(create_het(residue, m["atoms"], model))
        for res1, res2 in zip(residue_objects[:-1], residue_objects[1:]):
            res1.next = res2
        for ligand in ligands:
            ligand_objects.append(create_het(ligand, m["atoms"], model))
        rep = "".join([CODES.get(res, "X") for res in chain["full_sequence"]])
        try:
            Chain(*(residue_objects + ligand_objects), id=chain["id"], rep=rep)
        except: exit()
    for a in model.atoms():
        if not a.residue and not a.ligand and a.name:
            model.remove(a)
    bond_atoms(model, m["connections"])
    return model


def atom_dict_to_atom(atom_dict):
    """Converts an atom ``dict`` to a :py:class:`.Atom`.

    :param dict atom_dict: The atom dictionary to load.
    :rtype: :py:class:`.Atom`"""

    return Atom(
     atom_dict["element"], atom_dict["x"], atom_dict["y"], atom_dict["z"],
     id=atom_dict["id"], name=atom_dict["name"],
     charge=atom_dict["charge"], anisotropy=atom_dict["anisotropy"],
     bfactor=atom_dict["bfactor"] if atom_dict["bfactor"] else 0
    )


def create_het(het, atom_dicts, model):
    """Creates either a :py:class:`.Residue` or :py:class:`.Ligand` from the
    relevant dictionary, as well as a source of atom dictionaries and a
    :py:class:`.Model` to update.

    :param dict het: The residue/ligand dictionary to read.
    :param list atom_dicts: The list of atom dictionaries to look through.
    :param Model model: The :py:class:`.Model` to update.
    :rtype: ``Ligand`` or ``Residue``"""

    atom_dicts = [a for a in atom_dicts if a["full_res_id"] == het["id"]]
    alt_loc = None
    if any([atom["occupancy"] < 1 for atom in atom_dicts]):
        if any([atom["alt_loc"] for atom in atom_dicts]):
            alt_loc = sorted([atom["alt_loc"]
             for atom in atom_dicts if atom["alt_loc"]])[0]
    atoms = [model.atom(a["id"]) for a in atom_dicts if a["occupancy"] == 1
     or a["alt_loc"] is None or a["alt_loc"] == alt_loc]
    Het = Residue if atom_dicts[0]["polymer"] else Ligand
    return Het(*atoms, id=het["id"], name=het["name"])


def bond_atoms(model, connections):
    """Bonds the atoms of a :py:class:`.Model` in a sensible way.

    :param Model model: The ``Model`` to be connected up.
    :param list connections: The list of connections to use."""

    make_intra_residue_bonds(model.residues(), BONDS)
    make_inter_residue_bonds(model.residues())
    make_connections_bonds(model, connections)


def make_intra_residue_bonds(residues, bonds):
    """Takes some :py:class:`.Residue` objects and bonds together its atoms
    internally, using a ``dict`` and the residue names as a reference.

    :param residues: A collection of Residues.
    :param dict d: The reference ``dict``"""

    for residue in residues:
        res_ref = bonds.get(residue.name)
        if res_ref:
            for atom in residue.atoms():
                atom_ref = res_ref.get(atom.name)
                if atom_ref:
                    for atom2 in residue.atoms():
                        if atom2.name in atom_ref:
                            atom.bond_to(atom2)


def make_inter_residue_bonds(residues):
    """Takes some :py:class:`.Residue` objects and bonds them together with
    peptide bonds. If the relevant atoms are more than 5 Angstroms apart, no
    bond will be made.

    :param residues: A collection of Residues.`"""

    for residue in residues:
        if residue.next:
            c = residue.atom(name="C")
            n = residue.next.atom(name="N")
            if c and n and c.distance_to(n) < 5:
                c.bond_to(n)


def make_connections_bonds(model, connections):
    """Takes a :py:class:`.Model` and a connections ``list`` and connects the
    atoms according to the specifications in the ``list``.

    :param Model model: A Model to connect up.
    :param list connections: The connections list from a data dictionary"""

    for connection in connections:
        try:
            connection["bond_to"].remove(connection["atom"])
        except: pass
        atom = model.atom(id=connection["atom"])
        if atom:
            for other in connection["bond_to"]:
                other_atom = model.atom(id=other)
                if other_atom:
                    atom.bond_to(other_atom)



DATA_DICT = {
 "description": {
  "code": None,
  "title": None,
  "deposition_date": None,
  "classification": None,
  "keywords": [],
  "authors": []
 }, "experiment": {
  "technique": None,
  "source_organism": None,
  "expression_system": None
 }, "quality": {
  "resolution": None,
  "rvalue": None,
  "rfree": None
 }, "geometry": {
  "assemblies": []
 }, "models": []
}

MODEL_DICT = {
 "chains": [], "residues": [], "ligands": [], "atoms": [], "connections": []
}

CHAIN_DICT = {
 "id": None,
 "full_sequence": []
}

RESIDUE_DICT = {
 "id": None,
 "name": None,
}

ATOM_DICT = {
 "id": 0,
 "element": None,
 "name": None,
 "x": None, "y": None, "z": None,
 "bfactor": None,
 "charge": 0,
 "residue_id": None,
 "residue_name": None,
 "residue_insert": "",
 "chain_id": None,
 "occupancy": 1,
 "alt_loc": None,
 "anisotropy": [],
 "polymer": False,
 "full_res_id": None
}
