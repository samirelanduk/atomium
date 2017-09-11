"""Contains the functions for creating Models from PdbDataFiles."""

from ..structures.models import Model
from ..structures.chains import Chain
from ..structures.molecules import Residue, Molecule
from ..structures.atoms import Atom

def pdb_data_file_to_models(data_file):
    """Converts a :py:class:`.PdbDataFile` to a list of :py:class:`.Model`
    objects.

    :param PdbDataFile data_file: The PdbDataFile.
    :rtype: ``list``"""

    model_numbers = sorted(set([a["model"] for a in data_file.atoms]))
    models = []
    for number in model_numbers:
        model = Model()
        load_chains(data_file.atoms, model, number)
        load_molecules(data_file.heteroatoms, model, number)
        bond_atoms(data_file.connections, model)
        models.append(model)
    return models


def load_chains(atom_dicts, model, model_number):
    """Adds :py:class:`.Chain` objects to a :py:class:`.Model` from a list of
    atom ``dict`` objects, using only those relevant to the specific model
    number.

    :param list atom_dicts: The atom ``dict`` objects.
    :param Model model: The Model to update.
    :param int model_number: The model number to use."""

    atoms = [
     atom_dict_to_atom(d) for d in atom_dicts if d["model"] == model_number
    ]
    residues = atoms_to_residues(atoms)
    chains = residues_to_chains(residues)
    for chain in chains:
        model.add_chain(chain)


def load_molecules(atom_dicts, model, model_number):
    """Adds :py:class:`.Molecule` objects to a :py:class:`.Model` from a list of
    atom ``dict`` objects, using only those relevant to the specific model
    number.

    :param list atom_dicts: The atom ``dict`` objects.
    :param Model model: The Model to update.
    :param int model_number: The model number to use."""

    atoms = [
     atom_dict_to_atom(d) for d in atom_dicts if d["model"] == model_number
    ]
    molecules = atoms_to_residues(atoms, molecule=True)
    for mol in molecules:
        model.add_molecule(mol)


def bond_atoms(connections, model):
    """Bonds the atoms of a :py:class:`.Model` in a sensible way.

    :param Model model: The Model to be connected up."""

    from ..structures.reference import bonds
    make_intra_residue_bonds(model.residues(), bonds)
    make_inter_residue_bonds(model.residues())
    make_connections_bonds(model, connections)


def atom_dict_to_atom(d):
    """Turns an atom ``dict`` into an actual :py:class:`.Atom`.

    :param dict d: The ``dict`` to convert.
    :rtype: ``Atom``"""

    atom = Atom(
     d["element"] if d["element"] else "X",
     d["x"], d["y"], d["z"],
     atom_id=d["atom_id"], name=d["atom_name"],
     charge=d["charge"] if d["charge"] else 0,
     bfactor=d["temperature_factor"] if d["temperature_factor"] else 0
    )
    atom.temp_chain_id = d["chain_id"]
    atom.temp_residue_id = d["chain_id"] + str(d["residue_id"])
    if d["insert_code"]: atom.temp_residue_id += d["insert_code"]
    atom.temp_residue_name = d["residue_name"]
    return atom


def atoms_to_residues(atoms, molecule=False):
    """Takes a list of atoms and creates a list of :py:class:`.Residue` objects
    from them, with the order preserved. Alternatively, simple
    :py:class:`.Molecule` objects can be produced instead.

    :param atoms: A collection of :py:class:`.Atom` objects.
    :param bool molecule: If ``True`` Molecules will be returned rather than\
    Residues.
    :returns: ``list``"""

    residue_ids = [atom.temp_residue_id for atom in atoms]
    unique_residue_ids = list(sorted(set(residue_ids), key=residue_ids.index))
    residues = [{"res": res_id, "atoms": []} for res_id in unique_residue_ids]
    for atom in atoms:
        for residue in residues:
            if residue["res"] == atom.temp_residue_id:
                residue["atoms"].append(atom)
                break
    real_residues = []
    for r in residues:
        if molecule:
            real_residues.append(Molecule(
             *r["atoms"], molecule_id=r["res"], name=r["atoms"][0].temp_residue_name
            ))
        else:
            real_residues.append(Residue(
             *r["atoms"], residue_id=r["res"], name=r["atoms"][0].temp_residue_name
            ))
            real_residues[-1].temp_chain_id = r["atoms"][0].temp_chain_id
    for atom in atoms:
        del atom.temp_residue_id
        del atom.temp_residue_name
        del atom.temp_chain_id
    return real_residues


def residues_to_chains(residues):
    """Takes a list of residues and creates a list of :py:class:`.Chain` objects
    from them, with the order preserved.

    :param residues: A collection of :py:class:`.Residue` objects.
    :returns: ``list`` of ``Chain`` objects"""

    chain_ids = [res.temp_chain_id for res in residues]
    unique_chain_ids = list(sorted(set(chain_ids), key=chain_ids.index))
    chains = [{"chain": chain, "residues": []} for chain in unique_chain_ids]
    for residue in residues:
        for chain in chains:
            if chain["chain"] == residue.temp_chain_id:
                chain["residues"].append(residue)
                break
    real_chains = []
    for chain in chains:
        for index, residue in enumerate(chain["residues"][:-1]):
            residue.next(chain["residues"][index + 1])
        real_chains.append(Chain(*chain["residues"], chain_id=chain["chain"]))
    for residue in residues:
        del residue.temp_chain_id
    return real_chains


def make_intra_residue_bonds(residues, d):
    """Takes some :py:class:`.Residue` objects and bonds together its atoms
    internally, using a ``dict`` and the residue names as a reference.

    :param residues: A collection of Residues.
    :param dict d: The reference ``dict``"""

    for residue in residues:
        res_ref = d.get(residue.name())
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
    :param list connections: The connections list from a\
    :py:class:`.PdbDataFile`"""

    for connection in connections:
        atom = model.atom(atom_id=connection["atom"])
        if atom:
            for other in connection["bond_to"]:
                other_atom = model.atom(atom_id=other)
                if other_atom:
                    atom.bond(other_atom)
