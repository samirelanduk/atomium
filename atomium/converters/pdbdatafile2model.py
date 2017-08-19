"""Contains the function for creating Models from PdbDataFiles."""

from ..structures.models import Model
from ..structures.chains import Chain
from ..structures.molecules import Residue
from ..structures.atoms import Atom

def pdb_data_file_to_model(data_file):
    """Converts a :py:class:`.PdbDataFile` to a :py:class:`.Model`

    :param PdbDataFile data_file: The PdbDataFile.
    :rtype: ``Model``"""

    model = Model()
    load_chains(data_file.atoms, model)
    return model


def load_chains(atom_dicts, model):
    """Adds :py:class:`.Chain` objects to a :py:class:`.Model` from a list of
    atom ``dict``s.

    :param list atom_dicts: The atom ``dict`` objects.
    :param Model model: The Model to update."""

    atoms = [atom_dict_to_atom(d) for d in atom_dicts]
    residues = atoms_to_residues(atoms)
    chains = residues_to_chains(residues)
    for chain in chains:
        model.add_chain(chain)


def atom_dict_to_atom(d):
    """Turns an atom ``dict`` into an actual :py:class:`.Atom`.

    :param dict d: The ``dict`` to convert.
    :rtype: ``Atom``"""

    atom = Atom(
     d["element"] if d["element"] else "X",
     d["x"], d["y"], d["z"],
     atom_id=d["atom_id"], name=d["atom_name"],
     charge=d["charge"] if d["charge"] else 0
    )
    atom.temp_chain_id = d["chain_id"]
    atom.temp_residue_id = d["chain_id"] + str(d["residue_id"])
    if d["insert_code"]: atom.temp_residue_id += d["insert_code"]
    atom.temp_residue_name = d["residue_name"]
    return atom


def atoms_to_residues(atoms):
    """Takes a list of atoms and creates a list of :py:class:`.Residue` objects
    from them, with the order preserved.

    :param atoms: A collection of :py:class:`.Atom` objects.
    :returns: ``list`` of ``Residue``s"""

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
    :returns: ``list`` of ``Chain``s"""

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
