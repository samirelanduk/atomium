from ..structures import Model, Atom, GhostAtom, SmallMolecule, Residue, Chain
from ..structures import BindSite, AlphaHelix, BetaStrand, Complex
from ..pdb.pdbdatafile import PdbDataFile
from ..pdb import residues as residues_dict

def model_from_pdb_data_file(data_file, model_id=1):
    if not isinstance(data_file, PdbDataFile):
        raise TypeError("model_from_pdb_data_file can only convert PdbDataFiles")
    for model_dict in data_file.models():
        if model_dict["model_id"] == model_id:
            model = Model()
            model._source = data_file

            add_small_molecules_to_model(model, data_file, model_id)
            add_chains_to_model(model, data_file, model_id)

            return model
    raise ValueError("There is no model with ID %i" % model_id)


def add_small_molecules_to_model(model, data_file, model_id):
    heteroatoms = data_file.heteroatoms()
    molecule_names = set([a["residue_name"] for a in heteroatoms])
    molecule_ids = set([_mol_id_from_atom(a) for a in heteroatoms])
    for molecule_name in molecule_names:
        for molecule_id in molecule_ids:
            relevant_atoms = [a for a in heteroatoms if a["model_id"] == model_id
             and a["residue_name"] == molecule_name and _mol_id_from_atom(a) == molecule_id]
            if relevant_atoms:
                atoms = [Atom(
                 a["x"], a["y"], a["z"],
                 a["element"],
                 a["atom_id"],
                 a["atom_name"]
                ) for a in relevant_atoms]
                small_molecule = SmallMolecule(
                 molecule_id, molecule_name, *atoms
                )
                model.add_small_molecule(small_molecule)


def _get_top_atom_id(atoms=None, heteroatoms=None, missing_info=None):
    atom_id = max([atom["atom_id"] for atom in atoms]) if atoms else -1
    return atom_id


def add_chains_to_model(model, data_file, model_id):
    atoms = [a for a in data_file.atoms() if a["model_id"] == model_id]
    heteroatoms = [a for a in data_file.atoms() if a["model_id"] == model_id]
    atom_id = max([atom["atom_id"] for atom in atoms]) if atoms else -1
    het_id = max([atom["atom_id"] for atom in heteroatoms]) if heteroatoms else -1
    highest_id = max(atom_id, het_id)
    chain_ids = set([a["chain_id"] for a in atoms])
    for chain_id in chain_ids:
        residues = []
        _add_residues_to_chain(residues, chain_id, atoms)
        missing_residue_info = _get_missing_residue_info(data_file, chain_id)
        if missing_residue_info:
            _add_missing_residues_to_chain(residues, missing_residue_info, highest_id)
        chain = Chain(chain_id, *residues)
        model.add_chain(chain)


def _mol_id_from_atom(atom):
    return atom["chain_id"] + str(atom["residue_id"]) + atom["insert_code"]


def _add_residues_to_chain(residues, chain_id, atoms):
    relevant_atoms = [a for a in atoms if a["chain_id"] == chain_id]
    residue_ids = []
    for a in relevant_atoms:
        id_ = _mol_id_from_atom(a)
        if id_ not in residue_ids:
            residue_ids.append(id_)
    for residue_id in residue_ids:
        residue_atoms = [
         a for a in relevant_atoms if _mol_id_from_atom(a) == residue_id
        ]
        pdb_atoms = [Atom(
         a["x"], a["y"], a["z"],
         a["element"],
         a["atom_id"],
         a["atom_name"]
        ) for a in residue_atoms]
        residue = Residue(
         residue_id, residue_atoms[0]["residue_name"], *pdb_atoms
        )
        residues.append(residue)


def _get_missing_residue_info(data_file, chain_id):
    missing_residues_remark = data_file.get_remark_by_number(465)
    if missing_residues_remark:
        missing_residues = missing_residues_remark["content"].split("\n")[6:]
        missing_residues = [line.split() for line in missing_residues if line]
        missing_residues = [
         res for res in missing_residues if len(res) == 3 or int(res[0]) == model_id
        ]
        missing_residues = [
         res for res in missing_residues if res[-2] == chain_id
        ]
        missing_residue_ids = [(res[0], res[-2] + res[-1]) for res in missing_residues]

        # Remove duplicate residue IDs - I can't remember why this is needed...
        while len(set([res[1] for res in missing_residue_ids])) < len(missing_residue_ids):
            ids = [res[1] for res in missing_residue_ids]
            duplicates = [id_ for id_ in ids if ids.count(id_) > 1]
            id_to_remove = duplicates[0]
            for res in missing_residue_ids[::-1]:
                if res[1] == id_to_remove:
                    missing_residue_ids.remove(res)
                    break
        return missing_residue_ids


def _add_missing_residues_to_chain(residues, missing_residue_info, highest_id):
    missing_residues = []
    atom_id = highest_id + 1
    for missing_id in missing_residue_info:
        missing_residues.append(_create_missing_residue_from_id(missing_id, atom_id))
        atom_id += len(missing_residues[-1].atoms(atom_type="all"))

    for missing_residue in missing_residues:
        closest_smaller_id = None
        if not _residue_id_is_greater_than_residue_id(residues[0].residue_id(), missing_residue.residue_id()):
            closest_smaller_id = sorted(
             residues,
             key=lambda k: _residue_id_to_int(missing_residue.residue_id()) - _residue_id_to_int(
              k.residue_id()
             ) if _residue_id_to_int(missing_residue.residue_id()) > _residue_id_to_int(
              k.residue_id()
             ) else float("inf"))[0]
        residues.insert(
         residues.index(closest_smaller_id) + 1 if closest_smaller_id else 0,
         missing_residue
        )


def _create_missing_residue_from_id(missing_id, atom_id):
    lookup = residues_dict.connection_data.get(missing_id[0])
    residue = None
    if lookup:
        empty_atoms = []
        for atom_name in lookup:
            atom = GhostAtom(atom_name[0], atom_id, atom_name)
            atom_id += 1
            empty_atoms.append(atom)
        residue = Residue(missing_id[1], missing_id[0], *empty_atoms)
    else:
        empty_atoms = []
        for atom_name in ["C", "CA", "N"]:
            atom = GhostAtom(atom_name[0], atom_id, atom_name)
            atom_id += 1
            empty_atoms.append(atom)
        residue = Residue(missing_id[1], missing_id[0], *empty_atoms)
    return residue


def _residue_id_to_int(residue_id):
    int_component = int(
     "".join([char for char in residue_id if char in "-0123456789"])
    )
    char_component = residue_id[-1] if residue_id[-1] not in "-0123456789" else ""
    int_component *= 100
    return int_component + (ord(char_component) - 64 if char_component else 0)


def _residue_id_is_greater_than_residue_id(residue_id1, residue_id2):
    residues = []
    for residue_id in (residue_id1, residue_id2):
        chain_id = residue_id[0] if residue_id[0].isalpha() else ""
        number = int("".join([char for char in residue_id if char.isnumeric() or char == "-"]))
        insert = ord(residue_id[-1]) if residue_id[-1].isalpha() else 0
        residues.append((chain_id, number, insert))
    if residues[0][1] > residues[1][1]:
        return True
    elif residues[0][1] == residues[1][1] and residues[0][2] > residues[1][2]:
        return True
    else:
        return False
