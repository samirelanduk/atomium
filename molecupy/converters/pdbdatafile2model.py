from ..structures import Model, Atom, GhostAtom, SmallMolecule, Residue, Chain
from ..structures import BindSite, AlphaHelix, BetaStrand, Complex
from ..pdb.pdbdatafile import PdbDataFile

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


def add_chains_to_model(model, data_file, model_id):
    atoms = data_file.atoms()
    atom_id = max([atom["atom_id"] for atom in atoms]) * 100 if atoms else 1000
    chain_ids = set([a["chain_id"] for a in atoms])
    for chain_id in chain_ids:
        relevant_atoms = [a for a in atoms if a["model_id"] == model_id
         and a["chain_id"] == chain_id]
        residue_ids = []
        for a in relevant_atoms:
            id_ = _mol_id_from_atom(a)
            if id_ not in residue_ids:
                residue_ids.append(id_)
        residues = []
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
            while len(set([res[1] for res in missing_residue_ids])) < len(missing_residue_ids):
                ids = [res[1] for res in missing_residue_ids]
                duplicates = [id_ for id_ in ids if ids.count(id_) > 1]
                id_to_remove = duplicates[0]
                for res in missing_residue_ids[::-1]:
                    if res[1] == id_to_remove:
                        missing_residue_ids.remove(res)
                        break
            missing_residues = []
            for missing_id in missing_residue_ids:
                lookup = residues_dict.connection_data.get(missing_id[0])
                if lookup:
                    empty_atoms = []
                    for atom_name in lookup:
                        atom = GhostAtom(atom_name[0], atom_id, atom_name)
                        atom_id += 1
                        empty_atoms.append(atom)
                    residue = Residue(missing_id[1], missing_id[0], *empty_atoms)
                    missing_residues.append(residue)
                else:
                    empty_atoms = []
                    for atom_name in ["C", "CA", "N"]:
                        atom = Atom(atom_name[0], atom_id, atom_name)
                        atom_id += 1
                        empty_atoms.append(atom)
                    residue = Residue(missing_id[1], missing_id[0], *empty_atoms)
                    missing_residues.append(residue)
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
        chain = Chain(chain_id, *residues)
        model.add_chain(chain)


def _mol_id_from_atom(atom):
    return atom["chain_id"] + str(atom["residue_id"]) + atom["insert_code"]
