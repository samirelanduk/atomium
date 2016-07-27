from .structures import Model, PdbAtom, Atom, SmallMolecule, Residue, Chain
from . import residues as residues_dict

class Pdb:

    def __init__(self, data_file):
        self._data_file = data_file
        self._models = []
        for model_dict in self._data_file.models():
            model = Model()
            _give_model_small_molecules(model, data_file, model_dict["model_id"])
            _give_model_chains(model, data_file, model_dict["model_id"])
            _connect_atoms(model, data_file, model_dict["model_id"])
            _bond_residue_atoms(model, data_file, model_dict["model_id"])
            _bond_residues_together(model, data_file, model_dict["model_id"])
            self._models.append(model)


    def __repr__(self):
        return "<Pdb (%s)>" % (self.pdb_code() if self.pdb_code() else "????")


    def data_file(self):
        return self._data_file


    def classification(self):
    	return self._data_file.classification()


    def deposition_date(self):
    	return self._data_file.deposition_date()


    def pdb_code(self):
    	return self._data_file.pdb_code()


    def is_obsolete(self):
    	return self._data_file.is_obsolete()


    def obsolete_date(self):
    	return self._data_file.obsolete_date()


    def replacement_code(self):
    	return self._data_file.replacement_code()


    def title(self):
    	return self._data_file.title()


    def split_codes(self):
    	return self._data_file.split_codes()


    def caveat(self):
    	return self._data_file.caveat()


    def keywords(self):
    	return self._data_file.keywords()


    def experimental_techniques(self):
    	return self._data_file.experimental_techniques()


    def model_count(self):
    	return self._data_file.model_count()


    def model_annotations(self):
    	return self._data_file.model_annotations()


    def authors(self):
    	return self._data_file.authors()


    def revisions(self):
    	return self._data_file.revisions()


    def supercedes(self):
    	return self._data_file.supercedes()


    def supercede_date(self):
    	return self._data_file.supercede_date()


    def journal(self):
    	return self._data_file.journal()


    def models(self):
        return self._models


    def model(self):
        return self._models[0]



def _give_model_small_molecules(model, data_file, model_id):
    heteroatoms = data_file.heteroatoms()
    molecule_names = set([a["residue_name"] for a in heteroatoms])
    molecule_ids = set([a["chain_id"] + str(a["residue_id"]) + a["insert_code"]
     for a in heteroatoms])
    for molecule_name in molecule_names:
        for molecule_id in molecule_ids:
            relevant_atoms = [a for a in heteroatoms if a["model_id"] == model_id
             and a["residue_name"] == molecule_name and a["chain_id"] +
              str(a["residue_id"]) + a["insert_code"] == molecule_id]
            if relevant_atoms:
                atoms = [PdbAtom(
                 a["x"], a["y"], a["z"],
                 a["element"],
                 a["atom_id"],
                 a["atom_name"]
                ) for a in relevant_atoms]
                small_molecule = SmallMolecule(
                 molecule_id, molecule_name, *atoms
                )
                model.add_small_molecule(small_molecule)


def _give_model_chains(model, data_file, model_id):
    atoms = data_file.atoms()
    atom_id = max([atom["atom_id"] for atom in atoms]) * 100 if atoms else 1000
    chain_ids = set([a["chain_id"] for a in atoms])
    for chain_id in chain_ids:
        relevant_atoms = [a for a in atoms if a["model_id"] == model_id
         and a["chain_id"] == chain_id]
        residue_ids = []
        for a in relevant_atoms:
            id_ = str(a["residue_id"]) + a["insert_code"]
            if id_ not in residue_ids:
                residue_ids.append(id_)
        residues = []
        for residue_id in residue_ids:
            residue_atoms = [
             a for a in relevant_atoms if str(a["residue_id"]) + a["insert_code"] == residue_id
            ]
            pdb_atoms = [PdbAtom(
             a["x"], a["y"], a["z"],
             a["element"],
             a["atom_id"],
             a["atom_name"]
            ) for a in residue_atoms]
            residue = Residue(
             chain_id + residue_id, residue_atoms[0]["residue_name"], *pdb_atoms
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
            missing_residues = []
            for missing_id in missing_residue_ids:
                lookup = residues_dict.connection_data.get(missing_id[0])
                if lookup:
                    empty_atoms = []
                    for atom_name in lookup:
                        atom = Atom(atom_name[0], atom_id, atom_name)
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


def _connect_atoms(model, data_file, model_id):
    for connection in data_file.connections():
        atom = model.get_atom_by_id(connection["atom_id"])
        if atom:
            for bonded_atom_id in connection["bonded_atoms"]:
                bonded_atom = model.get_atom_by_id(bonded_atom_id)
                if bonded_atom:
                    atom.bond_to(bonded_atom)


def _bond_residue_atoms(model, data_file, model_id):
    for chain in model.chains():
        for residue in chain.residues(include_missing=False):
            lookup = residues_dict.connection_data.get(residue.residue_name())
            if lookup:
                for atom in residue.atoms():
                    atom_lookup = lookup.get(atom.atom_name())
                    if atom_lookup:
                        for other_atom_name in atom_lookup:
                            other_atom = residue.get_atom_by_name(other_atom_name)
                            if other_atom:
                                atom.bond_to(other_atom)


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


def _bond_residues_together(model, data_file, model_id):
    for chain in model.chains():
        for index, residue in enumerate(chain.residues()[:-1]):
            next_residue = chain.residues()[index + 1]
            residue.connect_to(next_residue)
            carboxy_atom = residue.get_atom_by_name("C", atom_type="pdb")
            amino_nitrogen = next_residue.get_atom_by_name("N", atom_type="pdb")
            if carboxy_atom and amino_nitrogen:
                carboxy_atom.bond_to(amino_nitrogen)
