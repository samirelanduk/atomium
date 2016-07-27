from .structures import Model, PdbAtom, SmallMolecule, Residue, Chain

class Pdb:

    def __init__(self, data_file):
        self._data_file = data_file
        self._models = []
        for model_dict in self._data_file.models():
            model = Model()
            _give_model_small_molecules(model, data_file, model_dict["model_id"])
            _give_model_chains(model, data_file, model_dict["model_id"])
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
        chain = Chain(chain_id, *residues)
        model.add_chain(chain)
