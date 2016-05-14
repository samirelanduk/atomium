from ...molecules import Atom, Molecule
from ...macromolecules import MacroModel

class Pdb:

    def __init__(self, data_file):
        self.data_file = data_file

        transfer_attrs = [
         "classification",
         "deposition_date",
         "pdb_code",
         "is_obsolete",
         "obsolete_date",
         "replacement_code",
         "title",
         "split_codes",
         "caveat",
         "keywords",
         "experimental_techniques",
         "model_num",
         "model_annotations",
         "authors",
         "revisions",
         "supercedes",
         "supercede_date",
         "journal"
        ]
        for attr in transfer_attrs:
            self.__dict__[attr] = self.data_file.__dict__[attr]

        self.models = []
        data_models = self.data_file.models if self.data_file.models else [{"model_id": 1}]
        for model in data_models:
            macro_model = MacroModel()
            add_small_molecules(macro_model, self.data_file, model["model_id"])
            self.models.append(macro_model)
        self.model = self.models[0]


    def __repr__(self):
        return "<Pdb (%s)>" % (self.pdb_code if self.pdb_code else "????")



def add_small_molecules(model, data_file, model_id):
    atoms = [atom for atom in data_file.heteroatoms if atom["model_id"] == model_id]
    small_molecule_names = set([a["residue_name"] for a in atoms])
    small_molecule_ids = set([a["residue_id"] for a in atoms])
    for small_molecule_name in small_molecule_names:
        for small_molecule_id in small_molecule_ids:
            relevant_atoms = [
             a for a in atoms if a["residue_id"] == small_molecule_id
              and a["residue_name"] == small_molecule_name
            ]
            if relevant_atoms:
                relevant_atoms = set([Atom(
                 atom["x"], atom["y"], atom["z"],
                 atom["element"],
                 atom_id=atom["atom_id"],
                 atom_name=atom["atom_name"]
                ) for atom in relevant_atoms])
                for connection in data_file.connections:
                    for atom in relevant_atoms:
                        if connection["atom_id"] == atom.atom_id:
                            for other_atom_id in connection["bonded_atoms"]:
                                for other_atom in relevant_atoms:
                                    if other_atom_id == other_atom.atom_id:
                                        atom.covalent_bond_to(other_atom)
                model.add_small_molecule(Molecule(
                 *relevant_atoms,
                 molecule_id=small_molecule_id,
                 molecule_name=small_molecule_name
                ))
