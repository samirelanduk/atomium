from .structures import Model, Atom, Ligand, Chain, Residue

def parse_mmcif_dict(mmcif_dict):
    

    # What are the entities here?
    entities = get_entities(mmcif_dict)

    # What structures are here?
    asyms = sorted(set([(atom["label_asym_id"], atom["label_entity_id"], atom["auth_asym_id"]) for atom in mmcif_dict["atom_site"]]))
    structures = [{
        "id": asym[0],
        "auth_id": asym[2],
        "entity": entities[asym[1]]
    } for asym in asyms]

    # Organise atoms
    atoms = {}
    for atom in mmcif_dict["atom_site"]:
        try:
            atoms[atom["label_asym_id"]].append(atom)
        except KeyError:
            atoms[atom["label_asym_id"]] = [atom]

    chains, ligands, waters = [], [], []
    for structure in structures:
        if structure["entity"]["type"] == "non-polymer":
            ligands.append(make_ligand(atoms[structure["id"]]))
        if structure["entity"]["type"] == "polymer":
            chains.append(make_chain(atoms[structure["id"]]))
        if structure["entity"]["type"] == "water":
            waters += make_waters(atoms[structure["id"]])

    
    model = Model(chains=chains, ligands=ligands, waters=waters)

    f = File(
        name=mmcif_dict.get("entry", [{}])[0].get("id"),
        models = [model]
    )
    f.source = mmcif_dict
    return f


def get_entities(mmcif_dict):
    entities = {e["id"]: {
        "id": e["id"], "type": e["type"], "name": e["pdbx_description"]
    } for e in mmcif_dict["entity"]}
    for entity in entities.values():
        if entity["type"] == "polymer":
            annotate_polymer_entity(entity, mmcif_dict)
        else:
            annotate_non_polymer_entity(entity, mmcif_dict)
    return entities


def annotate_polymer_entity(entity, mmcif_dict):
    pass


def annotate_non_polymer_entity(entity, mmcif_dict):
    pass


def make_chain(atoms):
    residues = []
    res_id = None
    residue_atoms = []
    for atom in atoms:
        if (atom["label_seq_id"], atom["pdbx_PDB_ins_code"]) != res_id:
            if residue_atoms: residues.append(Residue(*residue_atoms))
            res_id = (atom["label_seq_id"], atom["pdbx_PDB_ins_code"])
        residue_atoms.append(Atom(
            atom["type_symbol"], atom["Cartn_x"], atom["Cartn_y"], atom["Cartn_z"],
            atom["id"], atom["label_atom_id"], charge=atom["pdbx_formal_charge"],
            bvalue=atom["B_iso_or_equiv"]
        ))
    return Chain(*residues)
    

def make_ligand(atoms):
    ligand_atoms = [Atom(
        atom["type_symbol"], atom["Cartn_x"], atom["Cartn_y"], atom["Cartn_z"],
        atom["id"], atom["label_atom_id"], charge=atom["pdbx_formal_charge"],
        bvalue=atom["B_iso_or_equiv"]
    ) for atom in atoms]
    return Ligand(*ligand_atoms, name=atoms[0]["label_comp_id"], id=make_id(atoms[0]))


def make_waters(atoms):
    waters = []
    for atom in atoms:
        waters.append(Ligand(Atom(
            atom["type_symbol"], atom["Cartn_x"], atom["Cartn_y"], atom["Cartn_z"],
            atom["id"], atom["label_atom_id"], charge=atom["pdbx_formal_charge"],
            bvalue=atom["B_iso_or_equiv"]
        ), name=atom["label_comp_id"], id=make_id(atom)))
    return waters


def make_id(atom):
    chain = atom["auth_asym_id"]
    number = atom["auth_seq_id"]
    insertion = atom["pdbx_PDB_ins_code"].replace("?", "") or ""
    return f"{chain}.{number}{insertion}"



class File:
    
    def __init__(self, name="", models=None):
        self.name = name
        self.models = models or []
    

    def __repr__(self):
        return f"<File: {self.name}>" if self.name else "<File>"


    def __getattr__(self, name):
        if name in self.source: return self.source[name]
        if "__" in name:
            category, data = name.split("__")[:2]
            if category in self.source and data in self.source[category][0]:
                return self.source[category][0][data]
        raise AttributeError(f"File has no attribute '{name}'")


