from datetime import datetime
from .structures import Entity, Model, Atom, Ligand, Chain, Carbohydrate, Residue
from .assemblies import extract_assemblies
from .data import CODES

def parse_mmcif_dict(mmcif_dict):
    

    # What are the entities here?
    entities = get_entities(mmcif_dict)

    # What is the secondary structure
    secondary_structure = get_secondary_structure(mmcif_dict)

    # What models are here?
    model_atoms = divide_atoms_into_models(mmcif_dict)
    
    models = []
    for model in model_atoms:
        # What structures are here?
        structures = get_structures(model, entities)
        
        # Organise atoms by asym ID
        atoms = organise_atoms_by_asym_id(model)
        
        # Make objects
        chains, carbs, ligands, waters = [], [], [], []
        for structure in structures:
            if structure["entity"].type == "non-polymer":
                ligands.append(make_comp(atoms[structure["id"]], entities, ligand=True))
            if structure["entity"].type == "branched":
                carbs.append(make_chain(atoms[structure["id"]], entities, carb=True))
            if structure["entity"].type == "polymer":
                chains.append(make_chain(
                    atoms[structure["id"]], entities, secondary_structure
                ))
            if structure["entity"].type == "water":
                waters += make_waters(atoms[structure["id"]], entities)

        models.append(Model(chains=chains, carbohydrates=carbs, ligands=ligands, waters=waters))

    f = File(
        name=mmcif_dict.get("entry", [{}])[0].get("id"),
        models=models
    )
    f.source = mmcif_dict
    return f


def get_entities(mmcif_dict):
    entities = {e["id"]: Entity(
        id=e["id"], type=e["type"], name=e["pdbx_description"], sequence=""
    ) for e in mmcif_dict["entity"]}
    for entity in entities.values():
        if entity.type == "polymer":
            annotate_polymer_entity(entity, mmcif_dict)
        else:
            annotate_non_polymer_entity(entity, mmcif_dict)
    return entities


def get_secondary_structure(mmcif_dict):
    """Creates a dictionary of helices and strands, with each having a list of
    start and end residues.

    :param mmcif_dict: the .mmcif dict to read.
    :rtype: ``dict``"""

    helices, strands = [], []
    for helix in mmcif_dict.get("struct_conf", []):
        helices.append(["{}.{}{}".format(
            helix[f"{x}_auth_asym_id"], helix[f"{x}_auth_seq_id"],
            helix[f"pdbx_{x}_PDB_ins_code"].replace("?", ""),
        ) for x in ["beg", "end"]])
    for strand in mmcif_dict.get("struct_sheet_range", []):
        strands.append(["{}.{}{}".format(
            strand[f"{x}_auth_asym_id"], strand[f"{x}_auth_seq_id"],
            strand[f"pdbx_{x}_PDB_ins_code"].replace("?", ""),
        ) for x in ["beg", "end"]])
    return {"helices": helices, "strands": strands}


def annotate_polymer_entity(entity, mmcif_dict):
    for entity_poly in mmcif_dict.get("entity_poly"):
        if entity_poly["entity_id"] == entity.id:
            entity.sequence = entity_poly["pdbx_seq_one_letter_code"]
            break
    else: entity.sequence = ""


def annotate_non_polymer_entity(entity, mmcif_dict):
    pass


def divide_atoms_into_models(mmcif_dict):
    models = []
    model_id = None
    atoms = []
    for atom in mmcif_dict["atom_site"]:
        if atom["pdbx_PDB_model_num"] != model_id:
            if atoms: models.append(atoms)
            model_id = atom["pdbx_PDB_model_num"]
            atoms = []
        atoms.append(atom)
    if atoms: models.append(atoms)
    return models


def get_structures(atoms, entities):
    asyms = sorted(set([(
        atom["label_asym_id"], atom["label_entity_id"], atom["auth_asym_id"]
    ) for atom in atoms]))
    return [{
        "id": asym[0], "auth_id": asym[2], "entity": entities[asym[1]]
    } for asym in asyms]


def organise_atoms_by_asym_id(atoms):
    d = {}
    for atom in atoms:
        try:
            d[atom["label_asym_id"]].append(atom)
        except KeyError:
            d[atom["label_asym_id"]] = [atom]
    return d


def make_chain(atoms, entities, secondary_structure=None, carb=False):
    residues = []
    atoms_ = []
    res_id = None
    for atom in atoms:
        if (atom["auth_seq_id"], atom["pdbx_PDB_ins_code"]) != res_id:
            if atoms_:
                residues.append(atoms_)
                atoms_ = []
            res_id = (atom["auth_seq_id"], atom["pdbx_PDB_ins_code"])
        atoms_.append(atom)
    if atoms_: residues.append(atoms_)
    residues = [make_comp(residue, entities) for residue in residues]
    asym_id, auth_asym_id = atoms[0]["label_asym_id"], atoms[0]["auth_asym_id"]
    entity = entities[atoms[0]["label_entity_id"]]
    if carb:
        return Carbohydrate(
            *residues, id=atoms[0]["label_asym_id"], entity=entity,
            asym_id=asym_id, auth_asym_id=auth_asym_id, 
        )
    else:
        helices = [h for h in secondary_structure["helices"]
            if h[0].split(".")[0] == auth_asym_id]
        strands = [s for s in secondary_structure["strands"]
            if s[0].split(".")[0] == auth_asym_id]
        return Chain(
            *residues, id=atoms[0]["label_asym_id"], sequence=entity.sequence,
            asym_id=asym_id, auth_asym_id=auth_asym_id, entity=entity,
            helices=helices, strands=strands
        )
        

def make_comp(atoms, entities, ligand=False):
    alt_loc = None
    if any([float(atom["occupancy"]) < 1 for atom in atoms]):
        if any([atom["label_alt_id"] for atom in atoms]):
            alt_loc = sorted([atom["label_alt_id"] for atom in atoms if atom["label_alt_id"] != "."])[0]
    ligand_atoms = [make_atom(atom) for atom in atoms if float(atom["occupancy"]) == 1 or atom["label_alt_id"] == "." or atom["label_alt_id"] == alt_loc]
    Comp = Ligand if ligand else Residue
    entity = entities[atoms[0]["label_entity_id"]]
    return Comp(
        *ligand_atoms, name=atoms[0]["label_comp_id"], id=make_id(atoms[0]),
        **({"entity": entity} if ligand else {}),
        asym_id=atoms[0]["label_asym_id"], auth_asym_id=atoms[0]["auth_asym_id"]
    )
    


def make_waters(atoms, entities):
    waters = []
    for atom in atoms:
        entity = entities[atom["label_entity_id"]]
        waters.append(Ligand(
            make_atom(atom), name=atom["label_comp_id"], id=make_id(atom),
            asym_id=atom["label_asym_id"], auth_asym_id=atom["auth_asym_id"],
            water=True, entity=entity
        ))
    return waters


def make_atom(atom):
    return Atom(
        atom["type_symbol"],
        float(atom["Cartn_x"]), float(atom["Cartn_y"]), float(atom["Cartn_z"]),
        int(atom["id"]), atom["label_atom_id"],
        charge=float(atom["pdbx_formal_charge"].replace("?", "") or "0"),
        bvalue=float(atom["B_iso_or_equiv"].replace("?", "") or "0")
    )


def make_id(atom):
    chain = atom["auth_asym_id"]
    number = atom["auth_seq_id"]
    insert = atom["PDB_ins_code"] if "PDB_ins_code" in atom else atom["pdbx_PDB_ins_code"]
    insert = insert.replace("?", "") or ""
    return f"{chain}.{number}{insert}"



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
    

    @property
    def title(self):
        return self.source.get("struct", [{}])[0].get("title", "")


    @property
    def deposition_date(self):
        datestring = self.source.get("pdbx_database_status", [{}])[0].get("recvd_initial_deposition_date", "")
        return datetime.strptime(datestring, "%Y-%m-%d").date()


    @property
    def classification(self):
        return self.source.get("struct_keywords", [{}])[0].get("pdbx_keywords", "")


    @property
    def keywords(self):
        return self.source.get("struct_keywords", [{}])[0].get("text", "").split(", ")


    @property
    def authors(self):
        return [author["name"] for author in self.source.get("citation_author", [])]


    @property
    def technique(self):
        return self.source.get("exptl", [{}])[0].get("method", "")


    @property
    def source_organism(self):
        return self.source.get("entity_src_gen", [{}])[0].get("pdbx_gene_src_scientific_name", "")


    @property
    def expression_system(self):
        return self.source.get("entity_src_gen", [{}])[0].get("pdbx_host_org_scientific_name", "")


    @property
    def resolution(self):
        value = self.source.get("refine", [{}])[0].get("ls_d_res_high", "")
        return float(value) if value else None

    @property
    def rvalue(self):
        value = self.source.get("refine", [{}])[0].get("ls_R_factor_obs", "")
        return float(value) if value else None


    @property
    def rfree(self):
        value = self.source.get("refine", [{}])[0].get("ls_R_factor_R_free", "")
        return float(value) if value else None
    

    @property
    def missing_residues(self):
        return [{
            "id": make_id(res),
            "name": res["label_comp_id"]
        } for res in self.source.get("pdbx_unobs_or_zero_occ_residues", [])]


    @property
    def assemblies(self):
        return extract_assemblies(self.source)
    

    @property
    def model(self):
        return self.models[0] if self.models else None


