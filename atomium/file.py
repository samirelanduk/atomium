from datetime import datetime
from .structures import Model, Atom, Ligand, Chain, Residue
from .assemblies import extract_assemblies

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


