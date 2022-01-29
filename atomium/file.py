from datetime import datetime
from .structures import Model, Polymer, BranchedPolymer, NonPolymer, Water, Residue, Atom
from .assemblies import extract_assemblies

class File:
    
    def __init__(self, mmcif_dict):
        self.source = mmcif_dict
        self.name = mmcif_dict.get("entry", [{}])[0].get("id")
        self.models = make_models(mmcif_dict, self)
    

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
        classification = self.source.get("struct_keywords", [{}])[0].get("pdbx_keywords", "")
        return None if classification == "?" else classification


    @property
    def keywords(self):
        return self.source.get("struct_keywords", [{}])[0].get("text", "").split(", ")


    @property
    def authors(self):
        return [author["name"] for author in self.source.get("citation_author", [])]


    @property
    def technique(self):
        technique = self.source.get("exptl", [{}])[0].get("method", "")
        return None if technique == "?" else technique


    @property
    def source_organism(self):
        source_organism = self.source.get("entity_src_gen", [{}])[0].get("pdbx_gene_src_scientific_name", "")
        return None if source_organism == "?" else source_organism


    @property
    def expression_system(self):
        system = self.source.get("entity_src_gen", [{}])[0].get("pdbx_host_org_scientific_name", "")
        return None if system == "?" else system


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




def make_models(mmcif_dict, file):
    entities = make_entities(mmcif_dict)
    secondary_structure = get_secondary_structure(mmcif_dict)
    aniso = make_aniso(mmcif_dict)
    atom_lists = divide_atoms_into_models(mmcif_dict, aniso)
    models = []
    for model_atoms in atom_lists:
        asyms = organise_atoms_by_asym_id(model_atoms)
        molecules = [make_molecule(
            asym_id, atoms, entities, secondary_structure
        ) for asym_id, atoms in asyms.items()]
        models.append(Model(molecules=molecules, file=file))
    return models


def make_entities(mmcif_dict):
    entities = {}
    for row in mmcif_dict["entity"]:
        Parent = {
            "polymer": Polymer, "branched": BranchedPolymer,
            "non-polymer": NonPolymer, "water": Water
        }[row["type"]]
        cls = type(row["pdbx_description"], (Parent,), {
            "entity_id": row["id"], "entity_name": row["pdbx_description"]
        })
        if row["type"] == "polymer":
            annotate_polymer_entity(cls, mmcif_dict)
        entities[row["id"]] = cls
    return entities


def annotate_polymer_entity(Polymer, mmcif_dict):
    for entity_poly in mmcif_dict.get("entity_poly"):
        if entity_poly["entity_id"] == Polymer.entity_id:
            Polymer.sequence = entity_poly["pdbx_seq_one_letter_code"]
            break
    else: Polymer.sequence = ""


def get_secondary_structure(mmcif_dict):
    """Creates a dictionary of helices and strands, with each having a list of
    start and end residues.

    :param mmcif_dict: the .mmcif dict to read.
    :rtype: ``dict``"""

    helices, strands = [], []
    for helix in mmcif_dict.get("struct_conf", []):
        helices.append(["{}.{}{}".format(
            helix[f"{x}_label_asym_id"], helix[f"{x}_label_seq_id"],
            helix[f"pdbx_{x}_PDB_ins_code"].replace("?", ""),
        ) for x in ["beg", "end"]])
    for strand in mmcif_dict.get("struct_sheet_range", []):
        strands.append(["{}.{}{}".format(
            strand[f"{x}_label_asym_id"], strand[f"{x}_label_seq_id"],
            strand[f"pdbx_{x}_PDB_ins_code"].replace("?", ""),
        ) for x in ["beg", "end"]])
    return {"helices": helices, "strands": strands}


def make_aniso(mmcif_dict):
    """Makes a mapping of atom IDs to anisotropy information.

    :param mmcif_dict: the .mmcif dict to read.
    :rtype: ``dict``"""

    return {a["id"]: [
        float(a["U[{}][{}]".format(x, y)]) for x, y in ["11", "22", "33", "12", "13", "23"]
    ] for a in mmcif_dict.get("atom_site_anisotrop", [])}


def divide_atoms_into_models(mmcif_dict, aniso):
    models = []
    model_id = None
    atoms = []
    for atom in mmcif_dict["atom_site"]:
        atom["anisotropy"] = aniso.get(atom["id"])
        if atom["pdbx_PDB_model_num"] != model_id:
            if atoms: models.append(atoms)
            model_id = atom["pdbx_PDB_model_num"]
            atoms = []
        atoms.append(atom)
    if atoms: models.append(atoms)
    return models


def organise_atoms_by_asym_id(atoms):
    d = {}
    for atom in atoms:
        try:
            d[atom["label_asym_id"]].append(atom)
        except KeyError:
            d[atom["label_asym_id"]] = [atom]
    return d


def make_molecule(asym_id, atoms, entities, secondary_structure):
    entity_id = atoms[0]["label_entity_id"]
    auth_id = atoms[0]["auth_asym_id"]
    Molecule = entities[entity_id]
    if Molecule.__bases__[0] in (Polymer, BranchedPolymer):
        residues = make_residues(atoms)
        if Molecule.__bases__[0] is Polymer:
            helices = [h for h in secondary_structure["helices"]
                if h[0].split(".")[0] == asym_id]
            strands = [s for s in secondary_structure["strands"]
                if s[0].split(".")[0] == asym_id]
            return Molecule(
                id=asym_id, auth_id=auth_id, residues=residues,
                helices=helices, strands=strands
            )
        else:
            return Molecule(id=asym_id, auth_id=auth_id, residues=residues)
    elif Molecule.__bases__[0] in (NonPolymer, Water):
        comp_id = atoms[0]["label_comp_id"]
        return Molecule(id=asym_id, auth_id=auth_id, name=comp_id, atoms=[
            make_atom(a) for a in remove_alt_loc_atoms(atoms)
        ])


def make_residues(atoms):
    residues = []
    res_id = None
    atoms_ = []
    for atom in atoms:
        if (atom["auth_seq_id"], atom["pdbx_PDB_ins_code"]) != res_id:
            if atoms_:
                residues.append(atoms_)
                atoms_ = []
            res_id = (atom["auth_seq_id"], atom["pdbx_PDB_ins_code"])
        atoms_.append(atom)
    if atoms_: residues.append(atoms_)
    return [Residue(
        id=make_id(res_atoms[0]),
        name=res_atoms[0]["label_comp_id"],
        number=int(
            res_atoms[0]["label_seq_id"] if res_atoms[0]["label_seq_id"] != "."
            else res_atoms[0]["auth_seq_id"]
        ),
        atoms=[make_atom(a) for a in remove_alt_loc_atoms(res_atoms)],
    ) for res_atoms in residues]
    

def remove_alt_loc_atoms(atoms):
    alt_loc = None
    if any([atom["occupancy"] != "." and float(atom["occupancy"]) < 1 for atom in atoms]):
        if any([atom["label_alt_id"] and atom["label_alt_id"] != "." for atom in atoms]):
            alt_loc = sorted([atom["label_alt_id"] for atom in atoms if atom["label_alt_id"] != "."])[0]
    return [atom for atom in atoms if atom["occupancy"] == "." or float(atom["occupancy"]) == 1 or atom["label_alt_id"] == "." or atom["label_alt_id"] == alt_loc]


def make_atom(atom):
    return Atom(
        atom["type_symbol"],
        float(atom["Cartn_x"]), float(atom["Cartn_y"]), float(atom["Cartn_z"]),
        int(atom["id"]), atom["label_atom_id"],
        charge=float(atom["pdbx_formal_charge"].replace("?", "") or "0"),
        bvalue=float(atom["B_iso_or_equiv"].replace("?", "") or "0"),
        anisotropy=atom["anisotropy"]
    )


def make_id(atom):
    chain = atom["label_asym_id"]
    number = atom["label_seq_id"] if atom["label_seq_id"] != "." else atom["auth_seq_id"]
    insert = atom["PDB_ins_code"] if "PDB_ins_code" in atom else atom["pdbx_PDB_ins_code"]
    insert = insert.replace("?", "") or ""
    return f"{chain}.{number}{insert}"


