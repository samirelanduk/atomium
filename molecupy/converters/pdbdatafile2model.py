"""This module handles the logic of converting a :py:class:`.PdbDataFile` to a
:py:class:`.Model`"""

from ..structures import Model, Atom, GhostAtom, SmallMolecule, Residue, Chain
from ..structures import BindSite, AlphaHelix, BetaStrand, Complex
from ..pdb.pdbdatafile import PdbDataFile
from ..pdb import residues as residues_dict

def model_from_pdb_data_file(data_file, model_id=1):
    """Takes a :py:class:`.PdbDataFile`, converts it to a :py:class:`.Model`,
    and returns it.

    :py:class:`.PdbDataFile` objects can contain multiple models. By default,
    model 1 will be used, but you can specify specific models with the
    ``model_id`` argument.

    :param PdbDataFile data_file: The :py:class:`.PdbDataFile` to convert.
    :param int model_id: The ID of the model in the data fileto be used for\
    conversion.
    :rtype: :py:class:`.Model`"""

    if not isinstance(data_file, PdbDataFile):
        raise TypeError("model_from_pdb_data_file can only convert PdbDataFiles")
    for model_dict in data_file.models():
        if model_dict["model_id"] == model_id:
            model = Model()
            model._source = data_file

            add_small_molecules_to_model(model, data_file, model_id)
            add_chains_to_model(model, data_file, model_id)
            connect_atoms(model, data_file, model_id)
            bond_residue_atoms(model, data_file, model_id)
            bond_residues_together(model, data_file, model_id)
            make_disulphide_bonds(model, data_file, model_id)
            make_link_bonds(model, data_file, model_id)
            give_model_sites(model, data_file, model_id)
            map_sites_to_ligands(model, data_file, model_id)
            give_model_alpha_helices(model, data_file, model_id)
            give_model_beta_strands(model, data_file, model_id)
            give_model_complexes(model, data_file, model_id)

            return model
    raise ValueError("There is no model with ID %i" % model_id)


def add_small_molecules_to_model(model, data_file, model_id):
    """Takes a :py:class:`.Model` and creates :py:class:`.SmallMolecule`
    objects in it based on the ``heteroatoms`` in the provided
    :py:class:`.PdbDataFile`.

    :param Model model: the model to update.
    :param PdbDataFile data_file: The source Pdb Data File
    :param int model_id: The ID of the model in the data fileto be used for\
    conversion."""

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
    """Takes a :py:class:`.Model` and creates :py:class:`.Chain`
    objects in it based on the ``atoms`` in the provided
    :py:class:`.PdbDataFile`.

    :param Model model: the model to update.
    :param PdbDataFile data_file: The source Pdb Data File
    :param int model_id: The ID of the model in the data fileto be used for\
    conversion."""

    atoms = [a for a in data_file.atoms() if a["model_id"] == model_id]
    heteroatoms = [a for a in data_file.heteroatoms() if a["model_id"] == model_id]
    highest_id = _get_top_atom_id(atoms, heteroatoms)
    chain_ids = set([a["chain_id"] for a in atoms])
    for chain_id in sorted(chain_ids):
        residues = []
        _add_residues_to_chain(residues, chain_id, atoms)
        missing_residue_info = _get_missing_residue_info(data_file, chain_id)
        if missing_residue_info:
            _add_missing_residues_to_chain(residues, missing_residue_info, highest_id)
        chain = Chain(chain_id, *residues)
        highest_id = max(
         [atom.atom_id() for atom in chain.atoms(atom_type="all")]
        )
        model.add_chain(chain)


def connect_atoms(model, data_file, model_id):
    """Takes a :py:class:`.Model` and creates :py:class:`.Bond`
    objects between atoms in it based on the ``connections`` in the provided
    :py:class:`.PdbDataFile`.

    :param Model model: the model to update.
    :param PdbDataFile data_file: The source Pdb Data File
    :param int model_id: The ID of the model in the data fileto be used for\
    conversion."""

    for connection in data_file.connections():
        atom = model.get_atom_by_id(connection["atom_id"])
        if atom:
            for bonded_atom_id in connection["bonded_atoms"]:
                bonded_atom = model.get_atom_by_id(bonded_atom_id)
                if bonded_atom and atom is not bonded_atom:
                    atom.bond_to(bonded_atom)


def bond_residue_atoms(model, data_file, model_id):
    """Takes a :py:class:`.Model` and creates :py:class:`.Bond`
    objects within the residues of the Model, based on a pre-defined
    dictionary of how residues are connected internally.

    :param Model model: the model to update.
    :param PdbDataFile data_file: The source Pdb Data File
    :param int model_id: The ID of the model in the data fileto be used for\
    conversion."""

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


def bond_residues_together(model, data_file, model_id):
    """Takes a :py:class:`.Model` and creates :py:class:`.Bond`
    objects between the residues of chains in the model.

    :param Model model: the model to update.
    :param PdbDataFile data_file: The source Pdb Data File
    :param int model_id: The ID of the model in the data fileto be used for\
    conversion."""

    for chain in model.chains():
        for index, residue in enumerate(chain.residues()[:-1]):
            next_residue = chain.residues()[index + 1]
            residue.connect_to(next_residue)
            carboxy_atom = residue.get_atom_by_name("C",)
            amino_nitrogen = next_residue.get_atom_by_name("N")
            if carboxy_atom and amino_nitrogen:
                carboxy_atom.bond_to(amino_nitrogen)


def make_disulphide_bonds(model, data_file, model_id):
    """Takes a :py:class:`.Model` and creates disulphide :py:class:`.Bond`
    objects in it based on the ``ss_bonds`` in the provided
    :py:class:`.PdbDataFile`.

    :param Model model: the model to update.
    :param PdbDataFile data_file: The source Pdb Data File
    :param int model_id: The ID of the model in the data fileto be used for\
    conversion."""

    for disulphide_bond in data_file.ss_bonds():
        chain1 = model.get_chain_by_id(disulphide_bond["chain_id_1"])
        chain2 = model.get_chain_by_id(disulphide_bond["chain_id_2"])
        if chain1 and chain2:
            residue1 = chain1.get_residue_by_id(
             disulphide_bond["chain_id_1"] +
             str(disulphide_bond["residue_id_1"]) +
             disulphide_bond["insert_code_1"]
            )
            residue2 = chain2.get_residue_by_id(
             disulphide_bond["chain_id_2"] +
             str(disulphide_bond["residue_id_2"]) +
             disulphide_bond["insert_code_2"]
            )
            if residue1 and residue2:
                atom1 = residue1.get_atom_by_element("S")
                atom2 = residue2.get_atom_by_element("S")
                if atom1 and atom2 and atom1 is not atom2:
                    atom1.bond_to(atom2)


def make_link_bonds(model, data_file, model_id):
    """Takes a :py:class:`.Model` and creates specified :py:class:`.Bond`
    objects in it based on the ``links`` in the provided
    :py:class:`.PdbDataFile`.

    :param Model model: the model to update.
    :param PdbDataFile data_file: The source Pdb Data File
    :param int model_id: The ID of the model in the data fileto be used for\
    conversion."""

    for link in data_file.links():
        chain1 = model.get_chain_by_id(link["chain_id_1"])
        chain2 = model.get_chain_by_id(link["chain_id_2"])
        if chain1 and chain2:
            molecule1, molecule2 = [chain1.get_residue_by_id(
             link["chain_id_%i" % n] +
             str(link["residue_id_%i" % n]) +
             link["insert_code_%i" % n]
            ) for n in (1, 2)]
            if not molecule1:
                molecule1 = model.get_small_molecule_by_id(
                 link["chain_id_1"] +
                 str(link["residue_id_1"]) +
                 link["insert_code_1"]
                )
            if not molecule2:
                molecule2 = model.get_small_molecule_by_id(
                 link["chain_id_2"] +
                 str(link["residue_id_2"]) +
                 link["insert_code_2"]
                )
            if molecule1 and molecule2:
                atom1 = molecule1.get_atom_by_name(link["atom_1"])
                atom2 = molecule2.get_atom_by_name(link["atom_2"])
                if atom1 and atom2:
                    atom1.bond_to(atom2)


def give_model_sites(model, data_file, model_id):
    """Takes a :py:class:`.Model` and creates :py:class:`.BindSite`
    objects in it based on the ``sites`` in the provided
    :py:class:`.PdbDataFile`.

    :param Model model: the model to update.
    :param PdbDataFile data_file: The source Pdb Data File
    :param int model_id: The ID of the model in the data fileto be used for\
    conversion."""

    for site in data_file.sites():
        residues = [model.get_chain_by_id(residue["chain_id"]).get_residue_by_id(
         str(residue["chain_id"]) + str(residue["residue_id"]) + residue["insert_code"]
        ) if residue["chain_id"] in [
         chain.chain_id() for chain in model.chains()
        ] else None for residue in site["residues"]]
        residues = [r for r in residues if r]
        if residues:
            site = BindSite(
             site["site_id"],
             *residues
            )
            model.add_bind_site(site)


def map_sites_to_ligands(model, data_file, model_id):
    """Takes a :py:class:`.Model` and assocated ligands and binding sites to
    each other based on 800-remarks in the provided
    :py:class:`.PdbDataFile`.

    :param Model model: the model to update.
    :param PdbDataFile data_file: The source Pdb Data File
    :param int model_id: The ID of the model in the data fileto be used for\
    conversion."""

    remark_800 = data_file.get_remark_by_number(800)
    if remark_800:
        remark_lines = [
         line for line in remark_800["content"].split("\n") if line != "SITE"
        ]
        for index, line in enumerate(remark_lines):
            if line.startswith("SITE_IDENTIFIER"):
                site_id = line.split(":")[1].strip() if ":" in line else None
                if site_id:
                    for trailing_line in remark_lines[index+1:]:
                        if trailing_line.startswith("SITE_IDENTIFIER"): break
                        if trailing_line.startswith("SITE_DESCRIPTION"):
                            site = model.get_bind_site_by_id(site_id)
                            if site:
                                ligand_id = trailing_line.split()[-1]
                                if not ligand_id[0].isalpha():
                                    ligand_id = trailing_line.split()[-2] + ligand_id
                                ligand = model.get_small_molecule_by_id(ligand_id)
                                site.ligand(ligand)


def give_model_alpha_helices(model, data_file, model_id):
    """Takes a :py:class:`.Model` and creates :py:class:`.AlphaHelix`
    objects in it based on the ``helices`` in the provided
    :py:class:`.PdbDataFile`.

    :param Model model: the model to update.
    :param PdbDataFile data_file: The source Pdb Data File
    :param int model_id: The ID of the model in the data fileto be used for\
    conversion."""

    for helix in data_file.helices():
        chain = model.get_chain_by_id(helix["start_residue_chain_id"])
        if chain:
            start_residue = chain.get_residue_by_id(
             helix["start_residue_chain_id"] +
             str(helix["start_residue_id"]) +
             helix["start_residue_insert"]
            )
            end_residue = chain.get_residue_by_id(
             helix["end_residue_chain_id"] +
             str(helix["end_residue_id"]) +
             helix["end_residue_insert"]
            )
            if start_residue and end_residue:
                start_index = chain.residues().index(start_residue)
                end_index = chain.residues().index(end_residue)
                if end_index > start_index:
                    AlphaHelix(
                     helix["helix_name"],
                     *chain.residues()[start_index:end_index + 1],
                     comment=helix["comment"],
                     helix_class=HELIX_CLASSES.get(helix["helix_class"], HELIX_CLASSES[1])
                    )


def give_model_beta_strands(model, data_file, model_id):
    """Takes a :py:class:`.Model` and creates :py:class:`.BetaStrand`
    objects in it based on the ``sheets`` in the provided
    :py:class:`.PdbDataFile`.

    :param Model model: the model to update.
    :param PdbDataFile data_file: The source Pdb Data File
    :param int model_id: The ID of the model in the data fileto be used for\
    conversion."""

    for sheet in data_file.sheets():
        for strand in sheet["strands"]:
            chain = model.get_chain_by_id(strand["start_residue_chain_id"])
            if chain:
                start_residue = chain.get_residue_by_id(
                 strand["start_residue_chain_id"] +
                 str(strand["start_residue_id"]) +
                 strand["start_residue_insert"]
                )
                end_residue = chain.get_residue_by_id(
                 strand["end_residue_chain_id"] +
                 str(strand["end_residue_id"]) +
                 strand["end_residue_insert"]
                )
                if start_residue and end_residue:
                    start_index = chain.residues().index(start_residue)
                    end_index = chain.residues().index(end_residue)
                    BetaStrand(
                     str(strand["strand_id"]),
                     strand["sense"],
                     *chain.residues()[start_index:end_index + 1]
                    )


def give_model_complexes(model, data_file, model_id):
    """Takes a :py:class:`.Model` and creates :py:class:`.Complex`
    objects in it based on the ``compounds`` in the provided
    :py:class:`.PdbDataFile`.

    :param Model model: the model to update.
    :param PdbDataFile data_file: The source Pdb Data File
    :param int model_id: The ID of the model in the data fileto be used for\
    conversion."""
    
    for compound in data_file.compounds():
        chains = []
        for chain in model.chains():
            if chain.chain_id() in compound["CHAIN"]:
                chains.append(chain)
        if chains:
            complex_ = Complex(str(compound["MOL_ID"]), compound["MOLECULE"], *chains)
            model.add_complex(complex_)


def _get_top_atom_id(atoms=None, heteroatoms=None):
    atom_id = max([atom["atom_id"] for atom in atoms]) if atoms else -1
    hetero_id = max([atom["atom_id"] for atom in heteroatoms]) if heteroatoms else -1
    return max((atom_id, hetero_id))


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


HELIX_CLASSES = {
 1: "Right-handed alpha",
 2: "Right-handed omega",
 3: "Right-handed pi",
 4: "Right-handed gamma",
 5: "Right-handed 3 - 10",
 6: "Left-handed alpha",
 7: "Left-handed omega",
 8: "Left-handed gamma",
 9: "2 - 7 ribbon/helix",
 10: "Polyproline",
}
