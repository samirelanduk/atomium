"""This module contains creates the final Pdb object itself, and processes the
data contained in the data file."""

from .structures import Model, PdbAtom, Atom, SmallMolecule, Residue, Chain
from .structures import BindSite, AlphaHelix, BetaStrand
from . import residues as residues_dict

class Pdb:
    """A representation of a PDB file and its contents, including the structure.

    :param PdbDataFile data_file: The PDB data file with the parsed values."""

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
            _make_disulphide_bonds(model, data_file, model_dict["model_id"])
            _make_link_bonds(model, data_file, model_dict["model_id"])
            _give_model_sites(model, data_file, model_dict["model_id"])
            _map_sites_to_ligands(model, data_file, model_dict["model_id"])
            _give_model_alpha_helices(model, data_file, model_dict["model_id"])
            _give_model_beta_strands(model, data_file, model_dict["model_id"])
            self._models.append(model)


    def __repr__(self):
        return "<Pdb (%s)>" % (self.pdb_code() if self.pdb_code() else "????")


    def data_file(self):
        """The :py:class:`.PdbDataFile` from which the object was created.

        :rtype: ``PdbDataFile``"""

        return self._data_file


    def classification(self):
        """The PDB classification.

        :rtype: ``str``"""

        return self._data_file.classification()


    def deposition_date(self):
        """The date the PDB was deposited.

        :rtype: ``datetime.Date``"""

        return self._data_file.deposition_date()


    def pdb_code(self):
        """The PDB four-letter code.

        :rtype: ``str``"""

        return self._data_file.pdb_code()


    def is_obsolete(self):
        """``True`` if the PDB has been made obsolete by a newer PDB.

        :rtype: ``bool``"""

        return self._data_file.is_obsolete()


    def obsolete_date(self):
        """The date the PDB was made obsolete.

        :rtype: ``datetime.Date``"""

        return self._data_file.obsolete_date()


    def replacement_code(self):
        """The PDB code of the replacing PDB.

        :rtype: ``str``"""

        return self._data_file.replacement_code()


    def title(self):
        """The title of the PDB.

        :rtype: ``str``"""

        return self._data_file.title()


    def split_codes(self):
        """The PDB codes which complete this structure.

        :rtype: ``list``"""

        return self._data_file.split_codes()


    def caveat(self):
        """Any caveats for this structure.

        :rtype: ``str``"""

        return self._data_file.caveat()


    def keywords(self):
        """Keywords for this PDB.

        :rtype: ``list``"""

        return self._data_file.keywords()


    def experimental_techniques(self):
        """The experimental techniques used to produce this PDB.

        :rtype: ``list``"""

        return self._data_file.experimental_techniques()


    def model_count(self):
        """The number of models in this PDB.

        :rtype: ``int``"""

        return self._data_file.model_count()


    def model_annotations(self):
        """Annotations for the PDB's models.

        :rtype: ``list``"""

        return self._data_file.model_annotations()


    def authors(self):
        """The PDB's authors.

        :rtype: ``list``"""

        return self._data_file.authors()


    def revisions(self):
        """Any changes made to the PDB file.

        :rtype: ``list``"""

        return self._data_file.revisions()


    def supercedes(self):
        """The PDB codes that this PDB replaces.

        :rtype: ``list``"""

        return self._data_file.supercedes()


    def supercede_date(self):
        """The date this PDB replaced another.

        :rtype: ``datetime.Date``"""

        return self._data_file.supercede_date()


    def journal(self):
        """The publication information for this PDB.

        :rtype: ``dict``"""

        return self._data_file.journal()


    def models(self):
        """The PDB's models.

        :rtype: ``list``"""

        return self._models


    def model(self):
        """The first Model in the PDB models.

        :rtype: ``Model``"""

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
                if bonded_atom and atom is not bonded_atom:
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


def _make_disulphide_bonds(model, data_file, model_id):
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


def _make_link_bonds(model, data_file, model_id):
    for link in data_file.links():
        chain1 = model.get_chain_by_id(link["chain_id_1"])
        chain2 = model.get_chain_by_id(link["chain_id_2"])
        if chain1 and chain2:
            molecule1 = chain1.get_residue_by_id(
             link["chain_id_1"] +
             str(link["residue_id_1"]) +
             link["insert_code_1"]
            )
            molecule2 = chain2.get_residue_by_id(
             link["chain_id_2"] +
             str(link["residue_id_2"]) +
             link["insert_code_2"]
            )
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


def _give_model_sites(model, data_file, model_id):
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


def _map_sites_to_ligands(model, data_file, model_id):
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


def _give_model_alpha_helices(model, data_file, model_id):
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


def _give_model_beta_strands(model, data_file, model_id):
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
