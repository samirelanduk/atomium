"""This module contains creates the final Pdb object itself, and processes the
data contained in the data file."""

from ..structures import Model, Atom, GhostAtom, SmallMolecule, Residue, Chain
from ..structures import BindSite, AlphaHelix, BetaStrand, Complex
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
            _give_model_complexes(model, data_file, model_dict["model_id"])
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


def _give_model_complexes(model, data_file, model_id):
    for compound in data_file.compounds():
        chains = []
        for chain in model.chains():
            if chain.chain_id() in compound["CHAIN"]:
                chains.append(chain)
        if chains:
            complex_ = Complex(str(compound["MOL_ID"]), compound["MOLECULE"], *chains)
            model.add_complex(complex_)


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
