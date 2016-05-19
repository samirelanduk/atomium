"""This module contains creates the final Pdb object itself, and processes the
data contained in the data file."""

from .structures import PdbModel, PdbAtom, PdbSmallMolecule, PdbResidue, PdbChain, PdbSite
from .exceptions import *
from pprint import pprint

class Pdb:
    """A representation of a PDB file and its contents, including the structure.

    :param PdbDataFile data_file: The PDB data file with the parsed values.

    .. py:attribute:: data_file:

        The :py:class:`.PdbDataFile` from which the object was created.

    .. py:attribute:: classification:

        The PDB classification.

    .. py:attribute:: deposition_date:

        The date the PDB file was submitted.

    .. py:attribute:: pdb_code:

        The PDB's four character identifier.

    .. py:attribute:: is_obsolete:

        Returns ``True`` if the PDB has been obsoleted by another.

    .. py:attribute:: obsolete_date:

        The date on which the PDB was made obsolete.

    .. py:attribute:: replacement_code:

        The PDB which made this PDB obsolete.

    .. py:attribute:: title:

        The PDB's title.

    .. py:attribute:: split_codes:

        The other PDB codes which complete this structure.

    .. py:attribute:: caveat:

        The caveats which the PDB might have.

    .. py:attribute:: keywords:

        Keywords for this PDB.

    .. py:attribute:: experimental_techniques:

        The experimental technique(s) used to generate the PDB's data.

    .. py:attribute:: model_count:

        The number of models the PDB has.

    .. py:attribute:: model_annotations:

        Annotations for the PDB's models.

    .. py:attribute:: authors:

        The PDB's authors.

    .. py:attribute:: revisions:

        Any revisions made to the PDB after submission.

    .. py:attribute:: supercedes:

        The PDB(s) which this PDB replaced.

    .. py:attribute:: supercede_date:

        The date on which this PDB replaced another.

    .. py:attribute:: journal:

        The publication corresponding to this PDB as a ``dict``.

    .. py:attribute:: models:

        A list of :py:class:`.PdbModel` objects for this PDB.

    .. py:attribute:: model:

        The first :py:class:`.PdbModel` of this PDB."""

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
         "model_count",
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
        for model_dict in self.data_file.models:
            model = PdbModel()
            _give_model_small_molecules(model, self.data_file, model_dict["model_id"])
            _give_model_chains(model, self.data_file, model_dict["model_id"])
            _give_model_sites(model, self.data_file, model_dict["model_id"])
            _map_sites_to_ligands(model, self.data_file, model_dict["model_id"])
            self.models.append(model)
        self.model = self.models[0]


    def __repr__(self):
        return "<Pdb (%s)>" % (self.pdb_code if self.pdb_code else "????")



def _give_model_small_molecules(model, data_file, model_id):
    small_molecule_names = set([a["residue_name"] for a in data_file.heteroatoms])
    small_molecule_ids = set([a["residue_id"] for a in data_file.heteroatoms])
    for molecule_name in small_molecule_names:
        for molecule_id in small_molecule_ids:
            relevant_atoms = [a for a in data_file.heteroatoms
             if a["model_id"] == model_id and a["residue_name"] == molecule_name
              and a["residue_id"] == molecule_id]
            if relevant_atoms:
                chain_id = relevant_atoms[0]["chain_id"]
                chain_id = chain_id if chain_id else ""
                atoms = [PdbAtom(
                 a["x"], a["y"], a["z"],
                 a["element"],
                 a["atom_id"],
                 a["atom_name"]
                ) for a in relevant_atoms]
                small_molecule = PdbSmallMolecule(
                 "%s%i" % (chain_id, molecule_id), molecule_name, *atoms
                )
                model.add_small_molecule(small_molecule)


def _give_model_chains(model, data_file, model_id):
    chain_ids = set([a["chain_id"] for a in data_file.atoms])
    for chain_id in chain_ids:
        relevant_atoms = [a for a in data_file.atoms
         if a["model_id"] == model_id and a["chain_id"] == chain_id]
        residue_ids = set([str(a["residue_id"]) + a["insert_code"] for a in relevant_atoms])
        residues = []
        for residue_id in residue_ids:
            residue_atoms = [
             a for a in relevant_atoms if str(a["residue_id"]) + a["insert_code"] == residue_id
            ]
            atoms = [PdbAtom(
             a["x"], a["y"], a["z"],
             a["element"],
             a["atom_id"],
             a["atom_name"]
            ) for a in residue_atoms]
            residue = PdbResidue(
             chain_id + residue_id, residue_atoms[0]["residue_name"],
             *atoms
            )
            residues.append(residue)
        chain = PdbChain(chain_id, *residues)
        model.add_chain(chain)


def _give_model_sites(model, data_file, model_id):
    for site in data_file.sites:
        residues = [model.get_chain_by_id(residue["chain_id"]).get_residue_by_id(
         str(residue["chain_id"]) + str(residue["residue_id"]) + residue["insert_code"]
        ) for residue in site["residues"]]
        residues = [r for r in residues if isinstance(r, PdbResidue)]
        if residues:
            site = PdbSite(
             site["site_id"],
             *[r for r in residues if r]
            )
            model.add_site(site)


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
                            site = model.get_site_by_id(site_id)
                            if site:
                                ligand_id = "".join(trailing_line.split()[-2:])
                                ligand = model.get_small_molecule_by_id(ligand_id)
                                site.ligand = ligand
