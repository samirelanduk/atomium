"""This module handles the logic of converting a :py:class:`.Model` to a
:py:class:`.PdbDataFile`"""

from ..pdb.pdbdatafile import PdbDataFile
from ..structures.molecules import Residue, SmallMolecule
from ..structures.models import Model

def pdb_data_file_from_model(model):
    """Takes a :py:class:`.Model`, converts it to a :py:class:`.PdbdataFile`,
    and returns it.

    :param Model model: The :py:class:`.Model` to convert.
    :rtype: :py:class:`.PdbDataFile`"""

    if not isinstance(model, Model):
        raise TypeError("pdb_data_file_from_model can only convert Models")
    data_file = PdbDataFile()
    data_file._source = model

    add_complexes_to_data_file(data_file, model)
    add_atoms_to_data_file(data_file, model)
    add_connections_to_data_file(data_file, model)

    return data_file


def add_complexes_to_data_file(data_file, model):
    """Takes a :py:class:`.PdbDataFile` and updates its compounds based on
    the complexes in the provided :py:class:`.Model`

    :param PdbDataFile data_file: the PDB Data File to update.
    :param PdbDataFile model: The source Model"""

    for complex_ in sorted(list(model.complexes()), key=lambda k: k.complex_id()):
        data_file.compounds().append({
         "MOL_ID": int(complex_.complex_id()),
         "MOLECULE": complex_.complex_name(),
         "CHAIN": sorted([chain.chain_id() for chain in complex_.chains()])
        })


def add_atoms_to_data_file(data_file, model):
    """Takes a :py:class:`.PdbDataFile` and updates its atoms and heteroatoms
    based on the atoms in the provided :py:class:`.Model`

    :param PdbDataFile data_file: the PDB Data File to update.
    :param PdbDataFile model: The source Model"""

    for atom in sorted(list(model.atoms()), key=lambda k: k.atom_id()):
        residue_name, chain_id, residue_id, insert = None, None, None, None
        if atom.molecule():
            if isinstance(atom.molecule(), Residue):
                residue_name = atom.molecule().residue_name()
                residue_id = atom.molecule().residue_id()[1:]
                chain_id = atom.molecule().residue_id()[0]
                if atom.molecule().residue_id()[-1].isalpha():
                    insert = atom.molecule().residue_id()[-1]
                    residue_id = atom.molecule().residue_id()[1:-1]
            elif isinstance(atom.molecule(), SmallMolecule):
                residue_name = atom.molecule().molecule_name()
                residue_id = atom.molecule().molecule_id()[1:]
                chain_id = atom.molecule().molecule_id()[0]
                if atom.molecule().molecule_id()[-1].isalpha():
                    insert = atom.molecule().molecule_id()[-1]
                    residue_id = atom.molecule().molecule_id()[1:-1]
            residue_id = int(residue_id)
        atom_dict = {
         "atom_id": atom.atom_id(),
         "atom_name": atom.atom_name(),
         "alt_loc": None,
         "residue_name": residue_name,
         "chain_id": chain_id,
         "residue_id": residue_id,
         "insert_code": insert,
         "x": atom.x(),
         "y": atom.y(),
         "z": atom.z(),
         "occupancy": 1.0,
         "temperature_factor": 0.0,
         "element": atom.element(),
         "charge": None,
         "model_id": 1
        }
        if atom.molecule() and isinstance(atom.molecule(), Residue):
            data_file.atoms().append(atom_dict)
        else:
            data_file.heteroatoms().append(atom_dict)


def add_connections_to_data_file(data_file, model):
    """Takes a :py:class:`.PdbDataFile` and updates its connections based on
    the bonds in the provided :py:class:`.Model`

    :param PdbDataFile data_file: the PDB Data File to update.
    :param PdbDataFile model: The source Model"""

    for molecule in sorted(list(model.small_molecules()), key=lambda k: k.molecule_id()):
        for atom in sorted(list(molecule.atoms()), key=lambda k: k.atom_id()):
            other_atoms = sorted(list(atom.bonded_atoms()), key=lambda k: k.atom_id())
            connection = {
             "atom_id": atom.atom_id(),
             "bonded_atoms": [atom.atom_id() for atom in other_atoms]
            }
            data_file.connections().append(connection)
