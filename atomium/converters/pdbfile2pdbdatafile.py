"""Contains the function for creating PdbDataFiles from PdbFiles."""

from atomium.files.pdbdatafile import PdbDataFile

def pdb_file_to_pdb_data_file(pdb_file):
    """Converts a :py:class:`.PdbFile` to a :py:class:`.PdbDataFile`

    :param PdbFile pdb_file: The PdbFile.
    :rtype: ``PdbDataFile``"""

    data_file = PdbDataFile()
    extract_structure(pdb_file, data_file)
    return data_file


def extract_structure(pdb_file, data_file):
    """Takes a :py:class:`.PdbFile` and a  :py:class:`.PdbDataFile`, and updates
    the latter with ATOM, HETATM and CONECT information from the former.

    :param PdbFile pdb_file: The source PdbFile.
    :param PdbDataFile data_file: The PdbDataFile to update."""

    model_records = pdb_file.records(name="model")
    atom_records = pdb_file.records(name="atom")
    hetatm_records = pdb_file.records(name="hetatm")
    conect_records = pdb_file.records(name="conect")
    data_file.atoms = [
     atom_record_to_dict(rec, model_records) for rec in atom_records
    ]
    data_file.heteroatoms = [
     atom_record_to_dict(rec, model_records) for rec in hetatm_records
    ]
    data_file.connections = conect_records_to_list(conect_records)


def atom_record_to_dict(record, model_records):
    """Converts an ATOM or HETATM :py:class:`.PdbRecord` to a ``dict``.

    :param PdbRecord record: The record to parse.
    :param list model_records: The model records to use as number reference.
    :rtype: ``dict``"""

    charge = record[78:80]
    if isinstance(charge, str): charge = int(charge[::-1])
    model_locations = [rec.number() for rec in model_records]
    model_number, this_number = 1, record.number()
    for index, loc in enumerate(model_locations, start=1):
        if this_number > loc: model_number = index
    return {
     "atom_id": record[6:11],
     "atom_name": record[12:16],
     "alt_loc": record[16],
     "residue_name": record[17:20],
     "chain_id": record[21],
     "residue_id": record[22:26],
     "insert_code": record[26],
     "x": record[30:38],
     "y": record[38:46],
     "z": record[46:54],
     "occupancy": record[54:60],
     "temperature_factor": record[60:66],
     "element": record[76:78],
     "charge": charge,
     "model": model_number
    }


def conect_records_to_list(records):
    """Converts a sequence of CONECT :py:class:`.PdbRecord` objects to a list
    of connections.

    :param records: The records to parse.
    :rtype: ``list``"""

    atom_ids = sorted(list(set([r[6:11] for r in records])))
    return [{
     "atom": num,
     "bond_to": [int(n) for n in merge_records([
      r for r in records if r[6:11] == num
     ], 11).split()]
    } for num in atom_ids]


def merge_records(records, start, join=" ", dont_condense=""):
    """Gets a single continuous string from a sequence of records.

    :param list records: The records to merge.
    :param int start: The start point in each record.
    :param str join: The string to join on.
    :param str dont_condense: By default any spaces after spaces, semi-colons, \
    colons, commas and dashes will be removed, unless listed here.
    :rtype: ``str``"""

    string = join.join(
     str(record[start:]
    ) if record[start:] else "" for record in records)
    condense = [char for char in " ;:,-" if char not in dont_condense]
    for char in condense:
        string = string.replace(char + " ", char)
    return string
