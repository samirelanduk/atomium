import datetime
from ..pdb.pdbfile import PdbFile
from ..pdb.pdbdatafile import PdbDataFile

def pdb_data_file_from_pdb_file(pdb_file):
    if not isinstance(pdb_file, PdbFile):
        raise TypeError("pdb_data_file_from_pdb_file can only convert PdbFiles")
    data_file = PdbDataFile()

    process_header_records(data_file, pdb_file)
    process_obslte_records(data_file, pdb_file)
    process_title_records(data_file, pdb_file)
    return data_file


def process_header_records(data_file, pdb_file):
    header = pdb_file.get_record_by_name("HEADER")
    if header:
        data_file._classification = header[10:50]
        data_file._deposition_date = date_from_string(header[50:59])
        data_file._pdb_code = header[62:66]
    else:
        data_file._classification = None
        data_file._deposition_date = None
        data_file._pdb_code = None


def process_obslte_records(data_file, pdb_file):
    obslte = pdb_file.get_record_by_name("OBSLTE")
    if obslte:
        data_file._is_obsolete = True
        data_file._obsolete_date = date_from_string(obslte[11:20])
        data_file._replacement_code = obslte[31:35]
    else:
        data_file._is_obsolete = False
        data_file._obsolete_date = None
        data_file._replacement_code = None


def process_title_records(data_file, pdb_file):
    titles = pdb_file.get_records_by_name("TITLE")
    if titles:
        data_file._title = merge_records(titles, 10, dont_condense=",;:-")
    else:
        data_file._title = None


def date_from_string(s):
    """Gets a Date object from a PDB formatted date string.

    :param str s: A date in the format DD-MM-YY.
    :rtype: ``datetime.Date``"""

    if s:
        return datetime.datetime.strptime(
         s, "%d-%b-%y"
        ).date()
    else:
        return None


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