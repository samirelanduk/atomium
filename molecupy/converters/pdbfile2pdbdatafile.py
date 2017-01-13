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
    process_split_records(data_file, pdb_file)
    process_caveat_records(data_file, pdb_file)
    process_compnd_records(data_file, pdb_file)
    process_source_records(data_file, pdb_file)
    process_keywd_records(data_file, pdb_file)
    process_expdta_records(data_file, pdb_file)
    process_nummdl_records(data_file, pdb_file)
    process_mdltyp_records(data_file, pdb_file)
    process_author_records(data_file, pdb_file)
    process_revdat_records(data_file, pdb_file)
    process_sprsde_records(data_file, pdb_file)
    process_jrnl_records(data_file, pdb_file)
    process_remark_records(data_file, pdb_file)

    process_dbref_records(data_file, pdb_file)
    
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


def process_split_records(data_file, pdb_file):
    splits = pdb_file.get_records_by_name("SPLIT")
    data_file._split_codes = " ".join([r[10:] for r in splits]).split()


def process_caveat_records(data_file, pdb_file):
    caveats = pdb_file.get_records_by_name("CAVEAT")
    if caveats:
        data_file._caveat = merge_records(caveats, 19)
    else:
        data_file._caveat = None


def process_compnd_records(data_file, pdb_file):
    records = pdb_file.get_records_by_name("COMPND")
    data_file._compounds = records_to_token_value_dicts(records)


def process_source_records(data_file, pdb_file):
    records = pdb_file.get_records_by_name("SOURCE")
    data_file._sources = records_to_token_value_dicts(records)


def process_keywd_records(data_file, pdb_file):
    keywords = pdb_file.get_records_by_name("KEYWDS")
    if keywords:
        keyword_text = merge_records(keywords, 10)
        data_file._keywords = keyword_text.split(",")
    else:
        data_file._keywords = []


def process_expdta_records(data_file, pdb_file):
    expdta = pdb_file.get_records_by_name("EXPDTA")
    if expdta:
        expdta_text = merge_records(expdta, 10)
        data_file._experimental_techniques = expdta_text.split(";")
    else:
        data_file._experimental_techniques = []


def process_nummdl_records(data_file, pdb_file):
    nummdl = pdb_file.get_record_by_name("NUMMDL")
    if nummdl:
        data_file._model_count = nummdl[10:14]
    else:
        data_file._model_count = 1


def process_mdltyp_records(data_file, pdb_file):
    mdltyps = pdb_file.get_records_by_name("MDLTYP")
    if mdltyps:
        mdltyp_text = merge_records(mdltyps, 10, dont_condense=",")
        data_file._model_annotations = [
         ann.strip() for ann in mdltyp_text.split(";") if ann.strip()
        ]
    else:
        data_file._model_annotations = []


def process_author_records(data_file, pdb_file):
    authors = pdb_file.get_records_by_name("AUTHOR")
    if authors:
        data_file._authors = merge_records(authors, 10).split(",")
    else:
        data_file._authors = []


def process_revdat_records(data_file, pdb_file):
    revdats = pdb_file.get_records_by_name("REVDAT")
    if revdats:
        numbers = sorted(list(set([r[7:10] for r in revdats])))
        revisions = []
        for number in numbers:
            records = [r for r in revdats if r[7:10] == number]
            rec_types = merge_records(records, 39).split()
            revisions.append({
             "number": number,
             "date": date_from_string(records[0][13:22]),
             "type": records[0][31],
             "records": [r for r in rec_types if r]
            })
        data_file._revisions = revisions
    else:
        data_file._revisions = []


def process_sprsde_records(data_file, pdb_file):
    sprsde = pdb_file.get_record_by_name("SPRSDE")
    if sprsde:
        data_file._supercedes = sprsde[31:75].split()
        data_file._supercede_date = date_from_string(sprsde[11:20])
    else:
        data_file._supercedes = []
        data_file._supercede_date = None


def process_jrnl_records(data_file, pdb_file):
    jrnls = pdb_file.get_records_by_name("JRNL")
    if jrnls:
        journal = {}
        auths = [r for r in jrnls if r[12:16] == "AUTH"]
        journal["authors"] = merge_records(auths, 19).split(",") if auths else []
        titls = [r for r in jrnls if r[12:16] == "TITL"]
        journal["title"] = merge_records(titls, 19) if titls else None
        edits = [r for r in jrnls if r[12:16] == "EDIT"]
        journal["editors"] = merge_records(edits, 19).split(",") if edits else []
        refs = [r for r in jrnls if r[12:16] == "REF"]
        journal["reference"] = {} if refs else None
        if refs and "TO BE PUBLISHED" in refs[0].text():
            journal["reference"] = {
             "published": False, "publication": None,
             "volume": None, "page": None, "year": None
            }
        elif refs:
            journal["reference"] = {
             "published": True,
             "publication": refs[0][19:47],
             "volume": refs[0][51:55],
             "page": refs[0][56:61],
             "year": refs[0][62:66]
            }
        publs = [r for r in jrnls if r[12:16] == "PUBL"]
        journal["publisher"] = merge_records(
         publs, 19, dont_condense=",:;"
        ) if publs else None
        refns = [r for r in jrnls if r[12:16] == "REFN"]
        journal["reference_number"] = {
         "type": refns[0][35:39],
         "value": refns[0][40:65]
        } if refns else None
        pmids = [r for r in jrnls if r[12:16] == "PMID"]
        journal["pubmed"] = pmids[0].get_as_string(19, 79) if pmids else None
        dois = [r for r in jrnls if r[12:16] == "DOI"]
        journal["doi"] = dois[0][19:79] if dois else None
        data_file._journal = journal
    else:
        data_file._journal = None


def process_remark_records(data_file, pdb_file):
    remark_records = pdb_file.get_records_by_name("REMARK")
    if remark_records:
        remark_numbers = sorted(list(set([r[7:10] for r in remark_records])))
        remarks = []
        for number in remark_numbers:
            recs = [r for r in remark_records if r[7:10] == number]
            remark = {
             "number": number,
             "content": merge_records(recs[1:], 11, join="\n", dont_condense=" ,:;")
            }
            remarks.append(remark)
        data_file._remarks = remarks
    else:
        data_file._remarks = []


def process_dbref_records(data_file, pdb_file):
    dbrefs = pdb_file.get_records_by_name("DBREF")
    dbref1s = pdb_file.get_records_by_name("DBREF1")
    dbref2s = pdb_file.get_records_by_name("DBREF2")
    if dbrefs or dbref1s or dbref2s:
        dbreferences = [{
         "chain_id": r[12],
         "sequence_begin": r[14:18],
         "insert_begin": r[18] if r[18] else "",
         "sequence_end": r[20:24],
         "insert_end": r[24] if r[24] else "",
         "database": r[26:32],
         "accession": r.get_as_string(33, 40),
         "db_id": r[42:54],
         "db_sequence_begin": r[55:60],
         "db_insert_begin": r[60],
         "db_sequence_end": r[62:67],
         "db_insert_end": r[67]
        } for r in dbrefs]
        ref_pairs = zip(dbref1s, dbref2s)
        dbreferences += [{
         "chain_id": pair[0][12],
         "sequence_begin": pair[0][14:18],
         "insert_begin": pair[0][18] if pair[0][18] else "",
         "sequence_end": pair[0][20:24],
         "insert_end": pair[0][24] if pair[0][24] else "",
         "database": pair[0][26:32],
         "accession": pair[1].get_as_string(18, 40),
         "db_id": pair[0][47:67],
         "db_sequence_begin": pair[1][45:55],
         "db_insert_begin": None,
         "db_sequence_end": pair[1][57:67],
         "db_insert_end": None
        } for pair in ref_pairs]
        dbreferences = sorted(dbreferences, key=lambda k: k["chain_id"])
        data_file._dbreferences = dbreferences
    else:
        data_file._dbreferences = []


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


def records_to_token_value_dicts(records):
    """Produces a list of ``dict`` objects from the key-value pairs used in \
    COMPND and SOURCE records.

    :param list records: The records to use.
    :rtype: ``list``"""

    string = merge_records(records, 10)
    pairs = list(filter(None, string.split(";")))
    for pair_offset in range(1, len(pairs))[::-1]:
        if pairs[pair_offset].count(":") == 0:
            pairs[pair_offset-1] += "; " + pairs[pair_offset]
    pairs = [pair for pair in pairs if pair.count(":") == 1]
    pairs = [pair.split(":") for pair in pairs if pair]
    entities = []
    entity = {}
    for pair in pairs:
        if pair[1] == "NO":
            pair[1] = False
        elif pair[1] == "YES":
            pair[1] = True
        elif pair[1].isnumeric():
            pair[1] = int(pair[1])
        elif pair[0] == "CHAIN" or pair[0] == "SYNONYM":
            pair[1] = pair[1].replace(", ", ",").split(",")
        if pair[0] == "MOL_ID":
            if entity: entities.append(entity)
            entity = {}
        entity[pair[0]] = pair[1]
    if entity: entities.append(entity)
    return entities
