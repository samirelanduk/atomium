import datetime

class PdbDataFile:

    def __init__(self, pdb_file):
        self.pdb_file = pdb_file

        process_header(self)
        process_obslte(self)
        process_title(self)
        process_split(self)
        process_caveat(self)
        process_compnd(self)
        process_source(self)
        process_keywds(self)
        process_expdta(self)
        process_nummdl(self)
        process_mdltyp(self)
        process_author(self)
        process_revdat(self)
        process_sprsde(self)


    def __repr__(self):
        return "<PdbDataFile (????)>"



def date_from_string(s):
    return datetime.datetime.strptime(
     s, "%d-%b-%y"
    ).date()


def merge_records(records, start, join=" ", dont_condense=""):
    string = join.join(str(record[start:]) if record[start:] else "" for record in records)
    condense = [char for char in " ;:,-" if char not in dont_condense]
    for char in condense:
        string = string.replace(char + " ", char)
    return string


def records_to_token_value_dicts(records):
    string = merge_records(records, 10)
    pairs = string.split(";")
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


def process_header(data_file):
    header = data_file.pdb_file.get_record_by_name("HEADER")
    data_file.classification = header[10:50] if header else None
    data_file.deposition_date = date_from_string(header[50:59]) if header else None
    data_file.pdb_code = header[62:66] if header else None


def process_obslte(data_file):
    obslte = data_file.pdb_file.get_record_by_name("OBSLTE")
    data_file.is_obsolete = bool(obslte)
    data_file.obsolete_date = date_from_string(obslte[11:20]) if obslte else None
    data_file.replacement_code = obslte[31:35] if obslte else None


def process_title(data_file):
    titles = data_file.pdb_file.get_records_by_name("TITLE")
    title = merge_records(titles, 10, dont_condense=",;:-")
    data_file.title = title if title else None


def process_split(data_file):
    splits = data_file.pdb_file.get_records_by_name("SPLIT")
    data_file.split_codes = " ".join([r[10:] for r in splits]).split()


def process_caveat(data_file):
    caveats = data_file.pdb_file.get_records_by_name("CAVEAT")
    caveat = merge_records(caveats, 19)
    data_file.caveat = caveat if caveat else None


def process_compnd(data_file):
    records = data_file.pdb_file.get_records_by_name("COMPND")
    data_file.compounds = records_to_token_value_dicts(records)


def process_source(data_file):
    records = data_file.pdb_file.get_records_by_name("SOURCE")
    data_file.sources = records_to_token_value_dicts(records)


def process_keywds(data_file):
    keywords = data_file.pdb_file.get_records_by_name("KEYWDS")
    keyword_text = merge_records(keywords, 10)
    data_file.keywords = keyword_text.split(",") if keyword_text else []


def process_expdta(data_file):
    expdta = data_file.pdb_file.get_records_by_name("EXPDTA")
    expdta_text = merge_records(expdta, 10)
    data_file.experimental_techniques = expdta_text.split(";") if expdta_text else []


def process_nummdl(data_file):
    nummdl = data_file.pdb_file.get_record_by_name("NUMMDL")
    data_file.model_count = nummdl[10:14] if nummdl else 1


def process_mdltyp(data_file):
    mdltyps = data_file.pdb_file.get_records_by_name("MDLTYP")
    mdltyp_text = merge_records(mdltyps, 10, dont_condense=",")
    data_file.model_annotations = [
     ann.strip() for ann in mdltyp_text.split(";") if ann.strip()
    ]


def process_author(data_file):
    authors = data_file.pdb_file.get_records_by_name("AUTHOR")
    data_file.authors = merge_records(authors, 10).split(",") if authors else []


def process_revdat(data_file):
    revdats = data_file.pdb_file.get_records_by_name("REVDAT")
    numbers = sorted(list(set([r[7:10] for r in revdats])))
    data_file.revisions = []
    for number in numbers:
        records = [r for r in revdats if r[7:10] == number]
        rec_types = merge_records(records, 39).split()
        data_file.revisions.append({
         "number": number,
         "date": date_from_string(records[0][13:22]),
         "type": records[0][31],
         "records": [r for r in rec_types if r]
        })


def process_sprsde(data_file):
    sprsde = data_file.pdb_file.get_record_by_name("SPRSDE")
    data_file.supercedes = sprsde[31:75].split() if sprsde else []
    data_file.supercede_date = date_from_string(sprsde[11:20]) if sprsde else None
