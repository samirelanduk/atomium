import datetime

class PdbDataFile:

    def __init__(self, pdb_file):
        self._pdb_file = pdb_file


    def __repr__(self):
        return "<PdbDataFile (????)>"


    def pdb_file(self):
        return self._pdb_file


    def classification(self):
        header = self.pdb_file().get_record_by_name("HEADER")
        return header[10:50] if header else None


    def deposition_date(self):
        header = self.pdb_file().get_record_by_name("HEADER")
        return date_from_string(header[50:59]) if header else None


    def pdb_code(self):
        header = self.pdb_file().get_record_by_name("HEADER")
        return header[62:66] if header else None


    def is_obsolete(self):
        obslte = self.pdb_file().get_record_by_name("OBSLTE")
        return bool(obslte)


    def obsolete_date(self):
        obslte = self.pdb_file().get_record_by_name("OBSLTE")
        return date_from_string(obslte[11:20]) if obslte else None


    def replacement_code(self):
        obslte = self.pdb_file().get_record_by_name("OBSLTE")
        return obslte[31:35] if obslte else None



def date_from_string(s):
    return datetime.datetime.strptime(
     s, "%d-%b-%y"
    ).date()


def merge_records(records, start, join=" ", dont_condense=""):
    string = join.join(
     str(record[start:]
    ) if record[start:] else "" for record in records)
    condense = [char for char in " ;:,-" if char not in dont_condense]
    for char in condense:
        string = string.replace(char + " ", char)
    return string


def records_to_token_value_dicts(records):
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
