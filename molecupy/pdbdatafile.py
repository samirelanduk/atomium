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


    def title(self):
        titles = self.pdb_file().get_records_by_name("TITLE")
        title = merge_records(titles, 10, dont_condense=",;:-")
        return title if title else None


    def split_codes(self):
        splits = self.pdb_file().get_records_by_name("SPLIT")
        return " ".join([r[10:] for r in splits]).split()


    def caveat(self):
        caveats = self.pdb_file().get_records_by_name("CAVEAT")
        caveat = merge_records(caveats, 19)
        return caveat if caveat else None


    def compounds(self):
        records = self.pdb_file().get_records_by_name("COMPND")
        return records_to_token_value_dicts(records)


    def sources(self):
        records = self.pdb_file().get_records_by_name("SOURCE")
        return records_to_token_value_dicts(records)


    def keywords(self):
        keywords = self.pdb_file().get_records_by_name("KEYWDS")
        keyword_text = merge_records(keywords, 10)
        return keyword_text.split(",") if keyword_text else []


    def experimental_techniques(self):
        expdta = self.pdb_file().get_records_by_name("EXPDTA")
        expdta_text = merge_records(expdta, 10)
        return expdta_text.split(";") if expdta_text else []


    def model_count(self):
        nummdl = self.pdb_file().get_record_by_name("NUMMDL")
        return nummdl[10:14] if nummdl else 1


    def model_annotations(self):
        mdltyps = self.pdb_file().get_records_by_name("MDLTYP")
        mdltyp_text = merge_records(mdltyps, 10, dont_condense=",")
        return [
         ann.strip() for ann in mdltyp_text.split(";") if ann.strip()
        ]


    def authors(self):
        authors = self.pdb_file().get_records_by_name("AUTHOR")
        return merge_records(authors, 10).split(",") if authors else []


    def revisions(self):
        revdats = self.pdb_file().get_records_by_name("REVDAT")
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
        return revisions


    def supercedes(self):
        sprsde = self.pdb_file().get_record_by_name("SPRSDE")
        return sprsde[31:75].split() if sprsde else []


    def supercede_date(self):
        sprsde = self.pdb_file().get_record_by_name("SPRSDE")
        return date_from_string(sprsde[11:20]) if sprsde else None


    def journal(self):
        jrnls = self.pdb_file().get_records_by_name("JRNL")
        if not jrnls:
            return None
        else:
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
            journal["publisher"] = merge_records(publs, 19, dont_condense=",:;") if publs else None
            refns = [r for r in jrnls if r[12:16] == "REFN"]
            journal["reference_number"] = {
             "type": refns[0][35:39],
             "value": refns[0][40:65]
            } if refns else None
            pmids = [r for r in jrnls if r[12:16] == "PMID"]
            journal["pubmed"] = pmids[0].get_as_string(19, 79) if pmids else None
            dois = [r for r in jrnls if r[12:16] == "DOI"]
            journal["doi"] = dois[0][19:79] if dois else None
            return journal


    def remarks(self):
        remark_records = self.pdb_file().get_records_by_name("REMARK")
        remark_numbers = sorted(list(set([r[7:10] for r in remark_records])))
        remarks = []
        for number in remark_numbers:
            recs = [r for r in remark_records if r[7:10] == number]
            remark = {
             "number": number,
             "content": merge_records(recs[1:], 11, join="\n", dont_condense=" ,:;")
            }
            remarks.append(remark)
        return remarks


    def get_remark_by_number(self, number):
        for remark in self.remarks():
            if remark["number"] == number:
                return remark


    def dbreferences(self):
        dbrefs = self.pdb_file().get_records_by_name("DBREF")
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
        dbref1s = self.pdb_file().get_records_by_name("DBREF1")
        dbref2s = self.pdb_file().get_records_by_name("DBREF2")
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
        return dbreferences


    def sequence_differences(self):
        seqadvs = self.pdb_file().get_records_by_name("SEQADV")
        return [{
         "residue_name": r[12:15],
         "chain_id": r[16],
         "residue_id": r[18:22],
         "insert_code": r[22] if r[22] else "",
         "database": r[24:28],
         "accession": r[29:38],
         "db_residue_name": r[39:42],
         "db_residue_id": r[43:48],
         "conflict": r[49:70]
        } for r in seqadvs]


    def residue_sequences(self):
        seqres = self.pdb_file().get_records_by_name("SEQRES")
        chains = sorted(list(set([r[11] for r in seqres])))
        residue_sequences = []
        for chain in chains:
            records = [r for r in seqres if r[11] == chain]
            residue_sequences.append({
             "chain_id": chain,
             "length": records[0][13:17],
             "residues": merge_records(records, 19).split()
            })
        return residue_sequences


    def modified_residues(self):
        modres = self.pdb_file().get_records_by_name("MODRES")
        return [{
         "residue_name": r[12:15],
         "chain_id": r[16],
         "residue_id": r[18:22],
         "insert_code": r[22] if r[22] else "",
         "standard_resisdue_name": r[24:27],
         "comment": r[29:70]
        } for r in modres]


    def hets(self):
        hets = self.pdb_file().get_records_by_name("HET")
        return [{
         "het_name": r[7:10],
         "chain_id": r[12],
         "het_id": r[13:17],
         "insert_code": r[17] if r[17] else "",
         "atom_num": r[20:25],
         "description": r[30:70]
        } for r in hets]


    def het_names(self):
        hetnams = self.pdb_file().get_records_by_name("HETNAM")
        ids = list(set([r[11:14] for r in hetnams]))
        return {
         het_id: merge_records(
          [r for r in hetnams if r[11:14] == het_id], 15, dont_condense=":;"
         ) for het_id in ids
        }


    def het_synonyms(self):
        hetsyns = self.pdb_file().get_records_by_name("HETSYN")
        ids = list(set([r[11:14] for r in hetsyns]))
        return {
         het_id: merge_records(
          [r for r in hetsyns if r[11:14] == het_id], 15
         ).split(";") for het_id in ids
        }


    def formulae(self):
        formuls = self.pdb_file().get_records_by_name("FORMUL")
        ids = list(set([r[12:15] for r in formuls]))
        return {
         het_id: {
          "component_number": [r for r in formuls if r[12:15] == het_id][0][8:10],
          "is_water": [r for r in formuls if r[12:15] == het_id][0][18] == "*",
          "formula": merge_records(
           [r for r in formuls if r[12:15] == het_id], 19
          )
         } for het_id in ids
        }


    def helices(self):
        helix = self.pdb_file().get_records_by_name("HELIX")
        return [{
         "helix_id": r[7:10],
         "helix_name": r.get_as_string(11, 14),
         "start_residue_name": r[15:18],
         "start_residue_chain_id": r[19],
         "start_residue_id": r[21:25],
         "start_residue_insert": r[25] if r[25] else "",
         "end_residue_name": r[27:30],
         "end_residue_chain_id": r[31],
         "end_residue_id": r[33:37],
         "end_residue_insert": r[37] if r[37] else "",
         "helix_class": r[38:40],
         "comment": r[40:70],
         "length": r[71:76]
        } for r in helix]


    def sheets(self):
        sheet_records = self.pdb_file().get_records_by_name("SHEET")
        sheet_names = sorted(list(set([r[11:14] for r in sheet_records])))
        sheets = []
        for sheet_name in sheet_names:
            strands = [r for r in sheet_records if r[11:14] == sheet_name]
            sheets.append({
             "sheet_id": sheet_name,
             "strand_count": sheet_records[0][14:16],
             "strands": [{
              "strand_id": r[7:10],
              "start_residue_name": r[17:20],
              "start_residue_chain_id": r[21],
              "start_residue_id": r[22:26],
              "start_residue_insert": r[26] if r[26] else "",
              "end_residue_name": r[28:31],
              "end_residue_chain_id": r[32],
              "end_residue_id": r[33:37],
              "end_residue_insert": r[37] if r[37] else "",
              "sense": r[38:40] if r[38:40] else 0,
              "current_atom": r[41:45],
              "current_residue_name": r[45:48],
              "current_chain_id": r[49],
              "current_residue_id": r[50:54],
              "current_insert": r[54] if r[54] else "",
              "previous_atom": r[56:60],
              "previous_residue_name": r[60:63],
              "previous_chain_id": r[64],
              "previous_residue_id": r[65:69],
              "previous_insert": r[69] if r[69] else ""
             } for r in strands]
            })
        return sheets


    def ss_bonds(self):
        ssbonds = self.pdb_file().get_records_by_name("SSBOND")
        return [{
         "serial_num": r[7:10],
         "residue_name_1": r[11:14],
         "chain_id_1": r[15],
         "residue_id_1": r[17:21],
         "insert_code_1": r[21] if r[21] else "",
         "residue_name_2": r[25:28],
         "chain_id_2": r[29],
         "residue_id_2": r[31:35],
         "insert_code_2": r[35] if r[35] else "",
         "symmetry_1": r.get_as_string(59, 65),
         "symmetry_2": r.get_as_string(66, 72),
         "length": r[73:78]
        } for r in ssbonds]


    def links(self):
        links = self.pdb_file().get_records_by_name("LINK")
        return [{
         "atom_1": r[12:16],
         "alt_loc_1": r[16],
         "residue_name_1": r[17:20],
         "chain_id_1": r[21],
         "residue_id_1": r[22:26],
         "insert_code_1": r[26] if r[26] else "",
         "atom_2": r[42:46],
         "alt_loc_2": r[46],
         "residue_name_2": r[47:50],
         "chain_id_2": r[51],
         "residue_id_2": r[52:56],
         "insert_code_2": r[56] if r[56] else "",
         "symmetry_1": r.get_as_string(59, 65),
         "symmetry_2": r.get_as_string(66, 72),
         "length": r[73:78]
        } for r in links]


    def cis_peptides(self):
        cispeps = self.pdb_file().get_records_by_name("CISPEP")
        return [{
         "serial_num": r[7:10],
         "residue_name_1": r[11:14],
         "chain_id_1": r[15],
         "residue_id_1": r[17:21],
         "insert_1": r[21] if r[21] else "",
         "residue_name_2": r[25:28],
         "chain_id_2": r[29],
         "residue_id_2": r[31:35],
         "insert_2": r[35] if r[35] else "",
         "model_number": r[43:46],
         "angle": r[54:59]
        } for r in cispeps]


    def sites(self):
        site_records = self.pdb_file().get_records_by_name("SITE")
        site_names = sorted(list(set([r[11:14] for r in site_records])))
        sites = []
        for site_name in site_names:
            records = [r for r in site_records if r[11:14] == site_name]
            residues = []
            for r in records:
                for i in range(1, 5):
                    if r[(i * 11) + 7: (i * 11) + 17]:
                        residues.append({
                         "residue_name": r[(i * 11) + 7: (i * 11) + 10],
                         "chain_id": r[(i * 11) + 11],
                         "residue_id": r[(i * 11) + 12: (i * 11) + 16],
                         "insert_code":  r[(i * 11) + 16] if r[(i * 11) + 16] else ""
                        })
            sites.append({
             "site_id": site_name,
             "residue_count": records[0][15:17],
             "residues": residues
            })
        return sites



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
