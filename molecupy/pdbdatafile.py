"""This module performs the actual parsing of the PDB file, though it does not
process the values that it extracts."""

import datetime

class PdbDataFile:
    """This object is essentially a list of values extracted from a PDB file.
    Each PDB record type has a function which identifies it and extracts the
    information from it.

    :param PdbFile pdb_file: The PDB file to extract information from.

    .. py:attribute:: pdb_file:

        The :py:class:`.PdbFile` from which the object was created."""

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
        process_jrnl(self)
        process_remark(self)

        process_dbref(self)
        process_seqadv(self)
        process_seqres(self)
        process_modres(self)

        process_het(self)
        process_hetnam(self)
        process_hetsyn(self)
        process_formul(self)

        process_helix(self)
        process_sheet(self)

        process_ssbond(self)
        process_link(self)
        process_cispep(self)

        process_site(self)

        process_crystal(self)
        process_origx(self)
        process_scale(self)
        process_mtrix(self)

        process_model(self)
        process_atom(self)
        process_anisou(self)
        process_ter(self)
        process_hetatm(self)

        process_conect(self)

        process_master(self)


    def __repr__(self):
        return "<PdbDataFile (????)>"


    def get_remark_by_number(self, number):
        for remark in self.remarks:
            if remark["number"] == number:
                return remark



def date_from_string(s):
    """Gets a Date object from a PDB formatted date string.

    :param str s: A date in the format DD-MM-YY.
    :rtype: ``datetime.Date``"""

    return datetime.datetime.strptime(
     s, "%d-%b-%y"
    ).date()


def merge_records(records, start, join=" ", dont_condense=""):
    """Gets a single continuous string from a sequence of records.

    :param list records: The records to merge.
    :param int start: The start point in each record.
    :param str join: The string to join on.
    :param str dont_condense: By default any spaces after spaces, semi-colons, \
    colons, commas and dashes will be removed, unless listed here.
    :rtype: ``str``"""

    string = join.join(str(record[start:]) if record[start:] else "" for record in records)
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


def process_header(data_file):
    """Processes HEADER records.

    :param ``PdbDataFile`` data_file: The data file."""

    header = data_file.pdb_file.get_record_by_name("HEADER")
    data_file.classification = header[10:50] if header else None
    data_file.deposition_date = date_from_string(header[50:59]) if header else None
    data_file.pdb_code = header[62:66] if header else None


def process_obslte(data_file):
    """Processes OBSLTE records.

    :param ``PdbDataFile`` data_file: The data file."""

    obslte = data_file.pdb_file.get_record_by_name("OBSLTE")
    data_file.is_obsolete = bool(obslte)
    data_file.obsolete_date = date_from_string(obslte[11:20]) if obslte else None
    data_file.replacement_code = obslte[31:35] if obslte else None


def process_title(data_file):
    """Processes TITLE records.

    :param ``PdbDataFile`` data_file: The data file."""

    titles = data_file.pdb_file.get_records_by_name("TITLE")
    title = merge_records(titles, 10, dont_condense=",;:-")
    data_file.title = title if title else None


def process_split(data_file):
    """Processes SPLIT records.

    :param ``PdbDataFile`` data_file: The data file."""

    splits = data_file.pdb_file.get_records_by_name("SPLIT")
    data_file.split_codes = " ".join([r[10:] for r in splits]).split()


def process_caveat(data_file):
    """Processes CAVEAT records.

    :param ``PdbDataFile`` data_file: The data file."""

    caveats = data_file.pdb_file.get_records_by_name("CAVEAT")
    caveat = merge_records(caveats, 19)
    data_file.caveat = caveat if caveat else None


def process_compnd(data_file):
    """Processes COMPND records.

    :param ``PdbDataFile`` data_file: The data file."""

    records = data_file.pdb_file.get_records_by_name("COMPND")
    data_file.compounds = records_to_token_value_dicts(records)


def process_source(data_file):
    """Processes SOURCE records.

    :param ``PdbDataFile`` data_file: The data file."""

    records = data_file.pdb_file.get_records_by_name("SOURCE")
    data_file.sources = records_to_token_value_dicts(records)


def process_keywds(data_file):
    """Processes KEYWDS records.

    :param ``PdbDataFile`` data_file: The data file."""

    keywords = data_file.pdb_file.get_records_by_name("KEYWDS")
    keyword_text = merge_records(keywords, 10)
    data_file.keywords = keyword_text.split(",") if keyword_text else []


def process_expdta(data_file):
    """Processes EXPDTA records.

    :param ``PdbDataFile`` data_file: The data file."""

    expdta = data_file.pdb_file.get_records_by_name("EXPDTA")
    expdta_text = merge_records(expdta, 10)
    data_file.experimental_techniques = expdta_text.split(";") if expdta_text else []


def process_nummdl(data_file):
    """Processes NUMMDL records.

    :param ``PdbDataFile`` data_file: The data file."""

    nummdl = data_file.pdb_file.get_record_by_name("NUMMDL")
    data_file.model_count = nummdl[10:14] if nummdl else 1


def process_mdltyp(data_file):
    """Processes MDLTYP records.

    :param ``PdbDataFile`` data_file: The data file."""

    mdltyps = data_file.pdb_file.get_records_by_name("MDLTYP")
    mdltyp_text = merge_records(mdltyps, 10, dont_condense=",")
    data_file.model_annotations = [
     ann.strip() for ann in mdltyp_text.split(";") if ann.strip()
    ]


def process_author(data_file):
    """Processes AUTHOR records.

    :param ``PdbDataFile`` data_file: The data file."""

    authors = data_file.pdb_file.get_records_by_name("AUTHOR")
    data_file.authors = merge_records(authors, 10).split(",") if authors else []


def process_revdat(data_file):
    """Processes REVDAT records.

    :param ``PdbDataFile`` data_file: The data file."""

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
    """Processes SPRSDE records.

    :param ``PdbDataFile`` data_file: The data file."""

    sprsde = data_file.pdb_file.get_record_by_name("SPRSDE")
    data_file.supercedes = sprsde[31:75].split() if sprsde else []
    data_file.supercede_date = date_from_string(sprsde[11:20]) if sprsde else None


def process_jrnl(data_file):
    """Processes JRNL records.

    :param ``PdbDataFile`` data_file: The data file."""

    jrnls = data_file.pdb_file.get_records_by_name("JRNL")
    if not jrnls:
        data_file.journal = None
    else:
        data_file.journal = {}
        auths = [r for r in jrnls if r[12:16] == "AUTH"]
        data_file.journal["authors"] = merge_records(auths, 19).split(",") if auths else []
        titls = [r for r in jrnls if r[12:16] == "TITL"]
        data_file.journal["title"] = merge_records(titls, 19) if titls else None
        edits = [r for r in jrnls if r[12:16] == "EDIT"]
        data_file.journal["editors"] = merge_records(edits, 19).split(",") if edits else []
        refs = [r for r in jrnls if r[12:16] == "REF"]
        data_file.journal["reference"] = {} if refs else None
        if refs and "TO BE PUBLISHED" in refs[0].text:
            data_file.journal["reference"] = {
             "published": False, "publication": None,
             "volume": None, "page": None, "year": None
            }
        elif refs:
            data_file.journal["reference"] = {
             "published": True,
             "publication": refs[0][19:47],
             "volume": refs[0][51:55],
             "page": refs[0][56:61],
             "year": refs[0][62:66]
            }
        publs = [r for r in jrnls if r[12:16] == "PUBL"]
        data_file.journal["publisher"] = merge_records(publs, 19, dont_condense=",:;") if publs else None
        refns = [r for r in jrnls if r[12:16] == "REFN"]
        data_file.journal["reference_number"] = {
         "type": refns[0][35:39],
         "value": refns[0][40:65]
        } if refns else None
        pmids = [r for r in jrnls if r[12:16] == "PMID"]
        data_file.journal["pubmed"] = pmids[0].get_as_string(19, 79) if pmids else None
        dois = [r for r in jrnls if r[12:16] == "DOI"]
        data_file.journal["doi"] = dois[0][19:79] if dois else None


def process_remark(data_file):
    """Processes REMARK records.

    :param ``PdbDataFile`` data_file: The data file."""

    remarks = data_file.pdb_file.get_records_by_name("REMARK")
    remark_numbers = sorted(list(set([r[7:10] for r in remarks])))
    data_file.remarks = []
    for number in remark_numbers:
        recs = [r for r in remarks if r[7:10] == number]
        remark = {
         "number": number,
         "content": merge_records(recs[1:], 11, join="\n", dont_condense=" ,:;")
        }
        data_file.remarks.append(remark)


def process_dbref(data_file):
    """Processes DBREF, DBREF1 and DBREF2 records.

    :param ``PdbDataFile`` data_file: The data file."""

    dbrefs = data_file.pdb_file.get_records_by_name("DBREF")
    data_file.dbreferences = [{
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
    dbref1s = data_file.pdb_file.get_records_by_name("DBREF1")
    dbref2s = data_file.pdb_file.get_records_by_name("DBREF2")
    ref_pairs = zip(dbref1s, dbref2s)
    data_file.dbreferences += [{
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
    data_file.dbreferences = sorted(data_file.dbreferences, key=lambda k: k["chain_id"])


def process_seqadv(data_file):
    """Processes SEQADV records.

    :param ``PdbDataFile`` data_file: The data file."""

    seqadvs = data_file.pdb_file.get_records_by_name("SEQADV")
    data_file.sequence_differences = [{
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


def process_seqres(data_file):
    """Processes SEQRES records.

    :param ``PdbDataFile`` data_file: The data file."""

    seqres = data_file.pdb_file.get_records_by_name("SEQRES")
    chains = sorted(list(set([r[11] for r in seqres])))
    data_file.residue_sequences = []
    for chain in chains:
        records = [r for r in seqres if r[11] == chain]
        data_file.residue_sequences.append({
         "chain_id": chain,
         "length": records[0][13:17],
         "residues": merge_records(records, 19).split()
        })


def process_modres(data_file):
    """Processes MODRES records.

    :param ``PdbDataFile`` data_file: The data file."""

    modres = data_file.pdb_file.get_records_by_name("MODRES")
    data_file.modified_residues = [{
     "residue_name": r[12:15],
     "chain_id": r[16],
     "residue_id": r[18:22],
     "insert_code": r[22] if r[22] else "",
     "standard_resisdue_name": r[24:27],
     "comment": r[29:70]
    } for r in modres]


def process_het(data_file):
    """Processes HET records.

    :param ``PdbDataFile`` data_file: The data file."""

    hets = data_file.pdb_file.get_records_by_name("HET")
    data_file.hets = [{
     "het_name": r[7:10],
     "chain_id": r[12],
     "het_id": r[13:17],
     "insert_code": r[17] if r[17] else "",
     "atom_num": r[20:25],
     "description": r[30:70]
    } for r in hets]


def process_hetnam(data_file):
    """Processes HETNAM records.

    :param ``PdbDataFile`` data_file: The data file."""

    hetnams = data_file.pdb_file.get_records_by_name("HETNAM")
    ids = list(set([r[11:14] for r in hetnams]))
    data_file.het_names = {
     het_id: merge_records(
      [r for r in hetnams if r[11:14] == het_id], 15, dont_condense=":;"
     ) for het_id in ids
    }


def process_hetsyn(data_file):
    """Processes HETSYN records.

    :param ``PdbDataFile`` data_file: The data file."""

    hetsyns = data_file.pdb_file.get_records_by_name("HETSYN")
    ids = list(set([r[11:14] for r in hetsyns]))
    data_file.het_synonyms = {
     het_id: merge_records(
      [r for r in hetsyns if r[11:14] == het_id], 15
     ).split(";") for het_id in ids
    }


def process_formul(data_file):
    """Processes FORMUL records.

    :param ``PdbDataFile`` data_file: The data file."""

    formuls = data_file.pdb_file.get_records_by_name("FORMUL")
    ids = list(set([r[12:15] for r in formuls]))
    data_file.formulae = {
     het_id: {
      "component_number": [r for r in formuls if r[12:15] == het_id][0][8:10],
      "is_water": [r for r in formuls if r[12:15] == het_id][0][18] == "*",
      "formula": merge_records(
       [r for r in formuls if r[12:15] == het_id], 19
      )
     } for het_id in ids
    }


def process_helix(data_file):
    """Processes HELIX records.

    :param ``PdbDataFile`` data_file: The data file."""

    helix = data_file.pdb_file.get_records_by_name("HELIX")
    data_file.helices = [{
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


def process_sheet(data_file):
    """Processes SHEET records.

    :param ``PdbDataFile`` data_file: The data file."""

    sheets = data_file.pdb_file.get_records_by_name("SHEET")
    sheet_names = sorted(list(set([r[11:14] for r in sheets])))
    data_file.sheets = []
    for sheet_name in sheet_names:
        strands = [r for r in sheets if r[11:14] == sheet_name]
        data_file.sheets.append({
         "sheet_id": sheet_name,
         "strand_count": strands[0][14:16],
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
          "sense": r[38:40],
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


def process_ssbond(data_file):
    """Processes SSBOND records.

    :param ``PdbDataFile`` data_file: The data file."""

    ssbonds = data_file.pdb_file.get_records_by_name("SSBOND")
    data_file.ss_bonds = [{
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


def process_link(data_file):
    """Processes LINK records.

    :param ``PdbDataFile`` data_file: The data file."""

    links = data_file.pdb_file.get_records_by_name("LINK")
    data_file.links = [{
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


def process_cispep(data_file):
    """Processes CISPEP records.

    :param ``PdbDataFile`` data_file: The data file."""

    cispeps = data_file.pdb_file.get_records_by_name("CISPEP")
    data_file.cis_peptides = [{
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


def process_site(data_file):
    """Processes SITE records.

    :param ``PdbDataFile`` data_file: The data file."""

    sites = data_file.pdb_file.get_records_by_name("SITE")
    site_names = sorted(list(set([r[11:14] for r in sites])))
    data_file.sites = []
    for site_name in site_names:
        records = [r for r in sites if r[11:14] == site_name]
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
        data_file.sites.append({
         "site_id": site_name,
         "residue_count": records[0][15:17],
         "residues": residues
        })


def process_crystal(data_file):
    """Processes CRYST1 records.

    :param ``PdbDataFile`` data_file: The data file."""

    crystal = data_file.pdb_file.get_record_by_name("CRYST1")
    data_file.crystal_a = crystal[6:15] if crystal else None
    data_file.crystal_b = crystal[15:24] if crystal else None
    data_file.crystal_c = crystal[24:33] if crystal else None
    data_file.crystal_alpha = crystal[33:40] if crystal else None
    data_file.crystal_beta = crystal[40:47] if crystal else None
    data_file.crystal_gamma = crystal[47:54] if crystal else None
    data_file.crystal_s_group = crystal[55:66] if crystal else None
    data_file.crystal_z = crystal[66:70] if crystal else None


def process_origx(data_file):
    """Processes ORIGX1, ORIGX2 and ORIGX3 records.

    :param ``PdbDataFile`` data_file: The data file."""

    origx1 = data_file.pdb_file.get_record_by_name("ORIGX1")
    data_file.crystal_o11 = origx1[10:20] if origx1 else None
    data_file.crystal_o12 = origx1[20:30] if origx1 else None
    data_file.crystal_o13 = origx1[30:40] if origx1 else None
    data_file.crystal_t1 = origx1[45:55] if origx1 else None
    origx2 = data_file.pdb_file.get_record_by_name("ORIGX2")
    data_file.crystal_o21 = origx2[10:20] if origx2 else None
    data_file.crystal_o22 = origx2[20:30] if origx2 else None
    data_file.crystal_o23 = origx2[30:40] if origx2 else None
    data_file.crystal_t2 = origx2[45:55] if origx2 else None
    origx3 = data_file.pdb_file.get_record_by_name("ORIGX3")
    data_file.crystal_o31 = origx3[10:20] if origx3 else None
    data_file.crystal_o32 = origx3[20:30] if origx3 else None
    data_file.crystal_o33 = origx3[30:40] if origx3 else None
    data_file.crystal_t3 = origx3[45:55] if origx3 else None


def process_scale(data_file):
    """Processes SCALE1, SCALE2 and SCALE3 records.

    :param ``PdbDataFile`` data_file: The data file."""

    scale1 = data_file.pdb_file.get_record_by_name("SCALE1")
    data_file.crystal_s11 = scale1[10:20] if scale1 else None
    data_file.crystal_s12 = scale1[20:30] if scale1 else None
    data_file.crystal_s13 = scale1[30:40] if scale1 else None
    data_file.crystal_u1 = scale1[45:55] if scale1 else None
    scale2 = data_file.pdb_file.get_record_by_name("SCALE2")
    data_file.crystal_s21 = scale2[10:20] if scale2 else None
    data_file.crystal_s22 = scale2[20:30] if scale2 else None
    data_file.crystal_s23 = scale2[30:40] if scale2 else None
    data_file.crystal_u2 = scale2[45:55] if scale2 else None
    scale3 = data_file.pdb_file.get_record_by_name("SCALE3")
    data_file.crystal_s31 = scale3[10:20] if scale3 else None
    data_file.crystal_s32 = scale3[20:30] if scale3 else None
    data_file.crystal_s33 = scale3[30:40] if scale3 else None
    data_file.crystal_u3 = scale3[45:55] if scale3 else None


def process_mtrix(data_file):
    """Processes MTRIX1, MTRIX2 and MTRIX3 records.

    :param ``PdbDataFile`` data_file: The data file."""

    mtrix1 = data_file.pdb_file.get_record_by_name("MTRIX1")
    data_file.crystal_serial_1 = mtrix1[7:10] if mtrix1 else None
    data_file.crystal_m11 = mtrix1[10:20] if mtrix1 else None
    data_file.crystal_m12 = mtrix1[20:30] if mtrix1 else None
    data_file.crystal_m13 = mtrix1[30:40] if mtrix1 else None
    data_file.crystal_v1 = mtrix1[45:55] if mtrix1 else None
    data_file.crystal_i_given_1 = mtrix1[59] == 1 if mtrix1 else False
    mtrix2 = data_file.pdb_file.get_record_by_name("MTRIX2")
    data_file.crystal_serial_2 = mtrix2[7:10] if mtrix2 else None
    data_file.crystal_m21 = mtrix2[10:20] if mtrix2 else None
    data_file.crystal_m22 = mtrix2[20:30] if mtrix2 else None
    data_file.crystal_m23 = mtrix2[30:40] if mtrix2 else None
    data_file.crystal_v2 = mtrix2[45:55] if mtrix2 else None
    data_file.crystal_i_given_2 = mtrix2[59] == 1 if mtrix2 else False
    mtrix3 = data_file.pdb_file.get_record_by_name("MTRIX3")
    data_file.crystal_serial_3 = mtrix3[7:10] if mtrix3 else None
    data_file.crystal_m31 = mtrix3[10:20] if mtrix3 else None
    data_file.crystal_m32 = mtrix3[20:30] if mtrix3 else None
    data_file.crystal_m33 = mtrix3[30:40] if mtrix3 else None
    data_file.crystal_v3 = mtrix3[45:55] if mtrix3 else None
    data_file.crystal_i_given_3 = mtrix3[59] == 1 if mtrix3 else False


def process_model(data_file):
    """Processes MODEL records.

    :param ``PdbDataFile`` data_file: The data file."""

    models = data_file.pdb_file.get_records_by_name("MODEL")
    endmdls = data_file.pdb_file.get_records_by_name("ENDMDL")
    pairs = list(zip(models, endmdls))
    data_file.models = [{
     "model_id": pair[0][10:14],
     "start_record": pair[0].number,
     "end_record": pair[1].number,
    } for pair in pairs]
    if not pairs:
        data_file.models = [{
         "model_id": 1,
         "start_record": 0,
         "end_record": len(data_file.pdb_file.records),
        }]


def process_atom(data_file):
    """Processes ATOM records.

    :param ``PdbDataFile`` data_file: The data file."""

    atoms = data_file.pdb_file.get_records_by_name("ATOM")
    data_file.atoms = [{
     "atom_id": r[6:11],
     "atom_name": r[12:16],
     "alt_loc": r[16],
     "residue_name": r[17:20],
     "chain_id": r[21],
     "residue_id": r[22:26],
     "insert_code": r[26] if r[26] else "",
     "x": r[30:38],
     "y": r[38:46],
     "z": r[46:54],
     "occupancy": r[54:60],
     "temperature_factor": r[60:66],
     "element": r[76:78],
     "charge": r[78:80],
     "model_id": [m for m in data_file.models
      if r.number >= m["start_record"]
       and r.number <= m["end_record"]][0]["model_id"]
    } for r in atoms]


def process_anisou(data_file):
    """Processes ANISOU records.

    :param ``PdbDataFile`` data_file: The data file."""

    anisou = data_file.pdb_file.get_records_by_name("ANISOU")
    data_file.anisou = [{
     "atom_id": r[6:11],
     "atom_name": r[12:16],
     "alt_loc": r[16],
     "residue_name": r[17:20],
     "chain_id": r[21],
     "residue_id": r[22:26],
     "insert_code": r[26] if r[26] else "",
     "u11": r[29:35],
     "u22": r[36:42],
     "u33": r[43:49],
     "u12": r[50:56],
     "u13": r[57:63],
     "u23": r[64:70],
     "element": r[76:78],
     "charge": r[78:80],
     "model_id": [m for m in data_file.models
      if r.number >= m["start_record"]
       and r.number <= m["end_record"]][0]["model_id"]
    } for r in anisou]


def process_ter(data_file):
    """Processes TER records.

    :param ``PdbDataFile`` data_file: The data file."""

    ters = data_file.pdb_file.get_records_by_name("TER")
    data_file.termini = [{
     "atom_id": r[6:11],
     "residue_name": r[17:20],
     "chain_id": r[21],
     "residue_id": r[22:26],
     "insert_code": r[26] if r[26] else "",
     "model_id": [m for m in data_file.models
      if r.number >= m["start_record"]
       and r.number <= m["end_record"]][0]["model_id"]
    } for r in ters]


def process_hetatm(data_file):
    """Processes HETATM records.

    :param ``PdbDataFile`` data_file: The data file."""

    hetatms = data_file.pdb_file.get_records_by_name("HETATM")
    data_file.heteroatoms = [{
     "atom_id": r[6:11],
     "atom_name": r[12:16],
     "alt_loc": r[16],
     "residue_name": r.get_as_string(17, 20),
     "chain_id": r[21],
     "residue_id": r[22:26],
     "insert_code": r[26] if r[26] else "",
     "x": r[30:38],
     "y": r[38:46],
     "z": r[46:54],
     "occupancy": r[54:60],
     "temperature_factor": r[60:66],
     "element": r[76:78],
     "charge": r[78:80],
     "model_id": [m for m in data_file.models
      if r.number >= m["start_record"]
       and r.number <= m["end_record"]][0]["model_id"]
    } for r in hetatms]


def process_conect(data_file):
    """Processes CONECT records.

    :param ``PdbDataFile`` data_file: The data file."""

    conects = data_file.pdb_file.get_records_by_name("CONECT")
    atom_ids = sorted(list(set([r[6:11] for r in conects])))
    data_file.connections = [{
     "atom_id": num,
     "bonded_atoms": [int(n) for n in merge_records([
      r for r in conects if r[6:11] == num
     ], 11).split()]
    } for num in atom_ids]


def process_master(data_file):
    """Processes MASTER records.

    :param ``PdbDataFile`` data_file: The data file."""

    master = data_file.pdb_file.get_record_by_name("MASTER")
    data_file.master = {
      "remark_num": master[10:15],
      "het_num": master[20:25],
      "helix_num": master[25:30],
      "sheet_num": master[30:35],
      "site_num": master[40:45],
      "crystal_num": master[45:50],
      "coordinate_num": master[50:55],
      "ter_num": master[55:60],
      "conect_num": master[60:65],
      "seqres_num": master[65:70]
    } if master else None
