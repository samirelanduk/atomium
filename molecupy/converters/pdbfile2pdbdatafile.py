"""This module handles the logic of converting a :py:class:`.PdbFile` to a
:py:class:`.PdbDataFile`"""

import datetime
from ..pdb.pdbfile import PdbFile
from ..pdb.pdbdatafile import PdbDataFile

def pdb_data_file_from_pdb_file(pdb_file):
    """Takes a :py:class:`.PdbFile`, converts it to a :py:class:`.PdbDataFile`,
    and returns it.

    :param PdbFile pdb_file: The :py:class:`.PdbFile` to convert.
    :rtype: :py:class:`.PdbDataFile`"""

    if not isinstance(pdb_file, PdbFile):
        raise TypeError("pdb_data_file_from_pdb_file can only convert PdbFiles")
    data_file = PdbDataFile()
    data_file._source = pdb_file

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
    process_seqadv_records(data_file, pdb_file)
    process_seqres_records(data_file, pdb_file)
    process_modres_records(data_file, pdb_file)

    process_het_records(data_file, pdb_file)
    process_hetnam_records(data_file, pdb_file)
    process_hetsyn_records(data_file, pdb_file)
    process_formul_records(data_file, pdb_file)

    process_helix_records(data_file, pdb_file)
    process_sheet_records(data_file, pdb_file)

    process_ssbond_records(data_file, pdb_file)
    process_link_records(data_file, pdb_file)
    process_cispep_records(data_file, pdb_file)

    process_site_records(data_file, pdb_file)

    process_cryst1_records(data_file, pdb_file)
    process_origx_records(data_file, pdb_file)
    process_scale_records(data_file, pdb_file)
    process_mtrix_records(data_file, pdb_file)

    process_model_records(data_file, pdb_file)
    process_atom_records(data_file, pdb_file)
    process_anisou_records(data_file, pdb_file)
    process_ter_records(data_file, pdb_file)
    process_hetatm_records(data_file, pdb_file)

    process_conect_records(data_file, pdb_file)

    process_master_records(data_file, pdb_file)

    return data_file


def process_header_records(data_file, pdb_file):
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    HEADER records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

    header = pdb_file.get_record_by_name("HEADER")
    if header:
        data_file._classification = header[10:50]
        data_file._deposition_date = _date_from_string(header[50:59])
        data_file._pdb_code = header[62:66]
    else:
        data_file._classification = None
        data_file._deposition_date = None
        data_file._pdb_code = None


def process_obslte_records(data_file, pdb_file):
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    OBSLTE records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

    obslte = pdb_file.get_record_by_name("OBSLTE")
    if obslte:
        data_file._is_obsolete = True
        data_file._obsolete_date = _date_from_string(obslte[11:20])
        data_file._replacement_code = obslte[31:35]
    else:
        data_file._is_obsolete = False
        data_file._obsolete_date = None
        data_file._replacement_code = None


def process_title_records(data_file, pdb_file):
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    TITLE records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

    titles = pdb_file.get_records_by_name("TITLE")
    if titles:
        data_file._title = _merge_records(titles, 10, dont_condense=",;:-")
    else:
        data_file._title = None


def process_split_records(data_file, pdb_file):
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    SPLIT records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

    splits = pdb_file.get_records_by_name("SPLIT")
    data_file._split_codes = " ".join([r[10:] for r in splits]).split()


def process_caveat_records(data_file, pdb_file):
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    CAVEAT records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

    caveats = pdb_file.get_records_by_name("CAVEAT")
    if caveats:
        data_file._caveat = _merge_records(caveats, 19)
    else:
        data_file._caveat = None


def process_compnd_records(data_file, pdb_file):
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    COMPND records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

    records = pdb_file.get_records_by_name("COMPND")
    data_file._compounds = _records_to_token_value_dicts(records)


def process_source_records(data_file, pdb_file):
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    SOURCE records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

    records = pdb_file.get_records_by_name("SOURCE")
    data_file._sources = _records_to_token_value_dicts(records)


def process_keywd_records(data_file, pdb_file):
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    KEYWD records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

    keywords = pdb_file.get_records_by_name("KEYWDS")
    if keywords:
        keyword_text = _merge_records(keywords, 10)
        data_file._keywords = keyword_text.split(",")
    else:
        data_file._keywords = []


def process_expdta_records(data_file, pdb_file):
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    EXPDTA records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

    expdta = pdb_file.get_records_by_name("EXPDTA")
    if expdta:
        expdta_text = _merge_records(expdta, 10)
        data_file._experimental_techniques = expdta_text.split(";")
    else:
        data_file._experimental_techniques = []


def process_nummdl_records(data_file, pdb_file):
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    NUMMDL records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

    nummdl = pdb_file.get_record_by_name("NUMMDL")
    if nummdl:
        data_file._model_count = nummdl[10:14]
    else:
        data_file._model_count = 1


def process_mdltyp_records(data_file, pdb_file):
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    MDLTYP records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

    mdltyps = pdb_file.get_records_by_name("MDLTYP")
    if mdltyps:
        mdltyp_text = _merge_records(mdltyps, 10, dont_condense=",")
        data_file._model_annotations = [
         ann.strip() for ann in mdltyp_text.split(";") if ann.strip()
        ]
    else:
        data_file._model_annotations = []


def process_author_records(data_file, pdb_file):
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    AUTHOR records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

    authors = pdb_file.get_records_by_name("AUTHOR")
    if authors:
        data_file._authors = _merge_records(authors, 10).split(",")
    else:
        data_file._authors = []


def process_revdat_records(data_file, pdb_file):
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    REVDAT records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

    revdats = pdb_file.get_records_by_name("REVDAT")
    if revdats:
        numbers = sorted(list(set([r[7:10] for r in revdats])))
        revisions = []
        for number in numbers:
            records = [r for r in revdats if r[7:10] == number]
            rec_types = _merge_records(records, 39).split()
            revisions.append({
             "number": number,
             "date": _date_from_string(records[0][13:22]),
             "type": records[0][31],
             "records": [r for r in rec_types if r]
            })
        data_file._revisions = revisions
    else:
        data_file._revisions = []


def process_sprsde_records(data_file, pdb_file):
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    SPRSDE records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

    sprsde = pdb_file.get_record_by_name("SPRSDE")
    if sprsde:
        data_file._supercedes = sprsde[31:75].split()
        data_file._supercede_date = _date_from_string(sprsde[11:20])
    else:
        data_file._supercedes = []
        data_file._supercede_date = None


def process_jrnl_records(data_file, pdb_file):
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    JRNL records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

    jrnls = pdb_file.get_records_by_name("JRNL")
    if jrnls:
        journal = {}
        auths = [r for r in jrnls if r[12:16] == "AUTH"]
        journal["authors"] = _merge_records(auths, 19).split(",") if auths else []
        titls = [r for r in jrnls if r[12:16] == "TITL"]
        journal["title"] = _merge_records(titls, 19) if titls else None
        edits = [r for r in jrnls if r[12:16] == "EDIT"]
        journal["editors"] = _merge_records(edits, 19).split(",") if edits else []
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
        journal["publisher"] = _merge_records(
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
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    REMARK records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

    remark_records = pdb_file.get_records_by_name("REMARK")
    if remark_records:
        remark_numbers = sorted(list(set([r[7:10] for r in remark_records])))
        remarks = []
        for number in remark_numbers:
            recs = [r for r in remark_records if r[7:10] == number]
            remark = {
             "number": number,
             "content": _merge_records(recs[1:], 11, join="\n", dont_condense=" ,:;")
            }
            remarks.append(remark)
        data_file._remarks = remarks
    else:
        data_file._remarks = []


def process_dbref_records(data_file, pdb_file):
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    DBREF records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

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


def process_seqadv_records(data_file, pdb_file):
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    SEQADV records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

    seqadvs = pdb_file.get_records_by_name("SEQADV")
    if seqadvs:
        data_file._sequence_differences = [{
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
    else:
        data_file._sequence_differences = []


def process_seqres_records(data_file, pdb_file):
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    SEQRES records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

    seqres = pdb_file.get_records_by_name("SEQRES")
    if seqres:
        chains = sorted(list(set([r[11] for r in seqres])))
        residue_sequences = []
        for chain in chains:
            records = [r for r in seqres if r[11] == chain]
            residue_sequences.append({
             "chain_id": chain,
             "length": records[0][13:17],
             "residues": _merge_records(records, 19).split()
            })
        data_file._residue_sequences = residue_sequences
    else:
        data_file._residue_sequences = []


def process_modres_records(data_file, pdb_file):
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    MODRES records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

    modres = pdb_file.get_records_by_name("MODRES")
    if modres:
        data_file._modified_residues = [{
         "residue_name": r[12:15],
         "chain_id": r[16],
         "residue_id": r[18:22],
         "insert_code": r[22] if r[22] else "",
         "standard_resisdue_name": r[24:27],
         "comment": r[29:70]
        } for r in modres]
    else:
        data_file._modified_residues = []


def process_het_records(data_file, pdb_file):
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    HET records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

    hets = pdb_file.get_records_by_name("HET")
    if hets:
        data_file._hets = [{
         "het_name": r[7:10],
         "chain_id": r[12],
         "het_id": r[13:17],
         "insert_code": r[17] if r[17] else "",
         "atom_num": r[20:25],
         "description": r[30:70]
        } for r in hets]
    else:
        data_file._hets = []


def process_hetnam_records(data_file, pdb_file):
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    HETNAM records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

    hetnams = pdb_file.get_records_by_name("HETNAM")
    if hetnams:
        ids = list(set([r[11:14] for r in hetnams]))
        data_file._het_names = {
         het_id: _merge_records(
          [r for r in hetnams if r[11:14] == het_id], 15, dont_condense=":;"
         ) for het_id in ids
        }
    else:
        data_file._het_names = {}


def process_hetsyn_records(data_file, pdb_file):
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    HETSYN records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

    hetsyns = pdb_file.get_records_by_name("HETSYN")
    if hetsyns:
        ids = list(set([r[11:14] for r in hetsyns]))
        data_file._het_synonyms = {
         het_id: _merge_records(
          [r for r in hetsyns if r[11:14] == het_id], 15
         ).split(";") for het_id in ids
        }
    else:
        data_file._het_synonyms = {}


def process_formul_records(data_file, pdb_file):
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    FORMUL records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

    formuls = pdb_file.get_records_by_name("FORMUL")
    if formuls:
        ids = list(set([r[12:15] for r in formuls]))
        data_file._formulae = {
         het_id: {
          "component_number": [r for r in formuls if r[12:15] == het_id][0][8:10],
          "is_water": [r for r in formuls if r[12:15] == het_id][0][18] == "*",
          "formula": _merge_records(
           [r for r in formuls if r[12:15] == het_id], 19
          )
         } for het_id in ids
        }
    else:
        data_file._formulae = {}


def process_helix_records(data_file, pdb_file):
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    HELIX records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

    helix = pdb_file.get_records_by_name("HELIX")
    if helix:
        data_file._helices = [{
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
    else:
        data_file._helices = []


def process_sheet_records(data_file, pdb_file):
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    SHEET records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

    sheet_records = pdb_file.get_records_by_name("SHEET")
    if sheet_records:
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
            data_file._sheets = sheets
    else:
        data_file._sheets = []


def process_ssbond_records(data_file, pdb_file):
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    SSBOND records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

    ssbonds = pdb_file.get_records_by_name("SSBOND")
    if ssbonds:
        data_file._ss_bonds = [{
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
    else:
        data_file._ss_bonds= []


def process_link_records(data_file, pdb_file):
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    LINK records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

    links = pdb_file.get_records_by_name("LINK")
    if links:
        data_file._links = [{
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
    else:
        data_file._links = []


def process_cispep_records(data_file, pdb_file):
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    CISPEP records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

    cispeps = pdb_file.get_records_by_name("CISPEP")
    if cispeps:
        data_file._cis_peptides = [{
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
    else:
        data_file._cis_peptides = []


def process_site_records(data_file, pdb_file):
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    SITE records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

    site_records = pdb_file.get_records_by_name("SITE")
    if site_records:
        site_names = sorted(list(set(
         [r.get_as_string(11, 14) for r in site_records]
        )))
        sites = []
        for site_name in site_names:
            records = [
             r for r in site_records if r.get_as_string(11, 14) == site_name
            ]
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
        data_file._sites = sites
    else:
        data_file._sites = []


def process_cryst1_records(data_file, pdb_file):
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    CRYST1 records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

    crystal = pdb_file.get_record_by_name("CRYST1")
    if crystal:
        data_file._crystal = {
         "a": crystal[6:15],
         "b": crystal[15:24],
         "c": crystal[24:33],
         "alpha": crystal[33:40],
         "beta": crystal[40:47],
         "gamma": crystal[47:54],
         "s_group": crystal[55:66],
         "z": crystal[66:70]
        }
    else:
        data_file._crystal = None


def process_origx_records(data_file, pdb_file):
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    ORIGX records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

    origx1 = pdb_file.get_record_by_name("ORIGX1")
    origx2 = pdb_file.get_record_by_name("ORIGX2")
    origx3 = pdb_file.get_record_by_name("ORIGX3")
    if origx1 or origx2 or origx3:
        data_file._origx = {
         "o11": origx1[10:20] if origx1 else None,
         "o12": origx1[20:30] if origx1 else None,
         "o13": origx1[30:40] if origx1 else None,
         "t1": origx1[45:55] if origx1 else None,
         "o21": origx2[10:20] if origx2 else None,
         "o22": origx2[20:30] if origx2 else None,
         "o23": origx2[30:40] if origx2 else None,
         "t2": origx2[45:55] if origx2 else None,
         "o31": origx3[10:20] if origx3 else None,
         "o32": origx3[20:30] if origx3 else None,
         "o33": origx3[30:40] if origx3 else None,
         "t3": origx3[45:55] if origx3 else None,
        }
    else:
        data_file._origix = None


def process_scale_records(data_file, pdb_file):
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    SCALE records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

    scale1 = pdb_file.get_record_by_name("SCALE1")
    scale2 = pdb_file.get_record_by_name("SCALE2")
    scale3 = pdb_file.get_record_by_name("SCALE3")
    if scale1 or scale2 or scale3:
        data_file._scale = {
         "s11": scale1[10:20] if scale1 else None,
         "s12": scale1[20:30] if scale1 else None,
         "s13": scale1[30:40] if scale1 else None,
         "u1": scale1[45:55] if scale1 else None,
         "s21": scale2[10:20] if scale2 else None,
         "s22": scale2[20:30] if scale2 else None,
         "s23": scale2[30:40] if scale2 else None,
         "u2": scale2[45:55] if scale2 else None,
         "s31": scale3[10:20] if scale3 else None,
         "s32": scale3[20:30] if scale3 else None,
         "s33": scale3[30:40] if scale3 else None,
         "u3": scale3[45:55] if scale3 else None,
        }
    else:
        data_file._scale = None


def process_mtrix_records(data_file, pdb_file):
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    MTRIX records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

    matrx1 = pdb_file.get_record_by_name("MTRIX1")
    matrx2 = pdb_file.get_record_by_name("MTRIX2")
    matrx3 = pdb_file.get_record_by_name("MTRIX3")
    if matrx1 or matrx2 or matrx3:
        data_file._matrix = {
         "serial_1": matrx1[7:10] if matrx1 else None,
         "m11": matrx1[10:20] if matrx1 else None,
         "m12": matrx1[20:30] if matrx1 else None,
         "m13": matrx1[30:40] if matrx1 else None,
         "v1": matrx1[45:55] if matrx1 else None,
         "i_given_1": (matrx1[59] == 1) if matrx1 else None,
         "serial_2": matrx2[7:10] if matrx2 else None,
         "m21": matrx2[10:20] if matrx2 else None,
         "m22": matrx2[20:30] if matrx2 else None,
         "m23": matrx2[30:40] if matrx2 else None,
         "v2": matrx2[45:55] if matrx2 else None,
         "i_given_2": (matrx2[59] == 1) if matrx2 else None,
         "serial_3": matrx3[7:10] if matrx3 else None,
         "m31": matrx3[10:20] if matrx3 else None,
         "m32": matrx3[20:30] if matrx3 else None,
         "m33": matrx3[30:40] if matrx3 else None,
         "v3": matrx3[45:55] if matrx3 else None,
         "i_given_3": (matrx3[59] == 1) if matrx3 else None,
        }
    else:
        data_file._matrix = None


def process_model_records(data_file, pdb_file):
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    MODEL records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

    model_records = pdb_file.get_records_by_name("MODEL")
    if model_records:
        endmdls = pdb_file.get_records_by_name("ENDMDL")
        pairs = list(zip(model_records, endmdls))
        data_file._models = [{
         "model_id": pair[0][10:14],
         "start_record": pair[0].number(),
         "end_record": pair[1].number(),
        } for pair in pairs]
    else:
        data_file._models = [{
         "model_id": 1,
         "start_record": -1,
         "end_record": -1
        }]


def process_atom_records(data_file, pdb_file):
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    ATOM records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

    atoms = pdb_file.get_records_by_name("ATOM")
    if atoms:
        data_file._atoms = [{
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
         "model_id": [m for m in data_file.models()
          if r.number() >= m["start_record"]
           and r.number() <= m["end_record"]][0]["model_id"] if len(
            data_file.models()
             ) > 1 else 1
        } for r in atoms]
    else:
        data_file._atoms = []


def process_anisou_records(data_file, pdb_file):
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    ANISOU records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

    anisou = pdb_file.get_records_by_name("ANISOU")
    if anisou:
        data_file._anisou = [{
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
         "model_id": [m for m in data_file.models()
          if r.number() >= m["start_record"]
           and r.number() <= m["end_record"]][0]["model_id"] if len(
            data_file.models()
             ) > 1 else 1
        } for r in anisou]
    else:
        data_file._anisou = []


def process_ter_records(data_file, pdb_file):
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    TER records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

    ters = pdb_file.get_records_by_name("TER")
    if ters:
        data_file._termini = [{
         "atom_id": r[6:11],
         "residue_name": r[17:20],
         "chain_id": r[21],
         "residue_id": r[22:26],
         "insert_code": r[26] if r[26] else "",
         "model_id": [m for m in data_file.models()
          if r.number() >= m["start_record"]
           and r.number() <= m["end_record"]][0]["model_id"] if len(
            data_file.models()
             ) > 1 else 1
        } for r in ters]
    else:
        data_file._termini = []


def process_hetatm_records(data_file, pdb_file):
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    HETATM records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

    hetatms = pdb_file.get_records_by_name("HETATM")
    if hetatms:
        data_file._heteroatoms = [{
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
         "model_id": [m for m in data_file.models()
          if r.number() >= m["start_record"]
           and r.number() <= m["end_record"]][0]["model_id"] if len(
            data_file.models()
             ) > 1 else 1
        } for r in hetatms]
    else:
        data_file._heteroatoms = []


def process_conect_records(data_file, pdb_file):
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    CONECT records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

    conects = pdb_file.get_records_by_name("CONECT")
    if conects:
        atom_ids = sorted(list(set([r[6:11] for r in conects])))
        data_file._connections = []
        for atom_id in atom_ids:
            records = [r for r in conects if r[6:11] == atom_id]
            bonded_atoms = []
            for record in records:
                if record[11:16]: bonded_atoms.append(record[11:16])
                if record[16:21]: bonded_atoms.append(record[16:21])
                if record[21:26]: bonded_atoms.append(record[21:26])
                if record[26:31]: bonded_atoms.append(record[26:31])
            data_file._connections.append({
             "atom_id": atom_id,
             "bonded_atoms": bonded_atoms
            })
    else:
        data_file._connections = []


def process_master_records(data_file, pdb_file):
    """Takes a :py:class:`.PdbDataFile` and updates it based on the
    MASTER records in the provided :py:class:`.PdbFile`

    :param PdbDataFile data_file: the Data File to update.
    :param PdbFile pdb_file: The source Pdb File"""

    master = pdb_file.get_record_by_name("MASTER")
    data_file._master = {
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


def _date_from_string(s):
    if s:
        return datetime.datetime.strptime(
         s, "%d-%b-%y"
        ).date()
    else:
        return None


def _merge_records(records, start, join=" ", dont_condense=""):
    string = join.join(
     str(record[start:]
    ) if record[start:] else "" for record in records)
    condense = [char for char in " ;:,-" if char not in dont_condense]
    for char in condense:
        string = string.replace(char + " ", char)
    return string


def _records_to_token_value_dicts(records):
    string = _merge_records(records, 10)
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
