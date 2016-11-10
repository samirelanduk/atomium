"""This module performs the actual parsing of the PDB file, though it does not
process the values that it extracts."""

import datetime

class PdbDataFile:
    """This object is essentially a list of values extracted from a PDB file. It
    functions as a data sheet.

    :param PdbFile pdb_file: The PDB file to extract information from."""

    def __init__(self, pdb_file=None):
        self._original_pdb_file = pdb_file
        _process_header_records(self)
        _process_obslte_records(self)
        _process_title_records(self)
        _process_split_records(self)
        _process_caveat_records(self)
        _process_compnd_records(self)
        _process_source_records(self)
        _process_keywd_records(self)
        _process_expdta_records(self)
        _process_nummdl_records(self)
        _process_mdltyp_records(self)
        _process_author_records(self)
        _process_revdat_records(self)
        _process_sprsde_records(self)
        _process_jrnl_records(self)
        _process_remark_records(self)
        _process_dbref_records(self)
        _process_seqadv_records(self)
        _process_seqres_records(self)
        _process_modres_records(self)
        _process_het_records(self)
        _process_hetnam_records(self)
        _process_hetsyn_records(self)
        _process_formul_records(self)
        _process_helix_records(self)
        _process_sheet_records(self)
        _process_ssbond_records(self)
        _process_link_records(self)
        _process_cispep_records(self)
        _process_site_records(self)
        '''model_records = self.pdb_file().get_records_by_name("MODEL")
        endmdls = self.pdb_file().get_records_by_name("ENDMDL")
        pairs = list(zip(model_records, endmdls))
        self._models = [{
         "model_id": pair[0][10:14],
         "start_record": pair[0].number(),
         "end_record": pair[1].number(),
        } for pair in pairs]
        if not pairs:
            self._models = [{
             "model_id": 1,
             "start_record": 0,
             "end_record": len(self.pdb_file().records()),
            }]'''


    def __repr__(self):
        return "<PdbDataFile (%s)>" % (self.pdb_code() if self.pdb_code() else "????")


    def original_pdb_file(self):
        """The :py:class:`.PdbFile` from which the object was created.

        :rtype: ``PdbFile``"""

        return self._original_pdb_file


    def classification(self, classification=None):
        if classification:
            if not isinstance(classification, str):
                raise TypeError(
                 "classification must be str, not '%s'" % str(classification)
                )
            if len(classification) > 40:
                raise ValueError(
                 "classification must be <40 chars, not '%s'" % classification
                )
            self._classification = classification
        else:
            return self._classification


    def deposition_date(self, deposition_date=None):
        if deposition_date:
            if not isinstance(deposition_date, datetime.date):
                raise TypeError(
                 "deposition_date must be date, not '%s'" % str(deposition_date)
                )
            self._deposition_date = deposition_date
        else:
            return self._deposition_date


    def pdb_code(self, pdb_code=None):
        if pdb_code:
            if not isinstance(pdb_code, str):
                raise TypeError(
                 "pdb_code must be str, not '%s'" % str(pdb_code)
                )
            if len(pdb_code) != 4:
                raise ValueError(
                 "pdb_code must be 4 chars, not '%s'" % pdb_code
                )
            self._pdb_code = pdb_code
        else:
            return self._pdb_code


    def is_obsolete(self, is_obsolete=None):
        if is_obsolete is not None:
            if not isinstance(is_obsolete, bool):
                raise TypeError(
                 "is_obsolete must be bool, not '%s'" % str(is_obsolete)
                )
            self._is_obsolete = is_obsolete
        else:
            return self._is_obsolete


    def obsolete_date(self, obsolete_date=None):
        if obsolete_date:
            if not isinstance(obsolete_date, datetime.date):
                raise TypeError(
                 "obsolete_date must be date, not '%s'" % str(obsolete_date)
                )
            self._obsolete_date = obsolete_date
        else:
            return self._obsolete_date


    def replacement_code(self, replacement_code=None):
        if replacement_code:
            if not isinstance(replacement_code, str):
                raise TypeError(
                 "replacement_code must be str, not '%s'" % str(replacement_code)
                )
            if len(replacement_code) != 4:
                raise ValueError(
                 "replacement_code must be 4 chars, not '%s'" % replacement_code
                )
            self._replacement_code = replacement_code
        else:
            return self._replacement_code


    def title(self, title=None):
        if title:
            if not isinstance(title, str):
                raise TypeError("title must be str, not '%s'" % str(title))
            self._title = title
        else:
            return self._title


    def split_codes(self):
        return self._split_codes


    def caveat(self, caveat=None):
        if caveat:
            if not isinstance(caveat, str):
                raise TypeError("caveat must be str, not '%s'" % str(caveat))
            self._caveat = caveat
        else:
            return self._caveat


    def compounds(self):
        return self._compounds


    def sources(self):
        return self._sources


    def keywords(self):
        return self._keywords


    def experimental_techniques(self):
        return self._experimental_techniques


    def model_count(self, model_count=None):
        if model_count:
            if not isinstance(model_count, int):
                raise TypeError(
                 "model_count must be int, not '%s'" % str(model_count)
                )
            self._model_count = model_count
        else:
            return self._model_count


    def model_annotations(self):
        return self._model_annotations


    def authors(self):
        return self._authors


    def revisions(self):
        return self._revisions


    def supercedes(self):
        return self._supercedes


    def supercede_date(self, supercede_date=None):
        if supercede_date:
            if not isinstance(supercede_date, datetime.date):
                raise TypeError(
                 "supercede_date must be date, not '%s'" % str(supercede_date)
                )
            self._supercede_date = supercede_date
        else:
            return self._supercede_date


    def journal(self, journal=None):
        if journal:
            if not isinstance(journal, dict):
                raise TypeError(
                 "journal must be dict, not '%s'" % str(journal)
                )
            self._journal = journal
        else:
            return self._journal


    def remarks(self):
        return self._remarks


    def get_remark_by_number(self, number):
        for remark in self.remarks():
            if remark["number"] == number:
                return remark


    def dbreferences(self):
        return self._dbreferences


    def sequence_differences(self):
        return self._sequence_differences


    def residue_sequences(self):
        return self._residue_sequences


    def modified_residues(self):
        return self._modified_residues


    def hets(self):
        return self._hets


    def het_names(self):
        return self._het_names


    def het_synonyms(self):
        return self._het_synonyms


    def formulae(self):
        return self._formulae


    def helices(self):
        return self._helices


    def sheets(self):
        return self._sheets


    def ss_bonds(self):
        return self._ss_bonds


    def links(self):
        return self._links


    def cis_peptides(self):
        return self._cis_peptides


    def sites(self):
        return self._sites



def _process_header_records(data_file):
    if data_file.original_pdb_file():
        header = data_file.original_pdb_file().get_record_by_name("HEADER")
        if header:
            data_file._classification = header[10:50]
            data_file._deposition_date = date_from_string(header[50:59])
            data_file._pdb_code = header[62:66]
            return
    data_file._classification = None
    data_file._deposition_date = None
    data_file._pdb_code = None


def _process_obslte_records(data_file):
    if data_file.original_pdb_file():
        obslte = data_file.original_pdb_file().get_record_by_name("OBSLTE")
        if obslte:
            data_file._is_obsolete = True
            data_file._obsolete_date = date_from_string(obslte[11:20])
            data_file._replacement_code = obslte[31:35]
            return
    data_file._is_obsolete = False
    data_file._obsolete_date = None
    data_file._replacement_code = None


def _process_title_records(data_file):
    if data_file.original_pdb_file():
        titles = data_file.original_pdb_file().get_records_by_name("TITLE")
        if titles:
            data_file._title = merge_records(titles, 10, dont_condense=",;:-")
            return
    data_file._title = None


def _process_split_records(data_file):
    if data_file.original_pdb_file():
        splits = data_file.original_pdb_file().get_records_by_name("SPLIT")
        data_file._split_codes = " ".join([r[10:] for r in splits]).split()
        return
    data_file._split_codes = []


def _process_caveat_records(data_file):
    if data_file.original_pdb_file():
        caveats = data_file.original_pdb_file().get_records_by_name("CAVEAT")
        if caveats:
            data_file._caveat = merge_records(caveats, 19)
            return
    data_file._caveat = None


def _process_compnd_records(data_file):
    if data_file.original_pdb_file():
        records = data_file.original_pdb_file().get_records_by_name("COMPND")
        data_file._compounds = records_to_token_value_dicts(records)
        return
    data_file._compounds = []


def _process_source_records(data_file):
    if data_file.original_pdb_file():
        records = data_file.original_pdb_file().get_records_by_name("SOURCE")
        data_file._sources = records_to_token_value_dicts(records)
        return
    data_file._sources = []


def _process_keywd_records(data_file):
    if data_file.original_pdb_file():
        keywords = data_file.original_pdb_file().get_records_by_name("KEYWDS")
        if keywords:
            keyword_text = merge_records(keywords, 10)
            data_file._keywords = keyword_text.split(",")
            return
    data_file._keywords = []


def _process_expdta_records(data_file):
    if data_file.original_pdb_file():
        expdta = data_file.original_pdb_file().get_records_by_name("EXPDTA")
        if expdta:
            expdta_text = merge_records(expdta, 10)
            data_file._experimental_techniques = expdta_text.split(";")
            return
    data_file._experimental_techniques = []


def _process_nummdl_records(data_file):
    if data_file.original_pdb_file():
        nummdl = data_file.original_pdb_file().get_record_by_name("NUMMDL")
        if nummdl:
            data_file._model_count = nummdl[10:14]
            return
    data_file._model_count = 1


def _process_mdltyp_records(data_file):
    if data_file.original_pdb_file():
        mdltyps = data_file.original_pdb_file().get_records_by_name("MDLTYP")
        if mdltyps:
            mdltyp_text = merge_records(mdltyps, 10, dont_condense=",")
            data_file._model_annotations = [
             ann.strip() for ann in mdltyp_text.split(";") if ann.strip()
            ]
            return
    data_file._model_annotations = []


def _process_author_records(data_file):
    if data_file.original_pdb_file():
        authors = data_file.original_pdb_file().get_records_by_name("AUTHOR")
        if authors:
            data_file._authors = merge_records(authors, 10).split(",")
            return
    data_file._authors = []


def _process_revdat_records(data_file):
    if data_file.original_pdb_file():
        revdats = data_file.original_pdb_file().get_records_by_name("REVDAT")
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
            return
    data_file._revisions = []


def _process_sprsde_records(data_file):
    if data_file.original_pdb_file():
        sprsde = data_file.original_pdb_file().get_record_by_name("SPRSDE")
        if sprsde:
            data_file._supercedes = sprsde[31:75].split()
            data_file._supercede_date = date_from_string(sprsde[11:20])
            return
    data_file._supercedes = []
    data_file._supercede_date = None


def _process_jrnl_records(data_file):
    if data_file.original_pdb_file():
        jrnls = data_file.original_pdb_file().get_records_by_name("JRNL")
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
            data_file._journal = journal
            return
    else:
        data_file._journal = None


def _process_remark_records(data_file):
    if data_file.original_pdb_file():
        remark_records = data_file.original_pdb_file().get_records_by_name("REMARK")
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
            return
    data_file._remarks = []


def _process_dbref_records(data_file):
    if data_file.original_pdb_file():
        dbrefs = data_file.original_pdb_file().get_records_by_name("DBREF")
        dbref1s = data_file.original_pdb_file().get_records_by_name("DBREF1")
        dbref2s = data_file.original_pdb_file().get_records_by_name("DBREF2")
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
            return
    data_file._dbreferences = []


def _process_seqadv_records(data_file):
    if data_file.original_pdb_file():
        seqadvs = data_file.original_pdb_file().get_records_by_name("SEQADV")
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
            return
    data_file._sequence_differences = []


def _process_seqres_records(data_file):
    if data_file.original_pdb_file():
        seqres = data_file.original_pdb_file().get_records_by_name("SEQRES")
        if seqres:
            chains = sorted(list(set([r[11] for r in seqres])))
            residue_sequences = []
            for chain in chains:
                records = [r for r in seqres if r[11] == chain]
                residue_sequences.append({
                 "chain_id": chain,
                 "length": records[0][13:17],
                 "residues": merge_records(records, 19).split()
                })
            data_file._residue_sequences = residue_sequences
            return
    data_file._residue_sequences = []


def _process_modres_records(data_file):
    if data_file.original_pdb_file():
        modres = data_file.original_pdb_file().get_records_by_name("MODRES")
        if modres:
            data_file._modified_residues = [{
             "residue_name": r[12:15],
             "chain_id": r[16],
             "residue_id": r[18:22],
             "insert_code": r[22] if r[22] else "",
             "standard_resisdue_name": r[24:27],
             "comment": r[29:70]
            } for r in modres]
            return
    data_file._modified_residues = []


def _process_het_records(data_file):
    if data_file.original_pdb_file():
        hets = data_file.original_pdb_file().get_records_by_name("HET")
        if hets:
            data_file._hets = [{
             "het_name": r[7:10],
             "chain_id": r[12],
             "het_id": r[13:17],
             "insert_code": r[17] if r[17] else "",
             "atom_num": r[20:25],
             "description": r[30:70]
            } for r in hets]
            return
    data_file._hets = []


def _process_hetnam_records(data_file):
    if data_file.original_pdb_file():
        hetnams = data_file.original_pdb_file().get_records_by_name("HETNAM")
        if hetnams:
            ids = list(set([r[11:14] for r in hetnams]))
            data_file._het_names = {
             het_id: merge_records(
              [r for r in hetnams if r[11:14] == het_id], 15, dont_condense=":;"
             ) for het_id in ids
            }
            return
    data_file._het_names = {}


def _process_hetsyn_records(data_file):
    if data_file.original_pdb_file():
        hetsyns = data_file.original_pdb_file().get_records_by_name("HETSYN")
        if hetsyns:
            ids = list(set([r[11:14] for r in hetsyns]))
            data_file._het_synonyms = {
             het_id: merge_records(
              [r for r in hetsyns if r[11:14] == het_id], 15
             ).split(";") for het_id in ids
            }
            return
    data_file._het_synonyms = {}


def _process_formul_records(data_file):
    if data_file.original_pdb_file():
        formuls = data_file.original_pdb_file().get_records_by_name("FORMUL")
        if formuls:
            ids = list(set([r[12:15] for r in formuls]))
            data_file._formulae = {
             het_id: {
              "component_number": [r for r in formuls if r[12:15] == het_id][0][8:10],
              "is_water": [r for r in formuls if r[12:15] == het_id][0][18] == "*",
              "formula": merge_records(
               [r for r in formuls if r[12:15] == het_id], 19
              )
             } for het_id in ids
            }
            return
    data_file._formulae = {}


def _process_helix_records(data_file):
    if data_file.original_pdb_file():
        helix = data_file.original_pdb_file().get_records_by_name("HELIX")
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
            return
    data_file._helices = []


def _process_sheet_records(data_file):
    if data_file.original_pdb_file():
        sheet_records = data_file.original_pdb_file().get_records_by_name("SHEET")
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
                return
    data_file._sheets = []


def _process_ssbond_records(data_file):
    if data_file.original_pdb_file():
        ssbonds = data_file.original_pdb_file().get_records_by_name("SSBOND")
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
            return
    data_file._ss_bonds= []


def _process_link_records(data_file):
    if data_file.original_pdb_file():
        links = data_file.original_pdb_file().get_records_by_name("LINK")
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
            return
    data_file._links = []


def _process_cispep_records(data_file):
    if data_file.original_pdb_file():
        cispeps = data_file.original_pdb_file().get_records_by_name("CISPEP")
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
            return
    data_file._cis_peptides = []


def _process_site_records(data_file):
    if data_file.original_pdb_file():
        site_records = data_file.original_pdb_file().get_records_by_name("SITE")
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
            return
    data_file._sites = []


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


'''
    def sites(self):
        site_records = self.pdb_file().get_records_by_name("SITE")

        return sites


    def crystal_a(self):
        crystal = self.pdb_file().get_record_by_name("CRYST1")
        return crystal[6:15] if crystal else None


    def crystal_b(self):
        crystal = self.pdb_file().get_record_by_name("CRYST1")
        return crystal[15:24] if crystal else None


    def crystal_c(self):
        crystal = self.pdb_file().get_record_by_name("CRYST1")
        return crystal[24:33] if crystal else None


    def crystal_alpha(self):
        crystal = self.pdb_file().get_record_by_name("CRYST1")
        return crystal[33:40] if crystal else None


    def crystal_beta(self):
        crystal = self.pdb_file().get_record_by_name("CRYST1")
        return crystal[40:47] if crystal else None


    def crystal_gamma(self):
        crystal = self.pdb_file().get_record_by_name("CRYST1")
        return crystal[47:54] if crystal else None


    def crystal_s_group(self):
        crystal = self.pdb_file().get_record_by_name("CRYST1")
        return crystal[55:66] if crystal else None


    def crystal_z(self):
        crystal = self.pdb_file().get_record_by_name("CRYST1")
        return crystal[66:70] if crystal else None


    def origx(self, o):
        if o[1] == "1":
            origx1 = self.pdb_file().get_record_by_name("ORIGX1")
            if origx1:
                return {
                 "o11": origx1[10:20],
                 "o12": origx1[20:30],
                 "o13": origx1[30:40],
                 "t1": origx1[45:55]
                }.get(o)
            else:
                return None
        elif o[1] == "2":
            origx2 = self.pdb_file().get_record_by_name("ORIGX2")
            if origx2:
                return {
                 "o21": origx2[10:20],
                 "o22": origx2[20:30],
                 "o23": origx2[30:40],
                 "t2": origx2[45:55]
                }.get(o)
            else:
                return None
        elif o[1] == "3":
            origx3 = self.pdb_file().get_record_by_name("ORIGX3")
            if origx3:
                return {
                 "o31": origx3[10:20],
                 "o32": origx3[20:30],
                 "o33": origx3[30:40],
                 "t3": origx3[45:55]
                }.get(o)
            else:
                return None


    def scale(self, s):
        if s[1] == "1":
            scale1 = self.pdb_file().get_record_by_name("SCALE1")
            if scale1:
                return {
                 "s11": scale1[10:20],
                 "s12": scale1[20:30],
                 "s13": scale1[30:40],
                 "u1": scale1[45:55]
                }.get(s)
            else:
                return None
        elif s[1] == "2":
            scale2 = self.pdb_file().get_record_by_name("SCALE2")
            if scale2:
                return {
                 "s21": scale2[10:20],
                 "s22": scale2[20:30],
                 "s23": scale2[30:40],
                 "u2": scale2[45:55]
                }.get(s)
            else:
                return None
        elif s[1] == "3":
            scale3 = self.pdb_file().get_record_by_name("SCALE3")
            if scale3:
                return {
                 "s31": scale3[10:20],
                 "s32": scale3[20:30],
                 "s33": scale3[30:40],
                 "u3": scale3[45:55]
                }.get(s)
            else:
                return None


    def matrix(self, m):
        if m[1] == "1" or m[-2:] == "_1":
            matrix1 = self.pdb_file().get_record_by_name("MTRIX1")
            if matrix1:
                return {
                 "serial_1": matrix1[7:10],
                 "m11": matrix1[10:20],
                 "m12": matrix1[20:30],
                 "m13": matrix1[30:40],
                 "v1": matrix1[45:55],
                 "i_given_1": matrix1[59] == 1
                }.get(m)
            else:
                return None if "given" not in m else False
        elif m[1] == "2" or m[-2:] == "_2":
            matrix2 = self.pdb_file().get_record_by_name("MTRIX2")
            if matrix2:
                return {
                 "serial_2": matrix2[7:10],
                 "m21": matrix2[10:20],
                 "m22": matrix2[20:30],
                 "m23": matrix2[30:40],
                 "v2": matrix2[45:55],
                 "i_given_1": matrix2[59] == 1
                }.get(m)
            else:
                return None if "given" not in m else False
        elif m[1] == "3" or m[-2:] == "_3":
            matrix3 = self.pdb_file().get_record_by_name("MTRIX3")
            if matrix3:
                return {
                 "serial_3": matrix3[7:10],
                 "m31": matrix3[10:20],
                 "m32": matrix3[20:30],
                 "m33": matrix3[30:40],
                 "v3": matrix3[45:55],
                 "i_given_3": matrix3[59] == 1
                }.get(m)
            else:
                return None if "given" not in m else False


    def models(self):
        return self._models


    def atoms(self):
        atoms = self.pdb_file().get_records_by_name("ATOM")
        return [{
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
         "model_id": [m for m in self.models()
          if r.number() >= m["start_record"]
           and r.number() <= m["end_record"]][0]["model_id"]
        } for r in atoms]


    def anisou(self):
        anisou = self.pdb_file().get_records_by_name("ANISOU")
        return [{
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
         "model_id": [m for m in self.models()
          if r.number() >= m["start_record"]
           and r.number() <= m["end_record"]][0]["model_id"]
        } for r in anisou]


    def termini(self):
        ters = self.pdb_file().get_records_by_name("TER")
        return [{
         "atom_id": r[6:11],
         "residue_name": r[17:20],
         "chain_id": r[21],
         "residue_id": r[22:26],
         "insert_code": r[26] if r[26] else "",
         "model_id": [m for m in self.models()
          if r.number() >= m["start_record"]
           and r.number() <= m["end_record"]][0]["model_id"]
        } for r in ters]


    def heteroatoms(self):
        hetatms = self.pdb_file().get_records_by_name("HETATM")
        return [{
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
         "model_id": [m for m in self.models()
          if r.number() >= m["start_record"]
           and r.number() <= m["end_record"]][0]["model_id"]
        } for r in hetatms]


    def connections(self):
        conects = self.pdb_file().get_records_by_name("CONECT")
        atom_ids = sorted(list(set([r[6:11] for r in conects])))
        return [{
         "atom_id": num,
         "bonded_atoms": [int(n) for n in merge_records([
          r for r in conects if r[6:11] == num
         ], 11).split()]
        } for num in atom_ids]


    def master(self):
        master = self.pdb_file().get_record_by_name("MASTER")
        return {
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


'''
