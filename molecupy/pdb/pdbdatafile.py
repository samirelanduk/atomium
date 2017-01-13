"""This module performs the actual parsing of the PDB file, though it does not
process the values that it extracts."""

import datetime
import math
from .pdbfile import PdbRecord
from ..converters.pdbdatafile2pdbfile import pdb_file_from_pdb_data_file

class PdbDataFile:
    """This object is essentially a list of values extracted from a PDB file. It
    functions as a data sheet.

    :param PdbFile pdb_file: The PDB file to extract information from."""

    def __init__(self, pdb_file=None):
        '''self._original_pdb_file = pdb_file
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
        _process_cryst1_records(self)
        _process_origx_records(self)
        _process_scale_records(self)
        _process_matrix_records(self)
        _process_model_records(self)
        _process_atom_records(self)
        _process_anisou_records(self)
        _process_ter_records(self)
        _process_hetatm_records(self)
        _process_conect_records(self)
        _process_master_records(self)'''


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


    def crystal_a(self, crystal_a=None):
        if crystal_a is not None:
            if not isinstance(crystal_a, float):
                raise TypeError(
                 "crystal_a must be float, not '%s'" % str(crystal_a)
                )
            self._crystal_a = crystal_a
        else:
            return self._crystal_a


    def crystal_b(self, crystal_b=None):
        if crystal_b is not None:
            if not isinstance(crystal_b, float):
                raise TypeError(
                 "crystal_b must be float, not '%s'" % str(crystal_b)
                )
            self._crystal_b = crystal_b
        else:
            return self._crystal_b


    def crystal_c(self, crystal_c=None):
        if crystal_c is not None:
            if not isinstance(crystal_c, float):
                raise TypeError(
                 "crystal_c must be float, not '%s'" % str(crystal_c)
                )
            self._crystal_c = crystal_c
        else:
            return self._crystal_c


    def crystal_alpha(self, crystal_alpha=None):
        if crystal_alpha is not None:
            if not isinstance(crystal_alpha, float):
                raise TypeError(
                 "crystal_alpha must be float, not '%s'" % str(crystal_alpha)
                )
            self._crystal_alpha = crystal_alpha
        else:
            return self._crystal_alpha


    def crystal_beta(self, crystal_beta=None):
        if crystal_beta is not None:
            if not isinstance(crystal_beta, float):
                raise TypeError(
                 "crystal_beta must be float, not '%s'" % str(crystal_beta)
                )
            self._crystal_beta = crystal_beta
        else:
            return self._crystal_beta


    def crystal_gamma(self, crystal_gamma=None):
        if crystal_gamma is not None:
            if not isinstance(crystal_gamma, float):
                raise TypeError(
                 "crystal_gamma must be float, not '%s'" % str(crystal_gamma)
                )
            self._crystal_gamma = crystal_gamma
        else:
            return self._crystal_gamma


    def crystal_s_group(self, crystal_s_group=None):
        if crystal_s_group is not None:
            if not isinstance(crystal_s_group, str):
                raise TypeError(
                 "crystal_s_group must be str, not '%s'" % str(crystal_s_group)
                )
            self._crystal_s_group = crystal_s_group
        else:
            return self._crystal_s_group


    def crystal_z(self, crystal_z=None):
        if crystal_z is not None:
            if not isinstance(crystal_z, int):
                raise TypeError(
                 "crystal_z must be int, not '%s'" % str(crystal_z)
                )
            self._crystal_z = crystal_z
        else:
            return self._crystal_z


    def origx_o11(self, origx_o11=None):
        if origx_o11 is not None:
            if not isinstance(origx_o11, float):
                raise TypeError("origx_o11 must be float, not '%s'" % str(origx_o11))
            self._origx_o11 = origx_o11
        else:
            return self._origx_o11


    def origx_o12(self, origx_o12=None):
        if origx_o12 is not None:
            if not isinstance(origx_o12, float):
                raise TypeError("origx_o12 must be float, not '%s'" % str(origx_o12))
            self._origx_o12 = origx_o12
        else:
            return self._origx_o12


    def origx_o13(self, origx_o13=None):
        if origx_o13 is not None:
            if not isinstance(origx_o13, float):
                raise TypeError("origx_o13 must be float, not '%s'" % str(origx_o13))
            self._origx_o13 = origx_o13
        else:
            return self._origx_o13


    def origx_t1(self, origx_t1=None):
        if origx_t1 is not None:
            if not isinstance(origx_t1, float):
                raise TypeError("origx_t1 must be float, not '%s'" % str(origx_t1))
            self._origx_t1 = origx_t1
        else:
            return self._origx_t1


    def origx_o21(self, origx_o21=None):
        if origx_o21 is not None:
            if not isinstance(origx_o21, float):
                raise TypeError("origx_o21 must be float, not '%s'" % str(origx_o21))
            self._origx_o21 = origx_o21
        else:
            return self._origx_o21


    def origx_o22(self, origx_o22=None):
        if origx_o22 is not None:
            if not isinstance(origx_o22, float):
                raise TypeError("origx_o22 must be float, not '%s'" % str(origx_o22))
            self._origx_o22 = origx_o22
        else:
            return self._origx_o22


    def origx_o23(self, origx_o23=None):
        if origx_o23 is not None:
            if not isinstance(origx_o23, float):
                raise TypeError("origx_o23 must be float, not '%s'" % str(origx_o23))
            self._origx_o23 = origx_o23
        else:
            return self._origx_o23


    def origx_t2(self, origx_t2=None):
        if origx_t2 is not None:
            if not isinstance(origx_t2, float):
                raise TypeError("origx_t2 must be float, not '%s'" % str(origx_t2))
            self._origx_t2 = origx_t2
        else:
            return self._origx_t2


    def origx_o31(self, origx_o31=None):
        if origx_o31 is not None:
            if not isinstance(origx_o31, float):
                raise TypeError("origx_o31 must be float, not '%s'" % str(origx_o31))
            self._origx_o31 = origx_o31
        else:
            return self._origx_o31


    def origx_o32(self, origx_o32=None):
        if origx_o32 is not None:
            if not isinstance(origx_o32, float):
                raise TypeError("origx_o32 must be float, not '%s'" % str(origx_o32))
            self._origx_o32 = origx_o32
        else:
            return self._origx_o32


    def origx_o33(self, origx_o33=None):
        if origx_o33 is not None:
            if not isinstance(origx_o33, float):
                raise TypeError("origx_o33 must be float, not '%s'" % str(origx_o33))
            self._origx_o33 = origx_o33
        else:
            return self._origx_o33


    def origx_t3(self, origx_t3=None):
        if origx_t3 is not None:
            if not isinstance(origx_t3, float):
                raise TypeError("origx_t3 must be float, not '%s'" % str(origx_t3))
            self._origx_t3 = origx_t3
        else:
            return self._origx_t3


    def scale_s11(self, scale_s11=None):
        if scale_s11 is not None:
            if not isinstance(scale_s11, float):
                raise TypeError("scale_s11 must be float, not '%s'" % str(scale_s11))
            self._scale_s11 = scale_s11
        else:
            return self._scale_s11


    def scale_s12(self, scale_s12=None):
        if scale_s12 is not None:
            if not isinstance(scale_s12, float):
                raise TypeError("scale_s12 must be float, not '%s'" % str(scale_s12))
            self._scale_s12 = scale_s12
        else:
            return self._scale_s12


    def scale_s13(self, scale_s13=None):
        if scale_s13 is not None:
            if not isinstance(scale_s13, float):
                raise TypeError("scale_s13 must be float, not '%s'" % str(scale_s13))
            self._scale_s13 = scale_s13
        else:
            return self._scale_s13


    def scale_u1(self, scale_u1=None):
        if scale_u1 is not None:
            if not isinstance(scale_u1, float):
                raise TypeError("scale_u1 must be float, not '%s'" % str(scale_u1))
            self._scale_u1 = scale_u1
        else:
            return self._scale_u1


    def scale_s21(self, scale_s21=None):
        if scale_s21 is not None:
            if not isinstance(scale_s21, float):
                raise TypeError("scale_s21 must be float, not '%s'" % str(scale_s21))
            self._scale_s21 = scale_s21
        else:
            return self._scale_s21


    def scale_s22(self, scale_s22=None):
        if scale_s22 is not None:
            if not isinstance(scale_s22, float):
                raise TypeError("scale_s22 must be float, not '%s'" % str(scale_s22))
            self._scale_s22 = scale_s22
        else:
            return self._scale_s22


    def scale_s23(self, scale_s23=None):
        if scale_s23 is not None:
            if not isinstance(scale_s23, float):
                raise TypeError("scale_s23 must be float, not '%s'" % str(scale_s23))
            self._scale_s23 = scale_s23
        else:
            return self._scale_s23


    def scale_u2(self, scale_u2=None):
        if scale_u2 is not None:
            if not isinstance(scale_u2, float):
                raise TypeError("scale_u2 must be float, not '%s'" % str(scale_u2))
            self._scale_u2 = scale_u2
        else:
            return self._scale_u2


    def scale_s31(self, scale_s31=None):
        if scale_s31 is not None:
            if not isinstance(scale_s31, float):
                raise TypeError("scale_s31 must be float, not '%s'" % str(scale_s31))
            self._scale_s31 = scale_s31
        else:
            return self._scale_s31


    def scale_s32(self, scale_s32=None):
        if scale_s32 is not None:
            if not isinstance(scale_s32, float):
                raise TypeError("scale_s32 must be float, not '%s'" % str(scale_s32))
            self._scale_s32 = scale_s32
        else:
            return self._scale_s32


    def scale_s33(self, scale_s33=None):
        if scale_s33 is not None:
            if not isinstance(scale_s33, float):
                raise TypeError("scale_s33 must be float, not '%s'" % str(scale_s33))
            self._scale_s33 = scale_s33
        else:
            return self._scale_s33


    def scale_u3(self, scale_u3=None):
        if scale_u3 is not None:
            if not isinstance(scale_u3, float):
                raise TypeError("scale_u3 must be float, not '%s'" % str(scale_u3))
            self._scale_u3 = scale_u3
        else:
            return self._scale_u3


    def matrix_serial_1(self, matrix_serial_1=None):
        if matrix_serial_1 is not None:
            if not isinstance(matrix_serial_1, int):
                raise TypeError("matrix_serial_1 must be int, not '%s'" % str(matrix_serial_1))
            self._matrix_serial_1 = matrix_serial_1
        else:
            return self._matrix_serial_1


    def matrix_m11(self, matrix_m11=None):
        if matrix_m11 is not None:
            if not isinstance(matrix_m11, float):
                raise TypeError("matrix_m11 must be float, not '%s'" % str(matrix_m11))
            self._matrix_m11 = matrix_m11
        else:
            return self._matrix_m11


    def matrix_m12(self, matrix_m12=None):
        if matrix_m12 is not None:
            if not isinstance(matrix_m12, float):
                raise TypeError("matrix_m12 must be float, not '%s'" % str(matrix_m12))
            self._matrix_m12 = matrix_m12
        else:
            return self._matrix_m12


    def matrix_m13(self, matrix_m13=None):
        if matrix_m13 is not None:
            if not isinstance(matrix_m13, float):
                raise TypeError("matrix_m13 must be float, not '%s'" % str(matrix_m13))
            self._matrix_m13 = matrix_m13
        else:
            return self._matrix_m13


    def matrix_v1(self, matrix_v1=None):
        if matrix_v1 is not None:
            if not isinstance(matrix_v1, float):
                raise TypeError("matrix_v1 must be float, not '%s'" % str(matrix_v1))
            self._matrix_v1 = matrix_v1
        else:
            return self._matrix_v1


    def matrix_i_given_1(self, matrix_i_given_1=None):
        if matrix_i_given_1 is not None:
            if not isinstance(matrix_i_given_1, bool):
                raise TypeError("matrix_i_given_1 must be bool, not '%s'" % str(matrix_i_given_1))
            self._matrix_i_given_1 = matrix_i_given_1
        else:
            return self._matrix_i_given_1


    def matrix_serial_2(self, matrix_serial_2=None):
        if matrix_serial_2 is not None:
            if not isinstance(matrix_serial_2, int):
                raise TypeError("matrix_serial_2 must be int, not '%s'" % str(matrix_serial_2))
            self._matrix_serial_2 = matrix_serial_2
        else:
            return self._matrix_serial_2


    def matrix_m21(self, matrix_m21=None):
        if matrix_m21 is not None:
            if not isinstance(matrix_m21, float):
                raise TypeError("matrix_m21 must be float, not '%s'" % str(matrix_m21))
            self._matrix_m21 = matrix_m21
        else:
            return self._matrix_m21


    def matrix_m22(self, matrix_m22=None):
        if matrix_m22 is not None:
            if not isinstance(matrix_m22, float):
                raise TypeError("matrix_m22 must be float, not '%s'" % str(matrix_m22))
            self._matrix_m22 = matrix_m22
        else:
            return self._matrix_m22


    def matrix_m23(self, matrix_m23=None):
        if matrix_m23 is not None:
            if not isinstance(matrix_m23, float):
                raise TypeError("matrix_m23 must be float, not '%s'" % str(matrix_m23))
            self._matrix_m23 = matrix_m23
        else:
            return self._matrix_m23


    def matrix_v2(self, matrix_v2=None):
        if matrix_v2 is not None:
            if not isinstance(matrix_v2, float):
                raise TypeError("matrix_v2 must be float, not '%s'" % str(matrix_v2))
            self._matrix_v2 = matrix_v2
        else:
            return self._matrix_v2


    def matrix_i_given_2(self, matrix_i_given_2=None):
        if matrix_i_given_2 is not None:
            if not isinstance(matrix_i_given_2, bool):
                raise TypeError("matrix_i_given_2 must be bool, not '%s'" % str(matrix_i_given_2))
            self._matrix_i_given_2 = matrix_i_given_2
        else:
            return self._matrix_i_given_2


    def matrix_serial_3(self, matrix_serial_3=None):
        if matrix_serial_3 is not None:
            if not isinstance(matrix_serial_3, int):
                raise TypeError("matrix_serial_3 must be int, not '%s'" % str(matrix_serial_3))
            self._matrix_serial_3 = matrix_serial_3
        else:
            return self._matrix_serial_3


    def matrix_m31(self, matrix_m31=None):
        if matrix_m31 is not None:
            if not isinstance(matrix_m31, float):
                raise TypeError("matrix_m31 must be float, not '%s'" % str(matrix_m31))
            self._matrix_m31 = matrix_m31
        else:
            return self._matrix_m31


    def matrix_m32(self, matrix_m32=None):
        if matrix_m32 is not None:
            if not isinstance(matrix_m32, float):
                raise TypeError("matrix_m32 must be float, not '%s'" % str(matrix_m32))
            self._matrix_m32 = matrix_m32
        else:
            return self._matrix_m32


    def matrix_m33(self, matrix_m33=None):
        if matrix_m33 is not None:
            if not isinstance(matrix_m33, float):
                raise TypeError("matrix_m33 must be float, not '%s'" % str(matrix_m33))
            self._matrix_m33 = matrix_m33
        else:
            return self._matrix_m33


    def matrix_v3(self, matrix_v3=None):
        if matrix_v3 is not None:
            if not isinstance(matrix_v3, float):
                raise TypeError("matrix_v3 must be float, not '%s'" % str(matrix_v3))
            self._matrix_v3 = matrix_v3
        else:
            return self._matrix_v3


    def matrix_i_given_3(self, matrix_i_given_3=None):
        if matrix_i_given_3 is not None:
            if not isinstance(matrix_i_given_3, bool):
                raise TypeError("matrix_i_given_3 must be bool, not '%s'" % str(matrix_i_given_3))
            self._matrix_i_given_3 = matrix_i_given_3
        else:
            return self._matrix_i_given_3


    def models(self):
        return self._models


    def atoms(self):
        return self._atoms


    def anisou(self):
        return self._anisou


    def termini(self):
        return self._termini


    def heteroatoms(self):
        return self._heteroatoms


    def connections(self):
        return self._connections


    def master(self, master=None):
        if master:
            if not isinstance(master, dict):
                raise TypeError(
                 "master must be dict, not '%s'" % str(master)
                )
            self._master = master
        else:
            return self._master


    generate_pdb_file = pdb_file_from_pdb_data_file


'''

def process_revdat_records(data_file, pdb_file):
    if pdb_file:
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
            return
    data_file._revisions = []


def process_sprsde_records(data_file, pdb_file):
    if pdb_file:
        sprsde = pdb_file.get_record_by_name("SPRSDE")
        if sprsde:
            data_file._supercedes = sprsde[31:75].split()
            data_file._supercede_date = date_from_string(sprsde[11:20])
            return
    data_file._supercedes = []
    data_file._supercede_date = None


def process_jrnl_records(data_file, pdb_file):
    if pdb_file:
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


def process_remark_records(data_file, pdb_file):
    if pdb_file:
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
            return
    data_file._remarks = []


def process_dbref_records(data_file, pdb_file):
    if pdb_file:
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
            return
    data_file._dbreferences = []


def process_seqadv_records(data_file, pdb_file):
    if pdb_file:
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
            return
    data_file._sequence_differences = []


def process_seqres_records(data_file, pdb_file):
    if pdb_file:
        seqres = pdb_file.get_records_by_name("SEQRES")
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


def process_modres_records(data_file, pdb_file):
    if pdb_file:
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
            return
    data_file._modified_residues = []


def process_het_records(data_file, pdb_file):
    if pdb_file:
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
            return
    data_file._hets = []


def process_hetnam_records(data_file, pdb_file):
    if pdb_file:
        hetnams = pdb_file.get_records_by_name("HETNAM")
        if hetnams:
            ids = list(set([r[11:14] for r in hetnams]))
            data_file._het_names = {
             het_id: merge_records(
              [r for r in hetnams if r[11:14] == het_id], 15, dont_condense=":;"
             ) for het_id in ids
            }
            return
    data_file._het_names = {}


def process_hetsyn_records(data_file, pdb_file):
    if pdb_file:
        hetsyns = pdb_file.get_records_by_name("HETSYN")
        if hetsyns:
            ids = list(set([r[11:14] for r in hetsyns]))
            data_file._het_synonyms = {
             het_id: merge_records(
              [r for r in hetsyns if r[11:14] == het_id], 15
             ).split(";") for het_id in ids
            }
            return
    data_file._het_synonyms = {}


def process_formul_records(data_file, pdb_file):
    if pdb_file:
        formuls = pdb_file.get_records_by_name("FORMUL")
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


def process_helix_records(data_file, pdb_file):
    if pdb_file:
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
            return
    data_file._helices = []


def process_sheet_records(data_file, pdb_file):
    if pdb_file:
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
                return
    data_file._sheets = []


def process_ssbond_records(data_file, pdb_file):
    if pdb_file:
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
            return
    data_file._ss_bonds= []


def process_link_records(data_file, pdb_file):
    if pdb_file:
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
            return
    data_file._links = []


def process_cispep_records(data_file, pdb_file):
    if pdb_file:
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
            return
    data_file._cis_peptides = []


def process_site_records(data_file, pdb_file):
    if pdb_file:
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
            return
    data_file._sites = []


def process_cryst1_records(data_file, pdb_file):
    if pdb_file:
        crystal = pdb_file.get_record_by_name("CRYST1")
        if crystal:
            data_file._crystal_a = crystal[6:15]
            data_file._crystal_b = crystal[15:24]
            data_file._crystal_c = crystal[24:33]
            data_file._crystal_alpha = crystal[33:40]
            data_file._crystal_beta = crystal[40:47]
            data_file._crystal_gamma = crystal[47:54]
            data_file._crystal_s_group = crystal[55:66]
            data_file._crystal_z = crystal[66:70]
            return
    data_file._crystal_a = None
    data_file._crystal_b = None
    data_file._crystal_c = None
    data_file._crystal_alpha = None
    data_file._crystal_beta = None
    data_file._crystal_gamma = None
    data_file._crystal_s_group = None
    data_file._crystal_z = None


def process_origx_records(data_file, pdb_file):
    if pdb_file:
        origx1 = pdb_file.get_record_by_name("ORIGX1")
        origx2 = pdb_file.get_record_by_name("ORIGX2")
        origx3 = pdb_file.get_record_by_name("ORIGX3")
        if origx1 and origx2 and origx3:
            data_file._origx_o11 = origx1[10:20]
            data_file._origx_o12 = origx1[20:30]
            data_file._origx_o13 = origx1[30:40]
            data_file._origx_t1 = origx1[45:55]
            data_file._origx_o21 = origx2[10:20]
            data_file._origx_o22 = origx2[20:30]
            data_file._origx_o23 = origx2[30:40]
            data_file._origx_t2 = origx2[45:55]
            data_file._origx_o31 = origx3[10:20]
            data_file._origx_o32 = origx3[20:30]
            data_file._origx_o33 = origx3[30:40]
            data_file._origx_t3 = origx3[45:55]
            return
    data_file._origx_o11 = None
    data_file._origx_o12 = None
    data_file._origx_o13 = None
    data_file._origx_t1 = None
    data_file._origx_o21 = None
    data_file._origx_o22 = None
    data_file._origx_o23 = None
    data_file._origx_t2 = None
    data_file._origx_o31 = None
    data_file._origx_o32 = None
    data_file._origx_o33 = None
    data_file._origx_t3 = None


def process_scale_records(data_file, pdb_file):
    if pdb_file:
        scale1 = pdb_file.get_record_by_name("SCALE1")
        scale2 = pdb_file.get_record_by_name("SCALE2")
        scale3 = pdb_file.get_record_by_name("SCALE3")
        if scale1 and scale2 and scale3:
            data_file._scale_s11 = scale1[10:20]
            data_file._scale_s12 = scale1[20:30]
            data_file._scale_s13 = scale1[30:40]
            data_file._scale_u1 = scale1[45:55]
            data_file._scale_s21 = scale2[10:20]
            data_file._scale_s22 = scale2[20:30]
            data_file._scale_s23 = scale2[30:40]
            data_file._scale_u2 = scale2[45:55]
            data_file._scale_s31 = scale3[10:20]
            data_file._scale_s32 = scale3[20:30]
            data_file._scale_s33 = scale3[30:40]
            data_file._scale_u3 = scale3[45:55]
            return
    data_file._scale_s11 = None
    data_file._scale_s12 = None
    data_file._scale_s13 = None
    data_file._scale_u1 = None
    data_file._scale_s21 = None
    data_file._scale_s22 = None
    data_file._scale_s23 = None
    data_file._scale_u2 = None
    data_file._scale_s31 = None
    data_file._scale_s32 = None
    data_file._scale_s33 = None
    data_file._scale_u3 = None


def process_matrix_records(data_file, pdb_file):
    if pdb_file:
        matrx1 = pdb_file.get_record_by_name("MTRIX1")
        matrx2 = pdb_file.get_record_by_name("MTRIX2")
        matrx3 = pdb_file.get_record_by_name("MTRIX3")
        if matrx1 and matrx2 and matrx3:
            data_file._matrix_serial_1 = matrx1[7:10]
            data_file._matrix_m11 = matrx1[10:20]
            data_file._matrix_m12 = matrx1[20:30]
            data_file._matrix_m13 = matrx1[30:40]
            data_file._matrix_v1 = matrx1[45:55]
            data_file._matrix_i_given_1 = matrx1[59] == 1
            data_file._matrix_serial_2 = matrx2[7:10]
            data_file._matrix_m21 = matrx2[10:20]
            data_file._matrix_m22 = matrx2[20:30]
            data_file._matrix_m23 = matrx2[30:40]
            data_file._matrix_v2 = matrx2[45:55]
            data_file._matrix_i_given_2 = matrx2[59] == 1
            data_file._matrix_serial_3 = matrx3[7:10]
            data_file._matrix_m31 = matrx3[10:20]
            data_file._matrix_m32 = matrx3[20:30]
            data_file._matrix_m33 = matrx3[30:40]
            data_file._matrix_v3 = matrx3[45:55]
            data_file._matrix_i_given_3 = matrx3[59] == 1
            return
    data_file._matrix_serial_1 = None
    data_file._matrix_m11 = None
    data_file._matrix_m12 = None
    data_file._matrix_m13 = None
    data_file._matrix_v1 = None
    data_file._matrix_i_given_1 = None
    data_file._matrix_serial_2 = None
    data_file._matrix_m21 = None
    data_file._matrix_m22 = None
    data_file._matrix_m23 = None
    data_file._matrix_v2 = None
    data_file._matrix_i_given_2 = None
    data_file._matrix_serial_3 = None
    data_file._matrix_m31 = None
    data_file._matrix_m32 = None
    data_file._matrix_m33 = None
    data_file._matrix_v3 = None
    data_file._matrix_i_given_3 = None


def process_model_records(data_file, pdb_file):
    if pdb_file:
        model_records = pdb_file.get_records_by_name("MODEL")
        if model_records:
            endmdls = pdb_file.get_records_by_name("ENDMDL")
            pairs = list(zip(model_records, endmdls))
            data_file._models = [{
             "model_id": pair[0][10:14],
             "start_record": pair[0].number(),
             "end_record": pair[1].number(),
            } for pair in pairs]
            return
    data_file._models = [{
     "model_id": 1,
     "start_record": -1,
     "end_record": -1
    }]


def process_atom_records(data_file, pdb_file):
    if pdb_file:
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
            return
    data_file._atoms = []


def process_anisou_records(data_file, pdb_file):
    if pdb_file:
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
            return
    data_file._anisou = []


def process_ter_records(data_file, pdb_file):
    if pdb_file:
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
            return
    data_file._termini = []


def process_hetatm_records(data_file, pdb_file):
    if pdb_file:
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
            return
    data_file._heteroatoms = []


def process_conect_records(data_file, pdb_file):
    if pdb_file:
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
            return
    data_file._connections = []


def process_master_records(data_file, pdb_file):
    if pdb_file:
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
        return
    data_file._master = None'''
