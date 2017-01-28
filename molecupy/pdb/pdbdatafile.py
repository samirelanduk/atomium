"""This module performs the actual parsing of the PDB file, though it does not
process the values that it extracts."""

import datetime
import math
from .pdbfile import PdbRecord, PdbFile

class PdbDataFile:
    """This object is essentially a list of values extracted from a PDB file. It
    functions as a data sheet."""

    def __init__(self):
        self._source = None

        self._classification = None
        self._deposition_date = None
        self._pdb_code = None
        self._is_obsolete = False
        self._obsolete_date = None
        self._replacement_code = None
        self._title = None
        self._split_codes = []
        self._caveat = None
        self._compounds = []
        self._sources = []
        self._keywords = []
        self._experimental_techniques = []
        self._model_count = 1
        self._model_annotations = []
        self._authors = []
        self._revisions = []
        self._supercedes = []
        self._supercede_date = None
        self._journal = None
        self._remarks = []

        self._dbreferences = []
        self._sequence_differences = []
        self._residue_sequences = []
        self._modified_residues = []

        self._hets = []
        self._het_names = {}
        self._het_synonyms = {}
        self._formulae = {}

        self._helices = []
        self._sheets = []

        self._ss_bonds = []
        self._links = []
        self._cis_peptides = []

        self._sites = []

        self._crystal = None
        self._origix = None
        self._scale = None
        self._matrix = None

        self._models = [{"model_id": 1, "start_record": -1, "end_record": -1}]
        self._atoms = []
        self._anisou = []
        self._termini = []
        self._heteroatoms = []

        self._connections = []

        self._master = None


    def __repr__(self):
        return "<PdbDataFile (%s)>" % (self.pdb_code() if self.pdb_code() else "????")


    def source(self):
        """The object from which this PdbDataFile was created."""

        return self._source


    def to_pdb_file(self):
        """Converts the PdbDataFile to a :py:class:`.PdbFile`."""

        from ..converters.pdbdatafile2pdbfile import pdb_file_from_pdb_data_file
        return pdb_file_from_pdb_data_file(self)


    def classification(self, classification=None):
        """The classification of the PDB.

        :param str classification: if given, the classifcation will be set to\
        this.
        :rtype: ``str``"""

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


    def crystal(self, crystal=None):
        if crystal:
            if not isinstance(crystal, dict):
                raise TypeError(
                 "crystal must be dict, not '%s'" % str(crystal)
                )
            self._crystal = crystal
        else:
            return self._crystal


    def origx(self, origx=None):
        if origx:
            if not isinstance(origx, dict):
                raise TypeError(
                 "origx must be dict, not '%s'" % str(origx)
                )
            self._origx = origx
        else:
            return self._origx


    def scale(self, scale=None):
        if scale:
            if not isinstance(scale, dict):
                raise TypeError(
                 "scale must be dict, not '%s'" % str(scale)
                )
            self._scale = scale
        else:
            return self._scale


    def matrix(self, matrix=None):
        if matrix:
            if not isinstance(matrix, dict):
                raise TypeError(
                 "matrix must be dict, not '%s'" % str(matrix)
                )
            self._matrix = matrix
        else:
            return self._matrix


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
