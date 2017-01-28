"""This module contains creates the final Pdb object itself, and processes the
data contained in the data file."""

class Pdb:
    """A representation of a PDB file and its contents, including the structure.

    :param PdbDataFile data_file: The PDB data file with the parsed values."""

    def __init__(self, data_file):
        from ..converters.pdbdatafile2model import model_from_pdb_data_file
        self._data_file = data_file
        self._models = []
        for model_dict in self._data_file.models():
            model = model_from_pdb_data_file(data_file, model_dict["model_id"])
            self._models.append(model)


    def __repr__(self):
        return "<Pdb (%s)>" % (self.pdb_code() if self.pdb_code() else "????")


    def data_file(self):
        """The :py:class:`.PdbDataFile` from which the object was created.

        :rtype: ``PdbDataFile``"""

        return self._data_file


    def classification(self):
        """The PDB classification.

        :rtype: ``str``"""

        return self._data_file.classification()


    def deposition_date(self):
        """The date the PDB was deposited.

        :rtype: ``datetime.Date``"""

        return self._data_file.deposition_date()


    def pdb_code(self):
        """The PDB four-letter code.

        :rtype: ``str``"""

        return self._data_file.pdb_code()


    def is_obsolete(self):
        """``True`` if the PDB has been made obsolete by a newer PDB.

        :rtype: ``bool``"""

        return self._data_file.is_obsolete()


    def obsolete_date(self):
        """The date the PDB was made obsolete.

        :rtype: ``datetime.Date``"""

        return self._data_file.obsolete_date()


    def replacement_code(self):
        """The PDB code of the replacing PDB.

        :rtype: ``str``"""

        return self._data_file.replacement_code()


    def title(self):
        """The title of the PDB.

        :rtype: ``str``"""

        return self._data_file.title()


    def split_codes(self):
        """The PDB codes which complete this structure.

        :rtype: ``list``"""

        return self._data_file.split_codes()


    def caveat(self):
        """Any caveats for this structure.

        :rtype: ``str``"""

        return self._data_file.caveat()


    def keywords(self):
        """Keywords for this PDB.

        :rtype: ``list``"""

        return self._data_file.keywords()


    def experimental_techniques(self):
        """The experimental techniques used to produce this PDB.

        :rtype: ``list``"""

        return self._data_file.experimental_techniques()


    def model_count(self):
        """The number of models in this PDB.

        :rtype: ``int``"""

        return self._data_file.model_count()


    def model_annotations(self):
        """Annotations for the PDB's models.

        :rtype: ``list``"""

        return self._data_file.model_annotations()


    def authors(self):
        """The PDB's authors.

        :rtype: ``list``"""

        return self._data_file.authors()


    def revisions(self):
        """Any changes made to the PDB file.

        :rtype: ``list``"""

        return self._data_file.revisions()


    def supercedes(self):
        """The PDB codes that this PDB replaces.

        :rtype: ``list``"""

        return self._data_file.supercedes()


    def supercede_date(self):
        """The date this PDB replaced another.

        :rtype: ``datetime.Date``"""

        return self._data_file.supercede_date()


    def journal(self):
        """The publication information for this PDB.

        :rtype: ``dict``"""

        return self._data_file.journal()


    def models(self):
        """The PDB's models.

        :rtype: ``list``"""

        return list(self._models)


    def model(self):
        """The first :py:class:`.Model` in the PDB models.

        :rtype: ``Model``"""

        return self._models[0]
