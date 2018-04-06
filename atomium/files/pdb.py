"""Contains the Pdb class and functions for opening them."""

import datetime
from ..structures.models import Model

class Pdb:
    """A Pdb is used to represent a fully processed PDB file."""

    def __init__(self):
        self._models = []
        self._code, self._deposition_date = None, None
        self._title = None
        self._resolution = None
        self._organism = None
        self._expression_system = None
        self._technique = None
        self._classification = None
        self._rfactor = None


    def __repr__(self):
        num = len(self._models)
        return "<Pdb {}({} model{})>".format(
         self._code + " " if self._code else "", num, "" if num == 1 else "s"
        )


    @property
    def models(self):
        """Returns the :py:class:`.Model` objects that the Pdb contains.

        :rtype: ``tuple``"""

        return tuple(self._models)


    @property
    def model(self):
        """Returns the first :py:class:`.Model` that the Pdb file contains."""

        return self._models[0] if self._models else None


    @property
    def code(self):
        """The Pdb's 4-letter code.

        :raises TypeError: if a non-str code is set.
        :raises ValueError: if the PDB code given is invalid.
        :rtype: ``str``"""

        return self._code


    @code.setter
    def code(self, code):
        if not isinstance(code, str):
            raise TypeError("PDB code {} is not str".format(code))
        if len(code) != 4:
            raise ValueError("PDB code {} is not valid".format(code))
        self._code = code


    @property
    def deposition_date(self):
        """The Pdb's desposition date.

        :raises TypeError: if a non-date is given.
        :rtype: ``datetime.date``"""

        return self._deposition_date


    @deposition_date.setter
    def deposition_date(self, deposition_date):
        if not isinstance(deposition_date, datetime.date):
            raise TypeError("{} is not a Date".format(deposition_date))
        self._deposition_date = deposition_date


    @property
    def title(self):
        """The Pdb's title.

        :raises TypeError: if a non-str title is given.
        :rtype: ``str``"""

        return self._title


    @title.setter
    def title(self, title):
        if not isinstance(title, str):
            raise TypeError("PDB title {} is not str".format(title))
        self._title = title



    @property
    def resolution(self):
        """The Pdb's resolution.

        :raises TypeError: if a non-numeric resolution is given.
        :rtype: ``float``"""

        return self._resolution


    @resolution.setter
    def resolution(self, resolution):
        if not isinstance(resolution, (int, float)):
            raise TypeError("Resolution {} isn't numeric".format(resolution))
        self._resolution = resolution



    @property
    def rfactor(self):
        """The Pdb's R-factor.

        :raises TypeError: if a non-numeric rfactor is given.
        :rtype: ``float``"""

        return self._rfactor


    @rfactor.setter
    def rfactor(self, rfactor):
        if not isinstance(rfactor, (int, float)):
            raise TypeError("R-factor {} isn't numeric".format(rfactor))
        self._rfactor = rfactor


    @property
    def organism(self):
        """The Pdb's source organism.

        :raises TypeError: if a non-str organism is given.
        :rtype: ``str``"""

        return self._organism


    @organism.setter
    def organism(self, organism):
        if not isinstance(organism, str):
            raise TypeError("PDB organism {} is not str".format(organism))
        self._organism = organism


    @property
    def expression_system(self):
        """The Pdb's expression organism.

        :raises TypeError: if a non-str organism is given.
        :rtype: ``str``"""

        return self._expression_system


    @expression_system.setter
    def expression_system(self, expression_system):
        if not isinstance(expression_system, str):
            raise TypeError(
             "PDB expression organism {} is not str".format(expression_system)
            )
        self._expression_system = expression_system


    @property
    def technique(self):
        """The Pdb's experimental technique.

        :raises TypeError: if a non-str technique is given.
        :rtype: ``str``"""

        return self._technique


    @technique.setter
    def technique(self, technique):
        if not isinstance(technique, str):
            raise TypeError("Technique {} is not str".format(technique))
        self._technique = technique


    @property
    def classification(self):
        """The Pdb's classification.

        :raises TypeError: if a non-str classification is given.
        :rtype: ``str``"""

        return self._classification


    @classification.setter
    def classification(self, classification):
        if not isinstance(classification, str):
            raise TypeError(
             "Classification {} is not str".format(classification)
            )
        self._classification = classification


    def to_file_string(self):
        """Returns the file text that represents this Pdb.

        :rtype: ``str``"""

        from ..files.pdb2pdbdict import pdb_to_pdb_dict
        from ..files.pdbdict2pdbstring import pdb_dict_to_pdb_string
        pdb_dict = pdb_to_pdb_dict(self)
        return pdb_dict_to_pdb_string(pdb_dict)


    def save(self, path):
        """Saves the Pdb as a .pdb file.

        :param str path: The path to save to."""

        from ..files.utilities import string_to_file
        string_to_file(self.to_file_string(), path)
