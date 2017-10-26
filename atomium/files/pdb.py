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


    def __repr__(self):
        num = len(self._models)
        return "<Pdb {}({} model{})>".format(
         self._code + " " if self._code else "", num, "" if num == 1 else "s"
        )


    def models(self):
        """Returns the :py:class:`.Model` objects that the Pdb contains.

        :rtype: ``tuple``"""

        return tuple(self._models)


    def model(self, model=None):
        """Returns the first :py:class:`.Model` that the Pdb file contains."""

        return self._models[0] if self._models else None


    def code(self, code=None):
        """Returns or sets the Pdb's 4-letter code. If a value is given, the
        code will be updated.

        :param str code: if given, the code will be updated.
        :raises TypeError: if a non-str code is given.
        :raises ValueError: if the PDB code given is invalid.
        :rtype: ``str``"""

        if code:
            if not isinstance(code, str):
                raise TypeError("PDB code {} is not str".format(code))
            if len(code) != 4:
                raise ValueError("PDB code {} is not valid".format(code))
            self._code = code
        else:
            return self._code


    def deposition_date(self, deposition_date=None):
        """Returns or sets the Pdb's desposition date. If a value is given, the
        date will be updated.

        :param date deposition_date: if given, the date will be updated.
        :raises TypeError: if a non-date is given.
        :rtype: ``datetime.date``"""

        if deposition_date:
            if not isinstance(deposition_date, datetime.date):
                raise TypeError("{} is not a Date".format(deposition_date))
            self._deposition_date = deposition_date
        else:
            return self._deposition_date


    def title(self, title=None):
        """Returns or sets the Pdb's title. If a value is given, the title will
        be updated.

        :param str title: if given, the title will be updated.
        :raises TypeError: if a non-str title is given.
        :rtype: ``str``"""

        if title:
            if not isinstance(title, str):
                raise TypeError("PDB title {} is not str".format(title))
            self._title = title
        else:
            return self._title


    def resolution(self, resolution=None):
        """Returns or sets the Pdb's resolution. If a value is given, the
        resolution will be updated.

        :param str resolution: if given, the resolution will be updated.
        :raises TypeError: if a non-numeric resolution is given.
        :rtype: ``float``"""

        if resolution:
            if not isinstance(resolution, (int, float)):
                raise TypeError("Resolution {} isn't numeric".format(resolution))
            self._resolution = resolution
        else:
            return self._resolution


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
