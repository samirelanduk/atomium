"""Contains the Xyz class."""

from ..models.molecules import Model

class Xyz:
    """Represents .xyz files and the model they contain.

    :param str title: The title at the head of the file."""

    def __init__(self, title=""):
        self._title = title
        self._model = None


    def __repr__(self):
        return "<Xyz{}>".format(" (" + self._title + ")" if self._title else "")


    @property
    def title(self):
        """Returns the title associated with the .xyz file.

        :rtype: ``str``"""

        return self._title


    @title.setter
    def title(self, title):
        self._title = title


    @property
    def model(self):
        """Returns the :py:class:`.Model` that the .xyz file contains.

        :rtype: ``Model``"""

        return self._model


    def to_file_string(self):
        """Returns the file text that represents this Xyz.

        :rtype: ``str``"""

        from ..files.xyz2xyzdict import xyz_to_xyz_dict
        from ..files.xyzdict2xyzstring import xyz_dict_to_xyz_string
        xyz_dict = xyz_to_xyz_dict(self)
        return xyz_dict_to_xyz_string(xyz_dict)


    def save(self, path):
        """Saves the Xyz as a .xyz file.

        :param str path: The path to save to."""

        from ..files.utilities import string_to_file
        string_to_file(self.to_file_string(), path)
