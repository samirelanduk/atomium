"""Contains the Xyz class."""

from ..structures.models import Model

class Xyz:
    """Represents .xyz files and the model they contain.

    :pram str title: The title at the head of the file.
    :raises TypeError: if a title is given which is not a string."""

    def __init__(self, title=""):
        if not isinstance(title, str):
            raise TypeError("title must be str, not '{}'".format(title))
        self._title = title
        self._model = None


    def __repr__(self):
        return "<Xyz ({})>".format(self._title)


    def title(self, title=None):
        """Returns the title associated with the .xyz file. If a string is
        given, the title will be updated.

        :param str title: If given, the commnent will be updated to this.
        :raises TypeError: if the title given is not a string."""

        if title is None:
            return self._title
        else:
            if not isinstance(title, str):
                raise TypeError("title must be str, not '{}'".format(title))
            self._title = title


    def model(self, model=None):
        """Returns the :py:class:`.Model` that the .xyz file contains. If a
        model is given, the model will be changed to the new model.

        :param Model model: If given, the model will be updated to this.
        :raises TypeError: if the model given is not a :py:class:`.Model`."""

        if model is None:
            return self._model
        else:
            if not isinstance(model, Model):
                raise TypeError("model must be Model, not '{}'".format(model))
            self._model = model


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
