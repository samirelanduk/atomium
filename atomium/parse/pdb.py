"""Contains the Pdb class."""

from ..structures.models import Model

class Pdb:
    """A PdbDataFile is used to represent a fully processed PDB file."""

    def __init__(self):
        self._model = None


    def __repr__(self):
        return "<Pdb>"


    def model(self, model=None):
        """Returns the :py:class:`.Model` that the .pdb file contains. If a
        model is given, the model will be changed to the new model.

        :param Model model: If given, the model will be updated to this.
        :raises TypeError: if the model given is not a :py:class:`.Model`."""

        if model is None:
            return self._model
        else:
            if not isinstance(model, Model):
                raise TypeError("model must be Model, not '{}'".format(model))
            self._model = model
