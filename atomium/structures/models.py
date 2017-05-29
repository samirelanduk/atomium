"""Contains the Model class."""

from .molecules import AtomicStructure

class Model(AtomicStructure):
    """Base class: :py:class:`.AtomicStructure`

    Represents molecular systems.

    :param \*atoms: The atoms that make up the model."""

    def __init__(self, *atoms):
        AtomicStructure.__init__(self, *atoms)


    def to_file_string(self, file_format, description=""):
        """Converts a model to a filestring. Currently supported file formats
        are: .xyz.

        :param str file_format: The file format to use, in lowercase.
        :param str description: A model description to put in the file.
        :raises ValueError: if an unsopported file format is given."""
        
        if file_format == "xyz":
            from ..converters.model2xyzstring import model_to_xyz_string
            return model_to_xyz_string(self, description)
        else:
            raise ValueError("{} is not a valid file type".format(file_format))
