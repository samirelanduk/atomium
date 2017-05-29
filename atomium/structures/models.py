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


    def save(self, path, *args, **kwargs):
        """Saves the Model to file, in the format implied by the extension of
        the path you provide (i.e. giving a path ``/path/to/file.xyz`` will
        save as .xyz).

        :param str path: The path to save to. The extension you provide here is\
        important as atomium will use that to determine what file format to\
        save as.
        :param str description: A model description to put in the file."""

        file_format = path.split(".")[-1].lower()
        s = self.to_file_string(file_format, *args, **kwargs)
        from ..converters.strings import string_to_file
        string_to_file(s, path)
