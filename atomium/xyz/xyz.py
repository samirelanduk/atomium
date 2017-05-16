"""Contains the Xyz class."""

class Xyz:
    """Represents .xyz files and the model they contain.

    :pram str comment: The comment at the head of the file.
    :raises TypeError: if a comment is given which is not a string."""

    def __init__(self, comment=""):
        if not isinstance(comment, str):
            raise TypeError("comment must be str, not '{}'".format(comment))
        self._comment = comment
        self._model = None


    def __repr__(self):
        return "<Xyz ({})>".format(self._comment)
