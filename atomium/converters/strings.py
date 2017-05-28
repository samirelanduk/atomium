"""This module contains tools for processing strings into more manageable
forms, and back into strings."""

def string2lines(s):
    """Takes a string and breaks it into lines on line breaks. It will first
    convert any windows line breaks to unix line breaks.

    :param str s: The string to break up.
    :rtype: ``str``."""

    return s.replace("\r\n", "\n").split("\n")


def string_from_file(path):
    """Opens a file from the given path and returns the contents as a string.

    :param str path: The path to the file.
    :rtype: ``str``"""

    with open(path) as f:
        return f.read()



def string_to_file(string, path):
    """Saves a string to a given path as a file.

    :param str string: The string to save.
    :param str path: The file to save it in."""
    
    with open(path, "w") as f:
        f.write(string)
