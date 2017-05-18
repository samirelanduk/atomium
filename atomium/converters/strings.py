"""This module contains tools for processing strings into more manageable
forms, and back into strings."""

def string2lines(s):
    """Takes a string and breaks it into lines on line breaks. It will first
    convert any windows line breaks to unix line breaks.

    :param str s: The string to break up.
    :rtype: ``str``."""
    
    return s.replace("\r\n", "\n").split("\n")
