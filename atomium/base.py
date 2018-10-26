"""Decorators and metaclasses used by atomium structures."""

import re

def filter_objects(objects, key, value):
    """Takes a dictionary of objects, and filters them on object properties.

    :param dict objects: the dictionary of objects - the keys are unimportant.
    :param str key: the attribute to search. This can be an attribute of the\
    object, or attr__regex, or attr__gt etc.
    :param value: the value that the attribute must have.
    :rtype: ``dict``"""

    if "__" in key:
        attr, comp = key.split("__")
        if comp == "regex":
            objects = {i: o for i, o in objects.items()
             if re.match(value, getattr(o, attr))}
        else:
            comp = "__{}__".format(comp)
            objects = {i: o for i, o in objects.items()
             if getattr(getattr(o, attr), comp)(value)}
    else:
        objects = {i: o for i, o in objects.items()
         if getattr(o, key).__eq__(value)}
    return objects


def query(func, tuple_=False):
    """A decorator which can be applied to any function which returns a ``dict``
    and which takes no other paramters other than ``self``. It will query the
    returned objects by any keyword argument, or use a positional argument to
    search by ID.

    :param func: the function to modify.
    :param bool tuple_: if ``True``, objects will be returned in a tuple not a\
    set.
    :rtype: ``function``"""

    def structures(self, *args, **kwargs):
        objects = func(self)
        if len(args) == 1:
            return {objects[args[0]]} if args[0] in objects else set()
        for k, v in kwargs.items():
            objects = filter_objects(objects, k, v)
        t = tuple if tuple_ else set
        return t(objects.values())
    return structures


def getone(func):
    """A decorator which can be applied to any function which returns an
    iterable. It produces a function which just gets the first item in that
    iterable.

    :param func: the function to modify.
    :rtype: ``function``"""

    def structure(self, *args, **kwargs):
        for obj in func(self, *args, **kwargs): return obj
    return structure



class StructureClass(type):
    """A metaclass which can be applied to structure class. It will override
    the instantation behaviour so that all methods that belong to a preset
    list ('atoms', 'chains' etc.) will have the :py:func:`.query` decorator
    applied and a copy with the :py:func:`.getone` decorator applied."""

    METHODS = ["chains", "residues", "ligands", "waters", "molecules", "atoms"]

    def __new__(self, *args, **kwargs):
        cls = type.__new__(self, *args, **kwargs)
        for attribute in dir(cls):
            if attribute in cls.METHODS:
                setattr(cls, attribute, query(
                 getattr(cls, attribute),
                 tuple_=(attribute == "residues" and cls.__name__ == "Chain")
                ))
                setattr(cls, attribute[:-1], getone(getattr(cls, attribute)))
        return cls
