"""Decorators and metaclasses used by atomium structures."""

import re

def filter_objects(objects, key, value):
    """Takes a :py:class:`.StructureSet` of objects, and filters them on object
    properties.

    :param StructreSet objects: the dictionary of objects - the keys are\
    unimportant.
    :param str key: the attribute to search. This can be an attribute of the\
    object, or attr__regex, or attr__gt etc.
    :param value: the value that the attribute must have.
    :rtype: ``dict``"""

    if "__" in key:
        attr, comp = key.split("__")
        if comp == "regex":
            objects = StructureSet(*[
             s for s in objects.structures if re.match(value, getattr(s, attr))
            ])
        else:
            comp = "__{}__".format(comp)
            objects = StructureSet(*[s for s in objects.structures
             if getattr(getattr(s, attr), comp)(value)])
    else:
        objects = StructureSet(*[
         s for s in objects.structures if getattr(s, key).__eq__(value)
        ])
    return objects


def query(func, tuple_=False):
    """A decorator which can be applied to any function which returns a
    :py:class:`.StructureSet` and which takes no other paramters other than
    ``self``. It will query the returned objects by any keyword argument, or
    use a positional argument to search by ID.

    :param func: the function to modify.
    :param bool tuple_: if ``True``, objects will be returned in a tuple not a\
    set.
    :rtype: ``function``"""

    def structures(self, *args, **kwargs):
        objects = func(self)
        original = list(objects.structures)
        if len(args) == 1:
            return {objects.get(args[0])} if args[0] in objects.ids else set()
        for k, v in kwargs.items():
            objects = filter_objects(objects, k, v)
        if tuple_:
            return tuple(sorted(
             objects.structures, key=lambda s: original.index(s)
            ))
        else:
            return set(objects.structures)
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



class StructureSet:
    """A data structure for holding sub-structures. It stores them internally
    as a dictionary where they keys are IDs (to allow rapid lookup by ID) and
    the values are all structures with that ID (to allow for duplicate IDs).

    :param \* args: the structures that will make up the StructureSet."""

    def __init__(self, *args):
        self._d = {}
        for obj in args:
            if obj._id in self._d:
                self._d[obj._id].add(obj)
            else:
                self._d[obj._id] = {obj}


    def __add__(self, other):
        new = StructureSet()
        for s in (self, other):
            for key, value in s._d.items():
                if key in new._d:
                    new._d[key].update(value)
                else:
                    new._d[key] = value
        return new


    def __len__(self):
        return len(self.structures)


    def add(self, obj):
        """Adds a structure to the StructureSet.

        :param obj: the structure to add."""

        if obj._id in self._d:
            self._d[obj._id].add(obj)
        else:
            self._d[obj._id] = {obj}


    def remove(self, obj):
        """Removes a structure from the StructureSet.

        :param obj: the structure to remove."""

        self._d[obj._id].remove(obj)
        if not self._d[obj._id]:
            del self._d[obj._id]


    @property
    def ids(self):
        """Returns the IDs of the StructureSet.

        :rtype: ``set``"""

        return self._d.keys()


    @property
    def structures(self):
        """Returns the structures of the StructureSet.

        :rtype: ``list``"""

        structures = []
        for s in self._d.values(): structures += s
        return structures


    def get(self, id):
        """Gets a structure by ID.

        :returns: some structure."""

        matches = self._d.get(id, set())
        for match in matches: return match
