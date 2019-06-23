"""Decorators and metaclasses used by atomium structures."""

import re

def get_object_from_filter(obj, components):
    """Gets the object whose attributes are actually being queried, which may be
    a different object if there is a chain.

    :param obj: the intial object.
    :param list components: the components of the original key.
    :returns: the relevant object"""

    components = components[:]
    while len(components) > 2:
        obj = getattr(obj, components.pop(0))
    if len(components) == 2:
        if components[-1] != "regex":
            if not hasattr(obj, f"__{components[-1]}__"):
                obj = getattr(obj, components[0])
    return obj


def get_object_attribute_from_filter(obj, components):
    """Gets the object's value of some attribute based on a list of key
    components.

    :param obj: the object with attributes.
    :param list components: the components of the original key.
    :returns: the value"""

    try:
        return getattr(
         obj, components[-1] if hasattr(obj, components[-1]) else components[-2]
        )
    except: return None


def attribute_matches_value(attribute, value, components):
    """Checks if an attribute value matches a given value. The components given
    will determine whether an exact match is sought, or whether a more complex
    criterion is used.
    
    :param attribute: the value of an object's attribute.
    :param value: the value to match against.
    :param list components: the components of the original key.
    :rtype: ``bool``"""

    if components[-1] == "regex":
        return re.match(value, attribute)
    possible_magic = f"__{components[-1]}__"
    if hasattr(attribute, possible_magic):
        return getattr(attribute, possible_magic)(value)
    return getattr(attribute, "__eq__")(value)


def filter_objects(objects, key, value):
    """Takes a :py:class:`.StructureSet` of objects, and filters them on object
    properties.

    They key can be an attribute of the object, or a complex double-underscore
    separated chain of attributes.

    :param StructreSet objects: the dictionary of objects - the keys are\
    unimportant.
    :param str key: the attribute to search. This can be an attribute of the\
    object, or attr__regex, or attr__gt etc.
    :param value: the value that the attribute must have.
    :rtype: ``dict``"""

    components = key.split("__")
    matching_objects = []
    for structure in objects.structures:
        obj = get_object_from_filter(structure, components)
        attr = get_object_attribute_from_filter(obj, components)
        if attribute_matches_value(attr, value, components):
            matching_objects.append(structure)
    return StructureSet(*matching_objects)


def query(func, tuple_=False):
    """A decorator which can be applied to any function which returns a
    :py:class:`.StructureSet` and which takes no other parameters other than
    ``self``. It will query the returned objects by any keyword argument, or
    use a positional argument to search by ID.

    :param func: the function to modify.
    :param bool tuple_: if ``True``, objects will be returned in a tuple not a\
    set.
    :rtype: ``function``"""

    def structures(self, *args, **kwargs):
        objects = func(self)
        original = {s: n for n, s in enumerate(objects.structures)}
        if len(args) == 1:
            return {objects.get(args[0])} if args[0] in objects.ids else set()
        for k, v in kwargs.items():
            objects = filter_objects(objects, k, v)
        if tuple_:
            return tuple(sorted(
             objects.structures, key=lambda s: original[s]
            ))
        else:
            return set(objects.structures)
    return structures


def getone(func):
    """A decorator which can be applied to any function which returns an
    iterable. It produces a function which just gets the first item in that
    iterable.

    In atomium, various classes define methods like atoms, residues, etc. - this
    decorator can make a function like atom, residue which takes all the same
    params but just returns one object.

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
    """A data structure for holding structures. It stores them internally
    as a dictionary where they keys are IDs (to allow rapid lookup by ID) and
    the values are all structures with that ID (to allow for duplicate IDs).

    Two structure sets can be added together, but they are immutable - the
    structures they have when they are made is the structures they will always
    have.

    They're basically sets optimised to lookup things by ID.

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
        """Gets a structure by ID. If an ID points to multiple structures, just
        one will be returned.

        :returns: some structure."""

        matches = self._d.get(id, set())
        for match in matches: return match
