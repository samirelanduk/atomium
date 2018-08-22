"""Contains the File class."""

from math import inf
from copy import deepcopy
from ..models import Model

class File:
    """When a file is parsed, the result is a ``File``. It contains the
    structure of interest, as well as meta information."""

    def __init__(self):
        self._filetype = None
        self._code = None
        self._title = None
        self._deposition_date = None
        self._classification = None
        self._keywords = []
        self._authors = []
        self._technique = None
        self._source_organism = None
        self._expression_system = None
        self._resolution = None
        self._rvalue = None
        self._rfree = None
        self._assemblies = []
        self._models = []


    def __repr__(self):
        return "<{}{}File>".format(
         self._code if self._code else "",
         "." + self._filetype + " " if self._filetype else "",
        )


    @property
    def filetype(self):
        """The filetype that this File was created from, such as .pdb or
        .cif.

        :rtype: ``str``"""

        return self._filetype


    @filetype.setter
    def filetype(self, filetype):
        self._filetype = filetype


    @property
    def code(self):
        """The unique database identifer for this structure.

        :rtype: ``str``"""

        return self._code


    @code.setter
    def code(self, code):
        self._code = code


    @property
    def title(self):
        """The structure's text description.

        :rtype: ``str``"""

        return self._title


    @title.setter
    def title(self, title):
        self._title = title


    @property
    def deposition_date(self):
        """The date the structure was submitted for publication.

        :rtype: ``datetime.date``"""

        return self._deposition_date


    @deposition_date.setter
    def deposition_date(self, deposition_date):
        self._deposition_date = deposition_date


    @property
    def classification(self):
        """The structure's formal classification.

        :rtype: ``str``"""

        return self._classification


    @classification.setter
    def classification(self, classification):
        self._classification = classification


    @property
    def keywords(self):
        """The structure's keyword descriptors.

        :rtype: ``list``"""

        return self._keywords


    @property
    def authors(self):
        """The structure's authors.

        :rtype: ``list``"""

        return self._authors


    @property
    def technique(self):
        """The structure's experimental technique.

        :rtype: ``str``"""

        return self._technique


    @technique.setter
    def technique(self, technique):
        self._technique = technique


    @property
    def source_organism(self):
        """The structure's original organism.

        :rtype: ``float``"""

        return self._source_organism


    @source_organism.setter
    def source_organism(self, source_organism):
        self._source_organism = source_organism


    @property
    def expression_system(self):
        """The organism the structure was expressed in.

        :rtype: ``float``"""

        return self._expression_system


    @expression_system.setter
    def expression_system(self, expression_system):
        self._expression_system = expression_system


    @property
    def resolution(self):
        """The structure's resolution.

        :rtype: ``float``"""

        return self._resolution


    @resolution.setter
    def resolution(self, resolution):
        self._resolution = resolution


    @property
    def rvalue(self):
        """The structure's R-value.

        :rtype: ``float``"""

        return self._rvalue


    @rvalue.setter
    def rvalue(self, rvalue):
        self._rvalue = rvalue


    @property
    def rfree(self):
        """The structure's R-free value.

        :rtype: ``float``"""

        return self._rfree


    @rfree.setter
    def rfree(self, rfree):
        self._rfree = rfree


    @property
    def assemblies(self):
        """The structure's assembly instructions.

        :rtype: ``list``"""

        return self._assemblies


    @property
    def models(self):
        """The structure's models.

        :rtype: ``list``"""

        return self._models


    @property
    def model(self):
        """The structure's primary :py:class:`.Model`.

        :rtype: ``Model``"""

        if self._models: return self._models[0]


    def generate_assembly(self, id):
        """Creates a :py:class:`.Model` from the current model and the
        instructions contained in one of the Pdb's assemblies, which you
        specify.

        :param int id: The assembly to use to generate the assembly.
        :raises ValueError: if you give an ID which doesn't correspond to a\
        assembly.
        :rtype: ``Model``"""

        model = self._models[0]
        for assembly in self._assemblies:
            if assembly["id"] == id:
                break
        else:
            raise ValueError("No assembly with ID {}".format(id))
        new_chains = []
        for transformation in assembly["transformations"]:
            chains = [model.chain(id_) for id_ in transformation["chains"]]
            chains = [c for c in chains if c is not None]
            for chain in chains:
                new_chain = chain.copy()
                new_chain.transform(transformation["matrix"])
                new_chain.translate(transformation["vector"])
                new_chains.append(new_chain)
        return Model(*new_chains)


    @property
    def best_assembly(self):
        """Returns the 'best' biological assembly for this Pdb - the one with
        the lowest (most negative) delta energy.

        If there are no assemblies, ``None`` is returned.

        :rtype: ``Model``"""

        sorted_mol = sorted(
         self._assemblies,
         key=lambda b: inf if b["delta_energy"] is None else b["delta_energy"]
        )
        if sorted_mol:
            return sorted_mol[0]


    def generate_best_assembly(self):
        """Generates the 'best' biological assembly for this Pdb - the one with
        the lowest (most negative) delta energy.

        If there are no assemblies, the original model will be returned.

        :rtype: ``Model``"""

        best = self.best_assembly
        if best:
            return self.generate_assembly(best["id"])
        else:
            return self._models[0]


    def save(self, path):
        """Saves the File to the location specified. The file extension in the
        path given will be used to determine what file type to save as.

        Currently supported extensions are .pdb and .xyz (.cif coming soon).

        :param str path: The path to save to."""
        
        filestring = ""
        if path.endswith(".xyz"):
            from .xyz import file_to_xyz_string
            filestring = file_to_xyz_string(self)
        elif path.endswith(".pdb"):
            from .pdb import file_to_pdb_string
            filestring = file_to_pdb_string(self)
        else:
            raise ValueError("{} has an unknown file extension".format(path))
        with open(path, "w") as f:
            f.write(filestring)
