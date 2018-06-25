"""Contains the Pdb class and functions for opening them."""

import datetime
from math import inf
from ..models.molecules import Model

class Pdb:
    """A Pdb is used to represent a fully processed PDB file. It contains the
    :py:class:`.Model` and annotation information."""

    def __init__(self):
        self._models = []
        self._code, self._deposition_date = None, None
        self._title = None
        self._resolution = None
        self._organism = None
        self._expression_system = None
        self._technique = None
        self._classification = None
        self._rfactor = None
        self._rfree = None
        self._rcount = None
        self._keywords = []
        self._biomolecules = []


    def __repr__(self):
        num = len(self._models)
        return "<Pdb {}({} model{})>".format(
         self._code + " " if self._code else "", num, "" if num == 1 else "s"
        )


    @property
    def models(self):
        """Returns the :py:class:`.Model` objects that the Pdb contains.

        :rtype: ``tuple``"""

        return tuple(self._models)


    @property
    def model(self):
        """Returns the first :py:class:`.Model` that the Pdb file contains."""

        return self._models[0] if self._models else None


    @property
    def code(self):
        """The Pdb's 4-letter code - its unique identifier in the Protein Data
        Bank.

        :rtype: ``str``"""

        return self._code


    @code.setter
    def code(self, code):
        self._code = code


    @property
    def deposition_date(self):
        """The Pdb's desposition date - when it was submitted to the Protein
        Data Bank.

        :rtype: ``datetime.date``"""

        return self._deposition_date


    @deposition_date.setter
    def deposition_date(self, deposition_date):
        self._deposition_date = deposition_date


    @property
    def title(self):
        """The Pdb's title - a plain text description of its contents.

        :rtype: ``str``"""

        return self._title


    @title.setter
    def title(self, title):
        self._title = title



    @property
    def resolution(self):
        """The Pdb's resolution - a measure of the quality of the model(s)
        contained.

        :rtype: ``float``"""

        return self._resolution


    @resolution.setter
    def resolution(self, resolution):
        self._resolution = resolution



    @property
    def rfactor(self):
        """The Pdb's R-factor - a quality indicator which expresses the
        difference between the theoretical scattering pattern the model
        contained would produce, and the actual experimental data.

        :rtype: ``float``"""

        return self._rfactor


    @rfactor.setter
    def rfactor(self, rfactor):
        self._rfactor = rfactor


    @property
    def rfree(self):
        """The Pdb's Free R-factor - a less biased version of the R-factor.

        See `this explanation <https://bit.ly/2IaQ4ho>`_ for details.

        :rtype: ``float``"""

        return self._rfree


    @rfree.setter
    def rfree(self, rfree):
        self._rfree = rfree


    @property
    def rcount(self):
        """The Pdb's R-factor test set count - the number of observations used
        to calculate the free Rfactor.

        See `this explanation <https://bit.ly/2IaQ4ho>`_ for details.

        :rtype: ``float``"""

        return self._rcount


    @rcount.setter
    def rcount(self, rcount):
        self._rcount = rcount


    @property
    def organism(self):
        """The Pdb's source organism - where the molecule(s) contained come
        from.

        :rtype: ``str``"""

        return self._organism


    @organism.setter
    def organism(self, organism):
        self._organism = organism


    @property
    def expression_system(self):
        """The Pdb's expression organism - the organism the molecule(s)
        contained were actually expressed in for this experiment.

        :rtype: ``str``"""

        return self._expression_system


    @expression_system.setter
    def expression_system(self, expression_system):
        self._expression_system = expression_system


    @property
    def technique(self):
        """The Pdb's experimental technique, used to generate the model.

        :rtype: ``str``"""

        return self._technique


    @technique.setter
    def technique(self, technique):
        self._technique = technique


    @property
    def classification(self):
        """The Pdb's classification.

        :rtype: ``str``"""

        return self._classification


    @classification.setter
    def classification(self, classification):
        self._classification = classification


    @property
    def keywords(self):
        """The Pdb's keywords - helpful descriptive tags.

        :rtype: ``list``"""

        return self._keywords


    @property
    def biomolecules(self):
        """The Pdb's biomolecules - instructions for creating different
        biological assemblies using the Pdb's various chains.

        :rtype: ``list``"""

        return self._biomolecules


    def generate_assembly(self, id):
        """Creates a :py:class:`.Model` from the current model and the
        instructions contained in one of the Pdb's biomolecules, which you
        specify.

        :param int id: The biomolecule to use to generate the assembly.
        :raises ValueError: if you give an ID which doesn't correspond to a\
        biomolecule.
        :rtype: ``Model``"""

        model = self._models[0]
        for biomolecule in self._biomolecules:
            if biomolecule["id"] == id:
                break
        else:
            raise ValueError("No biomolecule with ID {}".format(id))
        new_chains = []
        for transformation in biomolecule["transformations"]:
            chains = [model.chain(id_) for id_ in transformation["chains"]]
            for chain in chains:
                new_chain = chain.copy()
                new_chain.transform(transformation["matrix"])
                new_chain.translate(transformation["vector"])
                new_chains.append(new_chain)
        return Model(*new_chains)


    @property
    def best_assembly(self):
        """Returms the 'best' biological assembly for this Pdb - the one with
        the lowest (most negative) delta energy.

        If there are no assemblies, ``None`` is returned.

        :rtype: ``Model``"""

        sorted_mol = sorted(
         self._biomolecules,
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


    def to_file_string(self):
        """Returns the file text that represents this Pdb.

        :rtype: ``str``"""

        from ..files.pdb2pdbdict import pdb_to_pdb_dict
        from ..files.pdbdict2pdbstring import pdb_dict_to_pdb_string
        pdb_dict = pdb_to_pdb_dict(self)
        return pdb_dict_to_pdb_string(pdb_dict)


    def save(self, path):
        """Saves the Pdb as a .pdb file.

        :param str path: The path to save to."""

        from ..files.utilities import string_to_file
        string_to_file(self.to_file_string(), path)
