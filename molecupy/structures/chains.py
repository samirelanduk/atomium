from .molecules import AtomicStructure, Residue, SmallMolecule
from . import matrix
from ..exceptions import NoResiduesError, BrokenHelixError, BrokenStrandError, DuplicateResiduesError

class ResiduicStructure(AtomicStructure):
    """Base class: :py:class:`.AtomicStructure`

    The base class for all structures which can be described as a set of
    residues.

    :param residues: A sequence of :py:class:`.Residue` objects in this\
    structure."""

    def __init__(self, *residues):
        if len(residues) == 0:
            raise NoResiduesError("Cannot make a ResiduicStructure with no residues")
        for residue in residues:
            if not isinstance(residue, Residue):
                raise TypeError(
                 "Can only make ResiduicStructure with Residues, not '%s'" % str(residue)
                )
        residue_ids = [residue.residue_id() for residue in residues]
        if len(residue_ids) != len(set(residue_ids)):
            duplicates = [id_ for id_ in residue_ids if residue_ids.count(id_) > 1]
            raise DuplicateResiduesError(
             "There are multiple residues with ID %s" % (duplicates[0])
            )
        self._residues = set(residues)


    def __repr__(self):
        return "<%s (%i residues)>" % (self.__class__.__name__, len(self._residues))


    def __getattr__(self, attribute):
        if attribute == "_atoms":
            atoms = set()
            for residue in self.residues():
                atoms.update(residue.atoms())
            return atoms
        else:
            return self.__getattribute__(attribute)


    def residues(self, include_missing=True):
        """Returns the residues in this structure as a ``set``.

        :param str include_missing: If ``False`` only residues present in the\
        PDB coordinates will be returned, and not missing ones.
        :rtype: ``set``"""

        if include_missing:
            return set(self._residues)
        else:
            return set([res for res in self._residues if not res.is_missing()])


    def add_residue(self, residue):
        """Adds a residue to the structure.

        :param Residue residue: The residue to add."""

        if not isinstance(residue, Residue):
            raise TypeError(
             "Can only add Residues to ResiduicStructures, not '%s'" % str(residue)
            )
        if residue.residue_id() in [residue.residue_id() for residue in self.residues()]:
            raise DuplicateResiduesError(
             "Cannot add residue with ID %i to %s as there is already an residue with that ID" % (
              residue.residue_id(), self
             )
            )
        self._residues.add(residue)


    def remove_residue(self, residue):
        """Removes a residue from the structure.

        :param Residue residue: The residue to remove."""

        self._residues.remove(residue)


    def get_residue_by_id(self, residue_id):
        """Returns the first residue that matches a given residue ID.

        :param str residue_id: The residue ID to search by.
        :rtype: :py:class:`.Residue` or ``None``"""

        if not isinstance(residue_id, str):
            raise TypeError(
             "Residue ID search must be by str, not '%s'" % str(residue_id)
            )
        for residue in self.residues():
            if residue.residue_id() == residue_id:
                return residue


    def get_residues_by_name(self, residue_name, include_missing=True):
        """Returns all the residues of a given name.

        :param str residue_name: The name to search by.
        :param str include_missing: If ``False`` only residues present in the\
        PDB coordinates will be returned, and not missing ones.
        :rtype: ``set`` of :py:class:`.Residue` objects."""

        if not isinstance(residue_name, str):
            raise TypeError(
             "Residue name search must be by str, not '%s'" % str(residue_name)
            )
        return set([
         residue for residue in self.residues(include_missing=include_missing)
          if residue.residue_name() == residue_name
        ])


    def get_residue_by_name(self, residue_name, include_missing=True):
        """Returns the first residue that matches a given name.

        :param str residue_name: The name to search by.
        :param str include_missing: If ``False`` only residues present in the\
        PDB coordinates will be returned, and not missing ones.
        :rtype: :py:class:`.Residue` or ``None``"""

        if not isinstance(residue_name, str):
            raise TypeError(
             "Residue name search must be by str, not '%s'" % str(residue_name)
            )
        for residue in self.residues(include_missing=include_missing):
            if residue.residue_name() == residue_name:
                return residue



class ResiduicSequence(ResiduicStructure):
    """Base class: :py:class:`ResiduicStructure`

    The base class for all structures which can be described as a sequence of
    residues.

    :param residues: A sequence of :py:class:`.Residue` objects in this\
    structure."""

    def __init__(self, *residues):
        ResiduicStructure.__init__(self, *residues)
        self._residues = list(residues)


    def residues(self, include_missing=True):
        """Returns the residues in this structure as a ``list``.

        :param str include_missing: If ``False`` only residues present in the\
        PDB coordinates will be returned, and not missing ones.
        :rtype: ``list``"""

        if include_missing:
            return list(self._residues)
        else:
            return [res for res in self._residues if not res.is_missing()]


    def add_residue(self, residue):
        """Adds a residue to the end of this sequence.

        :param Residue residue: The residue to add."""

        if not isinstance(residue, Residue):
            raise TypeError(
             "Can only add Residues to ResiduicSequences, not '%s'" % str(residue)
            )
        self._residues.append(residue)


    def sequence_string(self, include_missing=True):
        """Return the protein sequence of this chain as one letter codes.

        :param str include_missing: If ``False`` only residues present in the\
        PDB coordinates will be returned, and not missing ones.
        :rtype str: The protein sequence."""

        return "".join([RESIDUES.get(
         res.residue_name().upper(), "X"
        ) for res in self.residues(include_missing=include_missing)])



class Chain(ResiduicSequence):
    """Base class: :py:class:`ResiduicSequence`

    Represents chains - the polymeric units that make up most of PDB
    structures.

    :param chain_id: The chain's ID.
    :param residues: The residues in this chain."""

    def __init__(self, chain_id, *residues):
        if not isinstance(chain_id, str):
            raise TypeError("'%s' is not a valid chain_id" % str(chain_id))
        self._chain_id = chain_id
        ResiduicSequence.__init__(self, *residues)
        for residue in self._residues:
            residue._chain = self
        self._alpha_helices = set()
        self._beta_strands = set()
        self._model = None


    def __repr__(self):
        return "<Chain %s (%i residues)>" % (self._chain_id, len(self._residues))


    def chain_id(self):
        """Returns the chain's ID.

        :rtype: ``str``"""

        return self._chain_id


    def add_residue(self, residue):
        ResiduicSequence.add_residue(self, residue)
        residue._chain = self


    def remove_residue(self, residue):
        ResiduicSequence.remove_residue(self, residue)
        residue._chain = None


    def alpha_helices(self):
        """Returns the :py:class:`AlphaHelix` objects on this chain.

        :returns: ``set`` of ``AlphaHelix`` objects"""

        return set(self._alpha_helices)


    def beta_strands(self):
        """Returns the :py:class:`BetsStrand` objects on this chain.

        :returns: ``set`` of ``BetaStrand`` objects"""

        return set(self._beta_strands)


    def model(self):
        """Returns the :py:class:`.Model` that the chain inhabits.

        :rtype: ``Model``"""

        return self._model


    def get_helix_by_id(self, helix_id):
        """Retrurns the first alpha helix that matches a given helix ID.

        :param str helix_id: The helix ID to search by.
        :rtype: :py:class:`.AlphaHelix` or ``None``"""

        if not isinstance(helix_id, str):
            raise TypeError("Helix ID search must be by str, not '%s'" % str(helix_id))
        for helix in self.alpha_helices():
            if helix.helix_id() == helix_id:
                return helix


    def get_strand_by_id(self, strand_id):
        """Retrurns the first beta strand that matches a given strand ID.

        :param str strand_id: The strand ID to search by.
        :rtype: :py:class:`.BetsStrand` or ``None``"""

        if not isinstance(strand_id, str):
            raise TypeError("Strand ID search must be by str, not '%s'" % str(strand_id))
        for strand in self.beta_strands():
            if strand.strand_id() == strand_id:
                return strand


    generate_residue_distance_matrix = matrix.generate_residue_distance_matrix
    """Creates a 'distance matrix' as an
    `OmniCanvas <http://omnicanvas.readthedocs.io/>`_ canvas. The distance
    between any two residues are represented as gradients of colour.

    To output the resulting canvas to svg, the ``.save("filename.svg")``
    method can be used.

    :param int dimension: The width and height of the matrix in pixels
    :param int close_color: The colour of near residues (defaults to green)
    :param int far_color: The colour of far residues (defaults to red)
    :param int cutoff: The distance in Angstroms that the colour scale should end at (defaults to 40 Angstroms)
    :param subsequence: A sequence of two :py:class:`PdbResidue` objects that constitute the beginning and end of some subsequence, which will be marked on the matrix.
    :rtype: `OmniCanvas <http://omnicanvas.readthedocs.io/>`_ canvas"""



class BindSite(ResiduicStructure):
    """Base class: :py:class:`ResiduicStructure`

    Represents binding sites - the residue clusters that mediate ligand
    binding.

    :param site_id: The site's ID.
    :param residues: The residues in this chain."""

    def __init__(self, site_id, *residues):
        if not isinstance(site_id, str):
            raise TypeError("'%s' is not a valid site_id" % str(site_id))
        self._site_id = site_id
        self._ligand = None
        self._model = None
        ResiduicStructure.__init__(self, *residues)


    def __repr__(self):
        return "<BindSite %s (%i residues)>" % (self._site_id, len(self._residues))


    def site_id(self):
        """Returns the site's ID.

        :rtype: ``str``"""

        return self._site_id


    def ligand(self, ligand=None):
        """Returns or sets the site's :py:class:`.SmallMolecule` ligand.

        :param SmallMolecule ligand: If given, the ligand will be set to this.
        :rtype: ``SmallMolecule``"""

        if ligand is None:
            return self._ligand
        else:
            if not isinstance(ligand, SmallMolecule):
                raise TypeError(
                 "'%s' is not a valid ligand" % str(ligand)
                )
            self._ligand = ligand
            ligand._bind_site = self


    def model(self):
        """Returns the :py:class:`.Model` that the site inhabits.

        :rtype: ``Model``"""

        return self._model


    def continuous_sequence(self):
        """If the residues are on the same chain, this will return a continuous
        sequence that contains all residues in this site, otherwise ``None``.

        :rtype: ResiduicSequence"""

        if len(set([res.chain() for res in self.residues() if res.chain()])) == 1:
            chain = list(self.residues())[0].chain()
            min_index = min([chain.residues().index(res) for res in self.residues()])
            max_index = max([chain.residues().index(res) for res in self.residues()])
            return ResiduicSequence(*chain.residues()[min_index:max_index])
        else:
            return None



class AlphaHelix(ResiduicSequence):
    """Base class: :py:class:`ResiduicSequence`

    Represents alpha helices.

    :param str helix_id: The helix's ID.
    :param residues: The residues in this helix.
    :param str helix_class: The classification of the helix.
    :param str comment: Any comment associated with this helix."""

    def __init__(self, helix_id, *residues, helix_class=None, comment=None):
        if not isinstance(helix_id, str):
            raise TypeError("'%s' is not a valid helix_id" % str(helix_id))
        self._helix_id = helix_id
        if len(set([res.chain() for res in residues])) != 1:
            raise BrokenHelixError(
             "Cannot make helix %s with residues from multiple chains" % helix_id
            )
        ResiduicSequence.__init__(self, *residues)
        if helix_class is not None and not isinstance(helix_class, str):
            raise TypeError("'%s' is not a valid helix_class" % str(helix_class))
        self._helix_class = helix_class
        if comment is not None and not isinstance(comment, str):
            raise TypeError("'%s' is not a valid comment" % str(comment))
        self._comment = comment
        if self.chain():
            self.chain()._alpha_helices.add(self)


    def __repr__(self):
        return "<AlphaHelix %s (%i residues)>" % (self._helix_id, len(self._residues))


    def helix_id(self):
        """Returns the helix's ID.

        :rtype: ``str``"""

        return self._helix_id


    def helix_class(self, helix_class=None):
        """Returns or sets the helix's classification.

        :param str helix_class: If given, the class will be set to this.
        :rtype: ``str``"""

        if helix_class is None:
            return self._helix_class
        else:
            if not isinstance(helix_class, str):
                raise TypeError(
                 "'%s' is not a valid helix_class" % str(helix_class)
                )
            self._helix_class = helix_class


    def comment(self, comment=None):
        """Returns or sets the helix's comment.

        :param str comment: If given, the comment will be set to this.
        :rtype: ``str``"""

        if comment is None:
            return self._comment
        else:
            if not isinstance(comment, str):
                raise TypeError(
                 "'%s' is not a valid comment" % str(comment)
                )
            self._comment = comment


    def chain(self):
        """Returns the chain that this helix is on.

        :rtype: ``Chain``"""

        return self.residues()[0].chain()


    def add_residue(self, residue):
        if residue.chain() is not self.chain():
            raise BrokenHelixError(
             "Cannot add %s to %s as their chains don't match" % (str(residue), str(self))
            )
        ResiduicSequence.add_residue(self, residue)



class BetaStrand(ResiduicSequence):
    """Base class: :py:class:`ResiduicSequence`

    Represents beta strands.

    :param str strand_id: The strand's ID.
    :param residues: The residues in this strand.
    :param int sense: The sense of the strand with respect to the prior\
    strand."""

    def __init__(self, strand_id, sense, *residues):
        if not isinstance(strand_id, str):
            raise TypeError("'%s' is not a valid strand_id" % str(strand_id))
        self._strand_id = strand_id
        if not isinstance(sense, int):
            raise TypeError("'%s' is not a valid sense value" % str(sense))
        if not (-1 <= sense <= 1):
            raise ValueError("sense can only be -1, 0 or 1 - not %i" % sense)
        self._sense = sense
        if len(set([res.chain() for res in residues])) != 1:
            raise BrokenStrandError(
             "Cannot make strand %s with residues from multiple chains" % strand_id
            )
        ResiduicSequence.__init__(self, *residues)
        if self.chain():
            self.chain()._beta_strands.add(self)


    def __repr__(self):
        return "<BetaStrand %s (%i residues)>" % (self._strand_id, len(self._residues))


    def strand_id(self):
        """Returns the strand's ID.

        :rtype: ``str``"""

        return self._strand_id


    def sense(self, sense=None):
        """Returns or sets the strand's sense with respect to the previous
        strand.

        :param int sense: If given, the sense will be set to this.
        :rtype: ``int``"""

        if sense is None:
            return self._sense
        else:
            if not isinstance(sense, int):
                raise TypeError(
                 "'%s' is not a valid sense value" % str(sense)
                )
            if not (-1 <= sense <= 1):
                raise ValueError("sense can only be -1, 0 or 1 - not %i" % sense)
            self._sense = sense


    def chain(self):
        """Returns the chain that this strand is on.

        :rtype: ``Chain``"""

        return self.residues()[0].chain()


    def add_residue(self, residue):
        if residue.chain() is not self.chain():
            raise BrokenStrandError(
             "Cannot add %s to %s as their chains don't match" % (str(residue), str(self))
            )
        ResiduicSequence.add_residue(self, residue)



RESIDUES = {
 "GLY": "G", "ALA": "A", "LEU": "L", "MET": "M", "PHE": "F",
 "TRP": "W", "LYS": "K", "GLN": "Q", "GLU": "E", "SER": "S",
 "PRO": "P", "VAL": "V", "ILE": "I", "CYS": "C", "TYR": "Y",
 "HIS": "H", "ARG": "R", "ASN": "N", "ASP": "D", "THR": "T"
}
