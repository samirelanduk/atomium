"""This module contains the structures used to represent those found in PDB
files. This is where the bulk of the non-parsing work is done."""

import math
import omnicanvas
from collections import Counter
import warnings
from .exceptions import *
from pprint import pprint

class PdbAtom:
    """This class epresents atoms - the fundamental chemical building blocks.

    :param float x: The atom's x-coordinate.
    :param float y: The atom's y-coordinate.
    :param float z: The atom's z-coordinate.
    :param str element: The atom's element.
    :param int atom_id: The atom's id.
    :param str atom_name: The atom's name.

    .. py:attribute:: x:

        The atom's x-coordinate.

    .. py:attribute:: y:

        The atom's y-coordinate.

    .. py:attribute:: z:

        The atom's z-coordinate.

    .. py:attribute:: element:

        The atom's element.

    .. py:attribute:: atom_id:

        The atom's id.

    .. py:attribute:: atom_name:

        The atom's name.

    .. py:attribute:: covalent_bonds:

        A ``set`` of :py:class:`CovalentBond` objects attached to this atom.

    .. py:attribute:: molecule:

        The :py:class:`PdbResidue` or :py:class:`PdbSmallMolecule` the atom \
        belongs to."""

    def __init__(self, x, y, z, element, atom_id, atom_name):
        for coordinate in (x, y, z):
            if not isinstance(coordinate, float):
                raise TypeError("'%s' is not a valid coordinate" % str(coordinate))
        self.x = x
        self.y = y
        self.z = z

        if not isinstance(element, str):
            raise TypeError("'%s' is not a valid element" % str(element))
        if len(element) == 0:
            raise InvalidElementError("Atom's element can't be an empty string")
        elif len(element) > 2:
            raise InvalidElementError("'%s' is not a valid element" % element)
        self.element = element

        if not isinstance(atom_id, int):
            raise TypeError("'%s' is not a valid atom_id" % str(atom_id))
        self.atom_id = atom_id

        if not isinstance(atom_name, str):
            raise TypeError("'%s' is not a valid atom_name" % str(atom_name))
        self.atom_name = atom_name

        self.covalent_bonds = set()
        self.molecule = None


    def __repr__(self):
        return "<Atom %i (%s)>" % (self.atom_id, self.element)


    def get_mass(self):
        """Returns the atom's mass.

        :rtype: float"""

        return PERIODIC_TABLE.get(self.element.upper(), 0)


    def distance_to(self, other_atom):
        """Returns the distance to some other atom.

        :param PdbAtom other_atom: The other atom.
        :rtype: float"""

        x_sum = math.pow((other_atom.x - self.x), 2)
        y_sum = math.pow((other_atom.y - self.y), 2)
        z_sum = math.pow((other_atom.z - self.z), 2)
        distance = math.sqrt(x_sum + y_sum + z_sum)
        return distance


    def covalent_bond_to(self, other_atom):
        """Covalently bonds the atom to another atom.

        :param PdbAtom other_atom: The other atom."""

        CovalentBond(self, other_atom)


    def get_covalent_bonded_atoms(self):
        """Returns the atoms this atom is covalently bonded to.

        :returns: ``set`` of atoms bonded to this one."""

        covalent_bonded_atoms = set()
        for bond in self.covalent_bonds:
            for atom in bond.atoms:
                if atom is not self: covalent_bonded_atoms.add(atom)
        return covalent_bonded_atoms


    def get_covalent_accessible_atoms(self, already_checked=None):
        """Returns all the atoms reachable from this atom via covalent bonds.

        :returns: ``set`` of atoms."""

        already_checked = already_checked if already_checked else set()
        already_checked.add(self)
        while len(self.get_covalent_bonded_atoms().difference(already_checked)) > 0:
            picked = list(self.get_covalent_bonded_atoms().difference(already_checked))[0]
            picked.get_covalent_accessible_atoms(already_checked)
        already_checked_copy = already_checked.copy()
        already_checked_copy.discard(self)
        return already_checked_copy



class CovalentBond:
    """Represents a covalent bond between two atoms.

    :param PdbAtom atom1: The first atom.
    :param PdbAtom atom2: The second atom.

    .. py:attribute:: covalent_bonds:

        The ``set`` of :py:class:`PdbAtom` objects sharing this bond."""

    def __init__(self, atom1, atom2):
        if not isinstance(atom1, PdbAtom) or not isinstance(atom2, PdbAtom):
            raise TypeError(
             "Can only bond atoms, not %s to '%s'" % (str(atom1), str(atom2))
            )
        if atom1.distance_to(atom2) > 10:
            warning = "The bond between atom %s and atom %s is %.2f Angstroms" % (
             str(atom1.atom_id) if atom1.atom_id else atom1.element,
             str(atom2.atom_id) if atom2.atom_id else atom2.element,
             atom1.distance_to(atom2)
            )
            warnings.warn(warning, LongBondWarning)

        self.atoms = set((atom1, atom2))
        atom1.covalent_bonds.add(self)
        atom2.covalent_bonds.add(self)


    def __repr__(self):
        atom1_element, atom2_element = [a.element for a in self.atoms]
        return "<CovalentBond (%s-%s)>" % (atom1_element, atom2_element)


    def get_bond_length(self):
        """Returns the length of the covalent bond.

        :rtype: float"""

        atom1, atom2 = self.atoms
        return atom1.distance_to(atom2)



class AtomicStructure:
    """The base class for all structures which are composed of atoms.

    :param atoms: A sequence of :py:class:`PdbAtom` objects.

    .. py:attribute:: atoms:

        A ``set`` of :py:class:`PdbAtom` objects in this structure."""

    def __init__(self, *atoms):
        if not all(isinstance(atom, PdbAtom) for atom in atoms):
            non_atoms = [atom for atom in atoms if not isinstance(atom, PdbAtom)]
            raise TypeError("AtomicStructure needs atoms, not '%s'" % non_atoms[0])
        if not atoms:
            raise NoAtomsError("Cannot make AtomicStructure with zero atoms")
        atom_ids = [atom.atom_id for atom in atoms]
        if len(set(atom_ids)) < len(atom_ids):
            raise DuplicateAtomIdError("Cannot make AtomicStructure with duplicate atom_ids")
        self.atoms = set(atoms)


    def __repr__(self):
        return "<AtomicStructure (%i atoms)>" % len(self.atoms)


    def __contains__(self, atom):
        return atom in self.atoms


    def get_mass(self):
        """Returns the mass of the structure by summing the mass of all its
        atoms.

        :rtype: float"""

        return sum([atom.get_mass() for atom in self.atoms])


    def get_atom_by_name(self, atom_name):
        """Retrurns the first atom that matches a given atom name.

        :param str atom_name: The atom name to search by.
        :rtype: :py:class:`PdbAtom` or ``None``
        :raises TypeError: if the name given isn't a string."""

        if not isinstance(atom_name, str):
            raise TypeError("Atom name search must be by str, not '%s'" % str(atom_name))
        for atom in self.atoms:
            if atom.atom_name == atom_name:
                return atom


    def get_atoms_by_name(self, atom_name):
        """Retruns all the atoms that matches a given atom name.

        :param str atom_name: The atom name to search by.
        :rtype: ``set`` of :py:class:`PdbAtom` objects
        :raises TypeError: if the name given isn't a string."""

        if not isinstance(atom_name, str):
            raise TypeError("Atom name search must be by str, not '%s'" % str(atom_name))
        return set([atom for atom in self.atoms if atom.atom_name == atom_name])


    def get_atom_by_element(self, element):
        """Retrurns the first atom that matches a given element.

        :param str element: The element to search by.
        :rtype: :py:class:`PdbAtom` or ``None``
        :raises TypeError: if the element given isn't a string."""

        if not isinstance(element, str):
            raise TypeError("Atom element search must be by str, not '%s'" % str(element))
        for atom in self.atoms:
            if atom.element == element:
                return atom


    def get_atoms_by_element(self, element):
        """Retruns all the atoms a given element.

        :param str element: The element to search by.
        :rtype: ``set`` of :py:class:`PdbAtom` objects
        :raises TypeError: if the element given isn't a string."""

        if not isinstance(element, str):
            raise TypeError("Atom element search must be by str, not '%s'" % str(element))
        return set([atom for atom in self.atoms if atom.element == element])


    def get_atom_by_id(self, atom_id):
        """Retrurns the first atom that matches a given atom ID.

        :param str atom_id: The atom ID to search by.
        :rtype: :py:class:`PdbAtom` or ``None``
        :raises TypeError: if the name given isn't an integer."""

        if not isinstance(atom_id, int):
            raise TypeError("Atom ID search must be by int, not '%s'" % str(atom_id))
        for atom in self.atoms:
            if atom.atom_id == atom_id:
                return atom


    def get_formula(self):
        """Retrurns the formula (count of each atom) of the structure.

        :rtype: ``Counter``"""
        return Counter([
         atom.element for atom in self.atoms if atom.element != "H"
        ])


    def get_external_contacts_with(self, other_structure):
        """Returns the set of all 'contacts' with another atomic structure,
        where a contact is defined as any atom-atom pair with an inter-atomic
        distance less than or equal to four Angstroms.

        If the other atomic structure has atoms which are also in this atomic
        structure, those atoms will not be counted as part of the other
        structure.

        :param AtomicStructure other_structure: The other atomic\
        structure to compare to.
        :rtype: ``set`` of ``frozenset`` contacts."""

        contacts = set()
        our_atoms = list(self.atoms)
        their_atoms = [a for a in list(other_structure.atoms) if a not in our_atoms]
        for atom in our_atoms:
            for other_atom in their_atoms:
                if atom.distance_to(other_atom) <= 4:
                    contacts.add(frozenset((atom, other_atom)))
        return contacts


    def get_internal_contacts(self):
        """Returns the set of all atomic contacts within the atoms of an atomic
        structure, where a contact is defined as any atom-atom pair with an
        inter-atomic distance less than or equal to four Angstroms.

        Contacts between atoms covalently bonded to each other will be ignored,
        as will contacts between atoms separated by just two covalent bonds.

        :rtype: ``set`` of ``frozenset`` contacts."""

        contacts = set()
        atoms = list(self.atoms)
        for index, atom in enumerate(atoms[:-1]):
            too_close_atoms = set((atom,))
            for bonded_atom in atom.get_covalent_bonded_atoms():
                too_close_atoms.add(bonded_atom)
            second_tier_atoms = set()
            for close_atom in too_close_atoms:
                for bonded_atom in close_atom.get_covalent_bonded_atoms():
                    second_tier_atoms.add(bonded_atom)
            too_close_atoms.update(second_tier_atoms)
            for other_atom in atoms[index + 1:]:
                if other_atom not in too_close_atoms and atom.distance_to(other_atom) <= 4:
                    contacts.add(frozenset((atom, other_atom)))
        return contacts



class PdbSmallMolecule(AtomicStructure):
    """Base class: :py:class:`AtomicStructure`

    Represents the ligands, solvent molecules, and other non-polymeric
    molecules in a PDB structure.

    :param str molecule_id: The molecule's ID
    :param str molecule_name: The molecule's name
    :param atoms: The molecule's atoms

    .. py:attribute:: molecule_id:

        The molecule's ID.

    .. py:attribute:: molecule_name:

        The molecule's name.

   .. py:attribute:: model:

        The :py:class:`PdbModel` that the molecule belongs to."""

    def __init__(self, molecule_id, molecule_name, *atoms):
        if not isinstance(molecule_id, str):
            raise TypeError("'%s' is not a valid molecule_id" % str(molecule_id))
        self.molecule_id = molecule_id

        if not isinstance(molecule_name, str):
            raise TypeError("'%s' is not a valid molecule_name" % str(molecule_name))
        self.molecule_name = molecule_name

        AtomicStructure.__init__(self, *atoms)
        for atom in self.atoms:
            atom.molecule = self
        self.model = None


    def __repr__(self):
        return "<SmallMolecule (%s)>" % self.molecule_name


    def get_binding_site(self):
        """Returns the py:class:`PdbSite` that is listed as binding the
        molecule.

        :rtype: :py:class:`PdbSmallMolecule` or ``None``"""

        if self.model:
            for site in self.model.sites:
                if site.ligand is self:
                    return site


    def calculate_binding_site(self):
        """Attempts to predict the residues that the molecule binds to by using
        an atomic distance cutoff of 3.0 Angstroms.

        :rtype: :py:class:`PdbSmallMolecule` or ``None``"""

        if self.model:
            residues = set()
            for chain in self.model.chains:
                residues.update(chain.residues)
            close_residues = set()
            for atom in self.atoms:
                for residue in residues:
                    if any(r_atom.distance_to(atom) <= 3 for r_atom in residue.atoms):
                        close_residues.add(residue)
            if close_residues:
                site = PdbSite("calc", *close_residues)
                site.ligand = self
                return site



class PdbResidue(AtomicStructure):
    """Base class: :py:class:`AtomicStructure`

    A Residue on a chain.

    :param str residue_id: The residue's ID
    :param str residue_name: The residue's name
    :param atoms: The residue's atoms

    .. py:attribute:: residue_id:

        The residue's ID.

    .. py:attribute:: residue_name:

        The redidue's name - usually a standard three letter code.

    .. py:attribute:: downstream_residue:

        The next residue in the chain to which this residue is connected.

    .. py:attribute:: upstream_residue:

        The previous residue in the chain to which this residue is connected."""

    def __init__(self, residue_id, residue_name, *atoms):
        if not isinstance(residue_id, str):
            raise TypeError("'%s' is not a valid residue_id" % str(residue_id))
        self.residue_id = residue_id

        if not isinstance(residue_name, str):
            raise TypeError("'%s' is not a valid molecule_name" % str(residue_name))
        self.residue_name = residue_name

        AtomicStructure.__init__(self, *atoms)
        for atom in self.atoms:
            atom.molecule = self

        self.chain = None
        self.downstream_residue = None
        self.upstream_residue = None


    def __repr__(self):
        return "<Residue (%s)>" % self.residue_name


    def connect_to(self, downstream_residue, this_atom, their_atom):
        """Connects the residue to its downstream residue (in the case of
        proteins, proteins start at the N terminus and end at the C terminus).
        This creates the necessary :py:class`CovalentBond` and also lets the
        residue know what its neighbours are.

        :param PdbResidue downstream_residue: The next residue to connect to.
        :param PdbAtom this_atom: The atom in this residue facilitating the\
        connection.
        :param PdbAtom their_atom: The atom in the other residue facilitating\
        the connection.
        :raises ValueError: If the atoms are not in the correct residue."""

        if this_atom not in self:
            raise ValueError(
             "%s is not in %s - cannot connect." % (str(this_atom), str(self))
            )
        if their_atom not in downstream_residue:
            raise ValueError(
             "%s is not in %s - cannot connect." % (str(their_atom), str(downstream_residue))
            )
        this_atom.covalent_bond_to(their_atom)
        self.downstream_residue = downstream_residue
        downstream_residue.upstream_residue = self


    def get_alpha_carbon(self):
        """Returns the alpha carbon of this residue.

        :rtype: :py:class:`PdbAtom`"""

        atom = self.get_atom_by_name("CA")
        if not atom:
            atom = self.get_atom_by_element("C")
        if not atom:
            atom = list(self.atoms)[0]
        return atom



class ResiduicStructure(AtomicStructure):
    """Base class: :py:class:`AtomicStructure`

    The base class for all structures which can be described as a set of
    residues.

    :param residues: A ``set`` of :py:class:`PdbResidue` objects in this\
    structure.

    .. py:attribute:: residues:

         A ``set`` of :py:class:`PdbResidue` objects in this structure."""

    def __init__(self, *residues):
        if not all(isinstance(residue, PdbResidue) for residue in residues):
            non_residues = [
             residue for residue in residues if not isinstance(residue, PdbResidue)
            ]
            raise TypeError("ResiduicStructure needs residues, not '%s'" % non_residues[0])
        if not residues:
            raise NoResiduesError("Cannot make ResiduicStructure with zero residues")
        residue_ids = [residue.residue_id for residue in residues]
        if len(set(residue_ids)) < len(residue_ids):
            raise DuplicateResidueIdError("Cannot make ResiduicStructure with duplicate residue_ids")
        self.residues = set(residues)


    def __repr__(self):
        return "<ResiduicStructure (%i residues)>" % len(self.residues)


    def __getattr__(self, attribute):
        if attribute == "atoms":
            atoms = set()
            for residue in self.residues:
                atoms.update(residue.atoms)
            return atoms
        else:
            return self.__getattribute__(attribute)


    def __contains__(self, obj):
        return obj in self.residues or obj in self.atoms


    def get_residue_by_name(self, residue_name):
        """Retrurns the first residue that matches a given residue name.

        :param str residue_name: The residue name to search by.
        :rtype: :py:class:`PdbResidue` or ``None``
        :raises TypeError: if the name given isn't a string."""

        if not isinstance(residue_name, str):
            raise TypeError("Residue name search must be by str, not '%s'" % str(residue_name))
        for residue in self.residues:
            if residue.residue_name == residue_name:
                return residue


    def get_residues_by_name(self, residue_name):
        """Retruns all the residues that matches a given residue name.

        :param str residue_name: The residue name to search by.
        :rtype: ``set`` of :py:class:`PdbResidue` objects
        :raises TypeError: if the name given isn't a string."""

        if not isinstance(residue_name, str):
            raise TypeError("Residue name search must be by str, not '%s'" % str(residue_name))
        return set([residue for residue in self.residues if residue.residue_name == residue_name])


    def get_residue_by_id(self, residue_id):
        """Retrurns the first residue that matches a given residue ID.

        :param str residue_id: The residue ID to search by.
        :rtype: :py:class:`PdbResidue` or ``None``
        :raises TypeError: if the name given isn't an integer."""

        if not isinstance(residue_id, str):
            raise TypeError("Residue ID search must be by str, not '%s'" % str(residue_id))
        for residue in self.residues:
            if residue.residue_id == residue_id:
                return residue



class ResiduicSequence(ResiduicStructure):
    """Base class: :py:class:`ResiduicStructure`

    The base class for all structures which can be described as a sequence of
    residues.

    :param residues: A ``list`` of :py:class:`PdbResidue` objects in this\
    structure.

    .. py:attribute:: residues:

         A ``list`` of :py:class:`PdbResidue` objects in this structure."""

    def __init__(self, *residues):
        ResiduicStructure.__init__(self, *residues)
        self.residues = list(residues)


    def __repr__(self):
        return "<ResiduicSequence (%i residues)>" % len(self.residues)


    def get_sequence_string(self):
        """Return the protein sequence of this chain as one letter codes.

        :rtype str: The protein sequence."""

        return "".join(
         [RESIDUES.get(res.residue_name, "X") for res in self.residues]
        )



class PdbChain(ResiduicSequence):
    """Base class: :py:class:`ResiduicSequence`

    Represents PDB chains - the polymeric units that make up most of PDB
    structures.

    :param chain_id: The chain's ID.
    :param residues: The residues in this chain.

    .. py:attribute:: chain_id:

         A chain's ID.

    .. py:attribute:: model:

         The :py:class:`PdbModel` that the chain belongs to."""

    def __init__(self, chain_id, *residues):
        if not isinstance(chain_id, str):
            raise TypeError("'%s' is not a valid chain_id" % str(chain_id))
        self.chain_id = chain_id

        ResiduicSequence.__init__(self, *residues)
        for residue in self.residues:
            residue.chain = self
        self.model = None
        self.missing_residues = []
        self.alpha_helices = set()
        self.beta_strands = set()


    def __repr__(self):
        return "<Chain %s (%i residues)>" % (self.chain_id, len(self.residues))


    def get_alpha_helix_by_id(self, helix_id):
        """Retrurns the first alpha helix that matches a given molecule ID.

        :param str helix_id: The helix ID to search by.
        :rtype: :py:class:`PdbAlphaHelix` or ``None``
        :raises TypeError: if the name given isn't an string."""

        if not isinstance(helix_id, str):
            raise TypeError("Can only search alpha helix IDs by str")
        for alpha_helix in self.alpha_helices:
            if alpha_helix.helix_id == helix_id:
                return alpha_helix


    def get_beta_strand_by_id(self, strand_id):
        """Retrurns the first beta strand that matches a given molecule ID.

        :param str strand_id: The strand ID to search by.
        :rtype: :py:class:`PdbBetaStrand` or ``None``
        :raises TypeError: if the name given isn't an string."""
        
        if not isinstance(strand_id, int):
            raise TypeError("Can only search beta strand IDs by int")
        for beta_strand in self.beta_strands:
            if beta_strand.strand_id == strand_id:
                return beta_strand


    def get_all_residue_ids(self):
        """Returns a ``list`` of residue IDs for this chain, in order,
        `including` the residues of missing residues.

        :rtype: ``list`` of ``str`` objects."""

        residue_ids = [residue.residue_id for residue in self.residues]
        for missing_id in self.missing_residues:
            closest_smaller_id = None
            if not _residue_id_is_greater_than_residue_id(residue_ids[0], missing_id):
                closest_smaller_id = sorted(
                 residue_ids,
                 key=lambda k: _residue_id_to_int(missing_id) - _residue_id_to_int(
                  k
                 ) if _residue_id_to_int(missing_id) > _residue_id_to_int(
                  k
                 ) else float("inf"))[0]
            residue_ids.insert(
             residue_ids.index(closest_smaller_id) + 1 if closest_smaller_id else 0,
             missing_id
            )
        return residue_ids


    def generate_residue_distance_matrix(self, dimension=700, close_color=120,
     far_color=0, cutoff=40, subsequence=None):
        """Creates a 'distance matrix' as an
        `OmniCanvas <http://omnicanvas.readthedocs.io/>`_ canvas. The distance
        between any two residues are represented as gradients of colour.

        To output the resulting canvas to svg, the ``.save("filename.svg")``
        method can be used.

        :param int dimension: The width and height of the matrix in pixels
        :param int close_color: The colour of near residues (defaults to green)
        :param int far_color: The colour of far residues (defaults to red)
        :param int cutoff: The distance in Angstroms that the colour scale\
        should end at (defaults to 40 Angstroms)
        :param subsequence: A sequence of two :py:class:`PdbResidue` objects\
        that constitute the beginning and end of some subsequence, which will\
        be marked on the matrix.
        :rtype: `OmniCanvas <http://omnicanvas.readthedocs.io/>`_ canvas"""

        # Validation
        if not isinstance(close_color, int):
            raise TypeError("close_color must be int, not '%s'" % str(close_color))
        if not isinstance(far_color, int):
            raise TypeError("far_color must be int, not '%s'" % str(far_color))
        if not 0 <= close_color < 360:
            raise ValueError("close_color must be between 0 and 360, not %i" % close_color)
        if not 0 <= far_color < 360:
            raise ValueError("far_color must be between 0 and 360, not %i" % far_color)
        if not isinstance(cutoff, int) and not isinstance(cutoff, float):
            raise TypeError("cutoff must be numeric, not '%s'" % str(dimension))
        if not isinstance(dimension, int):
            raise TypeError("dimension must be int, not '%s'" % str(dimension))
        if subsequence:
            if not isinstance(subsequence, list) and not isinstance(subsequence, tuple):
                raise TypeError("subsequence must be a list or tuple")
            for residue in subsequence:
                if not isinstance(residue, PdbResidue):
                    raise TypeError("Only PdbResidues can be given for subsequence")
                if residue not in self.residues:
                    raise ValueError(
                     "%s does not have subsequence residue %s" % (str(self), str(residue))
                    )
            if len(subsequence) != 2:
                raise ValueError(
                 "Only two residues can be used to define a subsequence, not %s" % str(subsequence)
                )
        # Set up canvas
        matrix = omnicanvas.Canvas(dimension, dimension)
        residues = [self.get_residue_by_id(id_) for id_ in self.get_all_residue_ids()]

        # Set up parameters
        padding_proportion = 0.09
        padding = padding_proportion * dimension
        plot_dimension = dimension - (2 * padding)
        chain_length = len(self.get_all_residue_ids())
        cell_dimension = plot_dimension / chain_length
        plot_width = dimension - (2 * padding)
        bar_width = 4
        bar_left = (dimension / 2) - (bar_width / 2) + 5
        hypoteneuse = math.sqrt((plot_width ** 2) + (plot_width ** 2))
        bar_top = (dimension / 2) - (hypoteneuse / 2) + 5
        diagonal_chunk = hypoteneuse / len(residues)
        chain_color = 80
        helix_color = 325
        strand_color = 182
        tick = 0
        if len(residues) >= 10000:
            tick = 5000
        elif len(residues) >= 1000:
            tick = 500
        elif len(residues) >= 100:
            tick = 50
        elif len(residues) >= 10:
            tick = 5
        else:
            tick = 1

        # Calculate distances
        distances = []
        for index, residue2 in enumerate(residues[1:][::-1]):
            row = []
            for residue1 in residues[:0 - (index + 1)]:
                row.append({
                 "residue1": residue1,
                 "residue2": residue2,
                 "distance": residue1.get_alpha_carbon().distance_to(residue2.get_alpha_carbon()) if residue1 and residue2 else None
                })
            distances.append(row)

        # Add cells
        for row_index, row in enumerate(distances):
            for cell_index, cell in enumerate(row):
                color = "#FFFFFF"
                if cell["distance"] is not None:
                    fraction = cell["distance"] / cutoff if cell["distance"] <= cutoff else 1
                    if far_color >= close_color:
                        distance_from_start = fraction * (far_color - close_color)
                        color = close_color + distance_from_start
                    else:
                        distance_from_start = fraction * (close_color - far_color)
                        color = close_color - distance_from_start
                    color = omnicanvas.hsl_to_rgb(color, 100, 50)

                matrix.add_rectangle(
                 padding + (cell_index * cell_dimension),
                 padding + (row_index * cell_dimension),
                 cell_dimension + 0.5,
                 cell_dimension + 0.5,
                 fill_color=color,
                 line_width=0,
                 data={
                  "onmouseover": "cellHovered(this)",
                  "onmouseleave": "cellLeft(this)",
                  "data": "%s,%i,%s,%i,%.2f" % (
                   cell["residue1"].residue_name if cell["residue1"] else "???",
                   cell_index + 1,
                   cell["residue2"].residue_name if cell["residue2"] else "???",
                   len(residues) - row_index,
                   cell["distance"] if cell["distance"] is not None else 0.0
                  )
                 }
                )

        # Subsequence
        if subsequence:
            x = (residues.index(subsequence[0]) * cell_dimension) + padding
            y = dimension - (((residues.index(subsequence[1]) + 1) * cell_dimension) + padding)
            matrix.add_line(
             x, dimension - padding, x, y,
             line_width=1.5,
             line_style=".."
            )
            matrix.add_line(
             dimension - padding, y, x, y,
             line_width=1.5,
             line_style=".."
            )

        # Add gridlines, border and labels
        residue_number = 0
        while residue_number <= len(residues) - 1:
            x = padding + (residue_number * cell_dimension)
            y = dimension - x
            matrix.add_line(
             x, padding, x, dimension - padding
            )
            matrix.add_line(
             padding, y, dimension - padding, y
            )
            if residue_number != len(residues) - 1:
                matrix.add_text(
                 x + (0.5 * cell_dimension), padding * 0.75, str(residue_number + 1)
                )
            if residue_number != 0:
                matrix.add_text(
                 padding - 2, y - (0.5 * cell_dimension), str(residue_number + 1),
                 horizontal_align="left"
                )
            residue_number += tick
        matrix.add_text(dimension / 2, padding * 0.35, "Residue 1")
        matrix.add_text(
         padding * 0.35, dimension / 2, "Residue 2", rotation=(
          padding * 0.35, dimension / 2, 270
         )
        )
        matrix.add_polygon(
         dimension - padding, padding,
         dimension - padding, dimension - padding,
         padding, dimension - padding,
         line_width=2,
         line_color="#FFFFFF"
        )
        matrix.add_line(
         padding, padding, dimension - padding, padding, line_width=2
        )
        matrix.add_line(
         padding, padding, padding, dimension - padding, line_width=2
        )
        matrix.add_line(
         dimension - padding, padding,  padding, dimension - padding, line_width=2
        )

        # Add secondary structure
        matrix.add_rectangle(
         bar_left, bar_top, bar_width, hypoteneuse,
         rotation=((dimension / 2) + 5, (dimension / 2) + 5, 45),
         line_width=0,
         fill_color=omnicanvas.hsl_to_rgb(chain_color, 100, 50)
        )
        for helix in self.alpha_helices:
            start = (len(residues) - self.get_all_residue_ids().index(
             helix.residues[-1].residue_id
            )) - 1
            end = (len(residues) - self.get_all_residue_ids().index(
             helix.residues[0].residue_id
            )) - 1
            matrix.add_rectangle(
             bar_left - 1, bar_top + (diagonal_chunk * start),
             bar_width + 2, diagonal_chunk * ((end - start) + 1),
             fill_color=omnicanvas.hsl_to_rgb(helix_color, 100, 50),
             line_width=0,
             rotation=((dimension / 2) + 5, (dimension / 2) + 5, 45)
            )
        for strand in self.beta_strands:
            start = (len(residues) - self.get_all_residue_ids().index(
             strand.residues[-1].residue_id
            )) - 1
            end = (len(residues) - self.get_all_residue_ids().index(
             strand.residues[0].residue_id
            )) - 1
            matrix.add_rectangle(
             bar_left - 1, bar_top + (diagonal_chunk * start),
             bar_width + 2, diagonal_chunk * (end - start),
             fill_color=omnicanvas.hsl_to_rgb(strand_color, 100, 50),
             line_width=0,
             rotation=((dimension / 2) + 5, (dimension / 2) + 5, 45)
            )

        # Add legend
        legend_dimension = plot_width * 0.4
        legend_left = padding + (dimension / 2)
        legend_top = dimension - (padding + legend_dimension)
        scale_width = legend_dimension * 0.8
        scale_height = legend_dimension * 0.1
        scale_left = legend_left + (0.1 * legend_dimension)
        scale_right = legend_left + (0.9 * legend_dimension)
        scale_top = legend_top + (0.2 * legend_dimension)
        scale_bottom = legend_top + (0.3 * legend_dimension)
        scale_label_y = legend_top + (0.15 * legend_dimension)
        number_label_y = legend_top + (0.38 * legend_dimension)
        helix_top = legend_top + (0.5 * legend_dimension)
        helix_bottom = legend_top + (0.6 * legend_dimension)
        helix_width = legend_dimension * 0.4
        helix_left = scale_left
        helix_right = helix_left + helix_width
        strand_top = legend_top + (0.7 * legend_dimension)
        strand_bottom = legend_top + (0.8 * legend_dimension)
        x_pixels = range(math.floor(scale_left), math.ceil(scale_right))

        matrix.add_text(
         scale_left + (scale_width / 2),
         scale_label_y,
         "Distance (&#8491;ngstroms)",
         font_size=int(scale_width / 9)
        )
        for x_pixel in x_pixels:
            color = 0
            fraction = (x_pixel - scale_left) / (scale_right - scale_left)
            if far_color >= close_color:
                distance_from_start = fraction * (far_color - close_color)
                color = int(close_color + distance_from_start)
            else:
                distance_from_start = fraction * (close_color - far_color)
                color = int(close_color - distance_from_start)
            matrix.add_rectangle(
             x_pixel - 1, scale_top, 2, scale_bottom - scale_top,
             fill_color=omnicanvas.hsl_to_rgb(color, 100, 50),
             line_width=0
            )
        matrix.add_text(
         scale_left, number_label_y, "0",
         font_size=int(scale_width / 10)
        )
        matrix.add_text(
         scale_right, number_label_y, "%i+" % cutoff,
         font_size=int(scale_width / 10)
        )
        matrix.add_rectangle(
         helix_left, ((helix_bottom + helix_top) / 2) - ((bar_width / 2) + 0),
         helix_width, bar_width,
         fill_color=omnicanvas.hsl_to_rgb(chain_color, 100, 50),
         line_width=0
        )
        matrix.add_rectangle(
         helix_left + (0.1 * helix_width),
         ((helix_bottom + helix_top) / 2) - ((bar_width / 2) + 1),
         helix_width * 0.8,
         bar_width + 2,
         fill_color=omnicanvas.hsl_to_rgb(helix_color, 100, 50),
         line_width=0
        )
        matrix.add_rectangle(
         helix_left, ((strand_bottom + strand_top) / 2) - ((bar_width / 2) + 0),
         helix_width, bar_width,
         fill_color=omnicanvas.hsl_to_rgb(chain_color, 100, 50),
         line_width=0
        )
        matrix.add_rectangle(
         helix_left + (0.1 * helix_width),
         ((strand_bottom + strand_top) / 2) - ((bar_width / 2) + 1),
         helix_width * 0.8,
         bar_width + 2,
         fill_color=omnicanvas.hsl_to_rgb(strand_color, 100, 50),
         line_width=0
        )
        matrix.add_text(
         helix_right + (legend_dimension * 0.1),
         ((helix_bottom + helix_top) / 2),
         "&#945;-helix",
         font_size=int(scale_width / 10),
         horizontal_align="right"
        )
        matrix.add_text(
         helix_right + (legend_dimension * 0.1),
         ((strand_bottom + strand_top) / 2),
         "&#946;-helix",
         font_size=int(scale_width / 10),
         horizontal_align="right"
        )

        return matrix



class PdbSite(ResiduicStructure):
    """Base class: :py:class:`ResiduicStructure`

    Represents PDB binding sites - the residue clusters that mediate ligand
    binding.

    :param site_id: The site's ID.
    :param residues: The residues in this chain.

    .. py:attribute:: site_id:

         A site's ID.

    .. py:attribute:: ligand:

         The :py:class:`PdbSmallMolecule` that binds to this site.

    .. py:attribute:: model:

         The :py:class:`PdbModel` that the site belongs to."""

    def __init__(self, site_id, *residues):
        if not isinstance(site_id, str):
            raise TypeError("'%s' is not a valid site_id" % str(site_id))
        self.site_id = site_id

        ResiduicStructure.__init__(self, *residues)
        self.ligand = None
        self.model = None


    def __repr__(self):
        return "<Site %s (%i residues)>" % (self.site_id, len(self.residues))



class PdbAlphaHelix(ResiduicSequence):
    """Base class: :py:class:`ResiduicSequence`

    Represents alpha helices.

    :param str helix_id: The helix's ID.
    :param residues: The residues in this helix.
    :param str helix_class: The classification of the helix.
    :param str comment: Any comment associated with this helix.
    :raises: :class:`.BrokenHelixError` if the residues given are from different\
    chains.

    .. py:attribute:: helix_id:

         The helix's unique ID.

    .. py:attribute:: helix_class:

         The helix's classification, such as 'Right-handed pi' or\
         'polyproline'.

    .. py:attribute:: comment:

        Any comment associated with the helix."""

    def __init__(self, helix_id, *residues, helix_class=None, comment=None):
        if not isinstance(helix_id, str):
            raise TypeError("'%s' is not a valid helix_id" % str(helix_id))
        self.helix_id = helix_id

        if helix_class is not None and not isinstance(helix_class, str):
            raise TypeError("'%s' is not a valid helix_class" % str(helix_class))
        self.helix_class = helix_class

        if comment is not None and not isinstance(comment, str):
            raise TypeError("'%s' is not a valid comment" % str(comment))
        self.comment = comment

        ResiduicSequence.__init__(self, *residues)
        if len(set([res.chain for res in self.residues])) > 1:
            raise BrokenHelixError(
             "Cannot make an alpha helix with residues from different chains"
            )
        if self.get_chain():
            self.get_chain().alpha_helices.add(self)


    def __repr__(self):
        return "<AlphaHelix %s (%i residues)>" % (self.helix_id, len(self.residues))


    def get_chain(self):
        """Returns the chain that the helix exists on.

        :rtype: :py:class:`PdbChain` or ``None``"""

        return list(self.residues)[0].chain



class PdbBetaStrand(ResiduicSequence):
    """Base class: :py:class:`ResiduicSequence`

    Represents beta strands.

    :param str strand_id: The strand's ID.
    :param residues: The residues in this strand.
    :param int sense: The sense of the strand with respect to the prior\
    strand.
    :raises: :class:`.BrokenStrandError` if the residues given are from different\
    chains.

    .. py:attribute:: strand_id:

         The strand's unique ID.

    .. py:attribute:: sense:

         The sense of the strand with respect to the prior strand. -1 for\
         anti-parallel, 1 for parallel, and 0 if it is the first strand.."""

    def __init__(self, strand_id, sense, *residues):
        if not isinstance(strand_id, int):
            raise TypeError("'%s' is not a valid strand_id" % str(strand_id))
        self.strand_id = strand_id

        if not isinstance(sense, int):
            raise TypeError("'%s' is not a valid sense value" % str(sense))
        if not (-1 <= sense <= 1):
            raise ValueError("sense can only be -1, 0 or 1 - not %i" % sense)
        self.sense = sense

        ResiduicSequence.__init__(self, *residues)
        if len(set([res.chain for res in self.residues])) > 1:
            raise BrokenStrandError(
             "Cannot make a beta strand with residues from different chains"
            )
        if self.get_chain():
            self.get_chain().beta_strands.add(self)


    def __repr__(self):
        return "<BetaStrand %i (%i residues)>" % (self.strand_id, len(self.residues))


    def get_chain(self):
        """Returns the chain that the strand exists on.

        :rtype: :py:class:`PdbChain` or ``None``"""
        return list(self.residues)[0].chain



class PdbModel(AtomicStructure):
    """Base class: :py:class:`AtomicStructure`

    Represents the structural environment in which the other structures exist.

    .. py:attribute:: small_molecules:

         The ``set`` of :py:class:`PdbSmallMolecule` objects in this model.

    .. py:attribute:: chains:

         The ``set`` of :py:class:`PdbChain` objects in this model."""

    def __init__(self):
        self.small_molecules = set()
        self.chains = set()
        self.sites = set()


    def __repr__(self):
        return "<Model (%i atoms)>" % len(self.atoms)


    def __getattr__(self, attribute):
        if attribute == "atoms":
            atoms = set()
            for small_molecule in self.small_molecules:
                atoms.update(small_molecule.atoms)
            for chain in self.chains:
                atoms.update(chain.atoms)
            for site in self.sites:
                atoms.update(site.atoms)
            return atoms
        else:
            return self.__getattribute__(attribute)


    def add_small_molecule(self, small_molecule):
        """Validates and adds a :py:class:`PdbSmallMolecule` to the model.

        :param PdbSmallMolecule small_molecule: the molecule to add.
        :raises TypeError: if something other than a\
        :py:class:`PdbSmallMolecule` is added.
        :raises: :class:`.DuplicateSmallMoleculeError` if there is already a\
        :py:class:`PdbSmallMolecule` of the same ID."""

        if not(isinstance(small_molecule, PdbSmallMolecule)):
            raise TypeError("Can only add SmallMolecules with add_small_molecule()")
        existing_small_molecule_ids = [mol.molecule_id for mol in self.small_molecules]
        if small_molecule.molecule_id in existing_small_molecule_ids:
            raise DuplicateSmallMoleculeError(
             "Cannot have two small molecules with ID %s" % small_molecule.molecule_id
            )
        self.small_molecules.add(small_molecule)
        small_molecule.model = self


    def get_small_molecule_by_id(self, molecule_id):
        """Retrurns the first small molecule that matches a given molecule ID.

        :param str molecule_id: The molecule ID to search by.
        :rtype: :py:class:`PdbSmallMolecule` or ``None``
        :raises TypeError: if the name given isn't an string."""

        if not isinstance(molecule_id, str):
            raise TypeError("Can only search small molecule IDs by str")
        for small_molecule in self.small_molecules:
            if small_molecule.molecule_id == molecule_id:
                return small_molecule


    def get_small_molecule_by_name(self, molecule_name):
        """Retruns all the small molecules that matches a given molecule name.

        :param str molecule_name: The molecule name to search by.
        :rtype: ``set`` of :py:class:`PdbSmallMolecule` objects
        :raises TypeError: if the name given isn't a string."""

        if not isinstance(molecule_name, str):
            raise TypeError("Can only search small molecule names by string")
        for small_molecule in self.small_molecules:
            if small_molecule.molecule_name == molecule_name:
                return small_molecule


    def get_small_molecules_by_name(self, molecule_name):
        """Retruns all the small molecules that matches a given residue name.

        :param str molecule_name: The molecule name to search by.
        :rtype: ``set`` of :py:class:`PdbSmallMolecule` objects
        :raises TypeError: if the name given isn't a string."""

        if not isinstance(molecule_name, str):
            raise TypeError("Can only search small molecule names by string")
        return set([
         mol for mol in self.small_molecules if mol.molecule_name == molecule_name
        ])


    def add_chain(self, chain):
        """Validates and adds a :py:class:`PdbChain` to the model.

        :param PdbChain chain: the chain to add.
        :raises TypeError: if something other than a\
        :py:class:`PdbChain` is added.
        :raises: :class:`.DuplicateChainError` if there is already a\
        :py:class:`PdbChain` of the same ID."""

        if not(isinstance(chain, PdbChain)):
            raise TypeError("Can only add Chain with add_chain()")
        existing_chain_ids = [c.chain_id for c in self.chains]
        if chain.chain_id in existing_chain_ids:
            raise DuplicateChainError(
             "Cannot have two chains with ID %s" % chain.chain_id
            )
        self.chains.add(chain)
        chain.model = self


    def get_chain_by_id(self, chain_id):
        """Retrurns the first chain that matches a given chain ID.

        :param str chain_id: The chain ID to search by.
        :rtype: :py:class:`PdbChain` or ``None``
        :raises TypeError: if the name given isn't an string."""

        if not isinstance(chain_id, str):
            raise TypeError("Can only search chain IDs by str")
        for chain in self.chains:
            if chain.chain_id == chain_id:
                return chain


    def add_site(self, site):
        """Validates and adds a :py:class:`PdbSite` to the model.

        :param PdbSite site: the site to add.
        :raises TypeError: if something other than a\
        :py:class:`PdbSite` is added.
        :raises: :class:`.DuplicateSiteError` if there is already a\
        :py:class:`PdbSite` of the same ID."""

        if not(isinstance(site, PdbSite)):
            raise TypeError("Can only add Site with add_site()")
        existing_site_ids = [s.site_id for s in self.sites]
        if site.site_id in existing_site_ids:
            raise DuplicateSiteError(
             "Cannot have two sites with ID %s" % site.site_id
            )
        self.sites.add(site)
        site.model = self


    def get_site_by_id(self, site_id):
        """Retrurns the first site that matches a given site ID.

        :param str site_id: The site ID to search by.
        :rtype: :py:class:`PdbSite` or ``None``
        :raises TypeError: if the name given isn't an string."""

        if not isinstance(site_id, str):
            raise TypeError("Can only search site IDs by str")
        for site in self.sites:
            if site.site_id == site_id:
                return site



def _residue_id_to_int(residue_id):
    int_component = int(
     "".join([char for char in residue_id if char in "-0123456789"])
    )
    char_component = residue_id[-1] if residue_id[-1] not in "-0123456789" else ""
    int_component *= 100
    return int_component + (ord(char_component) - 64 if char_component else 0)


def _residue_id_is_greater_than_residue_id(residue_id1, residue_id2):
    residues = []
    for residue_id in (residue_id1, residue_id2):
        chain_id = residue_id[0] if residue_id[0].isalpha() else ""
        number = int("".join([char for char in residue_id if char.isnumeric() or char == "-"]))
        insert = ord(residue_id[-1]) if residue_id[-1].isalpha() else 0
        residues.append((chain_id, number, insert))
    if residues[0][1] > residues[1][1]:
        return True
    elif residues[0][1] == residues[1][1] and residues[0][2] > residues[1][2]:
        return True
    else:
        return False


PERIODIC_TABLE = {
 "H": 1.0079, "HE": 4.0026, "LI": 6.941, "BE": 9.0122, "B": 10.811,
  "C": 12.0107, "N": 14.0067, "O": 15.9994, "F": 18.9984, "NE": 20.1797,
   "NA": 22.9897, "MG": 24.305, "AL": 26.9815, "SI": 28.0855, "P": 30.9738,
    "S": 32.065, "CL": 35.453, "K": 39.0983, "AR": 39.948, "CA": 40.078,
     "SC": 44.9559, "TI": 47.867, "V": 50.9415, "CR": 51.9961, "MN": 54.938,
      "FE": 55.845, "NI": 58.6934, "CO": 58.9332, "CU": 63.546, "ZN": 65.39,
       "GA": 69.723, "GE": 72.64, "AS": 74.9216, "SE": 78.96, "BR": 79.904,
        "KR": 83.8, "RB": 85.4678, "SR": 87.62, "Y": 88.9059, "ZR": 91.224,
         "NB": 92.9064, "MO": 95.94, "TC": 98, "RU": 101.07, "RH": 102.9055,
          "PD": 106.42, "AG": 107.8682, "CD": 112.411, "IN": 114.818,
           "SN": 118.71, "SB": 121.76, "I": 126.9045, "TE": 127.6,
            "XE": 131.293, "CS": 132.9055, "BA": 137.327, "LA": 138.9055,
             "CE": 140.116, "PR": 140.9077, "ND": 144.24, "PM": 145,
              "SM": 150.36, "EU": 151.964, "GD": 157.25, "TB": 158.9253,
               "DY": 162.5, "HO": 164.9303, "ER": 167.259, "TM": 168.9342,
                "YB": 173.04, "LU": 174.967, "HF": 178.49, "TA": 180.9479,
                 "W": 183.84, "RE": 186.207, "OS": 190.23, "IR": 192.217,
                  "PT": 195.078, "AU": 196.9665, "HG": 200.59, "TL": 204.3833,
                   "PB": 207.2, "BI": 208.9804, "PO": 209, "AT": 210, "RN": 222,
                    "FR": 223, "RA": 226, "AC": 227, "PA": 231.0359,
                     "TH": 232.0381, "NP": 237, "U": 238.0289, "AM": 243,
                      "PU": 244, "CM": 247, "BK": 247, "CF": 251, "ES": 252,
                       "FM": 257, "MD": 258, "NO": 259, "RF": 261, "LR": 262,
                        "DB": 262, "BH": 264, "SG": 266, "MT": 268, "RG": 272,
                         "HS": 277
}

RESIDUES = {
 "GLY": "G", "ALA": "A", "LEU": "L", "MET": "M", "PHE": "F",
 "TRP": "W", "LYS": "K", "GLN": "Q", "GLU": "E", "SER": "S",
 "PRO": "P", "VAL": "V", "ILE": "I", "CYS": "C", "TYR": "Y",
 "HIS": "H", "ARG": "R", "ASN": "N", "ASP": "D", "THR": "T"
}
