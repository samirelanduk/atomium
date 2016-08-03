"""Contains classes for simple structures made of atoms."""

from collections import Counter
from .atoms import Atom, PdbAtom
from ..exceptions import NoAtomsError, MultipleResidueConnectionError, DuplicateAtomsError

class AtomicStructure:
    """The base class for all structures which are composed of atoms.

    :param atoms: A sequence of :py:class:`Atom` objects."""

    def __init__(self, *atoms):
        if len(atoms) == 0:
            raise NoAtomsError("Cannot make an AtomicStructure with no atoms")
        for atom in atoms:
            if not isinstance(atom, Atom):
                raise TypeError(
                 "Can only make AtomicStructures with Atoms, not '%s'" % str(atom)
                )
        atom_ids = [atom.atom_id() for atom in atoms]
        if len(atom_ids) != len(set(atom_ids)):
            raise DuplicateAtomsError(
             "Cannot make atomic structure with duplicate atom IDs"
            )
        self._atoms = set(atoms)


    def __repr__(self):
        return "<%s (%i atoms)>" % (self.__class__.__name__, len(self._atoms))


    def atoms(self, atom_type="all"):
        """Returns the atoms in this structure as a ``set``.

        :param str atom_type: The kind of atom to return. ``"all"`` will\
        return all atoms, ``"pdb"`` just PdbAtoms and ``"generic"`` will just\
        return generic non-PDB atoms.
        :rtype: ``set``"""

        if not isinstance(atom_type, str):
            raise TypeError("atom_type must be str, not '%s'" % str(atom_type))
        if atom_type == "pdb":
            return set(
             [atom for atom in self._atoms if isinstance(atom, PdbAtom)]
            )
        elif atom_type == "generic":
            return set(
             [atom for atom in self._atoms if not isinstance(atom, PdbAtom)]
            )
        elif atom_type == "all":
            return set(self._atoms)
        else:
            raise ValueError("'%s' is not a valid atom_type" % atom_type)


    def add_atom(self, atom):
        """Adds an atom to the structure.

        :param Atom atom: The atom to add."""

        if not isinstance(atom, Atom):
            raise TypeError(
             "Can only add Atoms to AtomicStructures, not '%s'" % str(atom)
            )
        if atom.atom_id() in [atom.atom_id() for atom in self.atoms()]:
            raise DuplicateAtomsError(
             "Cannot add atom with ID %i to %s as there is already an atom with that ID" % (
              atom.atom_id(), self
             )
            )
        self._atoms.add(atom)


    def remove_atom(self, atom):
        """Removes an atom from the structure.

        :param Atom atom: The atom to add."""

        self._atoms.remove(atom)


    def mass(self, atom_type="all"):
        """Returns the mass of the structure by summing the mass of all its
        atoms.

        :param str atom_type: The kind of atom to use. ``"all"`` will\
        use all atoms, ``"pdb"`` just PdbAtoms and ``"generic"`` will just\
        use generic non-PDB atoms.
        :rtype: float"""

        return sum([atom.mass() for atom in self.atoms(atom_type)])


    def formula(self, atom_type="all", include_hydrogens=False):
        """Retrurns the formula (count of each atom) of the structure.

        :param str atom_type: The kind of atom to use. ``"all"`` will\
        use all atoms, ``"pdb"`` just PdbAtoms and ``"generic"`` will just\
        use generic non-PDB atoms.
        :param bool include_hydrogens: determines whether hydrogen atoms should\
        be included.
        :rtype: ``Counter``"""

        return Counter([
         atom.element() for atom in self.atoms(atom_type=atom_type)
          if include_hydrogens or atom.element().upper() != "H"
        ])


    def contacts_with(self, other_atomic_structure, distance=4, include_hydrogens=True):
        """Returns the set of all 'contacts' with another atomic structure,
        where a contact is defined as any atom-atom pair with an inter-atomic
        distance less than or equal to some number of Angstroms.

        If the other atomic structure has atoms which are also in this atomic
        structure, those atoms will not be counted as part of the other
        structure.

        :param AtomicStructure other_structure: The other atomic\
        structure to compare to.
        :param distance: The distance to use (default is 4).
        :param bool include_hydrogens: determines whether hydrogen atoms should\
        be included.
        :rtype: ``set`` of ``frozenset`` contacts."""

        contacts = set()
        our_atoms = list(self.atoms(atom_type="pdb"))
        their_atoms = [
         a for a in list(other_atomic_structure.atoms(atom_type="pdb"))
          if a not in our_atoms
        ]
        if not include_hydrogens:
            our_atoms = [
             atom for atom in our_atoms if atom.element().upper() != "H"
            ]
            their_atoms = [
             atom for atom in their_atoms if atom.element().upper() != "H"
            ]
        for atom in our_atoms:
            for other_atom in their_atoms:
                if atom.distance_to(other_atom) <= distance:
                    contacts.add(frozenset((atom, other_atom)))
        return contacts


    def internal_contacts(self, distance=4, include_hydrogens=True):
        """Returns the set of all atomic contacts within the atoms of an atomic
        structure, where a contact is defined as any atom-atom pair with an
        inter-atomic distance less than or equal to four Angstroms.

        Contacts between atoms covalently bonded to each other will be ignored,
        as will contacts between atoms separated by just two covalent bonds.

        :param distance: The distance to use (default is 4).
        :param bool include_hydrogens: determines whether hydrogen atoms should\
        be included.
        :rtype: ``set`` of ``frozenset`` contacts."""

        contacts = set()
        atoms = list(self.atoms(atom_type="pdb"))
        if not include_hydrogens:
            atoms = [atom for atom in atoms if atom.element().upper() != "H"]
        for index, atom in enumerate(atoms[:-1]):
            too_close_atoms = set((atom,))
            for bonded_atom in atom.bonded_atoms():
                too_close_atoms.add(bonded_atom)
            second_tier_atoms = set()
            for close_atom in too_close_atoms:
                for bonded_atom in close_atom.bonded_atoms():
                    second_tier_atoms.add(bonded_atom)
            too_close_atoms.update(second_tier_atoms)
            for other_atom in atoms[index + 1:]:
                if other_atom not in too_close_atoms and atom.distance_to(other_atom) <= distance:
                    contacts.add(frozenset((atom, other_atom)))
        return contacts


    def predict_bind_site(self, distance=5, include_hydrogens=True):
        """Attempts to predict the residues that might make up the atomic
        structure's binding site by using atomic distances.

        :param distance: The distance to use (default is 5s).
        :param bool include_hydrogens: determines whether hydrogen atoms should\
        be included.
        :rtype: :py:class:`.BindSite` or ``None``"""

        from .chains import BindSite
        nearby_atoms = set()
        for atom in [atom for atom in self.atoms(atom_type="pdb")
         if (include_hydrogens or atom.element().upper() != "H")]:
            nearby_atoms.update(atom.local_atoms(
             distance=distance, include_hydrogens=include_hydrogens
            ))
        nearby_atoms = [atom for atom in nearby_atoms if atom not in self.atoms()]
        residues = set()
        for atom in nearby_atoms:
            if isinstance(atom.molecule(), Residue):
                residues.add(atom.molecule())
        return BindSite("CALC", *list(residues))


    def get_atom_by_id(self, atom_id):
        """Retrurns the first atom that matches a given atom ID.

        :param int atom_id: The atom ID to search by.
        :rtype: :py:class:`.Atom` or ``None``"""

        if not isinstance(atom_id, int):
            raise TypeError("Atom ID search must be by int, not '%s'" % str(atom_id))
        for atom in self.atoms():
            if atom.atom_id() == atom_id:
                return atom


    def get_atoms_by_element(self, element, atom_type="all"):
        """Retruns all the atoms a given element.

        :param str element: The element to search by.
        :param str atom_type: The kind of atom to use. ``"all"`` will\
        use all atoms, ``"pdb"`` just PdbAtoms and ``"generic"`` will just\
        use generic non-PDB atoms.
        :rtype: ``set`` of :py:class:`.Atom` objects."""

        if not isinstance(element, str):
            raise TypeError("Atom element search must be by str, not '%s'" % str(element))
        return set([
         atom for atom in self.atoms(atom_type=atom_type) if atom.element() == element
        ])


    def get_atom_by_element(self, element, atom_type="all"):
        """Retrurns the first atom that matches a given element.

        :param str element: The element to search by.
        :param str atom_type: The kind of atom to use. ``"all"`` will\
        use all atoms, ``"pdb"`` just PdbAtoms and ``"generic"`` will just\
        use generic non-PDB atoms.
        :rtype: :py:class:`.Atom` or ``None``"""

        if not isinstance(element, str):
            raise TypeError("Atom element search must be by str, not '%s'" % str(element))
        for atom in self.atoms(atom_type=atom_type):
            if atom.element() == element:
                return atom


    def get_atoms_by_name(self, atom_name, atom_type="all"):
        """Retruns all the atoms a given name.

        :param str atom_name: The name to search by.
        :param str atom_type: The kind of atom to use. ``"all"`` will\
        use all atoms, ``"pdb"`` just PdbAtoms and ``"generic"`` will just\
        use generic non-PDB atoms.
        :rtype: ``set`` of :py:class:`.Atom` objects."""

        if not isinstance(atom_name, str):
            raise TypeError("Atom name search must be by str, not '%s'" % str(atom_name))
        return set([
         atom for atom in self.atoms(atom_type=atom_type) if atom.atom_name() == atom_name
        ])


    def get_atom_by_name(self, atom_name, atom_type="all"):
        """Retrurns the first atom that matches a given name.

        :param str atom_name: The name to search by.
        :param str atom_type: The kind of atom to use. ``"all"`` will\
        use all atoms, ``"pdb"`` just PdbAtoms and ``"generic"`` will just\
        use generic non-PDB atoms.
        :rtype: :py:class:`.Atom` or ``None``"""

        if not isinstance(atom_name, str):
            raise TypeError("Atom name search must be by str, not '%s'" % str(atom_name))
        for atom in self.atoms(atom_type=atom_type):
            if atom.atom_name() == atom_name:
                return atom



class SmallMolecule(AtomicStructure):
    """Base class: :py:class:`AtomicStructure`

    Represents the ligands, solvent molecules, and other non-polymeric
    molecules in a PDB structure.

    :param str molecule_id: The molecule's ID.
    :param str molecule_name: The molecule's name.
    :param atoms: The molecule's atoms."""

    def __init__(self, molecule_id, molecule_name, *atoms):
        if not isinstance(molecule_id, str):
            raise TypeError("'%s' is not a valid molecule_id" % str(molecule_id))
        self._molecule_id = molecule_id
        if not isinstance(molecule_name, str):
            raise TypeError("'%s' is not a valid molecule_name" % str(molecule_name))
        self._molecule_name = molecule_name
        self._bind_site = None
        self._model = None
        AtomicStructure.__init__(self, *atoms)
        for atom in self._atoms:
            atom._molecule = self


    def __repr__(self):
        return "<SmallMolecule %s (%s)>" % (self._molecule_id, self._molecule_name)


    def molecule_id(self):
        """Returns the molecule's ID.

        :rtype: ``str``"""

        return self._molecule_id


    def molecule_name(self, molecule_name=None):
        """Returns or sets the molecule's name.

        :param str name: If given, the molecule's name will be set to this.
        :rtype: ``str``"""

        if molecule_name is None:
            return self._molecule_name
        else:
            if not isinstance(molecule_name, str):
                raise TypeError(
                 "'%s' is not a valid molecule_name" % str(molecule_name)
                )
            self._molecule_name = molecule_name


    def bind_site(self, bind_site=None):
        """Returns or sets the molecule's :py:class:`.BindSite`.

        :param BindSite bind_site: If given, the atom's bindsite will be set to this.
        :rtype: ``BindSite``"""

        if bind_site is None:
            return self._bind_site
        else:
            from .chains import BindSite
            if not isinstance(bind_site, BindSite):
                raise TypeError(
                 "'%s' is not a valid bind_site" % str(bind_site)
                )
            self._bind_site = bind_site
            bind_site._ligand = self


    def model(self):
        """Returns the :py:class:`.Model` that the molecule inhabits.

        :rtype: ``Model``"""

        return self._model


    def add_atom(self, atom):
        AtomicStructure.add_atom(self, atom)
        atom._molecule = self


    def remove_atom(self, atom):
        AtomicStructure.remove_atom(self, atom)
        atom._molecule = None



class Residue(AtomicStructure):
    """Base class: :py:class:`AtomicStructure`

    A Residue on a chain.

    :param str residue_id: The residue's ID.
    :param str residue_name: The residue's name.
    :param atoms: The residue's atoms."""

    def __init__(self, residue_id, residue_name, *atoms):
        if not isinstance(residue_id, str):
            raise TypeError("'%s' is not a valid residue_id" % str(residue_id))
        self._residue_id = residue_id
        if not isinstance(residue_name, str):
            raise TypeError("'%s' is not a valid residue_name" % str(residue_name))
        self._residue_name = residue_name
        AtomicStructure.__init__(self, *atoms)
        for atom in self._atoms:
            atom._molecule = self
        self._downstream_residue = None
        self._upstream_residue = None
        self._chain = None


    def __repr__(self):
        return "<Residue %s (%s)>" % (self._residue_id, self._residue_name)


    def residue_id(self):
        """Returns the residue's ID.

        :rtype: ``str``"""

        return self._residue_id


    def residue_name(self, residue_name=None):
        """Returns or sets the residue's name.

        :param str name: If given, the residue's name will be set to this.
        :rtype: ``str``"""

        if residue_name is None:
            return self._residue_name
        else:
            if not isinstance(residue_name, str):
                raise TypeError(
                 "'%s' is not a valid residue_name" % str(residue_name)
                )
            self._residue_name = residue_name


    def add_atom(self, atom):
        AtomicStructure.add_atom(self, atom)
        atom._molecule = self


    def remove_atom(self, atom):
        AtomicStructure.remove_atom(self, atom)
        atom._molecule = None


    def chain(self):
        """Returns the :py:class:`.Chain` that the residue is within.

        :rtype: ``Chain``"""

        return self._chain


    def is_missing(self):
        """Returns ``True`` if the residue was not observed in the experiment
        (and is therefore made up entirely of atoms with no coordinates).

        :rtype: ``bool``"""

        return not bool(self.atoms(atom_type="pdb"))


    def downstream_residue(self):
        """Returns the residue connected to this residue's carboxy end.

        :rtype: ``Residue``"""

        return self._downstream_residue


    def upstream_residue(self):
        """Returns the residue connected to this residue's amino end.

        :rtype: ``Residue``"""

        return self._upstream_residue


    def connect_to(self, downstream_residue):
        """Connects this residue to a downstream residue.

        :param Residue downstream_residue: The other residue."""

        if not isinstance(downstream_residue, Residue):
            raise TypeError(
             "Can only connect Residues to other Residues, not '%s'" % str(downstream_residue)
            )
        if self._downstream_residue is not None:
            raise MultipleResidueConnectionError(
             "%s already has a downstream residue (%s) and cannot connect to %s"
              % (self, self._downstream_residue, downstream_residue)
             )
        if downstream_residue._upstream_residue is not None:
            raise MultipleResidueConnectionError(
             "%s already has an upstream residue (%s) and so %s cannot connect to it"
              % (downstream_residue, downstream_residue._upstream_residue, self)
             )
        self._downstream_residue = downstream_residue
        downstream_residue._upstream_residue = self


    def disconnect_from(self, other_residue):
        """Breaks a connection with another residue.

        :param Residue other_residue: The other residue."""

        if self._downstream_residue is other_residue:
            self._downstream_residue = None
            other_residue._upstream_residue = None
        if self._upstream_residue is other_residue:
            self._upstream_residue = None
            other_residue._downstream_residue = None


    def alpha_carbon(self):
        """Attempts to retrieve the alpha carbon of the residue.

        :rtype: ``Atom``"""

        atom = self.get_atom_by_name("CA")
        if not atom:
            atom = self.get_atom_by_element("C")
        if not atom:
            atom = list(self.atoms())[0]
        return atom
