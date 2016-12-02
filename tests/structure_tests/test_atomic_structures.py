from collections import Counter
import warnings
from unittest import TestCase
import unittest.mock
from molecupy.structures import AtomicStructure, Atom, GhostAtom, BindSite, Residue
from molecupy import NoAtomsError, DuplicateAtomsError

class AtomicStructureTest(TestCase):

    def setUp(self):
        self.atoms = [unittest.mock.Mock(Atom) for _ in range(10)]
        self.ghost_atoms = [unittest.mock.Mock(GhostAtom) for _ in range(10)]
        self.all_atoms = self.atoms + self.ghost_atoms
        for index, atom in enumerate(self.all_atoms):
            atom.mass.return_value = index + 1
            atom.atom_id.return_value = index + 1



class AtomicStructureCreationTests(AtomicStructureTest):

    def test_can_create_atomic_structure_with_atoms(self):
        atomic_structure = AtomicStructure(*self.atoms)
        self.assertEqual(atomic_structure._atoms, set(self.atoms))


    def test_can_create_atomic_structure_with_ghost_atoms(self):
        atomic_structure = AtomicStructure(*self.ghost_atoms)
        self.assertEqual(atomic_structure._atoms, set(self.ghost_atoms))


    def test_can_create_atomic_structure_with_mixed_atoms(self):
        atomic_structure = AtomicStructure(*self.all_atoms)
        self.assertEqual(
         atomic_structure._atoms,
         set(self.ghost_atoms + self.atoms)
        )


    def test_can_only_create_atomic_structure_with_atoms(self):
        with self.assertRaises(TypeError):
            AtomicStructure("Atom1", "Atom2")


    def test_cannot_make_empty_atomic_structure(self):
        with self.assertRaises(NoAtomsError):
            AtomicStructure()


    def test_cannot_have_duplicate_atom_ids_in_atomic_structure(self):
        self.all_atoms[-1].atom_id.return_value = 19
        with self.assertRaises(DuplicateAtomsError):
            AtomicStructure(*self.all_atoms)


    def test_repr(self):
        atomic_structure = AtomicStructure(*self.atoms)
        self.assertEqual(str(atomic_structure), "<AtomicStructure (10 atoms)>")


    def test_atomic_structure_is_iterable_of_atoms(self):
        atomic_structure = AtomicStructure(*self.atoms)
        atoms_from_loop = set()
        for atom in atomic_structure:
            atoms_from_loop.add(atom)
        self.assertEqual(set(self.atoms), atoms_from_loop)



class AtomicStructurePropertyTests(AtomicStructureTest):

    def test_can_get_localised_atoms(self):
        atomic_structure = AtomicStructure(*self.all_atoms)
        self.assertEqual(
         atomic_structure.atoms(),
         set(self.atoms)
        )


    def test_can_get_ghost_atoms(self):
        atomic_structure = AtomicStructure(*self.all_atoms)
        self.assertEqual(
         atomic_structure.atoms(atom_type="ghost"),
         set(self.ghost_atoms)
        )


    def test_can_get_all_atoms(self):
        atomic_structure = AtomicStructure(*self.all_atoms)
        self.assertEqual(
         atomic_structure.atoms(atom_type="all"),
         set(self.all_atoms)
        )


    def test_atom_retrieval_must_be_of_valid_type(self):
        atomic_structure = AtomicStructure(*self.all_atoms)
        with self.assertRaises(TypeError):
            atomic_structure.atoms(atom_type=1)
        with self.assertRaises(ValueError):
            atomic_structure.atoms(atom_type="xyz")


    def test_atomic_structure_atoms_is_read_only(self):
        atom21 = unittest.mock.Mock(Atom)
        atomic_structure = AtomicStructure(*self.all_atoms)
        self.assertEqual(len(atomic_structure.atoms(atom_type="all")), 20)
        atomic_structure.atoms(atom_type="all").add(atom21)
        self.assertEqual(len(atomic_structure.atoms(atom_type="all")), 20)


    def test_can_add_atom(self):
        atom21 = unittest.mock.Mock(Atom)
        atomic_structure = AtomicStructure(*self.all_atoms)
        atomic_structure.add_atom(atom21)
        self.assertEqual(len(atomic_structure.atoms(atom_type="all")), 21)
        self.assertIn(atom21, atomic_structure.atoms(atom_type="all"))


    def test_can_only_add_atoms(self):
        atomic_structure = AtomicStructure(*self.all_atoms)
        with self.assertRaises(TypeError):
            atomic_structure.add_atom("atom21")


    def test_can_only_add_atoms_if_id_is_unique(self):
        atomic_structure = AtomicStructure(*self.all_atoms[:-1])
        self.all_atoms[-1].atom_id.return_value = 19
        with self.assertRaises(DuplicateAtomsError):
            atomic_structure.add_atom(self.all_atoms[-1])


    def test_can_remove_atoms(self):
        atomic_structure = AtomicStructure(*self.all_atoms)
        atomic_structure.remove_atom(self.atoms[5])
        self.assertEqual(len(atomic_structure.atoms(atom_type="all")), 19)
        self.assertNotIn(self.atoms[5], atomic_structure.atoms(atom_type="all"))
        atomic_structure.remove_atom(self.ghost_atoms[5])
        self.assertEqual(len(atomic_structure.atoms(atom_type="all")), 18)
        self.assertNotIn(
         self.ghost_atoms[5],
         atomic_structure.atoms(atom_type="all")
        )



class AtomicStructureMassTests(AtomicStructureTest):

    def test_can_get_atomic_structure_mass_all(self):
        atomic_structure = AtomicStructure(*self.all_atoms)
        self.assertEqual(atomic_structure.mass(atom_type="all"), 210)


    def test_can_get_atomic_structure_mass_localised(self):
        atomic_structure = AtomicStructure(*self.all_atoms)
        self.assertEqual(atomic_structure.mass(), 55)


    def test_can_get_atomic_structure_mass_ghost(self):
        atomic_structure = AtomicStructure(*self.all_atoms)
        self.assertEqual(atomic_structure.mass(atom_type="ghost"), 155)



class AtomicStructureFormulaTests(AtomicStructureTest):

    def setUp(self):
        AtomicStructureTest.setUp(self)
        for index, element in enumerate(["H", "H", "C", "C", "N", "C", "N", "F", "H", "H"]):
            self.atoms[index].element.return_value = element
        for index, element in enumerate(["H", "P", "C", "N", "N", "C", "N", "F", "H", "H"]):
            self.ghost_atoms[index].element.return_value = element
        self.atomic_structure = AtomicStructure(*self.all_atoms)


    def test_can_get_all_atom_formula_with_hydrogens(self):
        self.assertEqual(
         self.atomic_structure.formula(atom_type="all", include_hydrogens=True),
         Counter({"H": 7, "C": 5, "N": 5, "F": 2, "P": 1})
        )


    def test_can_get_all_atom_formula_without_hydrogens(self):
        self.assertEqual(
         self.atomic_structure.formula(atom_type="all", include_hydrogens=False),
         Counter({"C": 5, "N": 5, "F": 2, "P": 1})
        )


    def test_can_get_localised_atom_formula_with_hydrogens(self):
        self.assertEqual(
         self.atomic_structure.formula(include_hydrogens=True),
         Counter({"H": 4, "C": 3, "N": 2, "F": 1})
        )


    def test_can_get_localised_atom_formula_without_hydrogens(self):
        self.assertEqual(
         self.atomic_structure.formula(include_hydrogens=False),
         Counter({"C": 3, "N": 2, "F": 1})
        )


    def test_can_get_ghost_atom_formula_with_hydrogens(self):
        self.assertEqual(
         self.atomic_structure.formula(atom_type="ghost", include_hydrogens=True),
         Counter({"H": 3, "C": 2, "N": 3, "F": 1, "P": 1})
        )


    def test_can_get_ghost_atom_formula_without_hydrogens(self):
        self.assertEqual(
         self.atomic_structure.formula(atom_type="ghost", include_hydrogens=False),
         Counter({"C": 2, "N": 3, "F": 1, "P": 1})
        )


    def test_default_formula_is_localised_atoms_no_hydrogens(self):
        self.assertEqual(
         self.atomic_structure.formula(),
         self.atomic_structure.formula(include_hydrogens=False)
        )



class AtomRetrievalTests(AtomicStructureTest):

    def setUp(self):
        AtomicStructureTest.setUp(self)
        for index, atom in enumerate(self.all_atoms):
            atom.element.return_value = chr((index % 5) + 65)
            atom.atom_name.return_value = atom.element() + "X"
        self.atomic_structure = AtomicStructure(*self.all_atoms)


    def test_can_get_atom_by_id(self):
        self.assertIs(self.atomic_structure.get_atom_by_id(1), self.all_atoms[0])
        self.assertIs(self.atomic_structure.get_atom_by_id(5), self.all_atoms[4])
        self.assertIs(self.atomic_structure.get_atom_by_id(13), self.all_atoms[12])
        self.assertIs(self.atomic_structure.get_atom_by_id(24), None)


    def test_can_only_search_by_numeric_id(self):
        with self.assertRaises(TypeError):
            self.atomic_structure.get_atom_by_id("98")


    def test_can_get_all_atoms_by_element(self):
        self.assertEqual(
         self.atomic_structure.get_atoms_by_element("B", atom_type="all"),
         set((self.all_atoms[1], self.all_atoms[6], self.all_atoms[11], self.all_atoms[16]))
        )


    def test_can_get_localised_atoms_by_element(self):
        self.assertEqual(
         self.atomic_structure.get_atoms_by_element("B"),
         set((self.all_atoms[1], self.all_atoms[6]))
        )


    def test_can_get_ghost_atoms_by_element(self):
        self.assertEqual(
         self.atomic_structure.get_atoms_by_element("B", atom_type="ghost"),
         set((self.all_atoms[11], self.all_atoms[16]))
        )


    def test_failed_element_search_returns_empty_set(self):
        self.assertEqual(
         self.atomic_structure.get_atoms_by_element("F"),
         set()
        )


    def test_can_get_single_all_atom_by_element(self):
        self.assertIn(
         self.atomic_structure.get_atom_by_element("B", atom_type="all"),
         (self.all_atoms[1], self.all_atoms[6], self.all_atoms[11], self.all_atoms[16])
        )


    def test_can_get_single_localised_atom_by_element(self):
        self.assertIn(
         self.atomic_structure.get_atom_by_element("B"),
         (self.all_atoms[1], self.all_atoms[6])
        )


    def test_can_get_single_ghost_atom_by_element(self):
        self.assertIn(
         self.atomic_structure.get_atom_by_element("B", atom_type="ghost"),
         (self.all_atoms[11], self.all_atoms[16])
        )


    def test_failed_single_element_search_returns_none(self):
        self.assertEqual(
         self.atomic_structure.get_atom_by_element("F"),
         None
        )


    def test_can_only_search_by_string_element(self):
        with self.assertRaises(TypeError):
            self.atomic_structure.get_atom_by_element(None)
        with self.assertRaises(TypeError):
            self.atomic_structure.get_atoms_by_element(None)


    def test_can_get_all_atoms_by_name(self):
        self.assertEqual(
         self.atomic_structure.get_atoms_by_name("BX", atom_type="all"),
         set((self.all_atoms[1], self.all_atoms[6], self.all_atoms[11], self.all_atoms[16]))
        )


    def test_can_get_localised_atoms_by_name(self):
        self.assertEqual(
         self.atomic_structure.get_atoms_by_name("BX"),
         set((self.all_atoms[1], self.all_atoms[6]))
        )


    def test_can_get_ghost_atoms_by_name(self):
        self.assertEqual(
         self.atomic_structure.get_atoms_by_name("BX", atom_type="ghost"),
         set((self.all_atoms[11], self.all_atoms[16]))
        )


    def test_failed_name_search_returns_empty_set(self):
        self.assertEqual(
         self.atomic_structure.get_atoms_by_name("FX"),
         set()
        )


    def test_can_get_single_all_atom_by_name(self):
        self.assertIn(
         self.atomic_structure.get_atom_by_name("BX", atom_type="all"),
         (self.all_atoms[1], self.all_atoms[6], self.all_atoms[11], self.all_atoms[16])
        )


    def test_can_get_single_localised_atom_by_name(self):
        self.assertIn(
         self.atomic_structure.get_atom_by_name("BX"),
         (self.all_atoms[1], self.all_atoms[6])
        )


    def test_can_get_single_generic_atom_by_name(self):
        self.assertIn(
         self.atomic_structure.get_atom_by_name("BX", atom_type="ghost"),
         (self.all_atoms[11], self.all_atoms[16])
        )


    def test_failed_single_name_search_returns_none(self):
        self.assertEqual(
         self.atomic_structure.get_atom_by_name("FX"),
         None
        )


    def test_can_only_search_by_string_name(self):
        with self.assertRaises(TypeError):
            self.atomic_structure.get_atom_by_name(None)
        with self.assertRaises(TypeError):
            self.atomic_structure.get_atoms_by_name(None)



class AtomicStructureContactsTests(AtomicStructureTest):

    def setUp(self):
        x_values = [10.0, 20.0, 30.0, 40.0, 50.0, 80.0, 80.0, 80.0, 80.0, 80.0]
        y_values = [30.0, 30.0, 30.0, 30.0, 30.0, 10.0, 20.0, 30.0, 40.0, 50.0]
        self.atoms = [
         Atom(x_values[i], y_values[i], 10.0, "C", i + 1, "CX") for i in range(10)
        ]

    def test_can_get_contacts_between_atomic_structures(self):
        structure1 = AtomicStructure(*self.atoms[:5])
        structure2 = AtomicStructure(*self.atoms[5:])
        self.assertEqual(
         structure1.contacts_with(structure2, distance=30),
         set([frozenset([self.atoms[4], self.atoms[7]])])
        )
        self.assertEqual(
         structure1.contacts_with(structure2, distance=35),
         set([
          frozenset([self.atoms[4], self.atoms[7]]),
          frozenset([self.atoms[4], self.atoms[6]]),
          frozenset([self.atoms[4], self.atoms[8]])
         ])
        )
        self.assertEqual(
         structure1.contacts_with(structure2, distance=40),
         set([
          frozenset([self.atoms[4], self.atoms[7]]),
          frozenset([self.atoms[4], self.atoms[6]]),
          frozenset([self.atoms[4], self.atoms[8]]),
          frozenset([self.atoms[3], self.atoms[7]]),
          frozenset([self.atoms[4], self.atoms[5]]),
          frozenset([self.atoms[4], self.atoms[9]]),
         ])
        )


    def test_external_contacts_can_ignore_hydrogens(self):
        self.atoms[4].element("H")
        self.atoms[5].element("H")
        self.atoms[9].element("H")
        structure1 = AtomicStructure(*self.atoms[:5])
        structure2 = AtomicStructure(*self.atoms[5:])
        self.assertEqual(
         structure1.contacts_with(structure2, distance=30, include_hydrogens=False),
         set()
        )
        self.assertEqual(
         structure1.contacts_with(structure2, distance=35, include_hydrogens=False),
         set()
        )
        self.assertEqual(
         structure1.contacts_with(structure2, distance=40, include_hydrogens=False),
         set([frozenset([self.atoms[3], self.atoms[7]])])
        )


    def test_can_get_external_contacts_when_one_structure_is_part_of_the_other(self):
        structure1 = AtomicStructure(*self.atoms[:5])
        structure2 = AtomicStructure(*self.atoms)
        self.assertEqual(
         structure1.contacts_with(structure2, distance=30),
         set([frozenset([self.atoms[4], self.atoms[7]])])
        )
        self.assertEqual(
         structure1.contacts_with(structure2, distance=35),
         set([
          frozenset([self.atoms[4], self.atoms[7]]),
          frozenset([self.atoms[4], self.atoms[6]]),
          frozenset([self.atoms[4], self.atoms[8]])
         ])
        )
        self.assertEqual(
         structure1.contacts_with(structure2, distance=40),
         set([
          frozenset([self.atoms[4], self.atoms[7]]),
          frozenset([self.atoms[4], self.atoms[6]]),
          frozenset([self.atoms[4], self.atoms[8]]),
          frozenset([self.atoms[3], self.atoms[7]]),
          frozenset([self.atoms[4], self.atoms[5]]),
          frozenset([self.atoms[4], self.atoms[9]]),
         ])
        )


    def test_can_get_internal_contacts(self):
        warnings.simplefilter("ignore")
        for index, atom in enumerate(self.atoms[:-1]):
            atom.bond_to(self.atoms[index + 1])
        structure = AtomicStructure(*self.atoms)
        self.assertEqual(
         structure.internal_contacts(distance=31),
         set([
          frozenset([self.atoms[0], self.atoms[3]]),
          frozenset([self.atoms[1], self.atoms[4]]),
          frozenset([self.atoms[5], self.atoms[8]]),
          frozenset([self.atoms[6], self.atoms[9]]),
          frozenset([self.atoms[4], self.atoms[7]])
         ])
        )


    def test_internal_contacts_can_ignore_hydrogens(self):
        warnings.simplefilter("ignore")
        for index, atom in enumerate(self.atoms[:-1]):
            atom.bond_to(self.atoms[index + 1])
        self.atoms[4].element("H")
        self.atoms[5].element("H")
        self.atoms[9].element("H")
        structure = AtomicStructure(*self.atoms)
        self.assertEqual(
         structure.internal_contacts(distance=31, include_hydrogens=False),
         set([
          frozenset([self.atoms[0], self.atoms[3]])
         ])
        )



class BindSiteGenerationTests(AtomicStructureTest):

    def setUp(self):
        AtomicStructureTest.setUp(self)
        self.atomic_structure = AtomicStructure(*self.all_atoms[-13:])
        self.residue1 = unittest.mock.Mock(spec=Residue)
        self.residue2 = unittest.mock.Mock(spec=Residue)
        molecule1 = unittest.mock.Mock()
        self.atoms[-1].local_atoms.return_value = set(self.atoms[:2])
        self.atoms[-2].local_atoms.return_value = set(self.atoms[2:4])
        self.atoms[-3].local_atoms.return_value = set(self.atoms[4:6])
        self.atoms[0].molecule.return_value = self.residue1
        self.atoms[1].molecule.return_value = None
        self.atoms[2].molecule.return_value = None
        self.atoms[3].molecule.return_value = self.residue2
        self.atoms[3].element("H")
        self.atoms[4].molecule.return_value = None
        self.atoms[5].molecule.return_value = molecule1


    def test_can_find_bind_site(self):
        site = self.atomic_structure.predict_bind_site()
        self.assertIsInstance(site, BindSite)
        self.assertEqual(site.residues(), set([self.residue1, self.residue2]))


    def test_can_exclude_hydrogens(self):
        self.atoms[-2].local_atoms.return_value = set([self.atoms[2]])
        site = self.atomic_structure.predict_bind_site(include_hydrogens=False)
        self.assertEqual(site.residues(), set([self.residue1]))



class TransformationTest(TestCase):

    def setUp(self):
        self.atom1 = Atom(20.0, 20.0, 20.0, "X", 1, "X")
        self.atom2 = Atom(25.0, 20.0, 20.0, "X", 2, "X")
        self.atom3 = Atom(20.0, 20.0, 25.0, "X", 3, "X")
        self.atom4 = Atom(25.0, 20.0, 25.0, "X", 4, "X")
        self.atom5 = Atom(22.5, 25.0, 22.5, "X", 5, "X")
        self.atoms = [self.atom1, self.atom2, self.atom3, self.atom4, self.atom5]
        self.atom6 = GhostAtom("X", 6, "X")
        self.structure = AtomicStructure(
         self.atom1, self.atom2, self.atom3, self.atom4, self.atom5, self.atom6
        )



class TranslationTests(TransformationTest):

    def test_can_translate_structure_in_x_direction(self):
        self.structure.translate(10, 0, 0)
        self.assertEqual(
         (self.atom1.x(), self.atom1.y(), self.atom1.z()),
         (30, 20, 20)
        )
        self.assertEqual(
         (self.atom2.x(), self.atom2.y(), self.atom2.z()),
         (35, 20, 20)
        )
        self.assertEqual(
         (self.atom3.x(), self.atom3.y(), self.atom3.z()),
         (30, 20, 25)
        )
        self.assertEqual(
         (self.atom4.x(), self.atom4.y(), self.atom4.z()),
         (35, 20, 25)
        )
        self.assertEqual(
         (self.atom5.x(), self.atom5.y(), self.atom5.z()),
         (32.5, 25, 22.5)
        )
        self.structure.translate(-120, 0, 0)
        self.assertEqual(
         (self.atom1.x(), self.atom1.y(), self.atom1.z()),
         (-90, 20, 20)
        )
        self.assertEqual(
         (self.atom2.x(), self.atom2.y(), self.atom2.z()),
         (-85, 20, 20)
        )
        self.assertEqual(
         (self.atom3.x(), self.atom3.y(), self.atom3.z()),
         (-90, 20, 25)
        )
        self.assertEqual(
         (self.atom4.x(), self.atom4.y(), self.atom4.z()),
         (-85, 20, 25)
        )
        self.assertEqual(
         (self.atom5.x(), self.atom5.y(), self.atom5.z()),
         (-87.5, 25, 22.5)
        )


    def test_can_translate_structure_in_y_direction(self):
        self.structure.translate(0, 50, 0)
        self.assertEqual(
         (self.atom1.x(), self.atom1.y(), self.atom1.z()),
         (20, 70, 20)
        )
        self.assertEqual(
         (self.atom2.x(), self.atom2.y(), self.atom2.z()),
         (25, 70, 20)
        )
        self.assertEqual(
         (self.atom3.x(), self.atom3.y(), self.atom3.z()),
         (20, 70, 25)
        )
        self.assertEqual(
         (self.atom4.x(), self.atom4.y(), self.atom4.z()),
         (25, 70, 25)
        )
        self.assertEqual(
         (self.atom5.x(), self.atom5.y(), self.atom5.z()),
         (22.5, 75, 22.5)
        )
        self.structure.translate(0, -20, 0)
        self.assertEqual(
         (self.atom1.x(), self.atom1.y(), self.atom1.z()),
         (20, 50, 20)
        )
        self.assertEqual(
         (self.atom2.x(), self.atom2.y(), self.atom2.z()),
         (25, 50, 20)
        )
        self.assertEqual(
         (self.atom3.x(), self.atom3.y(), self.atom3.z()),
         (20, 50, 25)
        )
        self.assertEqual(
         (self.atom4.x(), self.atom4.y(), self.atom4.z()),
         (25, 50, 25)
        )
        self.assertEqual(
         (self.atom5.x(), self.atom5.y(), self.atom5.z()),
         (22.5, 55, 22.5)
        )


    def test_can_translate_structure_in_z_direction(self):
        self.structure.translate(0, 0, 100)
        self.assertEqual(
         (self.atom1.x(), self.atom1.y(), self.atom1.z()),
         (20, 20, 120)
        )
        self.assertEqual(
         (self.atom2.x(), self.atom2.y(), self.atom2.z()),
         (25, 20, 120)
        )
        self.assertEqual(
         (self.atom3.x(), self.atom3.y(), self.atom3.z()),
         (20, 20, 125)
        )
        self.assertEqual(
         (self.atom4.x(), self.atom4.y(), self.atom4.z()),
         (25, 20, 125)
        )
        self.assertEqual(
         (self.atom5.x(), self.atom5.y(), self.atom5.z()),
         (22.5, 25, 122.5)
        )
        self.structure.translate(0, 0, -1)
        self.assertEqual(
         (self.atom1.x(), self.atom1.y(), self.atom1.z()),
         (20, 20, 119)
        )
        self.assertEqual(
         (self.atom2.x(), self.atom2.y(), self.atom2.z()),
         (25, 20, 119)
        )
        self.assertEqual(
         (self.atom3.x(), self.atom3.y(), self.atom3.z()),
         (20, 20, 124)
        )
        self.assertEqual(
         (self.atom4.x(), self.atom4.y(), self.atom4.z()),
         (25, 20, 124)
        )
        self.assertEqual(
         (self.atom5.x(), self.atom5.y(), self.atom5.z()),
         (22.5, 25, 121.5)
        )


    def test_can_translate_structure_in_multiple_directions(self):
        self.structure.translate(19, -10432, 0.25)
        self.assertEqual(
         (self.atom1.x(), self.atom1.y(), self.atom1.z()),
         (39, -10412, 20.25)
        )
        self.assertEqual(
         (self.atom2.x(), self.atom2.y(), self.atom2.z()),
         (44, -10412, 20.25)
        )
        self.assertEqual(
         (self.atom3.x(), self.atom3.y(), self.atom3.z()),
         (39, -10412, 25.25)
        )
        self.assertEqual(
         (self.atom4.x(), self.atom4.y(), self.atom4.z()),
         (44, -10412, 25.25)
        )
        self.assertEqual(
         (self.atom5.x(), self.atom5.y(), self.atom5.z()),
         (41.5, -10407, 22.75)
        )



class RotationTests(TransformationTest):

    def test_can_rotate_around_x_axis_by_180(self):
        self.structure.rotate("x", 180)
        self.assertAlmostEqual(self.atom1.x(), 20, delta=0.005)
        self.assertAlmostEqual(self.atom1.y(), -20, delta=0.005)
        self.assertAlmostEqual(self.atom1.z(), -20, delta=0.005)
        self.assertAlmostEqual(self.atom2.x(), 25, delta=0.005)
        self.assertAlmostEqual(self.atom2.y(), -20, delta=0.005)
        self.assertAlmostEqual(self.atom2.z(), -20, delta=0.005)
        self.assertAlmostEqual(self.atom3.x(), 20, delta=0.005)
        self.assertAlmostEqual(self.atom3.y(), -20, delta=0.005)
        self.assertAlmostEqual(self.atom3.z(), -25, delta=0.005)
        self.assertAlmostEqual(self.atom4.x(), 25, delta=0.005)
        self.assertAlmostEqual(self.atom4.y(), -20, delta=0.005)
        self.assertAlmostEqual(self.atom4.z(), -25, delta=0.005)
        self.assertAlmostEqual(self.atom5.x(), 22.5, delta=0.005)
        self.assertAlmostEqual(self.atom5.y(), -25, delta=0.005)
        self.assertAlmostEqual(self.atom5.z(), -22.5, delta=0.005)


    def test_can_rotate_around_x_axis_by_90(self):
        self.structure.rotate("x", 90)
        self.assertAlmostEqual(self.atom1.x(), 20, delta=0.005)
        self.assertAlmostEqual(self.atom1.y(), -20, delta=0.005)
        self.assertAlmostEqual(self.atom1.z(), 20, delta=0.005)
        self.assertAlmostEqual(self.atom2.x(), 25, delta=0.005)
        self.assertAlmostEqual(self.atom2.y(), -20, delta=0.005)
        self.assertAlmostEqual(self.atom2.z(), 20, delta=0.005)
        self.assertAlmostEqual(self.atom3.x(), 20, delta=0.005)
        self.assertAlmostEqual(self.atom3.y(), -25, delta=0.005)
        self.assertAlmostEqual(self.atom3.z(), 20, delta=0.005)
        self.assertAlmostEqual(self.atom4.x(), 25, delta=0.005)
        self.assertAlmostEqual(self.atom4.y(), -25, delta=0.005)
        self.assertAlmostEqual(self.atom4.z(), 20, delta=0.005)
        self.assertAlmostEqual(self.atom5.x(), 22.5, delta=0.005)
        self.assertAlmostEqual(self.atom5.y(), -22.5, delta=0.005)
        self.assertAlmostEqual(self.atom5.z(), 25, delta=0.005)


    def test_can_rotate_around_x_axis_by_45(self):
        self.structure.rotate("x", 45)
        self.assertAlmostEqual(self.atom1.x(), 20, delta=0.005)
        self.assertAlmostEqual(self.atom1.y(), 0, delta=0.005)
        self.assertAlmostEqual(self.atom1.z(), 28.284271, delta=0.005)
        self.assertAlmostEqual(self.atom2.x(), 25, delta=0.005)
        self.assertAlmostEqual(self.atom2.y(), 0, delta=0.005)
        self.assertAlmostEqual(self.atom2.z(), 28.284271, delta=0.005)
        self.assertAlmostEqual(self.atom3.x(), 20, delta=0.005)
        self.assertAlmostEqual(self.atom3.y(), -3.535534, delta=0.005)
        self.assertAlmostEqual(self.atom3.z(), 31.819805, delta=0.005)
        self.assertAlmostEqual(self.atom4.x(), 25, delta=0.005)
        self.assertAlmostEqual(self.atom4.y(), -3.535534, delta=0.005)
        self.assertAlmostEqual(self.atom4.z(), 31.819805, delta=0.005)
        self.assertAlmostEqual(self.atom5.x(), 22.5, delta=0.005)
        self.assertAlmostEqual(self.atom5.y(), 1.767767, delta=0.005)
        self.assertAlmostEqual(self.atom5.z(), 33.587572, delta=0.005)


    def test_can_rotate_around_x_axis_outside_of_360_range(self):
        xyz_values_1 = []
        xyz_values_2 = []
        xyz_values_3 = []
        xyz_values = [xyz_values_1, xyz_values_2, xyz_values_3]
        for n in range(3):
            atom1 = Atom(20.0, 20.0, 20.0, "X", 1, "X")
            atom2 = Atom(25.0, 20.0, 20.0, "X", 2, "X")
            atom3 = Atom(20.0, 20.0, 25.0, "X", 3, "X")
            atom4 = Atom(25.0, 20.0, 25.0, "X", 4, "X")
            atom5 = Atom(22.5, 25.0, 22.5, "X", 5, "X")
            atoms = [atom1, atom2, atom3, atom4, atom5]
            structure = AtomicStructure(*atoms)
            structure.rotate("x", -90 + (360 * n))
            for atom in atoms:
                xyz_values[n].append(atom.x())
                xyz_values[n].append(atom.y())
                xyz_values[n].append(atom.z())
        for a in range(15):
            self.assertAlmostEqual(xyz_values_1[a], xyz_values_2[a], delta=0.005)
            self.assertAlmostEqual(xyz_values_3[a], xyz_values_2[a], delta=0.005)


    def test_can_rotate_around_y_axis_by_180(self):
        self.structure.rotate("y", 180)
        self.assertAlmostEqual(self.atom1.x(), -20, delta=0.005)
        self.assertAlmostEqual(self.atom1.y(), 20, delta=0.005)
        self.assertAlmostEqual(self.atom1.z(), -20, delta=0.005)
        self.assertAlmostEqual(self.atom2.x(), -25, delta=0.005)
        self.assertAlmostEqual(self.atom2.y(), 20, delta=0.005)
        self.assertAlmostEqual(self.atom2.z(), -20, delta=0.005)
        self.assertAlmostEqual(self.atom3.x(), -20, delta=0.005)
        self.assertAlmostEqual(self.atom3.y(), 20, delta=0.005)
        self.assertAlmostEqual(self.atom3.z(), -25, delta=0.005)
        self.assertAlmostEqual(self.atom4.x(), -25, delta=0.005)
        self.assertAlmostEqual(self.atom4.y(), 20, delta=0.005)
        self.assertAlmostEqual(self.atom4.z(), -25, delta=0.005)
        self.assertAlmostEqual(self.atom5.x(), -22.5, delta=0.005)
        self.assertAlmostEqual(self.atom5.y(), 25, delta=0.005)
        self.assertAlmostEqual(self.atom5.z(), -22.5, delta=0.005)


    def test_can_rotate_around_y_axis_by_90(self):
        self.structure.rotate("y", 90)
        self.assertAlmostEqual(self.atom1.x(), 20, delta=0.005)
        self.assertAlmostEqual(self.atom1.y(), 20, delta=0.005)
        self.assertAlmostEqual(self.atom1.z(), -20, delta=0.005)
        self.assertAlmostEqual(self.atom2.x(), 20, delta=0.005)
        self.assertAlmostEqual(self.atom2.y(), 20, delta=0.005)
        self.assertAlmostEqual(self.atom2.z(), -25, delta=0.005)
        self.assertAlmostEqual(self.atom3.x(), 25, delta=0.005)
        self.assertAlmostEqual(self.atom3.y(), 20, delta=0.005)
        self.assertAlmostEqual(self.atom3.z(), -20, delta=0.005)
        self.assertAlmostEqual(self.atom4.x(), 25, delta=0.005)
        self.assertAlmostEqual(self.atom4.y(), 20, delta=0.005)
        self.assertAlmostEqual(self.atom4.z(), -25, delta=0.005)
        self.assertAlmostEqual(self.atom5.x(), 22.5, delta=0.005)
        self.assertAlmostEqual(self.atom5.y(), 25, delta=0.005)
        self.assertAlmostEqual(self.atom5.z(), -22.5, delta=0.005)


    def test_can_rotate_around_y_axis_by_45(self):
        self.structure.rotate("y", 45)
        self.assertAlmostEqual(self.atom1.x(), 28.284271, delta=0.005)
        self.assertAlmostEqual(self.atom1.y(), 20, delta=0.005)
        self.assertAlmostEqual(self.atom1.z(), 0, delta=0.005)
        self.assertAlmostEqual(self.atom2.x(), 31.819805, delta=0.005)
        self.assertAlmostEqual(self.atom2.y(), 20, delta=0.005)
        self.assertAlmostEqual(self.atom2.z(), -3.535534, delta=0.005)
        self.assertAlmostEqual(self.atom3.x(), 31.819805, delta=0.005)
        self.assertAlmostEqual(self.atom3.y(), 20, delta=0.005)
        self.assertAlmostEqual(self.atom3.z(), 3.535534, delta=0.005)
        self.assertAlmostEqual(self.atom4.x(), 35.355339, delta=0.005)
        self.assertAlmostEqual(self.atom4.y(), 20, delta=0.005)
        self.assertAlmostEqual(self.atom4.z(), 0, delta=0.005)
        self.assertAlmostEqual(self.atom5.x(), 31.819805, delta=0.005)
        self.assertAlmostEqual(self.atom5.y(), 25, delta=0.005)
        self.assertAlmostEqual(self.atom5.z(), 0, delta=0.005)


    def test_can_rotate_around_y_axis_outside_of_360_range(self):
        xyz_values_1 = []
        xyz_values_2 = []
        xyz_values_3 = []
        xyz_values = [xyz_values_1, xyz_values_2, xyz_values_3]
        for n in range(3):
            atom1 = Atom(20.0, 20.0, 20.0, "X", 1, "X")
            atom2 = Atom(25.0, 20.0, 20.0, "X", 2, "X")
            atom3 = Atom(20.0, 20.0, 25.0, "X", 3, "X")
            atom4 = Atom(25.0, 20.0, 25.0, "X", 4, "X")
            atom5 = Atom(22.5, 25.0, 22.5, "X", 5, "X")
            atoms = [atom1, atom2, atom3, atom4, atom5]
            structure = AtomicStructure(*atoms)
            structure.rotate("y", -90 + (360 * n))
            for atom in atoms:
                xyz_values[n].append(atom.x())
                xyz_values[n].append(atom.y())
                xyz_values[n].append(atom.z())
        for a in range(15):
            self.assertAlmostEqual(xyz_values_1[a], xyz_values_2[a], delta=0.005)
            self.assertAlmostEqual(xyz_values_3[a], xyz_values_2[a], delta=0.005)



    def test_can_rotate_around_z_axis_by_180(self):
        self.structure.rotate("z", 180)
        self.assertAlmostEqual(self.atom1.x(), -20, delta=0.005)
        self.assertAlmostEqual(self.atom1.y(), -20, delta=0.005)
        self.assertAlmostEqual(self.atom1.z(), 20, delta=0.005)
        self.assertAlmostEqual(self.atom2.x(), -25, delta=0.005)
        self.assertAlmostEqual(self.atom2.y(), -20, delta=0.005)
        self.assertAlmostEqual(self.atom2.z(), 20, delta=0.005)
        self.assertAlmostEqual(self.atom3.x(), -20, delta=0.005)
        self.assertAlmostEqual(self.atom3.y(), -20, delta=0.005)
        self.assertAlmostEqual(self.atom3.z(), 25, delta=0.005)
        self.assertAlmostEqual(self.atom4.x(), -25, delta=0.005)
        self.assertAlmostEqual(self.atom4.y(), -20, delta=0.005)
        self.assertAlmostEqual(self.atom4.z(), 25, delta=0.005)
        self.assertAlmostEqual(self.atom5.x(), -22.5, delta=0.005)
        self.assertAlmostEqual(self.atom5.y(), -25, delta=0.005)
        self.assertAlmostEqual(self.atom5.z(), 22.5, delta=0.005)


    def test_can_rotate_around_z_axis_by_90(self):
        self.structure.rotate("z", 90)
        self.assertAlmostEqual(self.atom1.x(), -20, delta=0.005)
        self.assertAlmostEqual(self.atom1.y(), 20, delta=0.005)
        self.assertAlmostEqual(self.atom1.z(), 20, delta=0.005)
        self.assertAlmostEqual(self.atom2.x(), -20, delta=0.005)
        self.assertAlmostEqual(self.atom2.y(), 25, delta=0.005)
        self.assertAlmostEqual(self.atom2.z(), 20, delta=0.005)
        self.assertAlmostEqual(self.atom3.x(), -20, delta=0.005)
        self.assertAlmostEqual(self.atom3.y(), 20, delta=0.005)
        self.assertAlmostEqual(self.atom3.z(), 25, delta=0.005)
        self.assertAlmostEqual(self.atom4.x(), -20, delta=0.005)
        self.assertAlmostEqual(self.atom4.y(), 25, delta=0.005)
        self.assertAlmostEqual(self.atom4.z(), 25, delta=0.005)
        self.assertAlmostEqual(self.atom5.x(), -25, delta=0.005)
        self.assertAlmostEqual(self.atom5.y(), 22.5, delta=0.005)
        self.assertAlmostEqual(self.atom5.z(), 22.5, delta=0.005)


    def test_can_rotate_around_z_axis_by_45(self):
        self.structure.rotate("z", 45)
        self.assertAlmostEqual(self.atom1.x(), 0, delta=0.005)
        self.assertAlmostEqual(self.atom1.y(), 28.284271, delta=0.005)
        self.assertAlmostEqual(self.atom1.z(), 20, delta=0.005)
        self.assertAlmostEqual(self.atom2.x(), 3.535534, delta=0.005)
        self.assertAlmostEqual(self.atom2.y(), 31.819805, delta=0.005)
        self.assertAlmostEqual(self.atom2.z(), 20, delta=0.005)
        self.assertAlmostEqual(self.atom3.x(), 0, delta=0.005)
        self.assertAlmostEqual(self.atom3.y(), 28.284271, delta=0.005)
        self.assertAlmostEqual(self.atom3.z(), 25, delta=0.005)
        self.assertAlmostEqual(self.atom4.x(), 3.535534, delta=0.005)
        self.assertAlmostEqual(self.atom4.y(), 31.819805, delta=0.005)
        self.assertAlmostEqual(self.atom4.z(), 25, delta=0.005)
        self.assertAlmostEqual(self.atom5.x(), -1.767767, delta=0.005)
        self.assertAlmostEqual(self.atom5.y(), 33.587572, delta=0.005)
        self.assertAlmostEqual(self.atom5.z(), 22.5, delta=0.005)


    def test_can_rotate_around_z_axis_outside_of_360_range(self):
        xyz_values_1 = []
        xyz_values_2 = []
        xyz_values_3 = []
        xyz_values = [xyz_values_1, xyz_values_2, xyz_values_3]
        for n in range(3):
            atom1 = Atom(20.0, 20.0, 20.0, "X", 1, "X")
            atom2 = Atom(25.0, 20.0, 20.0, "X", 2, "X")
            atom3 = Atom(20.0, 20.0, 25.0, "X", 3, "X")
            atom4 = Atom(25.0, 20.0, 25.0, "X", 4, "X")
            atom5 = Atom(22.5, 25.0, 22.5, "X", 5, "X")
            atoms = [atom1, atom2, atom3, atom4, atom5]
            structure = AtomicStructure(*atoms)
            structure.rotate("z", -90 + (360 * n))
            for atom in atoms:
                xyz_values[n].append(atom.x())
                xyz_values[n].append(atom.y())
                xyz_values[n].append(atom.z())
        for a in range(15):
            self.assertAlmostEqual(xyz_values_1[a], xyz_values_2[a], delta=0.005)
            self.assertAlmostEqual(xyz_values_3[a], xyz_values_2[a], delta=0.005)



    def test_diagonal_structure_rotation(self):
        atoms = [Atom(float(i), float(i), float(i), "X", i, "X") for i in range(1, 6)]
        structure = AtomicStructure(*atoms)
        structure.rotate("x", -45)
        z = 1.414214
        self.assertAlmostEqual(atoms[0].x(), 1, delta=0.005)
        self.assertAlmostEqual(atoms[0].y(), z, delta=0.005)
        self.assertAlmostEqual(atoms[0].z(), 0, delta=0.005)
        self.assertAlmostEqual(atoms[1].x(), 2, delta=0.005)
        self.assertAlmostEqual(atoms[1].y(), 2 * z, delta=0.005)
        self.assertAlmostEqual(atoms[1].z(), 0, delta=0.005)
        self.assertAlmostEqual(atoms[2].x(), 3, delta=0.005)
        self.assertAlmostEqual(atoms[2].y(), 3 * z, delta=0.005)
        self.assertAlmostEqual(atoms[2].z(), 0, delta=0.005)
        self.assertAlmostEqual(atoms[3].x(), 4, delta=0.005)
        self.assertAlmostEqual(atoms[3].y(), 4 * z, delta=0.005)
        self.assertAlmostEqual(atoms[3].z(), 0, delta=0.005)
        self.assertAlmostEqual(atoms[4].x(), 5, delta=0.005)
        self.assertAlmostEqual(atoms[4].y(), 5 * z, delta=0.005)
        self.assertAlmostEqual(atoms[4].z(), 0, delta=0.005)


    def test_rotation_axis_must_be_string(self):
        with self.assertRaises(TypeError):
            self.structure.rotate(100, 90)


    def test_rotation_axis_must_be_valid(self):
        with self.assertRaises(ValueError):
            self.structure.rotate("zz", 90)
        with self.assertRaises(ValueError):
            self.structure.rotate("a", 90)


    def test_rotation_angle_must_be_numeric(self):
        with self.assertRaises(TypeError):
            self.structure.rotate("x", "90")
        self.structure.rotate("x", 90.7)



class MassCenterTests(AtomicStructureTest):

    def test_can_get_mass_centre_between_two_identical_atoms(self):
        atom1 = Atom(1.0, 1.0, 1.0, "H", 1, "H")
        atom2 = Atom(2.0, 2.0, 2.0, "H", 2, "H")
        structure = AtomicStructure(atom1, atom2)
        self.assertEqual(
         structure.center_of_mass(),
         (1.5, 1.5, 1.5)
        )


    def test_can_get_mass_center_between_multiple_identical_atoms_on_plane(self):
        atom1 = Atom(1.0, 1.0, 0.0, "H", 1, "H")
        atom2 = Atom(7.0, 1.0, 0.0, "H", 2, "H")
        atom3 = Atom(1.0, 7.0, 0.0, "H", 3, "H")
        atom4 = Atom(7.0, 7.0, 0.0, "H", 4, "H")
        structure = AtomicStructure(atom1, atom2, atom3, atom4)
        self.assertEqual(
         structure.center_of_mass(),
         (4.0, 4.0, 0.0)
        )


    def test_can_get_mass_center_between_multiple_identical_atoms_in_3d(self):
        atom1 = Atom(1.0, 1.0, 0.0, "H", 1, "H")
        atom2 = Atom(7.0, 1.0, 0.0, "H", 2, "H")
        atom3 = Atom(1.0, 7.0, 0.0, "H", 3, "H")
        atom4 = Atom(7.0, 7.0, 0.0, "H", 4, "H")
        atom5 = Atom(1.0, 1.0, 10.0, "H", 5, "H")
        atom6 = Atom(7.0, 1.0, 10.0, "H", 6, "H")
        atom7 = Atom(1.0, 7.0, 10.0, "H", 7, "H")
        atom8 = Atom(7.0, 7.0, 10.0, "H", 8, "H")
        structure = AtomicStructure(
         atom1, atom2, atom3, atom4, atom5, atom6, atom7, atom8
        )
        self.assertEqual(
         structure.center_of_mass(),
         (4.0, 4.0, 5.0)
        )


    def test_can_get_mass_centre_between_two_asymmetrical_atoms(self):
        atom1 = unittest.mock.Mock(Atom)
        atom1.x.return_value = 6.0
        atom1.y.return_value = 0.0
        atom1.z.return_value = 0.0
        atom1.mass.return_value = 1.0
        atom2 = unittest.mock.Mock(Atom)
        atom2.x.return_value = 10.0
        atom2.y.return_value = 0.0
        atom2.z.return_value = 0.0
        atom2.mass.return_value = 3.0

        structure = AtomicStructure(atom1, atom2)
        self.assertEqual(
         structure.center_of_mass(),
         (9.0, 0.0, 0.0)
        )


    def test_can_get_mass_centre_between_four_asymmetrical_atoms(self):
        atom1 = unittest.mock.Mock(Atom)
        atom1.x.return_value = 0.0
        atom1.y.return_value = 0.0
        atom1.z.return_value = -1.0
        atom1.mass.return_value = 9.0
        atom2 = unittest.mock.Mock(Atom)
        atom2.x.return_value = 0.0
        atom2.y.return_value = 0.0
        atom2.z.return_value = 2.0
        atom2.mass.return_value = 1.0
        atom3 = unittest.mock.Mock(Atom)
        atom3.x.return_value = 0.0
        atom3.y.return_value = 0.0
        atom3.z.return_value = 6.0
        atom3.mass.return_value = 5.0
        atom4 = unittest.mock.Mock(Atom)
        atom4.x.return_value = 0.0
        atom4.y.return_value = 0.0
        atom4.z.return_value = 8.0
        atom4.mass.return_value = 2.0
        structure = AtomicStructure(atom1, atom2, atom3, atom4)
        self.assertAlmostEqual(structure.center_of_mass()[-1], 2.29, delta=0.005)


    def test_can_get_mass_centre_between_asymmetrical_atoms_in_plane(self):
        atom1 = unittest.mock.Mock(Atom)
        atom1.x.return_value = -1.0
        atom1.y.return_value = 1.0
        atom1.z.return_value = 0.0
        atom1.mass.return_value = 5.0
        atom2 = unittest.mock.Mock(Atom)
        atom2.x.return_value = 4.0
        atom2.y.return_value = 3.0
        atom2.z.return_value = 0.0
        atom2.mass.return_value = 3.0
        atom3 = unittest.mock.Mock(Atom)
        atom3.x.return_value = 8.0
        atom3.y.return_value = 7.0
        atom3.z.return_value = 0.0
        atom3.mass.return_value = 9.0
        structure = AtomicStructure(atom1, atom2, atom3)
        self.assertAlmostEqual(structure.center_of_mass()[0], 4.65, delta=0.005)
        self.assertAlmostEqual(structure.center_of_mass()[1], 4.53, delta=0.005)



class RadiusOfGyrationTests(AtomicStructureTest):

    def test_can_get_Rg_with_two_identical_atoms(self):
        atom1 = Atom(0.0, 0.0, 0.0, "H", 1, "H")
        atom2 = Atom(10.0, 0.0, 0.0, "H", 2, "H")
        structure = AtomicStructure(atom1, atom2)
        self.assertEqual(
         structure.radius_of_gyration(),
         5
        )
