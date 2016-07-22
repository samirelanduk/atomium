from collections import Counter
import warnings
from unittest import TestCase
import unittest.mock
from molecupy.structures import AtomicStructure, PdbAtom, Atom
from molecupy import NoAtomsError

class AtomicStructureTest(TestCase):

    def setUp(self):
        self.pdb_atoms = [unittest.mock.Mock(spec=PdbAtom) for _ in range(10)]
        self.generic_atoms = [unittest.mock.Mock(spec=Atom) for _ in range(10)]
        self.all_atoms = self.pdb_atoms + self.generic_atoms
        for index, atom in enumerate(self.all_atoms):
            atom.mass.return_value = index + 1



class AtomicStructureCreationTests(AtomicStructureTest):

    def test_can_create_atomic_structure_with_pdb_atoms(self):
        atomic_structure = AtomicStructure(*self.pdb_atoms)
        self.assertEqual(atomic_structure._atoms, set(self.pdb_atoms))


    def test_can_create_atomic_structure_with_generic_atoms(self):
        atomic_structure = AtomicStructure(*self.generic_atoms)
        self.assertEqual(atomic_structure._atoms, set(self.generic_atoms))


    def test_can_create_atomic_structure_with_mixed_atoms(self):
        atomic_structure = AtomicStructure(*self.all_atoms)
        self.assertEqual(
         atomic_structure._atoms,
         set(self.generic_atoms + self.pdb_atoms)
        )


    def test_can_only_create_atomic_structure_with_atoms(self):
        with self.assertRaises(TypeError):
            AtomicStructure("Atom1", "Atom2")


    def test_cannot_make_empty_atomic_structure(self):
        with self.assertRaises(NoAtomsError):
            AtomicStructure()


    def test_repr(self):
        atomic_structure = AtomicStructure(*self.pdb_atoms)
        self.assertEqual(str(atomic_structure), "<AtomicStructure (10 atoms)>")



class AtomicStructurePropertyTests(AtomicStructureTest):

    def test_can_get_pdb_atoms(self):
        atomic_structure = AtomicStructure(*self.all_atoms)
        self.assertEqual(
         atomic_structure.atoms(atom_type="pdb"),
         set(self.pdb_atoms)
        )


    def test_can_get_generic_atoms(self):
        atomic_structure = AtomicStructure(*self.all_atoms)
        self.assertEqual(
         atomic_structure.atoms(atom_type="generic"),
         set(self.generic_atoms)
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


    def test_default_atom_retrieval_is_all(self):
        atomic_structure = AtomicStructure(*self.all_atoms)
        self.assertEqual(
         atomic_structure.atoms(),
         set(self.all_atoms)
        )


    def test_atomic_structure_atoms_is_read_only(self):
        atom21 = unittest.mock.Mock(spec=PdbAtom)
        atomic_structure = AtomicStructure(*self.all_atoms)
        self.assertEqual(len(atomic_structure.atoms()), 20)
        atomic_structure.atoms(atom_type="all").add(atom21)
        self.assertEqual(len(atomic_structure.atoms()), 20)


    def test_can_add_atom(self):
        atom21 = unittest.mock.Mock(spec=PdbAtom)
        atomic_structure = AtomicStructure(*self.all_atoms)
        atomic_structure.add_atom(atom21)
        self.assertEqual(len(atomic_structure.atoms()), 21)
        self.assertIn(atom21, atomic_structure.atoms())


    def test_can_only_add_atoms(self):
        atomic_structure = AtomicStructure(*self.all_atoms)
        with self.assertRaises(TypeError):
            atomic_structure.add_atom("atom21")


    def test_can_remove_atoms(self):
        atomic_structure = AtomicStructure(*self.all_atoms)
        atomic_structure.remove_atom(self.pdb_atoms[5])
        self.assertEqual(len(atomic_structure.atoms()), 19)
        self.assertNotIn(self.pdb_atoms[5], atomic_structure.atoms(atom_type="all"))
        atomic_structure.remove_atom(self.generic_atoms[5])
        self.assertEqual(len(atomic_structure.atoms()), 18)
        self.assertNotIn(self.generic_atoms[5], atomic_structure.atoms(atom_type="all"))



class AtomicStructureMassTests(AtomicStructureTest):

    def test_can_get_atomic_structure_mass_all(self):
        atomic_structure = AtomicStructure(*self.all_atoms)
        self.assertEqual(atomic_structure.mass(atom_type="all"), 210)


    def test_can_get_atomic_structure_mass_pdb(self):
        atomic_structure = AtomicStructure(*self.all_atoms)
        self.assertEqual(atomic_structure.mass(atom_type="pdb"), 55)


    def test_can_get_atomic_structure_mass_generic(self):
        atomic_structure = AtomicStructure(*self.all_atoms)
        self.assertEqual(atomic_structure.mass(atom_type="generic"), 155)


    def test_default_atomic_mass_is_all(self):
        atomic_structure = AtomicStructure(*self.all_atoms)
        self.assertEqual(atomic_structure.mass(), 210)



class AtomicStructureFormulaTests(AtomicStructureTest):

    def setUp(self):
        AtomicStructureTest.setUp(self)
        for index, element in enumerate(["H", "H", "C", "C", "N", "C", "N", "F", "H", "H"]):
            self.pdb_atoms[index].element.return_value = element
        for index, element in enumerate(["H", "P", "C", "N", "N", "C", "N", "F", "H", "H"]):
            self.generic_atoms[index].element.return_value = element
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


    def test_can_get_pdb_atom_formula_with_hydrogens(self):
        self.assertEqual(
         self.atomic_structure.formula(atom_type="pdb", include_hydrogens=True),
         Counter({"H": 4, "C": 3, "N": 2, "F": 1})
        )


    def test_can_get_pdb_atom_formula_without_hydrogens(self):
        self.assertEqual(
         self.atomic_structure.formula(atom_type="pdb", include_hydrogens=False),
         Counter({"C": 3, "N": 2, "F": 1})
        )


    def test_can_get_generic_atom_formula_with_hydrogens(self):
        self.assertEqual(
         self.atomic_structure.formula(atom_type="generic", include_hydrogens=True),
         Counter({"H": 3, "C": 2, "N": 3, "F": 1, "P": 1})
        )


    def test_can_get_generic_atom_formula_without_hydrogens(self):
        self.assertEqual(
         self.atomic_structure.formula(atom_type="generic", include_hydrogens=False),
         Counter({"C": 2, "N": 3, "F": 1, "P": 1})
        )


    def test_default_formula_is_all_atoms_no_hydrogens(self):
        self.assertEqual(
         self.atomic_structure.formula(),
         self.atomic_structure.formula(atom_type="all", include_hydrogens=False)
        )



class AtomRetrievalTests(AtomicStructureTest):

    def setUp(self):
        AtomicStructureTest.setUp(self)
        for index, atom in enumerate(self.all_atoms):
            atom.atom_id.return_value = index + 1
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


    def test_can_get_pdb_atoms_by_element(self):
        self.assertEqual(
         self.atomic_structure.get_atoms_by_element("B", atom_type="pdb"),
         set((self.all_atoms[1], self.all_atoms[6]))
        )


    def test_can_get_generic_atoms_by_element(self):
        self.assertEqual(
         self.atomic_structure.get_atoms_by_element("B", atom_type="generic"),
         set((self.all_atoms[11], self.all_atoms[16]))
        )


    def test_default_element_atom_retrieval_is_all(self):
        self.assertEqual(
         self.atomic_structure.get_atoms_by_element("B"),
         set((self.all_atoms[1], self.all_atoms[6], self.all_atoms[11], self.all_atoms[16]))
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


    def test_can_get_single_pdb_atom_by_element(self):
        self.assertIn(
         self.atomic_structure.get_atom_by_element("B", atom_type="pdb"),
         (self.all_atoms[1], self.all_atoms[6])
        )


    def test_can_get_single_generic_atom_by_element(self):
        self.assertIn(
         self.atomic_structure.get_atom_by_element("B", atom_type="generic"),
         (self.all_atoms[11], self.all_atoms[16])
        )


    def test_default_element_single_atom_retrieval_is_all(self):
        self.assertIn(
         self.atomic_structure.get_atom_by_element("B"),
         (self.all_atoms[1], self.all_atoms[6], self.all_atoms[11], self.all_atoms[16])
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


    def test_can_get_pdb_atoms_by_name(self):
        self.assertEqual(
         self.atomic_structure.get_atoms_by_name("BX", atom_type="pdb"),
         set((self.all_atoms[1], self.all_atoms[6]))
        )


    def test_can_get_generic_atoms_by_name(self):
        self.assertEqual(
         self.atomic_structure.get_atoms_by_name("BX", atom_type="generic"),
         set((self.all_atoms[11], self.all_atoms[16]))
        )


    def test_default_name_atom_retrieval_is_all(self):
        self.assertEqual(
         self.atomic_structure.get_atoms_by_name("BX"),
         set((self.all_atoms[1], self.all_atoms[6], self.all_atoms[11], self.all_atoms[16]))
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


    def test_can_get_single_pdb_atom_by_name(self):
        self.assertIn(
         self.atomic_structure.get_atom_by_name("BX", atom_type="pdb"),
         (self.all_atoms[1], self.all_atoms[6])
        )


    def test_can_get_single_generic_atom_by_name(self):
        self.assertIn(
         self.atomic_structure.get_atom_by_name("BX", atom_type="generic"),
         (self.all_atoms[11], self.all_atoms[16])
        )


    def test_default_name_single_atom_retrieval_is_all(self):
        self.assertIn(
         self.atomic_structure.get_atom_by_name("BX"),
         (self.all_atoms[1], self.all_atoms[6], self.all_atoms[11], self.all_atoms[16])
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
        self.pdb_atoms = [
         PdbAtom(x_values[i], y_values[i], 10.0, "C", i + 1, "CX") for i in range(10)
        ]

    def test_can_get_contacts_between_atomic_structures(self):
        structure1 = AtomicStructure(*self.pdb_atoms[:5])
        structure2 = AtomicStructure(*self.pdb_atoms[5:])
        self.assertEqual(
         structure1.contacts_with(structure2, distance=30),
         set([frozenset([self.pdb_atoms[4], self.pdb_atoms[7]])])
        )
        self.assertEqual(
         structure1.contacts_with(structure2, distance=35),
         set([
          frozenset([self.pdb_atoms[4], self.pdb_atoms[7]]),
          frozenset([self.pdb_atoms[4], self.pdb_atoms[6]]),
          frozenset([self.pdb_atoms[4], self.pdb_atoms[8]])
         ])
        )
        self.assertEqual(
         structure1.contacts_with(structure2, distance=40),
         set([
          frozenset([self.pdb_atoms[4], self.pdb_atoms[7]]),
          frozenset([self.pdb_atoms[4], self.pdb_atoms[6]]),
          frozenset([self.pdb_atoms[4], self.pdb_atoms[8]]),
          frozenset([self.pdb_atoms[3], self.pdb_atoms[7]]),
          frozenset([self.pdb_atoms[4], self.pdb_atoms[5]]),
          frozenset([self.pdb_atoms[4], self.pdb_atoms[9]]),
         ])
        )


    def test_external_contacts_can_ignore_hydrogens(self):
        self.pdb_atoms[4].element("H")
        self.pdb_atoms[5].element("H")
        self.pdb_atoms[9].element("H")
        structure1 = AtomicStructure(*self.pdb_atoms[:5])
        structure2 = AtomicStructure(*self.pdb_atoms[5:])
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
         set([frozenset([self.pdb_atoms[3], self.pdb_atoms[7]])])
        )


    def test_can_get_external_contacts_when_one_structure_is_part_of_the_other(self):
        structure1 = AtomicStructure(*self.pdb_atoms[:5])
        structure2 = AtomicStructure(*self.pdb_atoms)
        self.assertEqual(
         structure1.contacts_with(structure2, distance=30),
         set([frozenset([self.pdb_atoms[4], self.pdb_atoms[7]])])
        )
        self.assertEqual(
         structure1.contacts_with(structure2, distance=35),
         set([
          frozenset([self.pdb_atoms[4], self.pdb_atoms[7]]),
          frozenset([self.pdb_atoms[4], self.pdb_atoms[6]]),
          frozenset([self.pdb_atoms[4], self.pdb_atoms[8]])
         ])
        )
        self.assertEqual(
         structure1.contacts_with(structure2, distance=40),
         set([
          frozenset([self.pdb_atoms[4], self.pdb_atoms[7]]),
          frozenset([self.pdb_atoms[4], self.pdb_atoms[6]]),
          frozenset([self.pdb_atoms[4], self.pdb_atoms[8]]),
          frozenset([self.pdb_atoms[3], self.pdb_atoms[7]]),
          frozenset([self.pdb_atoms[4], self.pdb_atoms[5]]),
          frozenset([self.pdb_atoms[4], self.pdb_atoms[9]]),
         ])
        )


    def test_can_get_internal_contacts(self):
        warnings.simplefilter("ignore")
        for index, atom in enumerate(self.pdb_atoms[:-1]):
            atom.bond_to(self.pdb_atoms[index + 1])
        structure = AtomicStructure(*self.pdb_atoms)
        self.assertEqual(
         structure.internal_contacts(distance=31),
         set([
          frozenset([self.pdb_atoms[0], self.pdb_atoms[3]]),
          frozenset([self.pdb_atoms[1], self.pdb_atoms[4]]),
          frozenset([self.pdb_atoms[5], self.pdb_atoms[8]]),
          frozenset([self.pdb_atoms[6], self.pdb_atoms[9]]),
          frozenset([self.pdb_atoms[4], self.pdb_atoms[7]])
         ])
        )


    def test_internal_contacts_can_ignore_hydrogens(self):
        warnings.simplefilter("ignore")
        for index, atom in enumerate(self.pdb_atoms[:-1]):
            atom.bond_to(self.pdb_atoms[index + 1])
        self.pdb_atoms[4].element("H")
        self.pdb_atoms[5].element("H")
        self.pdb_atoms[9].element("H")
        structure = AtomicStructure(*self.pdb_atoms)
        self.assertEqual(
         structure.internal_contacts(distance=31, include_hydrogens=False),
         set([
          frozenset([self.pdb_atoms[0], self.pdb_atoms[3]])
         ])
        )
