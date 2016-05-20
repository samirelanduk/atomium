from collections import Counter
from unittest import TestCase
from molecupy import exceptions
from molecupy.structures import AtomicStructure, PdbAtom

class AtomicStructureTest(TestCase):

    def setUp(self):
        self.atom1 = PdbAtom(1.0, 1.0, 1.0, "H", 1, "H1")
        self.atom2 = PdbAtom(1.0, 1.0, 2.0, "C", 2, "CA")
        self.atom3 = PdbAtom(1.0, 1.0, 3.0, "O", 3, "OX1")


    def check_valid_atomic_structure(self, atomic_structure):
        self.assertIsInstance(atomic_structure, AtomicStructure)
        self.assertIsInstance(atomic_structure.atoms, set)
        for atom in atomic_structure.atoms:
            self.assertIsInstance(atom, PdbAtom)
        self.assertRegex(str(atomic_structure), r"<AtomicStructure \((\d+) atoms\)>")



class AtomicStructureCreationTests(AtomicStructureTest):

    def test_can_create_atomic_structure(self):
        atomic_structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.check_valid_atomic_structure(atomic_structure)


    def test_atomic_structure_needs_atoms(self):
        with self.assertRaises(TypeError):
            atomic_structure = AtomicStructure(self.atom1, self.atom2, "atom3")


    def test_atomic_structure_needs_at_least_one_atom(self):
        with self.assertRaises(exceptions.NoAtomsError):
            atomic_structure = AtomicStructure()


    def test_atoms_in_atomic_structure(self):
        atomic_structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(len(atomic_structure.atoms), 3)
        self.assertIn(self.atom1, atomic_structure)
        self.assertIn(self.atom2, atomic_structure)
        self.assertIn(self.atom3, atomic_structure)


    def test_atomic_structure_needs_unique_atom_ids(self):
        atom4 = PdbAtom(1.0, 1.0, 4.0, "H", 1, "H1")
        atom5 = PdbAtom(1.0, 1.0, 5.0, "H", 5, "H1")
        atom6 = PdbAtom(1.0, 1.0, 6.0, "H", 6, "H1")
        with self.assertRaises(exceptions.DuplicateAtomIdError):
            atomic_structure = AtomicStructure(
             self.atom1, self.atom2, self.atom3, atom4, atom5, atom6
            )
        atom4.atom_id = 4
        atomic_structure = AtomicStructure(
         self.atom1, self.atom2, self.atom3, atom4, atom5, atom6
        )



class AtomicStructureBehaviorTests(AtomicStructureTest):

    def test_can_get_atomic_structure_mass(self):
        atomic_structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertAlmostEqual(atomic_structure.get_mass(), 1 + 12 + 16, delta=0.5)



class AtomicStructureAtomRetrievalTests(AtomicStructureTest):

    def test_can_get_atom_by_name(self):
        atom4 = PdbAtom(1.0, 1.0, 4.0, "H", atom_id=4, atom_name="H1")
        atomic_structure = AtomicStructure(self.atom1, self.atom2, self.atom3, atom4)
        self.assertIn(atomic_structure.get_atom_by_name("H1"), (self.atom1, atom4))
        self.assertIs(atomic_structure.get_atom_by_name("xxx"), None)


    def test_can_get_multiple_atoms_by_name(self):
        atom4 = PdbAtom(1.0, 1.0, 4.0, "H", atom_id=4, atom_name="H1")
        atomic_structure = AtomicStructure(self.atom1, self.atom2, self.atom3, atom4)
        self.assertEqual(
         atomic_structure.get_atoms_by_name("H1"),
         set([self.atom1, atom4])
        )
        self.assertEqual(atomic_structure.get_atoms_by_name("XXX"), set())


    def test_can_only_search_by_string_name(self):
        atomic_structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        with self.assertRaises(TypeError):
            atomic_structure.get_atom_by_name(None)
        with self.assertRaises(TypeError):
            atomic_structure.get_atoms_by_name(None)


    def test_can_get_atom_by_element(self):
        atom4 = PdbAtom(1.0, 1.0, 4.0, "H", atom_id=4, atom_name="H1")
        atomic_structure = AtomicStructure(self.atom1, self.atom2, self.atom3, atom4)
        self.assertIn(atomic_structure.get_atom_by_element("H"), (self.atom1, atom4))
        self.assertIs(atomic_structure.get_atom_by_element("xxx"), None)


    def test_can_get_multiple_atoms_by_element(self):
        atom4 = PdbAtom(1.0, 1.0, 4.0, "H", atom_id=4, atom_name="H2")
        atomic_structure = AtomicStructure(self.atom1, self.atom2, self.atom3, atom4)
        self.assertEqual(
         atomic_structure.get_atoms_by_element("H"),
         set([self.atom1, atom4])
        )
        self.assertEqual(
         atomic_structure.get_atoms_by_element("C"),
         set([self.atom2])
        )
        self.assertEqual(atomic_structure.get_atoms_by_element("XXX"), set())


    def test_can_only_search_by_string_element(self):
        atomic_structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        with self.assertRaises(TypeError):
            atomic_structure.get_atom_by_element(None)
        with self.assertRaises(TypeError):
            atomic_structure.get_atoms_by_element(None)


    def test_can_get_atom_by_id(self):
        atomic_structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertIs(atomic_structure.get_atom_by_id(1), self.atom1)
        self.assertIs(atomic_structure.get_atom_by_id(2), self.atom2)
        self.assertIs(atomic_structure.get_atom_by_id(3), self.atom3)
        self.assertIs(atomic_structure.get_atom_by_id(4), None)


    def test_can_only_search_by_numeric_id(self):
        atomic_structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        with self.assertRaises(TypeError):
            atomic_structure.get_atom_by_id(None)




class FormulaTests(AtomicStructureTest):

    def test_can_get_formula(self):
        atom4 = PdbAtom(1.0, 1.0, 4.0, "O", 4, "O")
        atom5 = PdbAtom(1.0, 1.0, 5.0, "N", 5, "NA")
        atomic_structure = AtomicStructure(self.atom1, self.atom2, self.atom3, atom4, atom5)
        self.assertEqual(
         atomic_structure.get_formula(),
         Counter({"C": 1, "N": 1, "O": 2})
        )



class AtomicStructureContactsTest(AtomicStructureTest):

    def setUp(self):
        self.atom1 = PdbAtom(1.0, 10.0, 10.0, "C", 1, "C")
        self.atom2 = PdbAtom(1.0, 10.0, 12.0, "C", 2, "C")
        self.atom3 = PdbAtom(1.0, 10.0, 14.0, "C", 3, "C")
        self.atom4 = PdbAtom(1.0, 12.0, 14.0, "C", 4, "C")
        self.atom5 = PdbAtom(1.0, 14.0, 14.0, "C", 5, "C")
        self.atom6 = PdbAtom(1.0, 14.0, 12.0, "C", 6, "C")
        self.atom7 = PdbAtom(1.0, 14.0, 10.0, "C", 7, "C")
        self.atom8 = PdbAtom(1.0, 14.0, 8.0, "C", 8, "C")
        self.atom9 = PdbAtom(1.0, 14.0, 6.0, "C", 9, "C")
        self.atom10 = PdbAtom(1.0, 12.0, 6.0, "C", 10, "C")
        self.atom11 = PdbAtom(1.0, 10.0, 6.0, "C", 11, "C")
        self.atom12 = PdbAtom(1.0, 8.0, 6.0, "C", 12, "C")
        self.atom13 = PdbAtom(1.0, 6.0, 6.0, "C", 13, "C")
        self.atom14 = PdbAtom(1.0, 6.0, 8.0, "C", 14, "C")
        self.atom15 = PdbAtom(1.0, 6.0, 10.0, "C", 15, "C")
        self.atom16 = PdbAtom(1.0, 6.0, 12.0, "C", 16, "C")
        self.atom17 = PdbAtom(1.0, 6.0, 14.0, "C", 17, "C")
        self.atom18 = PdbAtom(1.0, 6.0, 16.0, "C", 18, "C")
        self.atom19 = PdbAtom(1.0, 6.0, 18.0, "C", 19, "C")
        self.atom20 = PdbAtom(1.0, 6.0, 20.0, "C", 20, "C")
        self.atom21 = PdbAtom(1.0, 6.0, 22.0, "C", 21, "C")
        self.structure1 = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.structure2 = AtomicStructure(
         self.atom1, self.atom2, self.atom3, self.atom4, self.atom5, self.atom6,
         self.atom7, self.atom8, self.atom9, self.atom10, self.atom11, self.atom12,
         self.atom13, self.atom14, self.atom15, self.atom16, self.atom17, self.atom18,
         self.atom19, self.atom20, self.atom21
        )
        for i in range(1, 21):
            self.structure2.get_atom_by_id(i).covalent_bond_to(
             self.structure2.get_atom_by_id(i + 1)
            )


    def test_can_get_contacts_between_atomic_structures(self):
        self.assertEqual(
         len(self.structure1.get_external_contacts_with(self.structure2)),
         9
        )
        self.assertEqual(
         self.structure1.get_external_contacts_with(self.structure2),
         set((
          frozenset((self.atom1, self.atom7)),
          frozenset((self.atom1, self.atom11)),
          frozenset((self.atom1, self.atom15)),
          frozenset((self.atom2, self.atom6)),
          frozenset((self.atom2, self.atom16)),
          frozenset((self.atom2, self.atom4)),
          frozenset((self.atom3, self.atom4)),
          frozenset((self.atom3, self.atom5)),
          frozenset((self.atom3, self.atom17)),
         ))
        )


    def test_can_internal_contacts(self):
        self.assertEqual(
         len(self.structure2.get_internal_contacts()),
         6
        )
        self.assertEqual(
         self.structure2.get_internal_contacts(),
         set((
          frozenset((self.atom1, self.atom7)),
          frozenset((self.atom1, self.atom11)),
          frozenset((self.atom1, self.atom15)),
          frozenset((self.atom2, self.atom6)),
          frozenset((self.atom2, self.atom16)),
          frozenset((self.atom3, self.atom17)),
         ))
        )
