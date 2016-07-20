from unittest import TestCase
import unittest.mock
from molecupy.structures import AtomicStructure, PdbAtom, Atom
from molecupy import NoAtomsError

class AtomicStructureTest(TestCase):

    def setUp(self):
        self.pdb_atoms = [unittest.mock.Mock(spec=PdbAtom) for _ in range(10)]
        self.generic_atoms = [unittest.mock.Mock(spec=Atom) for _ in range(10)]
        self.all_atoms = self.pdb_atoms + self.generic_atoms



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


    def test_default_atom_retrieval_is_pdb(self):
        atomic_structure = AtomicStructure(*self.all_atoms)
        self.assertEqual(
         atomic_structure.atoms(),
         set(self.pdb_atoms)
        )


    def test_atomic_structure_atoms_is_read_only(self):
        atom21 = unittest.mock.Mock(spec=PdbAtom)
        atomic_structure = AtomicStructure(*self.all_atoms)
        self.assertEqual(len(atomic_structure.atoms()), 10)
        atomic_structure.atoms(atom_type="all").add(atom21)
        self.assertEqual(len(atomic_structure.atoms()), 10)


    def test_can_add_atom(self):
        atom21 = unittest.mock.Mock(spec=PdbAtom)
        atomic_structure = AtomicStructure(*self.all_atoms)
        atomic_structure.add_atom(atom21)
        self.assertEqual(len(atomic_structure.atoms()), 11)
        self.assertIn(atom21, atomic_structure.atoms())


    def test_can_only_add_atoms(self):
        atomic_structure = AtomicStructure(*self.all_atoms)
        with self.assertRaises(TypeError):
            atomic_structure.add_atom("atom21")


    def test_can_remove_atoms(self):
        atomic_structure = AtomicStructure(*self.all_atoms)
        atomic_structure.remove_atom(self.pdb_atoms[5])
        self.assertEqual(len(atomic_structure.atoms(atom_type="all")), 19)
        self.assertNotIn(self.pdb_atoms[5], atomic_structure.atoms(atom_type="all"))
        atomic_structure.remove_atom(self.generic_atoms[5])
        self.assertEqual(len(atomic_structure.atoms(atom_type="all")), 18)
        self.assertNotIn(self.generic_atoms[5], atomic_structure.atoms(atom_type="all"))
