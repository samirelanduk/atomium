from unittest import TestCase
import unittest.mock
from molecupy.structures import AtomicStructure, PdbAtom, Atom
from molecupy import NoAtomsError

class AtomicStructureTest(TestCase):

    def setUp(self):
        self.pdb_atoms = [unittest.mock.Mock(spec=PdbAtom) for _ in range(10)]
        self.generic_atoms = [unittest.mock.Mock(spec=Atom) for _ in range(10)]



class AtomicStructureCreationTests(AtomicStructureTest):

    def test_can_create_atomic_structure_with_pdb_atoms(self):
        atomic_structure = AtomicStructure(*self.pdb_atoms)
        self.assertEqual(atomic_structure._atoms, set(self.pdb_atoms))


    def test_can_create_atomic_structure_with_generic_atoms(self):
        atomic_structure = AtomicStructure(*self.generic_atoms)
        self.assertEqual(atomic_structure._atoms, set(self.generic_atoms))


    def test_can_create_atomic_structure_with_mixed_atoms(self):
        all_atoms = self.pdb_atoms + self.generic_atoms
        self.assertEqual(len(all_atoms), 20)
        atomic_structure = AtomicStructure(*all_atoms)
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
