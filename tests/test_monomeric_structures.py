from unittest import TestCase
from molecupy import exceptions
from molecupy.molecules import Atom, AtomicStructure
from molecupy.polymers import Monomer, MonomericStructure

class MonomericStructureTest(TestCase):

    def setUp(self):
        self.atom1 = Atom(1.0, 1.0, 1.0, "H", atom_id=1, atom_name="H1")
        self.atom2 = Atom(1.0, 1.0, 2.0, "C", atom_id=2, atom_name="CA")
        self.atom3 = Atom(1.0, 1.0, 3.0, "O", atom_id=3, atom_name="OX1")
        self.atom2.covalent_bond_to(self.atom1)
        self.atom2.covalent_bond_to(self.atom3)
        self.monomer1 = Monomer(1, "MON1", self.atom1, self.atom2, self.atom3)
        self.atom4 = Atom(1.0, 1.0, 4.0, "H", atom_id=4, atom_name="H1")
        self.atom5 = Atom(1.0, 1.0, 5.0, "C", atom_id=5, atom_name="CA")
        self.atom6 = Atom(1.0, 1.0, 6.0, "O", atom_id=6, atom_name="OX1")
        self.atom5.covalent_bond_to(self.atom4)
        self.atom5.covalent_bond_to(self.atom6)
        self.monomer2 = Monomer(2, "MON2", self.atom4, self.atom5, self.atom6)
        self.atom7 = Atom(1.0, 1.0, 7.0, "H", atom_id=7, atom_name="H1")
        self.atom8 = Atom(1.0, 1.0, 8.0, "C", atom_id=8, atom_name="CA")
        self.atom9 = Atom(1.0, 1.0, 9.0, "O", atom_id=9, atom_name="OX1")
        self.atom8.covalent_bond_to(self.atom7)
        self.atom8.covalent_bond_to(self.atom9)
        self.monomer3 = Monomer(3, "MON3", self.atom7, self.atom8, self.atom9)


    def check_valid_monomeric_structure(self, monomeric_structure):
        self.assertIsInstance(monomeric_structure, MonomericStructure)
        self.assertIsInstance(monomeric_structure, AtomicStructure)
        self.assertIsInstance(monomeric_structure.monomers, set)
        for monomer in monomeric_structure.monomers:
            self.assertIsInstance(monomer, Monomer)
            for atom in monomer.atoms:
                self.assertIn(atom, monomeric_structure.atoms)
        self.assertRegex(
         str(monomeric_structure),
         r"<MonomericStructure \((\d+) monomers\)>"
        )



class MonomericStructureCreationTests(MonomericStructureTest):

    def test_can_create_monomeric_structure(self):
        monomeric_structure = MonomericStructure(
         self.monomer1,
         self.monomer2,
         self.monomer3
        )
        self.check_valid_monomeric_structure(monomeric_structure)


    def test_can_monomeric_structure_needs_monomers(self):
        with self.assertRaises(TypeError):
            monomeric_structure = MonomericStructure(
             self.monomer1,
             self.monomer2,
             self.atom9
            )


    def test_monomeric_structure_needs_at_least_one_monomer(self):
        with self.assertRaises(exceptions.NoMonomersError):
            monomeric_structure = MonomericStructure()


    def test_monomeric_structure_needs_unique_monomeric_ids(self):
        self.monomer3.monomer_id = 2
        with self.assertRaises(exceptions.DuplicateMonomerIdError):
            monomeric_structure = MonomericStructure(
             self.monomer1, self.monomer2, self.monomer3
            )
        self.monomer3.monomer_id = 3
        monomeric_structure = MonomericStructure(
         self.monomer1, self.monomer2, self.monomer3
        )


    def test_monomerics_in_monomeric_structure(self):
        monomeric_structure = MonomericStructure(
         self.monomer1,
         self.monomer2,
         self.monomer3
        )
        self.assertEqual(len(monomeric_structure.monomers), 3)
        self.assertIn(self.monomer1, monomeric_structure)
        self.assertIn(self.monomer2, monomeric_structure)
        self.assertIn(self.monomer3, monomeric_structure)
