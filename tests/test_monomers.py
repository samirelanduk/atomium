from unittest import TestCase
from molecupy import exceptions
from molecupy.atomic import Molecule, Atom
from molecupy.polymers import Monomer

class MonomerTest(TestCase):

    def setUp(self):
        self.atom1 = Atom(1.0, 1.0, 1.0, "H", atom_id=1, atom_name="H1")
        self.atom2 = Atom(1.0, 1.0, 2.0, "C", atom_id=2, atom_name="CA")
        self.atom3 = Atom(1.0, 1.0, 3.0, "O", atom_id=3, atom_name="OX1")
        self.atom2.covalent_bond_to(self.atom1)
        self.atom2.covalent_bond_to(self.atom3)
        self.atom4 = Atom(1.0, 1.0, 4.0, "H", atom_id=4, atom_name="H1")
        self.atom5 = Atom(1.0, 1.0, 5.0, "C", atom_id=5, atom_name="CA")
        self.atom6 = Atom(1.0, 1.0, 6.0, "O", atom_id=6, atom_name="OX1")
        self.atom5.covalent_bond_to(self.atom4)
        self.atom5.covalent_bond_to(self.atom6)
        self.atom7 = Atom(1.0, 1.0, 7.0, "H", atom_id=7, atom_name="H1")
        self.atom8 = Atom(1.0, 1.0, 8.0, "C", atom_id=8, atom_name="CA")
        self.atom9 = Atom(1.0, 1.0, 9.0, "O", atom_id=9, atom_name="OX1")
        self.atom8.covalent_bond_to(self.atom7)
        self.atom8.covalent_bond_to(self.atom9)


    def check_valid_monomer(self, monomer):
        self.assertIsInstance(monomer, Monomer)
        self.assertIsInstance(monomer, Molecule)
        for atom in monomer.atoms:
            self.assertEqual(atom.molecule, monomer)
        self.assertIsInstance(monomer.monomer_id, int)
        self.assertIsInstance(monomer.monomer_name, str)
        with self.assertRaises(AttributeError):
            monomer.molecule_id
        with self.assertRaises(AttributeError):
            monomer.molecule_name
        self.assertRegex(str(monomer), r"<Monomer \((.+)\)>")



class MonomerCreationTest(MonomerTest):

    def test_can_create_monomer(self):
        monomer = Monomer(1, "MON", self.atom1, self.atom2, self.atom3)
        self.check_valid_monomer(monomer)



class MonomerConnectionTest(MonomerTest):

    def setUp(self):
        MonomerTest.setUp(self)
        self.monomer1 = Monomer(1, "MON1", self.atom1, self.atom2, self.atom3)
        self.monomer2 = Monomer(2, "MON2", self.atom4, self.atom5, self.atom6)
        self.monomer3 = Monomer(3, "MON3", self.atom7, self.atom8, self.atom9)
        self.monomer1.connect_to(self.monomer2, self.atom3, self.atom4)
        self.monomer2.connect_to(self.monomer3, self.atom6, self.atom7)


    def test_can_connect_monomers(self):
        self.assertIn(self.atom4, self.atom3.get_covalent_bonded_atoms())
        self.assertIn(self.atom3, self.atom4.get_covalent_bonded_atoms())
        self.assertIn(self.atom7, self.atom6.get_covalent_bonded_atoms())
        self.assertIn(self.atom6, self.atom7.get_covalent_bonded_atoms())


    def test_monomer_connection_needs_correct_atoms(self):
        self.monomer1 = Monomer(1, "MON1", self.atom1, self.atom2, self.atom3)
        self.monomer2 = Monomer(2, "MON2", self.atom4, self.atom5, self.atom6)
        with self.assertRaises(exceptions.InvalidAtomError):
            self.monomer1.connect_to(self.monomer2, self.atom3, self.atom7)
        with self.assertRaises(exceptions.InvalidAtomError):
            self.monomer1.connect_to(self.monomer2, self.atom4, self.atom5)


    def test_can_get_connected_monomers(self):
        self.assertEqual(self.monomer1.downstream_monomer, self.monomer2)
        self.assertEqual(self.monomer2.downstream_monomer, self.monomer3)
        self.assertEqual(self.monomer3.upstream_monomer, self.monomer2)
        self.assertEqual(self.monomer2.upstream_monomer, self.monomer1)


    def test_can_get_all_connected_monomers(self):
        self.assertEqual(
         self.monomer1.get_downstream_monomers(),
         set([self.monomer2, self.monomer3])
        )
        self.assertEqual(
         self.monomer2.get_downstream_monomers(),
         set([self.monomer3])
        )
        self.assertEqual(
         self.monomer3.get_upstream_monomers(),
         set([self.monomer2, self.monomer1])
        )
        self.assertEqual(
         self.monomer2.get_upstream_monomers(),
         set([self.monomer1])
        )
        self.assertEqual(
         self.monomer1.get_accessible_monomers(),
         set([self.monomer3, self.monomer2])
        )
        self.assertEqual(
         self.monomer2.get_accessible_monomers(),
         set([self.monomer1, self.monomer3])
        )
        self.assertEqual(
         self.monomer3.get_accessible_monomers(),
         set([self.monomer1, self.monomer2])
        )
