from unittest import TestCase
from molecupy import exceptions
from molecupy.molecules import Molecule, Atom
from molecupy.macromolecules import Residue

class ResidueTest(TestCase):

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


    def check_valid_residue(self, residue):
        self.assertIsInstance(residue, Residue)
        self.assertIsInstance(residue, Molecule)
        for atom in residue.atoms:
            self.assertEqual(atom.molecule, residue)
        self.assertIsInstance(residue.residue_id, int)
        self.assertIsInstance(residue.residue_name, str)
        with self.assertRaises(AttributeError):
            residue.molecule_id
        with self.assertRaises(AttributeError):
            residue.molecule_name
        self.assertRegex(str(residue), r"<Residue \((.+)\)>")



class ResidueCreationTest(ResidueTest):

    def test_can_create_residue(self):
        residue = Residue(1, "MON", self.atom1, self.atom2, self.atom3)
        self.check_valid_residue(residue)



class ResidueConnectionTest(ResidueTest):

    def setUp(self):
        ResidueTest.setUp(self)
        self.residue1 = Residue(1, "MON1", self.atom1, self.atom2, self.atom3)
        self.residue2 = Residue(2, "MON2", self.atom4, self.atom5, self.atom6)
        self.residue3 = Residue(3, "MON3", self.atom7, self.atom8, self.atom9)
        self.residue1.connect_to(self.residue2, self.atom3, self.atom4)
        self.residue2.connect_to(self.residue3, self.atom6, self.atom7)


    def test_can_connect_residues(self):
        self.assertIn(self.atom4, self.atom3.get_covalent_bonded_atoms())
        self.assertIn(self.atom3, self.atom4.get_covalent_bonded_atoms())
        self.assertIn(self.atom7, self.atom6.get_covalent_bonded_atoms())
        self.assertIn(self.atom6, self.atom7.get_covalent_bonded_atoms())


    def test_residue_connection_needs_correct_atoms(self):
        self.residue1 = Residue(1, "MON1", self.atom1, self.atom2, self.atom3)
        self.residue2 = Residue(2, "MON2", self.atom4, self.atom5, self.atom6)
        with self.assertRaises(exceptions.InvalidAtomError):
            self.residue1.connect_to(self.residue2, self.atom3, self.atom7)
        with self.assertRaises(exceptions.InvalidAtomError):
            self.residue1.connect_to(self.residue2, self.atom4, self.atom5)


    def test_can_get_connected_residues(self):
        self.assertEqual(self.residue1.downstream_residue, self.residue2)
        self.assertEqual(self.residue2.downstream_residue, self.residue3)
        self.assertEqual(self.residue3.upstream_residue, self.residue2)
        self.assertEqual(self.residue2.upstream_residue, self.residue1)


    def test_can_get_all_connected_residues(self):
        self.assertEqual(
         self.residue1.get_downstream_residues(),
         (self.residue2, self.residue3)
        )
        self.assertEqual(
         self.residue2.get_downstream_residues(),
         (self.residue3,)
        )
        self.assertEqual(
         self.residue3.get_upstream_residues(),
         (self.residue2, self.residue1)
        )
        self.assertEqual(
         self.residue2.get_upstream_residues(),
         (self.residue1,)
        )
        self.assertEqual(
         self.residue1.get_accessible_residues(),
         set([self.residue3, self.residue2])
        )
        self.assertEqual(
         self.residue2.get_accessible_residues(),
         set([self.residue1, self.residue3])
        )
        self.assertEqual(
         self.residue3.get_accessible_residues(),
         set([self.residue1, self.residue2])
        )


    def test_can_deal_with_cyclic_residue_chains(self):
        atom10 = Atom(1.0, 2.0, 7.0, "H", atom_id=10, atom_name="H1")
        atom11 = Atom(1.0, 2.0, 5.0, "C", atom_id=11, atom_name="CA")
        atom12 = Atom(1.0, 2.0, 3.0, "O", atom_id=12, atom_name="OX1")
        atom10.covalent_bond_to(atom11)
        atom11.covalent_bond_to(atom12)
        residue4 = Residue(4, "MON4", atom10, atom11, atom12)
        self.residue3.connect_to(residue4, self.atom9, atom10)
        residue4.connect_to(self.residue1, atom12, self.atom1)
        self.assertEqual(
         self.residue1.get_downstream_residues(),
         (self.residue2, self.residue3, residue4)
        )
        self.assertEqual(
         self.residue3.get_upstream_residues(),
         (self.residue2, self.residue1, residue4)
        )
        self.assertEqual(
         self.residue1.get_accessible_residues(),
         set([self.residue3, self.residue2, residue4])
        )
