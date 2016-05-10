from unittest import TestCase
from molecupy import exceptions
from molecupy.molecules import Atom, Molecule
from molecupy.macromolecules import Residue, ResiduicStructure, Chain, MacroModel

class ChainTest(TestCase):

    def setUp(self):
        self.atom1 = Atom(1.0, 1.0, 1.0, "H", atom_id=1, atom_name="H1")
        self.atom2 = Atom(1.0, 1.0, 2.0, "C", atom_id=2, atom_name="CA")
        self.atom3 = Atom(1.0, 1.0, 3.0, "O", atom_id=3, atom_name="OX1")
        self.atom2.covalent_bond_to(self.atom1)
        self.atom2.covalent_bond_to(self.atom3)
        self.residue1 = Residue(1, "MON1", self.atom1, self.atom2, self.atom3)
        self.atom4 = Atom(1.0, 1.0, 4.0, "H", atom_id=4, atom_name="H1")
        self.atom5 = Atom(1.0, 1.0, 5.0, "C", atom_id=5, atom_name="CA")
        self.atom6 = Atom(1.0, 1.0, 6.0, "O", atom_id=6, atom_name="OX1")
        self.atom5.covalent_bond_to(self.atom4)
        self.atom5.covalent_bond_to(self.atom6)
        self.residue2 = Residue(2, "MON2", self.atom4, self.atom5, self.atom6)
        self.atom7 = Atom(1.0, 1.0, 7.0, "H", atom_id=7, atom_name="H1")
        self.atom8 = Atom(1.0, 1.0, 8.0, "C", atom_id=8, atom_name="CA")
        self.atom9 = Atom(1.0, 1.0, 9.0, "O", atom_id=9, atom_name="OX1")
        self.atom8.covalent_bond_to(self.atom7)
        self.atom8.covalent_bond_to(self.atom9)
        self.residue3 = Residue(3, "MON3", self.atom7, self.atom8, self.atom9)
        self.residue1.connect_to(self.residue2, self.atom3, self.atom4)
        self.residue2.connect_to(self.residue3, self.atom6, self.atom7)


    def check_valid_chain(self, chain, check_id=False):
        self.assertIsInstance(chain, Chain)
        self.assertIsInstance(chain, ResiduicStructure)
        self.assertIsInstance(chain, Molecule)
        self.assertIsInstance(chain.residues, tuple)
        if check_id:
            self.assertIsInstance(chain.chain_id, str)
        for residue in chain.residues:
            self.assertEqual(residue.chain, chain)
        if chain.model is not None:
            self.assertIsInstance(chain.model, MacroModel)
        self.assertRegex(
         str(chain),
         r"<Chain \((\d+) residues\)>"
        )



class ChainCreationTests(ChainTest):

    def test_can_create_chain(self):
        chain = Chain(
         self.residue1,
         self.residue2,
         self.residue3
        )
        self.check_valid_chain(chain)


    def test_can_create_chain_with_id(self):
        chain = Chain(
         self.residue1,
         self.residue2,
         self.residue3,
         chain_id="X"
        )
        self.check_valid_chain(chain, check_id=True)


    def test_chain_id_must_be_str(self):
        with self.assertRaises(TypeError):
            chain = Chain(
             self.residue1,
             self.residue2,
             self.residue3,
             chain_id=1
            )


    def test_chain_needs_connected_residues(self):
        atom10 = Atom(1.0, 1.0, 10.0, "H", atom_id=10, atom_name="H1")
        atom11 = Atom(1.0, 1.0, 11.0, "C", atom_id=11, atom_name="CA")
        atom12 = Atom(1.0, 1.0, 12.0, "O", atom_id=12, atom_name="OX1")
        self.atom9.covalent_bond_to(atom10)
        atom10.covalent_bond_to(atom11)
        atom11.covalent_bond_to(atom12)
        residue4 = Residue(4, "MON4", atom10, atom11, atom12)
        with self.assertRaises(exceptions.BrokenChainError):
            chain = Chain(
             self.residue1,
             self.residue2,
             self.residue3,
             residue4
            )



class ChainBehaviorTests(ChainTest):

    def test_chain_is_indexable(self):
        chain = Chain(
         self.residue1,
         self.residue2,
         self.residue3
        )
        self.assertEqual(chain[0], self.residue1)
        self.assertEqual(chain[1], self.residue2)
        self.assertEqual(chain[2], self.residue3)
        self.assertEqual(chain[0:2], (self.residue1, self.residue2))
        self.assertEqual(chain[1:3], (self.residue2, self.residue3))
        self.assertEqual(
         chain[::-1],
         (self.residue3, self.residue2, self.residue1)
        )
