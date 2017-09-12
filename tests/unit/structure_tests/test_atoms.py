from unittest import TestCase
from unittest.mock import patch, Mock
from atomium.structures.atoms import Atom, Bond

class AtomCreationTests(TestCase):

    def test_can_create_atom(self):
        atom = Atom("C", 2, 3, 5)
        self.assertEqual(atom._element, "C")
        self.assertEqual(atom._x, 2)
        self.assertEqual(atom._y, 3)
        self.assertEqual(atom._z, 5)
        self.assertEqual(atom._id, None)
        self.assertEqual(atom._name, None)
        self.assertEqual(atom._charge, 0)
        self.assertEqual(atom._bfactor, 0)
        self.assertEqual(atom._bonds, set())
        self.assertEqual(atom._residue, None)
        self.assertEqual(atom._chain, None)
        self.assertEqual(atom._molecule, None)
        self.assertEqual(atom._model, None)


    def test_atom_element_must_be_str(self):
        with self.assertRaises(TypeError):
            Atom(1, 2, 3, 5)


    def test_atom_element_must_be_1_or_2_chars(self):
        with self.assertRaises(ValueError):
            Atom("", 2, 3, 5)
        with self.assertRaises(ValueError):
            Atom("XXX", 2, 3, 5)
        Atom("XX", 2, 3, 5)
        Atom("X", 2, 3, 5)


    def test_atom_x_coord_must_be_number(self):
        with self.assertRaises(TypeError):
            Atom("C", "2", 3, 5)
        Atom("C", 2.5, 3, 5)


    def test_atom_y_coord_must_be_number(self):
        with self.assertRaises(TypeError):
            Atom("C", 2, "3", 5)
        Atom("C", 2, 3.5, 5)


    def test_atom_z_coord_must_be_number(self):
        with self.assertRaises(TypeError):
            Atom("C", 2, 3, "5")
        Atom("C", 2, 3, 5.5)


    def test_can_create_atom_with_id(self):
        atom = Atom("C", 2, 3, 5, atom_id=20)
        self.assertEqual(atom._id, 20)


    def test_atom_id_must_be_integer(self):
        with self.assertRaises(TypeError):
            Atom("C", 2, 3, 5, atom_id=20.5)


    def test_can_create_atom_with_name(self):
        atom = Atom("C", 2, 3, 5, name="CA")
        self.assertEqual(atom._name, "CA")


    def test_atom_name_must_be_str(self):
        with self.assertRaises(TypeError):
            Atom("C", 2, 3, 5, name=20.5)


    def test_can_create_atom_with_charge(self):
        atom = Atom("C", 2, 3, 5, charge=-2.5)
        self.assertEqual(atom._charge, -2.5)


    def test_atom_charge_must_be_number(self):
        with self.assertRaises(TypeError):
            Atom("C", 2, 3, 5, charge="20.5")
        Atom("C", 2, 3, 5, charge=10)


    def test_can_create_atom_with_bfactor(self):
        atom = Atom("C", 2, 3, 5, bfactor=2.5)
        self.assertEqual(atom._bfactor, 2.5)


    def test_atom_bfactor_must_be_number(self):
        with self.assertRaises(TypeError):
            Atom("C", 2, 3, 5, bfactor="20.5")
        Atom("C", 2, 3, 5, bfactor=10)


    def test_atom_bfactor_must_be_positive(self):
        with self.assertRaises(ValueError):
            Atom("C", 2, 3, 5, bfactor=-3)



class AtomReprTests(TestCase):

    def test_atom_repr(self):
        atom = Atom("C", 2, 3, 5)
        self.assertEqual(str(atom), "<C Atom at (2, 3, 5)>")


    def test_atom_repr_with_id(self):
        atom = Atom("C", 2, 3, 5, atom_id=15)
        self.assertEqual(str(atom), "<C Atom 15 at (2, 3, 5)>")



class AtomElementTests(TestCase):

    def test_element_property(self):
        atom = Atom("C", 2, 3, 5)
        self.assertIs(atom._element, atom.element())


    def test_can_update_element(self):
        atom = Atom("C", 2, 3, 5)
        atom.element("N")
        self.assertEqual(atom._element, "N")


    def test_atom_element_must_be_str(self):
        atom = Atom("C", 2, 3, 5)
        with self.assertRaises(TypeError):
            atom.element(1)


    def test_atom_element_must_be_1_or_2_chars(self):
        atom = Atom("C", 2, 3, 5)
        with self.assertRaises(ValueError):
            atom.element("")
        with self.assertRaises(ValueError):
            atom.element("XXX")
        atom.element("XX")



class AtomXTests(TestCase):

    def test_x_property(self):
        atom = Atom("C", 2, 3, 5)
        self.assertIs(atom._x, atom.x())


    def test_can_update_x(self):
        atom = Atom("C", 2, 3, 5)
        atom.x(6)
        self.assertEqual(atom._x, 6)


    def test_atom_x_must_be_numeric(self):
        atom = Atom("C", 2, 3, 5)
        with self.assertRaises(TypeError):
            atom.x("4")
        atom.x(4.5)



class AtomYTests(TestCase):

    def test_y_property(self):
        atom = Atom("C", 2, 3, 5)
        self.assertIs(atom._y, atom.y())


    def test_can_update_y(self):
        atom = Atom("C", 2, 3, 5)
        atom.y(6)
        self.assertEqual(atom._y, 6)


    def test_atom_y_must_be_numeric(self):
        atom = Atom("C", 2, 3, 5)
        with self.assertRaises(TypeError):
            atom.y("4")
        atom.y(4.5)



class AtomZTests(TestCase):

    def test_z_property(self):
        atom = Atom("C", 2, 3, 5)
        self.assertIs(atom._z, atom.z())


    def test_can_update_z(self):
        atom = Atom("C", 2, 3, 5)
        atom.z(6)
        self.assertEqual(atom._z, 6)


    def test_atom_z_must_be_numeric(self):
        atom = Atom("C", 2, 3, 5)
        with self.assertRaises(TypeError):
            atom.z("4")
        atom.z(4.5)



class AtomLocationTests(TestCase):

    def test_atom_location(self):
        atom = Atom("C", 20, 30, 50)
        self.assertEqual(atom.location(), (20, 30, 50))



class AtomIdTests(TestCase):

    def test_id_property(self):
        atom = Atom("C", 2, 3, 5, atom_id=100)
        self.assertIs(atom._id, atom.atom_id())



class AtomNameTests(TestCase):

    def test_name_property(self):
        atom = Atom("C", 2, 3, 5, name="CA")
        self.assertIs(atom._name, atom.name())


    def test_can_update_name(self):
        atom = Atom("C", 2, 3, 5, name="CA")
        atom.name("CB")
        self.assertEqual(atom._name, "CB")


    def test_atom_name_must_be_str(self):
        atom = Atom("C", 2, 3, 5, name="CA")
        with self.assertRaises(TypeError):
            atom.name(4)



class AtomResidueTests(TestCase):

    def test_residue_property(self):
        residue = Mock()
        atom = Atom("C", 2, 3, 5)
        atom._residue = residue
        self.assertIs(atom.residue(), residue)



class AtomChainTests(TestCase):

    def test_chain_property(self):
        chain = Mock()
        atom = Atom("C", 2, 3, 5)
        atom._chain = chain
        self.assertIs(atom.chain(), chain)



class AtomMoleculeTests(TestCase):

    def test_molecule_property(self):
        molecule = Mock()
        atom = Atom("C", 2, 3, 5)
        atom._molecule = molecule
        self.assertIs(atom.molecule(), molecule)



class AtomModelTests(TestCase):

    def test_model_property(self):
        model = Mock()
        atom = Atom("C", 2, 3, 5)
        atom._model = model
        self.assertIs(atom.model(), model)



class AtomChargeTests(TestCase):

    def test_charge_property(self):
        atom = Atom("C", 2, 3, 5, charge=0.5)
        self.assertIs(atom._charge, atom.charge())


    def test_can_update_charge(self):
        atom = Atom("C", 2, 3, 5)
        atom.charge(6)
        self.assertEqual(atom._charge, 6)


    def test_atom_charge_must_be_numeric(self):
        atom = Atom("C", 2, 3, 5)
        with self.assertRaises(TypeError):
            atom.charge("4")
        atom.charge(4.5)



class AtomBfactorTests(TestCase):

    def test_bfactor_property(self):
        atom = Atom("C", 2, 3, 5, bfactor=0.5)
        self.assertIs(atom._bfactor, atom.bfactor())


    def test_can_update_bfactor(self):
        atom = Atom("C", 2, 3, 5)
        atom.bfactor(6)
        self.assertEqual(atom._bfactor, 6)


    def test_atom_bfactor_must_be_numeric(self):
        atom = Atom("C", 2, 3, 5)
        with self.assertRaises(TypeError):
            atom.bfactor("4")
        atom.bfactor(4.5)


    def test_atom_bfactor_must_be_positive(self):
        atom = Atom("C", 2, 3, 5)
        with self.assertRaises(ValueError):
            atom.bfactor(-9)



class AtomBondsTests(TestCase):

    def test_bonds_property(self):
        atom = Atom("C", 2, 3, 5)
        atom._bonds = set(("bond1", "bond2"))
        self.assertEqual(atom.bonds(), atom._bonds)
        self.assertIsNot(atom.bonds(), atom._bonds)



class BondedAtomTests(TestCase):

    @patch("atomium.structures.atoms.Atom.bonds")
    def test_can_get_bonded_atoms(self, mock_bonds):
        atom = Atom("C", 2, 3, 5)
        bond1, bond2, bond3 = Mock(Bond), Mock(Bond), Mock(Bond)
        atom2, atom3, atom4 = Mock(Atom), Mock(Atom), Mock(Atom)
        bond1.atoms.return_value = set([atom, atom2])
        bond2.atoms.return_value = set([atom, atom3])
        bond3.atoms.return_value = set([atom, atom4])
        mock_bonds.return_value = set([bond1, bond2, bond3])
        self.assertEqual(atom.bonded_atoms(), set([atom2, atom3, atom4]))



class AtomBondingTests(TestCase):

    @patch("atomium.structures.atoms.Bond")
    def test_can_bond_other_atom(self, mock_bond):
        atom1 = Atom("C", 2, 3, 5)
        atom2 = Mock(Atom)
        atom1.bond(atom2)
        mock_bond.assert_called_with(atom1, atom2)


    @patch("atomium.structures.atoms.Bond")
    @patch("atomium.structures.atoms.Atom.bonded_atoms")
    def test_cant_make_second_bond(self, mock_atoms, mock_bond):
        atom1 = Atom("C", 2, 3, 5)
        atom2 = Mock(Atom)
        mock_atoms.return_value = set([atom2])
        atom1.bond(atom2)
        self.assertFalse(mock_bond.called)



class AtomBondBreakingTests(TestCase):

    @patch("atomium.structures.atoms.Atom.bonds")
    def test_can_break_bond_between_atoms(self, mock_bonds):
        atom1 = Atom("C", 2, 3, 5)
        atom2 = Mock(Atom)
        bond = Mock(Bond)
        bond.atoms.return_value = set([atom1, atom2])
        mock_bonds.return_value = set([bond])
        atom2.bonds.return_value = set([bond])
        atom1.unbond(atom2)
        bond.destroy.assert_called_with()


    def test_bond_breaking_needs_atoms(self):
        atom1 = Atom("C", 2, 3, 5)
        with self.assertRaises(TypeError):
            atom1.unbond("atom2")


    @patch("atomium.structures.atoms.Atom.bonds")
    def test_bond_breaking_needs_other_atom(self, mock_bonds):
        atom1 = Atom("C", 2, 3, 5)
        bond = Mock(Bond)
        bond.atoms.return_value = set([atom1, Mock()])
        mock_bonds.return_value = set([bond])
        with self.assertRaises(ValueError):
            atom1.unbond(atom1)


    @patch("atomium.structures.atoms.Atom.bonds")
    def test_bond_unbreaking_needs_actual_bond(self, mock_bonds):
        atom1 = Atom("C", 2, 3, 5)
        atom2 = Mock(Atom)
        bond = Mock(Bond)
        bond.atoms.return_value = set([atom1, Mock(Atom)])
        mock_bonds.return_value = set([bond])
        atom2.bonds.return_value = set([bond])
        with self.assertRaises(ValueError):
            atom1.unbond(atom2)



class AtomBondWithTests(TestCase):

    @patch("atomium.structures.atoms.Atom.bonds")
    def test_can_get_bond_with_atom(self, mock_bonds):
        atom1 = Atom("C", 2, 3, 5)
        atom2 = Mock(Atom)
        atom3 = Mock(Atom)
        atom4 = Mock(Atom)
        bonda = Mock(Bond)
        bondb = Mock(Bond)
        bonda.atoms.return_value = set([atom1, atom2])
        bondb.atoms.return_value = set([atom1, atom3])
        mock_bonds.return_value = set([bonda, bondb])
        self.assertIs(atom1.bond_with(atom2), bonda)
        self.assertIs(atom1.bond_with(atom3), bondb)
        self.assertIs(atom1.bond_with(atom4), None)
        self.assertIs(atom1.bond_with(atom1), None)


    def test_bond_with_needs_atom(self):
        atom1 = Atom("C", 2, 3, 5)
        with self.assertRaises(TypeError):
            atom1.bond_with("atom2")



class AtomMassTests(TestCase):

    def test_known_element_mass(self):
        atom = Atom("C", 2, 3, 5)
        self.assertAlmostEqual(atom.mass(), 12, delta=0.1)
        atom._element = "H"
        self.assertAlmostEqual(atom.mass(), 1, delta=0.1)


    def test_atom_mass_case_insensitive(self):
        atom = Atom("he", 2, 3, 5)
        self.assertAlmostEqual(atom.mass(), 4, delta=0.1)
        atom = Atom("He", 2, 3, 5)
        self.assertAlmostEqual(atom.mass(), 4, delta=0.1)
        atom = Atom("hE", 2, 3, 5)
        self.assertAlmostEqual(atom.mass(), 4, delta=0.1)
        atom = Atom("HE", 2, 3, 5)
        self.assertAlmostEqual(atom.mass(), 4, delta=0.1)


    def test_unknown_atom_mass(self):
        atom = Atom("XX", 2, 3, 5)
        self.assertEqual(atom.mass(), 0)



class AtomDistanceToTests(TestCase):

    def test_can_get_distance_between_atoms(self):
        atom1 = Atom("C", 4, 8, 3)
        atom2 = Mock(Atom)
        atom2.location.return_value = (2, 3, 5)
        self.assertAlmostEqual(atom1.distance_to(atom2), 5.744, delta=0.001)


    def test_atom_distance_can_be_zero(self):
        atom1 = Atom("C", 4, 8, 3)
        atom2 = Mock(Atom)
        atom2.location.return_value = (4, 8, 3)
        self.assertEqual(atom1.distance_to(atom2), 0)


    def test_other_atom_must_be_atom(self):
        atom1 = Atom("C", 4, 8, 3)
        atom2 = "atom"
        with self.assertRaises(TypeError):
            atom1.distance_to(atom2)


    def test_can_get_distance_to_xyz_tuple(self):
        atom1 = Atom("C", 4, 8, 3)
        atom2 = (2, 3, 5)
        self.assertAlmostEqual(atom1.distance_to(atom2), 5.744, delta=0.001)
