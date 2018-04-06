import math
from unittest import TestCase
from unittest.mock import patch, Mock, PropertyMock
from atomium.structures.atoms import Atom, Bond, atom_query

class AtomCreationTests(TestCase):

    def test_can_create_atom(self):
        atom = Atom("C")
        self.assertEqual(atom._element, "C")
        self.assertEqual(atom._x, 0)
        self.assertEqual(atom._y, 0)
        self.assertEqual(atom._z, 0)
        self.assertEqual(atom._id, 0)
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


    def test_can_create_atom_with_location(self):
        atom = Atom("C", 2, 3, 4)
        self.assertEqual(atom._x, 2)
        self.assertEqual(atom._y, 3)
        self.assertEqual(atom._z, 4)


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
        atom = Atom("C", 2, 3, 5, id=20)
        self.assertEqual(atom._id, 20)


    def test_id_must_be_integer(self):
        with self.assertRaises(TypeError):
            Atom("C", 2, 3, 5, id=20.5)


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
        Atom("C", 2, 3, 5, bfactor=-3)



class AtomReprTests(TestCase):

    def test_atom_repr(self):
        atom = Atom("C", 2, 3, 5, id=15)
        self.assertEqual(str(atom), "<C Atom 15 at (2, 3, 5)>")



class AtomElementTests(TestCase):

    def test_element_property(self):
        atom = Atom("C", 2, 3, 5)
        self.assertIs(atom._element, atom.element)


    def test_can_update_element(self):
        atom = Atom("C", 2, 3, 5)
        atom.element = "N"
        self.assertEqual(atom._element, "N")


    def test_atom_element_must_be_str(self):
        atom = Atom("C", 2, 3, 5)
        with self.assertRaises(TypeError):
            atom.element = 1



class AtomXTests(TestCase):

    def test_x_property(self):
        atom = Atom("C", 2, 3, 5)
        self.assertIs(atom._x, atom.x)


    def test_can_update_x(self):
        atom = Atom("C", 2, 3, 5)
        atom.x = 6
        self.assertEqual(atom._x, 6)


    def test_atom_x_must_be_numeric(self):
        atom = Atom("C", 2, 3, 5)
        with self.assertRaises(TypeError):
            atom.x = "4"
        atom.x = 4.5



class AtomYTests(TestCase):

    def test_y_property(self):
        atom = Atom("C", 2, 3, 5)
        self.assertIs(atom._y, atom.y)


    def test_can_update_y(self):
        atom = Atom("C", 2, 3, 5)
        atom.y = 6
        self.assertEqual(atom._y, 6)


    def test_atom_y_must_be_numeric(self):
        atom = Atom("C", 2, 3, 5)
        with self.assertRaises(TypeError):
            atom.y = "4"
        atom.y = 4.5



class AtomZTests(TestCase):

    def test_x_property(self):
        atom = Atom("C", 2, 3, 5)
        self.assertIs(atom._z, atom.z)


    def test_can_update_x(self):
        atom = Atom("C", 2, 3, 5)
        atom.z = 6
        self.assertEqual(atom._z, 6)


    def test_atom_x_must_be_numeric(self):
        atom = Atom("C", 2, 3, 5)
        with self.assertRaises(TypeError):
            atom.z = "4"
        atom.z = 4.5



class AtomIdTests(TestCase):

    def test_id_property(self):
        atom = Atom("C", id=100)
        self.assertIs(atom._id, atom.id)



class AtomNameTests(TestCase):

    def test_name_property(self):
        atom = Atom("C", name="CA")
        self.assertIs(atom._name, atom.name)


    def test_can_update_name(self):
        atom = Atom("C", name="CA")
        atom.name = "CB"
        self.assertEqual(atom._name, "CB")


    def test_atom_name_must_be_str(self):
        atom = Atom("C", name="CA")
        with self.assertRaises(TypeError):
            atom.name = 4



class AtomChargeTests(TestCase):

    def test_charge_property(self):
        atom = Atom("C", charge=0.5)
        self.assertIs(atom._charge, atom.charge)


    def test_can_update_charge(self):
        atom = Atom("C")
        atom.charge = 6
        self.assertEqual(atom._charge, 6)


    def test_atom_charge_must_be_numeric(self):
        atom = Atom("C")
        with self.assertRaises(TypeError):
            atom.charge = "4"
        atom.charge = 4.5
        atom.charge = -9



class AtomBfactorTests(TestCase):

    def test_bfactor_property(self):
        atom = Atom("C", bfactor=0.5)
        self.assertIs(atom._bfactor, atom.bfactor)


    def test_can_update_bfactor(self):
        atom = Atom("C")
        atom.bfactor = 6
        self.assertEqual(atom._bfactor, 6)


    def test_atom_bfactor_must_be_numeric(self):
        atom = Atom("C")
        with self.assertRaises(TypeError):
            atom.bfactor = "4"
        atom.bfactor = 4.5
        atom.bfactor = -9



class AtomLocationTests(TestCase):

    def test_atom_location(self):
        atom = Atom("C", 20, 30, 50)
        self.assertEqual(atom.location, (20, 30, 50))



class AtomRoundingTests(TestCase):

    def test_can_round_atom_location(self):
        atom = Atom("C", 0.1111, 0.4444, 0.8888)
        atom.trim(None)
        self.assertEqual(atom._x, 0.1111)
        self.assertEqual(atom._y, 0.4444)
        self.assertEqual(atom._z, 0.8888)
        atom.trim(3)
        self.assertEqual(atom._x, 0.111)
        self.assertEqual(atom._y, 0.444)
        self.assertEqual(atom._z, 0.889)
        atom.trim(2)
        self.assertEqual(atom._x, 0.11)
        self.assertEqual(atom._y, 0.44)
        self.assertEqual(atom._z, 0.89)
        atom.trim(1)
        self.assertEqual(atom._x, 0.1)
        self.assertEqual(atom._y, 0.4)
        self.assertEqual(atom._z, 0.9)
        atom.trim(0)
        self.assertEqual(atom._x, 0)
        self.assertEqual(atom._y, 0)
        self.assertEqual(atom._z, 1)



class AtomTranslationTests(TestCase):

    @patch("atomium.structures.atoms.Atom.trim")
    def test_can_translate_atoms(self, mock_trim):
        atom = Atom("C", 20, 30, 50)
        atom.translate(5, -4, 12)
        self.assertEqual(atom._x, 25)
        self.assertEqual(atom._y, 26)
        self.assertEqual(atom._z, 62)
        mock_trim.assert_called_with(12)


    @patch("atomium.structures.atoms.Atom.trim")
    def test_can_translate_atoms_with_iterable(self, mock_trim):
        atom = Atom("C", 20, 30, 50)
        atom.translate((5, -4, 12))
        self.assertEqual(atom._x, 25)
        self.assertEqual(atom._y, 26)
        self.assertEqual(atom._z, 62)
        mock_trim.assert_called_with(12)


    @patch("atomium.structures.atoms.Atom.trim")
    def test_can_translate_atoms_and_specify_trim(self, mock_trim):
        atom = Atom("C", 20, 30, 50)
        atom.translate(5, -4, 12, trim=18)
        mock_trim.assert_called_with(18)
        atom.translate((5, -4, 12), trim=None)
        mock_trim.assert_called_with(None)



class AtomMovingTests(TestCase):

    def test_can_move_atom(self):
        atom = Atom("C", 20, 30, 50)
        atom.move_to(1, 2, 3)
        self.assertEqual(atom._x, 1)
        self.assertEqual(atom._y, 2)
        self.assertEqual(atom._z, 3)


    def test_coordinates_must_be_numeric(self):
        atom = Atom("C", 20, 30, 50)
        with self.assertRaises(TypeError):
            atom.move_to(1, 2, "3")
        with self.assertRaises(TypeError):
            atom.move_to(1, "2", 3)
        with self.assertRaises(TypeError):
            atom.move_to("1", 2, 3)
        atom.move_to(1.1, 2.1, 3.1)



class AtomRotationMatrixGenerationTests(TestCase):

    def test_can_generate_rotation_matrix_x(self):
        atom = Atom("C")
        matrix = atom.generate_rotation_matrix(math.pi / 2, "x")
        self.assertEqual(
         [[round(val, 12) for val in row] for row in matrix],
         [[1, 0, 0], [0, 0, -1], [0, 1, 0]]
        )
        matrix = atom.generate_rotation_matrix(math.pi / -4, "x")
        self.assertEqual(
         [[round(val, 3) for val in row] for row in matrix],
         [[1, 0, 0], [0, 0.707, 0.707], [0, -0.707, 0.707]]
        )


    def test_can_generate_rotation_matrix_y(self):
        atom = Atom("C")
        matrix = atom.generate_rotation_matrix(math.pi / 2, "y")
        self.assertEqual(
         [[round(val, 12) for val in row] for row in matrix],
          [[0, 0, 1], [0, 1, 0], [-1, 0, 0]]
        )
        matrix = atom.generate_rotation_matrix(math.pi / 18, "y")
        self.assertEqual(
         [[round(val, 3) for val in row] for row in matrix],
         [[0.985, 0, 0.174], [0, 1, 0], [-0.174, 0, 0.985]]
        )


    def test_can_generate_rotation_matrix_z(self):
        atom = Atom("C")
        matrix = atom.generate_rotation_matrix(math.pi / 2, "z")
        self.assertEqual(
         [[round(val, 12) for val in row] for row in matrix],
          [[0, -1, 0], [1, 0, 0], [0, 0, 1]]
        )
        matrix = atom.generate_rotation_matrix(math.pi / -36, "z")
        self.assertEqual(
         [[round(val, 3) for val in row] for row in matrix],
         [[0.996, 0.087, 0], [-0.087, 0.996, 0], [0, 0, 1]]
        )



class AtomRotationTests(TestCase):

    @patch("atomium.structures.atoms.Atom.generate_rotation_matrix")
    @patch("atomium.structures.atoms.Atom.location", new_callable=PropertyMock)
    @patch("atomium.structures.atoms.Atom.trim")
    def test_can_rotate(self, mock_trim, mock_loc, mock_matrix):
        matrix = Mock()
        mock_matrix.return_value = matrix
        matrix.dot.return_value = [10, 20, 30]
        mock_loc.return_value = (0.1, 0.2, 0.3)
        atom = Atom("C", 1, 1, 1)
        atom.rotate(9, "x")
        mock_matrix.assert_called_with(9, "x")
        mock_loc.assert_called_with()
        matrix.dot.assert_called_with((0.1, 0.2, 0.3))
        self.assertEqual(atom._x, 10)
        self.assertEqual(atom._y, 20)
        self.assertEqual(atom._z, 30)
        mock_trim.assert_called_with(12)


    @patch("atomium.structures.atoms.Atom.generate_rotation_matrix")
    @patch("atomium.structures.atoms.Atom.location", new_callable=PropertyMock)
    @patch("atomium.structures.atoms.Atom.trim")
    def test_can_rotate_no_trim(self, mock_trim, mock_loc, mock_matrix):
        matrix = Mock()
        mock_matrix.return_value = matrix
        matrix.dot.return_value = [10, 20, 30]
        atom = Atom("C", 1, 1, 1)
        atom.rotate(9, "x", trim=None)
        mock_trim.assert_called_with(None)


    @patch("atomium.structures.atoms.Atom.generate_rotation_matrix")
    @patch("atomium.structures.atoms.Atom.location", new_callable=PropertyMock)
    @patch("atomium.structures.atoms.Atom.trim")
    def test_can_rotate_in_degrees(self, mock_trim, mock_loc, mock_matrix):
        matrix = Mock()
        mock_matrix.return_value = matrix
        matrix.dot.return_value = [10, 20, 30]
        atom = Atom("C", 1, 1, 1)
        atom.rotate(90, "x", degrees=True)
        mock_matrix.assert_called_with(math.pi / 2, "x")


    def test_axis_must_be_valid(self):
        atom = Atom("C", 1, 1, 1)
        with self.assertRaises(ValueError):
            atom.rotate(90, "s")



class AtomMassTests(TestCase):

    def test_known_element_mass(self):
        atom = Atom("C", 2, 3, 5)
        self.assertAlmostEqual(atom.mass, 12, delta=0.1)
        atom._element = "H"
        self.assertAlmostEqual(atom.mass, 1, delta=0.1)


    def test_atom_mass_case_insensitive(self):
        atom = Atom("he", 2, 3, 5)
        self.assertAlmostEqual(atom.mass, 4, delta=0.1)
        atom = Atom("He", 2, 3, 5)
        self.assertAlmostEqual(atom.mass, 4, delta=0.1)
        atom = Atom("hE", 2, 3, 5)
        self.assertAlmostEqual(atom.mass, 4, delta=0.1)
        atom = Atom("HE", 2, 3, 5)
        self.assertAlmostEqual(atom.mass, 4, delta=0.1)


    def test_unknown_atom_mass(self):
        atom = Atom("XX", 2, 3, 5)
        self.assertEqual(atom.mass, 0)



class AtomDistanceToTests(TestCase):

    def test_can_get_distance_between_atoms(self):
        atom1 = Atom("C", 4, 8, 3)
        atom2 = Mock(Atom)
        atom2.location = (2, 3, 5)
        self.assertAlmostEqual(atom1.distance_to(atom2), 5.744, delta=0.001)


    def test_atom_distance_can_be_zero(self):
        atom1 = Atom("C", 4, 8, 3)
        atom2 = Mock(Atom)
        atom2.location = (4, 8, 3)
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



class AtomResidueTests(TestCase):

    def test_residue_property(self):
        residue = Mock()
        atom = Atom("C", 2, 3, 5)
        atom._residue = residue
        self.assertIs(atom.residue, residue)



class AtomChainTests(TestCase):

    def test_chain_property(self):
        chain = Mock()
        atom = Atom("C", 2, 3, 5)
        atom._chain = chain
        self.assertIs(atom.chain, chain)



class AtomMoleculeTests(TestCase):

    def test_molecule_property(self):
        molecule = Mock()
        atom = Atom("C", 2, 3, 5)
        atom._molecule = molecule
        self.assertIs(atom.molecule, molecule)



class AtomModelTests(TestCase):

    def test_model_property(self):
        model = Mock()
        atom = Atom("C", 2, 3, 5)
        atom._model = model
        self.assertIs(atom.model, model)


class AtomBondsTests(TestCase):

    def test_bonds_property(self):
        atom = Atom("C", 2, 3, 5)
        atom._bonds = set(("bond1", "bond2"))
        self.assertEqual(atom.bonds, atom._bonds)
        self.assertIsNot(atom.bonds, atom._bonds)



class BondedAtomTests(TestCase):

    @patch("atomium.structures.atoms.Atom.bonds", new_callable=PropertyMock)
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
        atom1.bond_to(atom2)
        mock_bond.assert_called_with(atom1, atom2)


    @patch("atomium.structures.atoms.Bond")
    @patch("atomium.structures.atoms.Atom.bonded_atoms")
    def test_cant_make_second_bond(self, mock_atoms, mock_bond):
        atom1 = Atom("C", 2, 3, 5)
        atom2 = Mock(Atom)
        mock_atoms.return_value = set([atom2])
        atom1.bond_to(atom2)
        self.assertFalse(mock_bond.called)



class AtomBondBreakingTests(TestCase):

    @patch("atomium.structures.atoms.Atom.bonds", new_callable=PropertyMock)
    def test_can_break_bond_between_atoms(self, mock_bonds):
        atom1 = Atom("C", 2, 3, 5)
        atom2 = Mock(Atom)
        bond = Mock(Bond)
        bond.atoms.return_value = set([atom1, atom2])
        mock_bonds.return_value = set([bond])
        atom2.bonds.return_value = set([bond])
        atom1.unbond_from(atom2)
        bond.destroy.assert_called_with()


    def test_bond_breaking_needs_atoms(self):
        atom1 = Atom("C", 2, 3, 5)
        with self.assertRaises(TypeError):
            atom1.unbond_from("atom2")


    @patch("atomium.structures.atoms.Atom.bonds")
    def test_bond_breaking_needs_other_atom(self, mock_bonds):
        atom1 = Atom("C", 2, 3, 5)
        bond = Mock(Bond)
        bond.atoms.return_value = set([atom1, Mock()])
        mock_bonds.return_value = set([bond])
        with self.assertRaises(ValueError):
            atom1.unbond_from(atom1)


    @patch("atomium.structures.atoms.Atom.bonds")
    def test_bond_unbreaking_needs_actual_bond(self, mock_bonds):
        atom1 = Atom("C", 2, 3, 5)
        atom2 = Mock(Atom)
        bond = Mock(Bond)
        bond.atoms.return_value = set([atom1, Mock(Atom)])
        mock_bonds.return_value = set([bond])
        atom2.bonds.return_value = set([bond])
        with self.assertRaises(ValueError):
            atom1.unbond_from(atom2)



class AtomBondWithTests(TestCase):

    @patch("atomium.structures.atoms.Atom.bonds", new_callable=PropertyMock)
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



class NearbyAtomTests(TestCase):

    def test_atom_with_no_model_has_no_nearby_atoms(self):
        atom = Atom("C", 4, 8, 3)
        self.assertEqual(atom.nearby_atoms(cutoff=1), set())
        self.assertEqual(atom.nearby_atoms(cutoff=100), set())


    @patch("atomium.structures.atoms.Atom.location", new_callable=PropertyMock)
    def test_can_get_nearby_atoms(self, mock_loc):
        model = Mock()
        mock_loc.return_value = (1, 2, 3)
        atom = Atom("C", 4, 8, 3)
        atom._model = model
        model.atoms_in_sphere.return_value = [1, 2, atom, 4]
        atoms = atom.nearby_atoms(4, a=1, b=2)
        model.atoms_in_sphere.assert_called_with(1, 2, 3, 4, a=1, b=2)
        self.assertEqual(atoms, [1, 2, 4])



class AtomCopyingTests(TestCase):

    def test_can_make_atom_copy(self):
        atom = Atom("he", 2, 3, 5, 10, "H1", 0.5, 0.4)
        copy = atom.copy()
        self.assertIsInstance(atom, Atom)
        self.assertIsNot(atom, copy)
        self.assertEqual(copy._element, "he")
        self.assertEqual(copy._x, 2)
        self.assertEqual(copy._y, 3)
        self.assertEqual(copy._z, 5)
        self.assertEqual(copy._id, 10)
        self.assertEqual(copy._name, "H1")
        self.assertEqual(copy._charge, 0.5)
        self.assertEqual(copy._bfactor, 0.4)
        self.assertEqual(copy._bonds, set())
        self.assertEqual(copy._residue, None)
        self.assertEqual(copy._chain, None)
        self.assertEqual(copy._molecule, None)
        self.assertEqual(copy._model, None)



class AtomSelectorTests(TestCase):

    def setUp(self):
        self.atoms = [
         Atom("C", 2, 3, 5, name="CA", id=15, charge=1),
         Atom("C", 2, 3, 5, name="CA", id=16, charge=1),
         Atom("P", 2, 3, 5, name="PB", id=17, charge=1),
         Atom("H", 2, 3, 5, name="H1", id=18, charge=1),
         Atom("ZN", 2, 3, 5, name="ZN", id=19, charge=1)
        ]
        self.atoms[0]._residue, self.atoms[1]._residue = "H", "H"
        def func(a, b, c=20):
            return self.atoms
        self.func = atom_query(func)


    def test_atoms_decorator_can_do_nothing(self):
        self.assertEqual(self.func(1, 2, c=3), self.atoms)


    def test_atom_query_can_search_ids(self):
        self.assertEqual(self.func(1, 2, id=16, c=3), set([self.atoms[1]]))


    def test_atom_query_can_search_names(self):
        self.assertEqual(self.func(1, 2, name="CA", c=3), set(self.atoms[:2]))
        self.assertEqual(self.func(1, 2, name="PB", c=3), set(self.atoms[2:3]))


    def test_atom_query_can_search_elements(self):
        self.assertEqual(self.func(1, 2, element="C", c=3), set(self.atoms[:2]))
        self.assertEqual(self.func(1, 2, element="c", c=3), set(self.atoms[:2]))


    def test_can_exclude_hydrogen(self):
        self.assertEqual(
         self.func(1, 2, hydrogen=False, c=3),
         set(self.atoms) - set(self.atoms[3:4])
        )


    def test_can_exclude_heteroatoms(self):
        self.assertEqual(self.func(1, 2, het=False, c=3), set(self.atoms[:2]))


    def test_can_exclude_metals(self):
        self.assertEqual(self.func(1, 2, metal=False, c=3), set(self.atoms[:4]))


    def test_can_combine_queries(self):
        self.assertEqual(
         self.func(1, 2, hydrogen=False, element="C", c=3), set(self.atoms[:2])
        )
