import math
import numpy as np
from unittest import TestCase
from unittest.mock import Mock, patch, MagicMock, PropertyMock
from atomium.structures import Atom, Residue, Ligand

class AtomCreationTests(TestCase):

    def test_can_create_atoms(self):
        atom = Atom("C", 1, 2, 3, 10, "CA", -1, 0.76, [1, 2, 3, 4, 5, 6])
        self.assertEqual(atom._element, "C")
        self.assertEqual(atom._x, 1)
        self.assertEqual(atom._y, 2)
        self.assertEqual(atom._z, 3)
        self.assertEqual(atom._id, 10)
        self.assertEqual(atom._name, "CA")
        self.assertEqual(atom._charge, -1)
        self.assertEqual(atom._bvalue, 0.76)
        self.assertEqual(atom._anisotropy, [1, 2, 3, 4, 5, 6])
        self.assertIsNone(atom._structure)
        self.assertIsNone(atom._model)



class AtomReprTests(TestCase):

    def test_atom_repr(self):
        atom = Atom("C", 1, 2, 3, 10, "CA", -1, 0.76, [1, 2, 3, 4, 5, 6])
        self.assertEqual(repr(atom), "<Atom 10 (CA)>")



class MultiAtomTranslationTests(TestCase):

    def test_can_translate_atoms(self):
        atoms = [Mock(_x=1, _y=2, _z=3), Mock(_x=4, _y=5, _z=6)]
        Atom.translate_atoms([10, 20, -30], atoms[0], atoms[1])
        self.assertEqual(atoms[0]._x, 11)
        self.assertEqual(atoms[0]._y, 22)
        self.assertEqual(atoms[0]._z, -27)
        self.assertEqual(atoms[1]._x, 14)
        self.assertEqual(atoms[1]._y, 25)
        self.assertEqual(atoms[1]._z, -24)



class MultiAtomTransformationTests(TestCase):

    def test_can_transform_atoms(self):
        atoms = [Mock(location=[1, 2, 3]), Mock(location=[4, 5, 6])]
        Atom.transform_atoms([[1, 2, 3], [4, 5, 6], [7, 8, 9]], atoms[0], atoms[1])
        self.assertEqual(atoms[0]._x, 14)
        self.assertEqual(atoms[0]._y, 32)
        self.assertEqual(atoms[0]._z, 50)
        self.assertEqual(atoms[1]._x, 32)
        self.assertEqual(atoms[1]._y, 77)
        self.assertEqual(atoms[1]._z, 122)



class MultiAtomRotationTests(TestCase):

    @patch("atomium.structures.Atom.transform_atoms")
    def test_can_rotate_atoms(self, mock_trans):
        atoms = [Mock(), Mock()]
        Atom.rotate_atoms(math.pi / 2, "x", atoms[0], atoms[1])
        matrix_used = mock_trans.call_args_list[0][0][0]
        self.assertAlmostEqual(matrix_used[0][0], 1, delta=0.000001)
        self.assertAlmostEqual(matrix_used[0][1], 0, delta=0.000001)
        self.assertAlmostEqual(matrix_used[0][2], 0, delta=0.000001)
        self.assertAlmostEqual(matrix_used[1][0], 0, delta=0.000001)
        self.assertAlmostEqual(matrix_used[1][1], 0, delta=0.000001)
        self.assertAlmostEqual(matrix_used[1][2], -1, delta=0.000001)
        self.assertAlmostEqual(matrix_used[2][0], 0, delta=0.000001)
        self.assertAlmostEqual(matrix_used[2][1], 1, delta=0.000001)
        self.assertAlmostEqual(matrix_used[2][2], 0, delta=0.000001)


    @patch("atomium.structures.Atom.transform_atoms")
    def test_axis_must_be_valid(self, mock_trans):
        atoms = [Mock(), Mock()]
        with self.assertRaises(ValueError):
            Atom.rotate_atoms(math.pi / 2, "q", atoms[0], atoms[1])



class AtomElementPropertyTests(TestCase):

    def test_can_get_element(self):
        atom = Atom("C", 1, 2, 3, 10, "CA", -1, 0.76, [1, 2, 3, 4, 5, 6])
        self.assertIs(atom._element, atom.element)


    def test_can_set_element(self):
        atom = Atom("C", 1, 2, 3, 10, "CA", -1, 0.76, [1, 2, 3, 4, 5, 6])
        atom.element = "VALUE"
        self.assertEqual(atom._element, "VALUE")



class AtomXPropertyTests(TestCase):

    def test_can_get_x(self):
        atom = Atom("C", 1, 2, 3, 10, "CA", -1, 0.76, [1, 2, 3, 4, 5, 6])
        self.assertIs(atom._x, atom.x)


    def test_can_set_x(self):
        atom = Atom("C", 1, 2, 3, 10, "CA", -1, 0.76, [1, 2, 3, 4, 5, 6])
        atom.x = "VALUE"
        self.assertEqual(atom._x, "VALUE")



class AtomYPropertyTests(TestCase):

    def test_can_get_y(self):
        atom = Atom("C", 1, 2, 3, 10, "CA", -1, 0.76, [1, 2, 3, 4, 5, 6])
        self.assertIs(atom._y, atom.y)


    def test_can_set_y(self):
        atom = Atom("C", 1, 2, 3, 10, "CA", -1, 0.76, [1, 2, 3, 4, 5, 6])
        atom.y = "VALUE"
        self.assertEqual(atom._y, "VALUE")



class AtomZPropertyTests(TestCase):

    def test_can_get_z(self):
        atom = Atom("C", 1, 2, 3, 10, "CA", -1, 0.76, [1, 2, 3, 4, 5, 6])
        self.assertIs(atom._z, atom.z)


    def test_can_set_z(self):
        atom = Atom("C", 1, 2, 3, 10, "CA", -1, 0.76, [1, 2, 3, 4, 5, 6])
        atom.z = "VALUE"
        self.assertEqual(atom._z, "VALUE")



class AtomIdPropertyTests(TestCase):

    def test_can_get_id(self):
        atom = Atom("C", 1, 2, 3, 10, "CA", -1, 0.76, [1, 2, 3, 4, 5, 6])
        self.assertIs(atom._id, atom.id)



class AtomNamePropertyTests(TestCase):

    def test_can_get_name(self):
        atom = Atom("C", 1, 2, 3, 10, "CA", -1, 0.76, [1, 2, 3, 4, 5, 6])
        self.assertIs(atom._name, atom.name)


    def test_can_set_name(self):
        atom = Atom("C", 1, 2, 3, 10, "CA", -1, 0.76, [1, 2, 3, 4, 5, 6])
        atom.name = "VALUE"
        self.assertEqual(atom._name, "VALUE")



class AtomChargePropertyTests(TestCase):

    def test_can_get_charge(self):
        atom = Atom("C", 1, 2, 3, 10, "CA", -1, 0.76, [1, 2, 3, 4, 5, 6])
        self.assertIs(atom._charge, atom.charge)


    def test_can_set_charge(self):
        atom = Atom("C", 1, 2, 3, 10, "CA", -1, 0.76, [1, 2, 3, 4, 5, 6])
        atom.charge = "VALUE"
        self.assertEqual(atom._charge, "VALUE")



class AtomBvaluePropertyTests(TestCase):

    def test_can_get_bvalue(self):
        atom = Atom("C", 1, 2, 3, 10, "CA", -1, 0.76, [1, 2, 3, 4, 5, 6])
        self.assertIs(atom._bvalue, atom.bvalue)


    def test_can_set_bvalue(self):
        atom = Atom("C", 1, 2, 3, 10, "CA", -1, 0.76, [1, 2, 3, 4, 5, 6])
        atom.bvalue = "VALUE"
        self.assertEqual(atom._bvalue, "VALUE")



class AtomAnisotropyPropertyTests(TestCase):

    def test_can_get_anisotropy(self):
        atom = Atom("C", 1, 2, 3, 10, "CA", -1, 0.76, [1, 2, 3, 4, 5, 6])
        self.assertIs(atom._anisotropy, atom.anisotropy)



class AtomLocationTests(TestCase):

    def test_can_get_location(self):
        atom = Atom("C", 1, 2, 3, 10, "CA", -1, 0.76, [1, 2, 3, 4, 5, 6])
        self.assertEqual(atom.location, (1, 2, 3))



class AtomDistanceToTests(TestCase):

    def test_can_get_distance_between_atoms(self):
        atom1 = Atom("C", 4, 8, 3, 10, "CA", -1, 0.76, [1, 2, 3, 4, 5, 6])
        atom2 = Mock(Atom, location = (2, 3, 5))
        self.assertAlmostEqual(atom1.distance_to(atom2), 5.744, delta=0.001)


    def test_atom_distance_can_be_zero(self):
        atom1 = Atom("C", 4, 8, 3, 10, "CA", -1, 0.76, [1, 2, 3, 4, 5, 6])
        atom2 = Mock(Atom, location = (4, 8, 3))
        self.assertEqual(atom1.distance_to(atom2), 0)


    def test_can_get_distance_to_xyz_tuple(self):
        atom1 = Atom("C", 4, 8, 3, 10, "CA", -1, 0.76, [1, 2, 3, 4, 5, 6])
        atom2 = (2, 3, 5)
        self.assertAlmostEqual(atom1.distance_to(atom2), 5.744, delta=0.001)



class AtomMassTests(TestCase):

    def test_known_element_mass(self):
        atom = Atom("C", 4, 8, 3, 10, "CA", -1, 0.76, [1, 2, 3, 4, 5, 6])
        self.assertAlmostEqual(atom.mass, 12, delta=0.1)
        atom._element = "H"
        self.assertAlmostEqual(atom.mass, 1, delta=0.1)


    def test_atom_mass_case_insensitive(self):
        atom = Atom("he", 4, 8, 3, 10, -1, "H", 0.76, [1, 2, 3, 4, 5, 6])
        self.assertAlmostEqual(atom.mass, 4, delta=0.1)
        atom = Atom("He", 4, 8, 3, 10, -1, "H", 0.76, [1, 2, 3, 4, 5, 6])
        self.assertAlmostEqual(atom.mass, 4, delta=0.1)


    def test_unknown_atom_mass(self):
        atom = Atom("XX", 4, 8, 3, 10, "he", -1, 0.76, [1, 2, 3, 4, 5, 6])
        self.assertEqual(atom.mass, 0)



class AtomMetalTests(TestCase):

    def test_some_atoms_are_metals(self):
        self.assertFalse(Atom("C", 4, 8, 3, 10, "CA", -1, 0.76, [1,]).is_metal)
        self.assertFalse(Atom("c", 2, 3, 5, 10, "CA", -1, 0.76, [1,]).is_metal)
        self.assertTrue(Atom("FE", 2, 3, 5, 10, "CA", -1, 0.76, [1,]).is_metal)
        self.assertTrue(Atom("zn", 2, 3, 5, 10, "CA", -1, 0.76, [1,]).is_metal)



class AtomStructureTests(TestCase):

    def test_can_get_structure(self):
        atom = Atom("C", 4, 8, 3, 10, "CA", -1, 0.76, [1, 2, 3, 4, 5, 6])
        atom._structure = 1000
        self.assertIs(atom.structure, atom._structure)



class AtomChainTests(TestCase):

    def test_can_get_chain(self):
        atom = Atom("C", 4, 8, 3, 10, "CA", -1, 0.76, [1, 2, 3, 4, 5, 6])
        atom._structure = Mock(chain=1000)
        self.assertEqual(atom.chain, 1000)


    def test_can_get_chain_no_structure(self):
        atom = Atom("C", 4, 8, 3, 10, "CA", -1, 0.76, [1, 2, 3, 4, 5, 6])
        self.assertEqual(atom.chain, None)



class AtomModelTests(TestCase):

    @patch("atomium.structures.Atom.chain", new_callable=PropertyMock)
    def test_can_get_model(self, mock_chain):
        mock_chain.return_value = Mock(model=1000)
        atom = Atom("C", 4, 8, 3, 10, "CA", -1, 0.76, [1, 2, 3, 4, 5, 6])
        self.assertEqual(atom.model, 1000)


    @patch("atomium.structures.Atom.chain", new_callable=PropertyMock)
    def test_can_get_model_no_structure(self, mock_chain):
        mock_chain.return_value = None
        atom = Atom("C", 4, 8, 3, 10, "CA", -1, 0.76, [1, 2, 3, 4, 5, 6])
        self.assertEqual(atom.model, None)



class AtomAngleTests(TestCase):

    @patch("atomium.structures.Atom.location", new_callable=PropertyMock)
    def test_angle_between_atoms(self, mock_loc):
        mock_loc.return_value = (1, 1, 1)
        atom = Atom("C", 4, 8, 3, 10, "CA", -1, 0.76, [1, 2, 3, 4, 5, 6])
        atom1 = Mock(location=(2, 2, 2))
        atom2 = Mock(location=(-1, -1, -1))
        self.assertEqual(atom.angle(atom1, atom2), math.pi)


    @patch("atomium.structures.Atom.location", new_callable=PropertyMock)
    def test_angle_between_atoms_when_superimposed(self, mock_loc):
        mock_loc.return_value = (1, 1, 1)
        atom = Atom("C", 4, 8, 3, 10, "CA", -1, 0.76, [1, 2, 3, 4, 5, 6])
        atom1 = Mock(location=(2, 2, 2))
        atom2 = Mock(location=(1, 1, 1))
        self.assertEqual(atom.angle(atom1, atom2), 0)



class AtomNearbyAtomTests(TestCase):

    @patch("atomium.structures.Atom.model", new_callable=PropertyMock)
    @patch("atomium.structures.Atom.location", new_callable=PropertyMock)
    def test_can_get_nearby_atoms(self, mock_loc, mock_model):
        mock_loc.return_value = 90
        model = Mock()
        atom = Atom("C", 4, 8, 3, 10, "CA", -1, 0.76, [1, 2, 3, 4, 5, 6])
        mock_model.return_value = model
        model.atoms_in_sphere.return_value = {1, 5, 7, atom}
        self.assertEqual(atom.nearby_atoms(1, 4, b=2), {1, 5, 7})
        model.atoms_in_sphere.assert_called_with(90, 1, 4, b=2)


    @patch("atomium.structures.Atom.model", new_callable=PropertyMock)
    def test_can_get_nearby_atoms_no_model(self, mock_model):
        atom = Atom("C", 4, 8, 3, 10, "CA", -1, 0.76, [1, 2, 3, 4, 5, 6])
        mock_model.return_value = None
        self.assertEqual(atom.nearby_atoms(1, 4, b=2), set())



class AtomNearbyStructuresTests(TestCase):

    @patch("atomium.structures.Atom.nearby_atoms")
    @patch("atomium.structures.Atom.structure", new_callable=PropertyMock)
    def test_can_get_nearby_structures(self, mock_struct, mock_atoms):
        residues = [Mock(Residue), Mock(Residue)]
        ligands = [Mock(Ligand, water=False), Mock(Ligand, water=False)]
        waters = [Mock(Ligand, water=True), Mock(Ligand, water=True)]
        mock_struct.return_value = residues[0]
        mock_atoms.return_value = [
         Mock(structure=residues[0]), Mock(structure=residues[1]),
         Mock(structure=ligands[0]), Mock(structure=ligands[1]),
         Mock(structure=waters[0]), Mock(structure=waters[1]),
        ]
        atom = Atom("C", 4, 8, 3, 10, "CA", -1, 0.76, [1, 2, 3, 4, 5, 6])
        self.assertEqual(atom.nearby_structures(1, 2, a=9), {
         residues[1], ligands[0], ligands[1]
        })
        mock_atoms.assert_called_with(1, 2, a=9)


    @patch("atomium.structures.Atom.nearby_atoms")
    @patch("atomium.structures.Atom.structure", new_callable=PropertyMock)
    def test_can_get_nearby_structures_kwargs(self, mock_struct, mock_atoms):
        residues = [Mock(Residue), Mock(Residue)]
        ligands = [Mock(Ligand, water=False), Mock(Ligand, water=False)]
        waters = [Mock(Ligand, water=True), Mock(Ligand, water=True)]
        mock_struct.return_value = residues[0]
        mock_atoms.return_value = [
         Mock(structure=residues[0]), Mock(structure=residues[1]),
         Mock(structure=ligands[0]), Mock(structure=ligands[1]),
         Mock(structure=waters[0]), Mock(structure=waters[1]),
        ]
        atom = Atom("C", 4, 8, 3, 10, "CA", -1, 0.76, [1, 2, 3, 4, 5, 6])
        self.assertEqual(atom.nearby_structures(1, 2, a=9, residues=False, ligands=False, waters=True), {
         waters[0], waters[1]
        })
        mock_atoms.assert_called_with(1, 2, a=9)



class AtomTranslationTests(TestCase):

    @patch("atomium.structures.Atom.translate_atoms")
    @patch("atomium.structures.Atom.trim")
    def test_can_translate_atom(self, mock_trim, mock_trans):
        atom = Atom("C", 4, 8, 3, 10, "CA", -1, 0.76, [1, 2, 3, 4, 5, 6])
        atom.translate(1, 2, 3)
        mock_trans.assert_called_with((1, 2, 3), atom)
        mock_trim.assert_called_with(12)


    @patch("atomium.structures.Atom.translate_atoms")
    @patch("atomium.structures.Atom.trim")
    def test_can_translate_atom_with_vector(self, mock_trim, mock_trans):
        atom = Atom("C", 4, 8, 3, 10, "CA", -1, 0.76, [1, 2, 3, 4, 5, 6])
        atom.translate((1, 2, 3), trim=9)
        mock_trans.assert_called_with((1, 2, 3), atom)
        mock_trim.assert_called_with(9)



class AtomTransformationTests(TestCase):

    @patch("atomium.structures.Atom.transform_atoms")
    @patch("atomium.structures.Atom.trim")
    def test_can_transform_atom(self, mock_trim, mock_trans):
        atom = Atom("C", 4, 8, 3, 10, "CA", -1, 0.76, [1, 2, 3, 4, 5, 6])
        atom.transform([1], trim=9)
        mock_trans.assert_called_with([1], atom)
        mock_trim.assert_called_with(9)



class AtomTransformationTests(TestCase):

    @patch("atomium.structures.Atom.rotate_atoms")
    @patch("atomium.structures.Atom.trim")
    def test_can_rotate_atom(self, mock_trim, mock_rot):
        atom = Atom("C", 4, 8, 3, 10, "CA", -1, 0.76, [1, 2, 3, 4, 5, 6])
        atom.rotate(1, 2, trim=9)
        mock_rot.assert_called_with(1, 2, atom)
        mock_trim.assert_called_with(9)



class AtomMovingTests(TestCase):

    def test_can_move_atom(self):
        atom = Atom("C", 20, 30, 40, 10, "CA", -1, 0.76, [1,])
        atom.move_to(1, 2, 3)
        self.assertEqual(atom._x, 1)
        self.assertEqual(atom._y, 2)
        self.assertEqual(atom._z, 3)



class AtomTrimmingTests(TestCase):

    def test_can_not_round_atom_location(self):
        atom = Atom("C", 0.1111, 0.4444, 0.8888, 10, "CA", -1, 0.76, [1,])
        atom.trim(None)
        self.assertEqual((atom._x, atom._y, atom._z), (0.1111, 0.4444, 0.8888))


    def test_can_round_atom_location(self):
        atom = Atom("C", 0.1111, 0.4444, 0.8888, 10, "CA", -1, 0.76, [1,])
        atom.trim(3)
        self.assertEqual((atom._x, atom._y, atom._z), (0.111, 0.444, 0.889))
        atom.trim(2)
        self.assertEqual((atom._x, atom._y, atom._z), (0.11, 0.44, 0.89))
        atom.trim(1)
        self.assertEqual((atom._x, atom._y, atom._z), (0.1, 0.4, 0.9))
        atom.trim(0)
        self.assertEqual((atom._x, atom._y, atom._z), (0, 0, 1))



class AtomEquivalenceTests(TestCase):

    @patch("atomium.structures.Atom.location", new_callable=PropertyMock)
    def test_atoms_equivalent(self, mock_loc):
        mock_loc.return_value = (1, 2, 3)
        atom = Atom("C", 20, 30, 40, 10, "CA", -1, 0.76, [1,])
        other = Mock(Atom, _element="C", _id=10, _name="CA", location=(1, 2, 3),
         _charge=-1, _bvalue=0.76, _anisotropy=[1,])
        self.assertTrue(atom.equivalent_to(other))
        other._bvalue = 0.75
        self.assertFalse(atom.equivalent_to(other))



class AtomCopyingTests(TestCase):

    def test_can_copy_atom(self):
        atom = Atom("C", 4, 8, 3, 10, "CA", -1, 0.76, [1, 2, 3, 4, 5, 6])
        p = patch("atomium.structures.Atom")
        mock_atom = p.start()
        try:
            copy = atom.copy()
            mock_atom.assert_called_with(
             "C", 4, 8, 3, 10, "CA", -1, 0.76, [1, 2, 3, 4, 5, 6]
            )
        finally: p.stop()
