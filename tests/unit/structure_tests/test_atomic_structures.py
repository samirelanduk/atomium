from collections import Counter
from math import pi
from unittest import TestCase
from unittest.mock import Mock, patch
from atomium.structures.atoms import Atom
from atomium.structures.molecules import AtomicStructure, Residue

class AtomicStructureTest(TestCase):

    def setUp(self):
        self.atom1 = Mock(Atom)
        self.atom2 = Mock(Atom)
        self.atom3 = Mock(Atom)
        self.atom1.mass.return_value = 8
        self.atom2.mass.return_value = 5.1
        self.atom3.mass.return_value = 4
        self.atom1.element.return_value = "A"
        self.atom2.element.return_value = "B"
        self.atom3.element.return_value = "B"
        self.atom1.atom_id.return_value = 500
        self.atom2.atom_id.return_value = None
        self.atom3.atom_id.return_value = 700
        self.atom1.name.return_value = "CA"
        self.atom2.name.return_value = "CA"
        self.atom3.name.return_value = "NY"
        self.atoms = [self.atom1, self.atom2, self.atom3]



class AtomicStructureCreationTests(AtomicStructureTest):

    def test_can_create_atomic_structure_with_id_atoms(self):
        self.atom2.atom_id.return_value = 600
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(structure._id_atoms, {
         500: self.atom1, 600: self.atom2, 700: self.atom3
        })
        self.assertEqual(structure._atoms, set())


    def test_can_create_atomic_structure_with_no_id_atoms(self):
        self.atom1.atom_id.return_value = None
        self.atom3.atom_id.return_value = None
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(structure._id_atoms, {})
        self.assertEqual(structure._atoms, set(self.atoms))


    def test_can_create_atomic_structure_with_mixed_id_atoms(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(structure._atoms, set([self.atoms[1]]))
        self.assertEqual(structure._id_atoms, {500: self.atom1, 700: self.atom3})


    def test_atomic_structure_needs_atoms(self):
        with self.assertRaises(TypeError):
            AtomicStructure("self.atom1", self.atom2, self.atom3)


    def test_atomic_structure_will_accept_atomic_structures(self):
        structure = Mock(AtomicStructure)
        structure.atoms.return_value = set(self.atoms[1:])
        structure2 = AtomicStructure(self.atom1, structure)
        self.assertEqual(structure2._atoms, set([self.atoms[1]]))
        self.assertEqual(structure2._id_atoms, {500: self.atom1, 700: self.atom3})



class AtomicStructureReprTests(AtomicStructureTest):

    def test_atomic_structure_repr(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(str(structure), "<AtomicStructure (3 atoms)>")



class AtomicStructureContainerTests(AtomicStructureTest):

    def test_atomic_structure_is_container_of_its_atoms(self):
        structure = AtomicStructure(self.atom1, self.atom2)
        self.assertIn(self.atom1, structure)
        self.assertIn(self.atom2, structure)
        self.assertNotIn(self.atom3, structure)



class AtomicStructureAtomsTests(AtomicStructureTest):

    def test_can_get_all_atoms(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(structure.atoms(), set(self.atoms))


    def test_can_get_atoms_by_element(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(structure.atoms(element="A"), set(self.atoms[:1]))
        self.assertEqual(structure.atoms(element="b"), set(self.atoms[1:]))
        self.assertEqual(structure.atoms(element="C"), set())


    def test_can_exlcude_atoms_by_element(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(structure.atoms(exclude="A"), set(self.atoms[1:]))
        self.assertEqual(structure.atoms(exclude="B"), set(self.atoms[:1]))
        self.assertEqual(structure.atoms(exclude="C"), set(self.atoms))


    def test_can_get_atoms_by_id(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(structure.atoms(atom_id=500), set([self.atoms[0]]))
        self.assertEqual(structure.atoms(atom_id=700), set([self.atoms[2]]))
        self.assertEqual(structure.atoms(atom_id=300), set())


    def test_can_get_atoms_by_name(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(structure.atoms(name="CA"), set(self.atoms[:2]))
        self.assertEqual(structure.atoms(name="NY"), set([self.atoms[2]]))
        self.assertEqual(structure.atoms(name="CB"), set())



class AtomicStructureAtomTest(AtomicStructureTest):

    @patch("atomium.structures.molecules.AtomicStructure.atoms")
    def test_atom_calls_atoms(self, mock_atoms):
        mock_atoms.return_value = set([self.atoms[0]])
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        atom = structure.atom(fufu=500, element="A")
        mock_atoms.assert_called_with(fufu=500, element="A")
        self.assertIs(atom, self.atom1)


    @patch("atomium.structures.molecules.AtomicStructure.atoms")
    def test_atom_can_return_none(self, mock_atoms):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        mock_atoms.return_value = set()
        self.assertIs(structure.atom(element="C"), None)


    @patch("atomium.structures.molecules.AtomicStructure.atoms")
    def test_atom_can_get_atom_by_id(self, mock_atoms):
        mock_atoms.return_value = set([self.atoms[1]])
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        atom = structure.atom(atom_id=500)
        self.assertFalse(mock_atoms.called)
        self.assertIs(atom, self.atom1)



class AtomicStructureAtomAdditionTests(AtomicStructureTest):

    def test_can_add_atoms_to_structure(self):
        structure = AtomicStructure(self.atom1)
        self.assertEqual(structure._atoms, set())
        self.assertEqual(structure._id_atoms, {500: self.atom1})
        structure.add_atom(self.atom2)
        self.assertEqual(structure._atoms, set([self.atom2]))
        self.assertEqual(structure._id_atoms, {500: self.atom1})
        structure.add_atom(self.atom3)
        self.assertEqual(structure._atoms, set([self.atom2]))
        self.assertEqual(structure._id_atoms, {500: self.atom1, 700: self.atom3})


    def test_can_only_add_atoms(self):
        structure = AtomicStructure(self.atom1, self.atom2)
        with self.assertRaises(TypeError):
            structure.add_atom("atom")



class AtomicStructureAtomRemovalTests(AtomicStructureTest):

    def test_can_remove_atoms_from_structure(self):
        structure = AtomicStructure(self.atom1, self.atom2)
        self.assertEqual(structure._atoms, set([self.atom2]))
        self.assertEqual(structure._id_atoms, {500: self.atom1})
        structure.remove_atom(self.atom1)
        self.assertEqual(structure._atoms, set([self.atom2]))
        self.assertEqual(structure._id_atoms, {})
        structure.remove_atom(self.atom2)
        self.assertEqual(structure._atoms, set())
        self.assertEqual(structure._id_atoms, {})



class AtomicStructurePairwiseAtomTests(AtomicStructureTest):

    @patch("atomium.structures.molecules.AtomicStructure.atoms")
    def test_can_get_pairwise_atoms(self, mock_atoms):
        mock_atoms.return_value = set(self.atoms)
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        collected_atoms = []
        for pair in structure.pairwise_atoms():
            self.assertIsNot(pair[0], pair[1])
            collected_atoms.append(pair[0])
            collected_atoms.append(pair[1])
        for atom in self.atoms:
            self.assertEqual(collected_atoms.count(atom), 2)


    @patch("atomium.structures.molecules.AtomicStructure.atoms")
    def test_can_get_pairwise_atoms_one(self, mock_atoms):
        mock_atoms.return_value = set(self.atoms[:1])
        structure = AtomicStructure(self.atom1)
        self.assertEqual(len(list(structure.pairwise_atoms())), 0)



class AtomicStructureMassTests(AtomicStructureTest):

    @patch("atomium.structures.molecules.AtomicStructure.atoms")
    def test_structure_mass_is_sum_of_atom_masses(self, mock_atoms):
        mock_atoms.return_value = set(self.atoms)
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(structure.mass(), 17.1)



class AtomicStructureChargeTests(AtomicStructureTest):

    @patch("atomium.structures.molecules.AtomicStructure.atoms")
    def test_structure_charge_is_sum_of_atom_charges(self, mock_atoms):
        mock_atoms.return_value = set(self.atoms)
        self.atom1.charge.return_value = 0.2
        self.atom2.charge.return_value = -1.4
        self.atom3.charge.return_value = 0.6
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertAlmostEqual(structure.charge(), -0.6, delta=0.000005)



class AtomicStructureFormulaTests(AtomicStructureTest):

    @patch("atomium.structures.molecules.AtomicStructure.atoms")
    def test_can_get_formula(self, mock_atoms):
        mock_atoms.return_value = set(self.atoms)
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(structure.formula(), Counter({"A":1, "B":2}))



class AtomicStructureRoundingTests(AtomicStructureTest):

    @patch("atomium.structures.molecules.AtomicStructure.atoms")
    def test_can_round_atom_locations(self, mock_atoms):
        mock_atoms.return_value = set(self.atoms)
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        structure.round(10)
        self.atom1.round.assert_called_with(10)
        self.atom2.round.assert_called_with(10)
        self.atom3.round.assert_called_with(10)



class AtomicStructureOrientationTests(AtomicStructureTest):

    def setUp(self):
        AtomicStructureTest.setUp(self)
        self.patch1 = patch("atomium.structures.molecules.AtomicStructure.translate")
        self.patch2 = patch("atomium.structures.molecules.AtomicStructure.rotate")
        self.mock_translate = self.patch1.start()
        self.mock_rotate = self.patch2.start()
        self.atom1.location.return_value = (1, 2, 3)
        self.atom2.location.return_value = (1, 1, 1)
        self.atom3.location.return_value = (7, 8, 9)


    def tearDown(self):
        self.patch1.stop()
        self.patch2.stop()


    def test_can_orient_to_origin(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        structure.orient(self.atom1)
        self.mock_translate.assert_called_with(-1, -2, -3)


    def test_atom1_must_be_atom(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        with self.assertRaises(TypeError):
            structure.orient("atom")


    def test_can_orient_onto_x_axis(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        structure.orient(self.atom1, atom2=self.atom2, axis="x")
        self.mock_translate.assert_called_with(-1, -2, -3)
        first_rotation_call = self.mock_rotate.call_args_list[0][0]
        self.assertAlmostEqual(first_rotation_call[0], pi / 4, delta=0.00005)
        self.assertEqual(first_rotation_call[1], "y")
        second_rotation_call = self.mock_rotate.call_args_list[1][0]
        self.assertAlmostEqual(second_rotation_call[0], -0.9553166, delta=0.005)
        self.assertEqual(second_rotation_call[1], "z")


    def test_can_orient_onto_y_axis(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        structure.orient(self.atom1, atom2=self.atom2, axis="y")
        self.mock_translate.assert_called_with(-1, -2, -3)
        first_rotation_call = self.mock_rotate.call_args_list[0][0]
        self.assertAlmostEqual(first_rotation_call[0], pi / 4, delta=0.00005)
        self.assertEqual(first_rotation_call[1], "z")
        second_rotation_call = self.mock_rotate.call_args_list[1][0]
        self.assertAlmostEqual(second_rotation_call[0], -0.9553166, delta=0.005)
        self.assertEqual(second_rotation_call[1], "x")


    def test_can_orient_onto_z_axis(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        structure.orient(self.atom1, atom2=self.atom2, axis="z")
        self.mock_translate.assert_called_with(-1, -2, -3)
        first_rotation_call = self.mock_rotate.call_args_list[0][0]
        self.assertAlmostEqual(first_rotation_call[0], pi / 4, delta=0.00005)
        self.assertEqual(first_rotation_call[1], "x")
        second_rotation_call = self.mock_rotate.call_args_list[1][0]
        self.assertAlmostEqual(second_rotation_call[0], -0.9553166, delta=0.005)
        self.assertEqual(second_rotation_call[1], "y")


    def test_atom2_must_be_atom(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        with self.assertRaises(TypeError):
            structure.orient(self.atom2, atom2="atom")


    def test_axis_must_be_valid(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        with self.assertRaises(ValueError) as e:
            structure.orient(self.atom1, atom2=self.atom2, axis="7")
        self.assertIn("axis", str(e.exception))


    def test_can_rotate_onto_plane(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        structure.orient(
         self.atom1, atom2=self.atom2, axis="z", atom3=self.atom3, plane="xz"
        )
        self.mock_translate.assert_called_with(-1, -2, -3)
        first_rotation_call = self.mock_rotate.call_args_list[0][0]
        self.assertAlmostEqual(first_rotation_call[0], pi / 4, delta=0.00005)
        self.assertEqual(first_rotation_call[1], "x")
        second_rotation_call = self.mock_rotate.call_args_list[1][0]
        self.assertAlmostEqual(second_rotation_call[0], -0.9553166, delta=0.005)
        self.assertEqual(second_rotation_call[1], "y")
        third_rotation_call = self.mock_rotate.call_args_list[2][0]
        self.assertAlmostEqual(third_rotation_call[0], -0.8519663, delta=0.005)
        self.assertEqual(third_rotation_call[1], "z")


    def test_atom3_must_be_atom(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        with self.assertRaises(TypeError):
            structure.orient(self.atom1, atom2=self.atom2, atom3="atom3")


    def test_plane_must_be_valid(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        with self.assertRaises(ValueError) as e:
            structure.orient(self.atom1, atom2=self.atom2, atom3=self.atom3, plane="dtg")
        self.assertIn("plane", str(e.exception))
        with self.assertRaises(ValueError) as e:
            structure.orient(self.atom1, atom2=self.atom2, atom3=self.atom3, plane="xyz")
        self.assertIn("plane", str(e.exception))
        with self.assertRaises(ValueError) as e:
            structure.orient(self.atom1, atom2=self.atom2, atom3=self.atom3, plane="yz")
        self.assertIn("x-axis", str(e.exception))




class AtomicStructureTranslationTests(AtomicStructureTest):

    @patch("atomium.structures.molecules.AtomicStructure.atoms")
    def test_structure_translation(self, mock_atoms):
        mock_atoms.return_value = set(self.atoms)
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        structure.translate(5, 4, -2)
        for atom in self.atoms:
            atom.translate.assert_called_with(5, 4, -2)



class AtomicStructureRotationTests(AtomicStructureTest):

    @patch("atomium.structures.molecules.AtomicStructure.atoms")
    def test_structure_rotation(self, mock_atoms):
        mock_atoms.return_value = set(self.atoms)
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        structure.rotate(0.1, "z")
        for atom in self.atoms:
            atom.rotate.assert_called_with(0.1, "z")



class AtomicStructureCenterOfMassTests(AtomicStructureTest):

    @patch("atomium.structures.molecules.AtomicStructure.atoms")
    @patch("atomium.structures.molecules.AtomicStructure.mass")
    def test_can_get_center_of_mass_when_equal_mass(self, mock_mass, mock_atoms):
        mock_atoms.return_value = set(self.atoms[:2])
        mock_mass.return_value = 20
        self.atom1._x, self.atom1._y, self.atom1._z = 0, 0, 0
        self.atom2._x, self.atom2._y, self.atom2._z = 1, 1, 1
        self.atom1.mass.return_value = 10
        self.atom2.mass.return_value = 10
        structure = AtomicStructure(self.atom1, self.atom2)
        self.assertEqual(structure.center_of_mass(), (0.5, 0.5, 0.5))


    @patch("atomium.structures.molecules.AtomicStructure.atoms")
    @patch("atomium.structures.molecules.AtomicStructure.mass")
    def test_can_get_center_of_mass_when_unequal_mass(self, mock_mass, mock_atoms):
        mock_atoms.return_value = set(self.atoms[:2])
        mock_mass.return_value = 40
        self.atom1._x, self.atom1._y, self.atom1._z = 0, 0, 0
        self.atom2._x, self.atom2._y, self.atom2._z = 1, 1, 1
        self.atom1.mass.return_value = 10
        self.atom2.mass.return_value = 30
        structure = AtomicStructure(self.atom1, self.atom2)
        self.assertEqual(structure.center_of_mass(), (0.75, 0.75, 0.75))



class AtomicStructureRadiusOfGyrationTests(AtomicStructureTest):

    @patch("atomium.structures.molecules.AtomicStructure.atoms")
    @patch("atomium.structures.molecules.AtomicStructure.center_of_mass")
    def test_can_get_radius_of_gyration(self, mock_center, mock_atoms):
        mock_atoms.return_value = set(self.atoms[:2])
        mock_center.return_value = (5, 0, 0)
        self.atom1.distance_to.return_value = 5
        self.atom2.distance_to.return_value = 5
        self.atom1._x, self.atom1._y, self.atom1._z = 0, 0, 0
        self.atom2._x, self.atom2._y, self.atom2._z = 10, 0, 0
        structure = AtomicStructure(self.atom1, self.atom2)
        self.assertEqual(structure.radius_of_gyration(), 5)
        self.atom1.distance_to.assert_called_with((5, 0, 0))
        self.atom2.distance_to.assert_called_with((5, 0, 0))



class AtomicStructurePatternMatchingTests(AtomicStructureTest):

    def test_pattern_must_be_structure(self):
        structure = AtomicStructure(self.atom1, self.atom2)
        with self.assertRaises(TypeError):
            structure.find_pattern("pattern")



class AtomicStructureCopyTests(AtomicStructureTest):

    def test_can_create_copy_of_atomic_structure(self):
        new_atoms = [Mock(Atom), Mock(Atom)]
        self.atom1.copy.return_value, self.atom2.copy.return_value = new_atoms
        new_atoms[0].atom_id.return_value = None
        new_atoms[1].atom_id.return_value = None
        structure = AtomicStructure(self.atom1, self.atom2)
        copy = structure.copy()
        self.assertEqual(copy._atoms, set(new_atoms))



class AtomicStructureToStringTests(AtomicStructureTest):

    @patch("atomium.files.pdb2pdbdict.structure_to_pdb_dict")
    @patch("atomium.files.pdbdict2pdbstring.pdb_dict_to_pdb_string")
    def test_can_save_as_pdb_string_with_description(self, mock_string, mock_dict):
        structure = AtomicStructure(*self.atoms)
        pdb_dict = {}
        mock_string.return_value = "filecontents"
        mock_dict.return_value = pdb_dict
        s = structure.to_file_string("pdb", description="TTT")
        mock_dict.assert_called_with(structure)
        mock_string.assert_called_with({"title": "TTT"})
        self.assertEqual(s, "filecontents")


    @patch("atomium.files.xyz2xyzdict.structure_to_xyz_dict")
    @patch("atomium.files.xyzdict2xyzstring.xyz_dict_to_xyz_string")
    def test_can_save_as_xyz_string_with_description(self, mock_string, mock_dict):
        structure = AtomicStructure(*self.atoms)
        xyz_dict = {}
        mock_string.return_value = "filecontents"
        mock_dict.return_value = xyz_dict
        s = structure.to_file_string("xyz", description="DDD")
        mock_dict.assert_called_with(structure)
        mock_string.assert_called_with({"title": "DDD"})
        self.assertEqual(s, "filecontents")


    def test_invalid_file_format_is_error(self):
        structure = AtomicStructure(*self.atoms)
        with self.assertRaises(ValueError):
            structure.to_file_string("nosuchfile")



class AtomicStructureSavingTests(AtomicStructureTest):

    @patch("atomium.structures.molecules.AtomicStructure.to_file_string")
    @patch("atomium.files.utilities.string_to_file")
    def test_saving_uses_correct_functions(self, mock_save, mock_string):
        mock_string.return_value = "filestring"
        structure = AtomicStructure(*self.atoms)
        structure.save("file.xyz", "a description")
        mock_string.assert_called_with("xyz", "a description")
        mock_save.assert_called_with("filestring", "file.xyz")


    @patch("atomium.structures.molecules.AtomicStructure.to_file_string")
    @patch("atomium.files.utilities.string_to_file")
    def test_file_format_extracted(self, mock_save, mock_string):
        mock_string.return_value = "filestring"
        structure = AtomicStructure(*self.atoms)
        structure.save("path/to/file.dfghdfg", "a description")
        mock_string.assert_called_with("dfghdfg", "a description")
        mock_save.assert_called_with("filestring", "path/to/file.dfghdfg")
