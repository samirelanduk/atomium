from collections import Counter
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
        self.atom2.atom_id.return_value = 600
        self.atom3.atom_id.return_value = 700
        self.atom1.name.return_value = "CA"
        self.atom2.name.return_value = "CA"
        self.atom3.name.return_value = "NY"
        self.atoms = [self.atom1, self.atom2, self.atom3]



class AtomicStructureCreationTests(AtomicStructureTest):

    def test_can_create_atomic_structure(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(structure._atoms, set(self.atoms))


    def test_atomic_structure_needs_atoms(self):
        with self.assertRaises(TypeError):
            AtomicStructure("self.atom1", self.atom2, self.atom3)



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
        self.assertEqual(structure.atoms(), structure._atoms)
        self.assertIsNot(structure.atoms(), structure._atoms)


    def test_can_get_atoms_by_element(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(structure.atoms(element="A"), set(self.atoms[:1]))
        self.assertEqual(structure.atoms(element="B"), set(self.atoms[1:]))
        self.assertEqual(structure.atoms(element="C"), set())


    def test_can_get_atoms_by_id(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(structure.atoms(atom_id=500), set([self.atoms[0]]))
        self.assertEqual(structure.atoms(atom_id=600), set([self.atoms[1]]))
        self.assertEqual(structure.atoms(atom_id=300), set())


    def test_can_get_atoms_by_name(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(structure.atoms(name="CA"), set(self.atoms[:2]))
        self.assertEqual(structure.atoms(name="NY"), set([self.atoms[2]]))
        self.assertEqual(structure.atoms(name="CB"), set())



class AtomicStructureAtomTest(AtomicStructureTest):

    @patch("atomium.structures.molecules.AtomicStructure.atoms")
    def test_atom_calls_atoms(self, mock_atoms):
        mock_atoms.return_value = set(self.atoms[:1])
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        atom = structure.atom(element="A")
        mock_atoms.assert_called_with(element="A")
        self.assertIs(atom, self.atom1)


    @patch("atomium.structures.molecules.AtomicStructure.atoms")
    def test_atom_can_return_none(self, mock_atoms):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        mock_atoms.return_value = set()
        self.assertIs(structure.atom(element="C"), None)


    @patch("atomium.structures.molecules.AtomicStructure.atoms")
    def test_atom_can_get_atom_by_id_and_element(self, mock_atoms):
        mock_atoms.return_value = set([self.atoms[0]])
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        atom = structure.atom(atom_id=500, element="A")
        mock_atoms.assert_called_with(atom_id=500, element="A")
        self.assertIs(atom, self.atom1)



class AtomicStructureAtomAdditionTests(AtomicStructureTest):

    def test_can_add_atoms_to_structure(self):
        structure = AtomicStructure(self.atom1, self.atom2)
        structure.add_atom(self.atom3)
        self.assertEqual(structure._atoms, set(self.atoms))


    def test_can_only_add_atoms(self):
        structure = AtomicStructure(self.atom1, self.atom2)
        with self.assertRaises(TypeError):
            structure.add_atom("atom")



class AtomicStructureAtomRemovalTests(AtomicStructureTest):

    def test_can_remove_atoms_from_structure(self):
        structure = AtomicStructure(self.atom1, self.atom2)
        structure.remove_atom(self.atom1)
        self.assertEqual(structure._atoms, set(self.atoms[1:2]))



class AtomicStructureMassTests(AtomicStructureTest):

    def test_structure_mass_is_sum_of_atom_masses(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(structure.mass(), 17.1)



class AtomicStructureChargeTests(AtomicStructureTest):

    def test_structure_charge_is_sum_of_atom_charges(self):
        self.atom1.charge.return_value = 0.2
        self.atom2.charge.return_value = -1.4
        self.atom3.charge.return_value = 0.6
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertAlmostEqual(structure.charge(), -0.6, delta=0.000005)



class AtomicStructureFormulaTests(AtomicStructureTest):

    def test_can_get_formula(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(structure.formula(), Counter({"A":1, "B":2}))



class AtomicStructureTranslationTests(AtomicStructureTest):

    @patch("atomium.structures.molecules.translate")
    def test_translation_uses_geometrica(self, mock_translate):
        structure = AtomicStructure(self.atom1, self.atom2)
        mock_translate.return_value = [(6, 6, 1), (9, 9, 4)]
        self.atom1._x, self.atom1._y, self.atom1._z = 1, 2, 3
        self.atom2._x, self.atom2._y, self.atom2._z = 4, 5, 6
        structure.translate(5, 4, -2)
        mock_translate.assert_called_with(list(structure._atoms), 5, 4, -2)
        self.assertEqual(
         set([
          (self.atom1._x, self.atom1._y, self.atom1._z),
          (self.atom2._x, self.atom2._y, self.atom2._z)
         ]),
         set([(6, 6, 1), (9, 9, 4)])
        )



class AtomicStructureRotationTests(AtomicStructureTest):

    @patch("atomium.structures.molecules.rotate")
    def test_rotation_uses_geometrica(self, mock_rotate):
        structure = AtomicStructure(self.atom1, self.atom2)
        mock_rotate.return_value = [(6, 6, 1), (9, 9, 4)]
        self.atom1._x, self.atom1._y, self.atom1._z = 1, 2, 3
        self.atom2._x, self.atom2._y, self.atom2._z = 4, 5, 6
        structure.rotate("x", 90)
        mock_rotate.assert_called_with(list(structure._atoms), "x", 90)
        self.assertEqual(
         set([
          (self.atom1._x, self.atom1._y, self.atom1._z),
          (self.atom2._x, self.atom2._y, self.atom2._z)
         ]),
         set([(6, 6, 1), (9, 9, 4)])
        )



class AtomicStructureCenterOfMassTests(AtomicStructureTest):

    @patch("atomium.structures.molecules.AtomicStructure.mass")
    def test_can_get_center_of_mass_when_equal_mass(self, mock_mass):
        mock_mass.return_value = 20
        self.atom1._x, self.atom1._y, self.atom1._z = 0, 0, 0
        self.atom2._x, self.atom2._y, self.atom2._z = 1, 1, 1
        self.atom1.mass.return_value = 10
        self.atom2.mass.return_value = 10
        structure = AtomicStructure(self.atom1, self.atom2)
        self.assertEqual(structure.center_of_mass(), (0.5, 0.5, 0.5))


    @patch("atomium.structures.molecules.AtomicStructure.mass")
    def test_can_get_center_of_mass_when_unequal_mass(self, mock_mass):
        mock_mass.return_value = 40
        self.atom1._x, self.atom1._y, self.atom1._z = 0, 0, 0
        self.atom2._x, self.atom2._y, self.atom2._z = 1, 1, 1
        self.atom1.mass.return_value = 10
        self.atom2.mass.return_value = 30
        structure = AtomicStructure(self.atom1, self.atom2)
        self.assertEqual(structure.center_of_mass(), (0.75, 0.75, 0.75))



class AtomicStructureRadiusOfGyrationTests(AtomicStructureTest):

    @patch("atomium.structures.molecules.AtomicStructure.center_of_mass")
    def test_can_get_radius_of_gyration(self, mock_center):
        mock_center.return_value = (5, 0, 0)
        self.atom1.distance_to.return_value = 5
        self.atom2.distance_to.return_value = 5
        self.atom1._x, self.atom1._y, self.atom1._z = 0, 0, 0
        self.atom2._x, self.atom2._y, self.atom2._z = 10, 0, 0
        structure = AtomicStructure(self.atom1, self.atom2)
        self.assertEqual(structure.radius_of_gyration(), 5)
        self.atom1.distance_to.assert_called_with((5, 0, 0))
        self.atom2.distance_to.assert_called_with((5, 0, 0))



class AtomicStructureToStringTests(AtomicStructureTest):

    @patch("atomium.converters.structure2xyzstring.structure_to_xyz_string")
    def test_can_save_as_xyz_string(self, mock_convert):
        mock_convert.return_value = "filestring"
        structure = AtomicStructure(*self.atoms)
        s = structure.to_file_string("xyz")
        self.assertEqual(s, "filestring")
        mock_convert.assert_called_with(structure, "")


    @patch("atomium.converters.structure2xyzstring.structure_to_xyz_string")
    def test_can_save_as_xyz_string_with_comment(self, mock_convert):
        mock_convert.return_value = "filestring"
        structure = AtomicStructure(*self.atoms)
        s = structure.to_file_string("xyz", description="A description")
        self.assertEqual(s, "filestring")
        mock_convert.assert_called_with(structure, "A description")


    def test_invalid_file_format_is_error(self):
        structure = AtomicStructure(*self.atoms)
        with self.assertRaises(ValueError):
            structure.to_file_string("nosuchfile")



class AtomicStructureSavingTests(AtomicStructureTest):

    @patch("atomium.structures.molecules.AtomicStructure.to_file_string")
    @patch("atomium.converters.strings.string_to_file")
    def test_saving_uses_correct_functions(self, mock_save, mock_string):
        mock_string.return_value = "filestring"
        structure = AtomicStructure(*self.atoms)
        structure.save("file.xyz", "a description")
        mock_string.assert_called_with("xyz", "a description")
        mock_save.assert_called_with("filestring", "file.xyz")


    @patch("atomium.structures.molecules.AtomicStructure.to_file_string")
    @patch("atomium.converters.strings.string_to_file")
    def test_file_format_extracted(self, mock_save, mock_string):
        mock_string.return_value = "filestring"
        structure = AtomicStructure(*self.atoms)
        structure.save("path/to/file.dfghdfg", "a description")
        mock_string.assert_called_with("dfghdfg", "a description")
        mock_save.assert_called_with("filestring", "path/to/file.dfghdfg")
