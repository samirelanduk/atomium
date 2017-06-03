from unittest import TestCase
from unittest.mock import Mock, patch
from atomium.converters.structure2xyzstring import structure_to_xyz_string
from atomium.structures.molecules import AtomicStructure
from atomium.structures.atoms import Atom

class AtomicStructureToXyzTest(TestCase):

    def setUp(self):
        self.structure = Mock(AtomicStructure)
        self.atoms = [Mock(Atom) for _ in range(3)]
        for index, atom in enumerate(self.atoms):
            atom.element.return_value = chr(65 + index)
            atom.x.return_value = (index * 5) + ((1 + index) / 10)
            atom.y.return_value = (index * -5) + ((1 + index) / 100)
            atom.z.return_value = (index * -1000) + ((1 + index) / 1000)
        self.structure.atoms.return_value = set(self.atoms)


    def test_can_get_string_from_structure(self):
        self.assertEqual(structure_to_xyz_string(self.structure), (
         "3\n"
         "\n"
         "A      0.100      0.010      0.001\n"
         "B      5.200     -4.980   -999.998\n"
         "C     10.300     -9.970  -1999.997"
        ))


    def test_structure_needed(self):
        with self.assertRaises(TypeError):
            structure_to_xyz_string("structure")


    def test_can_take_comment(self):
        self.assertEqual(structure_to_xyz_string(self.structure, comment="comment"), (
         "3\n"
         "comment\n"
         "A      0.100      0.010      0.001\n"
         "B      5.200     -4.980   -999.998\n"
         "C     10.300     -9.970  -1999.997"
        ))


    def test_comment_must_be_str(self):
        with self.assertRaises(TypeError):
            structure_to_xyz_string(self.structure, comment=100)
