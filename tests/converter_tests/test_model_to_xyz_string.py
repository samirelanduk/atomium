from unittest import TestCase
from unittest.mock import Mock, patch
from atomium.converters.model2xyzstring import model_to_xyz_string
from atomium.structures.models import Model
from atomium.structures.atoms import Atom

class ModelToXyzTest(TestCase):

    def setUp(self):
        self.model = Mock(Model)
        self.atoms = [Mock(Atom) for _ in range(3)]
        for index, atom in enumerate(self.atoms):
            atom.element.return_value = chr(65 + index)
            atom.x.return_value = (index * 5) + ((1 + index) / 10)
            atom.y.return_value = (index * -5) + ((1 + index) / 100)
            atom.z.return_value = (index * -1000) + ((1 + index) / 1000)
        self.model.atoms.return_value = set(self.atoms)


    def test_can_get_string_from_model(self):
        self.assertEqual(model_to_xyz_string(self.model), (
         "3\n"
         "\n"
         "A      0.100      0.010      0.001\n"
         "B      5.200     -4.980   -999.998\n"
         "C     10.300     -9.970  -1999.997"
        ))
