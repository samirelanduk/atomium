from unittest import TestCase
from unittest.mock import patch, Mock, PropertyMock
from atomium.models.molecules import Complex
from atomium.models.structures import AtomStructure

class ComplexCreationTests(TestCase):

    def test_can_create_complex(self):
        complex = Complex()
        self.assertIsInstance(complex, AtomStructure)
