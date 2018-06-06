from unittest import TestCase
from unittest.mock import patch, Mock, PropertyMock
from atomium.models.molecules import Het
from atomium.models.structures import AtomStructure

class HetCreationTests(TestCase):

    def test_can_create_het(self):
        het = Het()
        self.assertIsInstance(het, AtomStructure)
