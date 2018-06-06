from unittest import TestCase
from unittest.mock import patch, Mock, PropertyMock
from atomium.models.molecules import Het, Ligand

class HetCreationTests(TestCase):

    def test_can_create_ligand(self):
        ligand = Ligand()
        self.assertIsInstance(ligand, Het)
