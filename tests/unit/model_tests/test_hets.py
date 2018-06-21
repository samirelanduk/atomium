from unittest import TestCase
from unittest.mock import patch, Mock, PropertyMock
from atomium.models.molecules import Het
from atomium.models.structures import AtomStructure

class HetCreationTests(TestCase):

    def test_can_create_het(self):
        het = Het()
        self.assertIsInstance(het, AtomStructure)



class HetMembershipPropertiesTests(TestCase):

    def test_hets_have_correct_properties(self):
        self.assertIn("model", Het.__dict__)
        self.assertIn("chain", Het.__dict__)
