from unittest import TestCase
from unittest.mock import patch, Mock, PropertyMock
from atomium.models.molecules import Model
from atomium.models.structures import AtomStructure

class ModelCreationTests(TestCase):

    def test_can_create_model(self):
        model = Model()
        self.assertIsInstance(model, AtomStructure)



class ModelMembershipPropertiesTests(TestCase):

    def test_models_have_correct_properties(self):
        for attr in ["ligand", "residue", "chain"]:
            self.assertIn(attr, Model.__dict__)
            self.assertIn(attr + "s", Model.__dict__)



class ModelCopyingTests(TestCase):

    @patch("atomium.models.molecules.Model.chains")
    def test_can_copy_model_chains(self, mock_chains):
        chains = [Mock(), Mock()]
        mock_chains.return_value = chains
        chains[0].copy.return_value, chains[1].copy.return_value = "AB"
        model = Model()
        patcher = patch("atomium.models.molecules.Model")
        mock_model = patcher.start()
        try:
            new_model = model.copy()
            mock_model.assert_called_with("A", "B")
        finally:
            patcher.stop()


    @patch("atomium.models.molecules.Model.chains")
    def test_can_copy_model_chains_and_atoms(self, mock_chains):
        chains = [Mock(), Mock()]
        mock_chains.return_value = chains
        chains[0].copy.return_value, chains[1].copy.return_value = "AB"
        model = Model()
        atoms = [Mock(), Mock(), Mock(_chain=None)]
        model._atoms = set(atoms)
        patcher = patch("atomium.models.molecules.Model")
        mock_model = patcher.start()
        mock_model.return_value = Mock(_atoms=set(atoms[:-1]))
        try:
            new_model = model.copy()
            mock_model.assert_called_with("A", "B")
            atoms[-1].copy.assert_called_with()
            self.assertFalse(atoms[0].called)
            self.assertFalse(atoms[0].called)
            self.assertEqual(new_model._atoms, set(atoms[:-1] + [atoms[-1].copy.return_value]))
        finally:
            patcher.stop()
