from unittest import TestCase
from unittest.mock import Mock, patch, PropertyMock, MagicMock
from atomium.structures import Model, Chain, Ligand, AtomStructure

class ModelTest(TestCase):

    def setUp(self):
        self.chains = [Mock(Chain, _id=1), Mock(Chain, _id=2)]
        self.ligands = [Mock(Ligand, _water=False, _id=3), Mock(Ligand, _water=False, _id=4)]
        self.waters = [Mock(Ligand, _water=True, _id=5), Mock(Ligand, _water=True, _id=6)]
        self.patch1 = patch("atomium.structures.StructureSet")
        self.mock_set = self.patch1.start()
        self.mock_set.return_value = Mock(structures=[])


    def tearDown(self):
        self.patch1.stop()


class ModelCreationTests(ModelTest):

    def test_can_create_empty_model(self):
        model = Model()
        self.assertEqual(model._chains, self.mock_set.return_value)
        self.assertEqual(model._ligands, self.mock_set.return_value)
        self.assertEqual(model._waters, self.mock_set.return_value)
        self.mock_set.assert_called_with()
        self.assertEqual(self.mock_set.call_count, 3)
        self.assertIsInstance(model, AtomStructure)


    def test_can_create_model_with_stuff(self):
        model = Model(*(self.chains + self.ligands + self.waters))
        self.assertEqual(model._chains, self.mock_set.return_value)
        self.assertEqual(model._ligands, self.mock_set.return_value)
        self.assertEqual(model._waters, self.mock_set.return_value)
        self.mock_set.assert_called_with()
        self.assertEqual(self.mock_set.call_count, 3)
        model._chains.add.assert_any_call(self.chains[0])
        model._chains.add.assert_any_call(self.chains[1])
        model._ligands.add.assert_any_call(self.ligands[0])
        model._ligands.add.assert_any_call(self.ligands[1])
        model._waters.add.assert_any_call(self.waters[0])
        model._waters.add.assert_any_call(self.waters[1])



class ModelReprTests(ModelTest):

    def test_empty_model_repr(self):
        model = Model()
        self.mock_set.return_value.__len__ = MagicMock()
        self.mock_set.return_value.__len__.return_value = 0
        self.assertEqual(repr(model), "<Model (0 chains, 0 ligands)>")


    def test_model_repr(self):
        model = Model()
        self.mock_set.return_value.__len__ = MagicMock()
        self.mock_set.return_value.__len__.side_effect = (1, 1, 2, 2)
        self.assertEqual(repr(model), "<Model (1 chain, 2 ligands)>")



class ModelContainerTests(ModelTest):

    @patch("atomium.structures.Model.molecules")
    @patch("atomium.structures.Model.residues")
    @patch("atomium.structures.Model.atoms")
    def test_model_is_container_of_molecules(self, mock_at, mock_res, mock_mol):
        model = Model()
        mock_mol.return_value = [1, 2]
        mock_res.return_value = [3, 4]
        mock_at.return_value = [5, 6]
        self.assertIn(1, model)
        self.assertIn(3, model)
        self.assertIn(5, model)
        self.assertNotIn(10, model)



class ModelChainsTests(ModelTest):
    """Impossible to unit test"""



class ModelLigandsTests(ModelTest):
    """Impossible to unit test"""



class ModelWatersTests(ModelTest):
    """Impossible to unit test"""



class ModelMoleculesTests(ModelTest):

    def test_can_get_molecules(self):
        model = Model(self.chains[0], self.ligands[0], self.waters[0])
        self.mock_set.return_value.__add__ = MagicMock()
        _value = self.mock_set.return_value
        self.assertEqual(model.molecules(), set())



class ModelResiduesTests(ModelTest):

    def test_can_get_residues(self):
        model = Model(*self.chains)
        model._chains.structures = self.chains
        residues = [Mock() for _ in range(6)]
        self.chains[0].residues.return_value = residues[:3]
        self.chains[1].residues.return_value = residues[3:]
        model.residues()
        self.mock_set.assert_called_with(*residues)



class ModelAtomsTests(ModelTest):

    @patch("atomium.structures.Model.molecules")
    def test_can_get_atoms(self, mock_mol):
        molecules = [Mock() for _ in range(3)]
        molecules[0]._atoms.structures = {2, 4}
        molecules[1]._atoms.structures = {6}
        molecules[2]._atoms.side_effect = Exception
        molecules[2]._residues.structures = {Mock(_atoms=Mock(structures={9})), Mock(_atoms=Mock(structures={11}))}
        mock_mol.return_value = molecules
        model = Model()
        model.atoms()
        self.mock_set.assert_called_with(2, 4, 6, 9, 11)



class ModelDehydrationTests(ModelTest):

    def test_can_dehydrate_model(self):
        model = Model(*(self.chains + self.ligands + self.waters))
        model._chains, model._ligands, model._waters = "CA", "LI", "WA"
        model.dehydrate()
        self.assertEqual(model._waters, self.mock_set.return_value)
        self.assertEqual(len(model._ligands), 2)
        self.assertEqual(len(model._chains), 2)



class ModelAddingTests(ModelTest):

    def test_can_add_molecule(self):
        model = Model()
        model.add(self.chains[0])
        model._chains.add.assert_called_with(self.chains[0])
        model.add(self.ligands[0])
        model._ligands.add.assert_called_with(self.ligands[0])
        model.add(self.waters[0])
        model._waters.add.assert_called_with(self.waters[0])



class ModelRemovingTests(ModelTest):

    def test_can_remove_molecule(self):
        model = Model(*(self.chains + self.ligands + self.waters))
        model.remove(self.chains[0])
        model._chains.remove.assert_called_with(self.chains[0])
        model.remove(self.ligands[0])
        model._ligands.remove.assert_called_with(self.ligands[0])
        model.remove(self.waters[0])
        model._waters.remove.assert_called_with(self.waters[0])
