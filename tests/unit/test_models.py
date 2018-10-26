from unittest import TestCase
from unittest.mock import Mock, patch, PropertyMock, MagicMock
from atomium.structures import Model, Chain, Ligand, AtomStructure

class ModelTest(TestCase):

    def setUp(self):
        self.chains = [Mock(Chain, _id=1), Mock(Chain, _id=2)]
        self.ligands = [Mock(Ligand, _water=False, _id=3), Mock(Ligand, _water=False, _id=4)]
        self.waters = [Mock(Ligand, _water=True, _id=5), Mock(Ligand, _water=True, _id=6)]


class ModelCreationTests(ModelTest):

    def test_can_create_empty_model(self):
        model = Model()
        self.assertEqual(model._chains, {})
        self.assertEqual(model._ligands, {})
        self.assertEqual(model._waters, {})
        self.assertIsInstance(model, AtomStructure)


    def test_can_create_model_with_stuff(self):
        model = Model(*(self.chains + self.ligands + self.waters))
        self.assertEqual(model._chains, {1: self.chains[0], 2: self.chains[1]})
        self.assertEqual(model._ligands, {3: self.ligands[0], 4: self.ligands[1]})
        self.assertEqual(model._waters, {5: self.waters[0], 6: self.waters[1]})
        for mol in self.chains: self.assertIs(mol._model, model)
        for mol in self.ligands: self.assertIs(mol._model, model)
        for mol in self.waters: self.assertIs(mol._model, model)



class ModelReprTests(ModelTest):

    def test_empty_model_repr(self):
        model = Model()
        self.assertEqual(repr(model), "<Model (0 chains, 0 ligands)>")


    def test_model_repr(self):
        model = Model()
        model._chains = [1]
        model._ligands = [2, 3]
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

    def test_can_get_chains(self):
        model = Model(self.chains[0], self.ligands[0], self.waters[0])
        self.assertEqual(model.chains(), {self.chains[0]})



class ModelLigandsTests(ModelTest):

    def test_can_get_ligands(self):
        model = Model(self.chains[0], self.ligands[0], self.waters[0])
        self.assertEqual(model.ligands(), {self.ligands[0]})



class ModelWatersTests(ModelTest):

    def test_can_get_waters(self):
        model = Model(self.chains[0], self.ligands[0], self.waters[0])
        self.assertEqual(model.waters(), {self.waters[0]})


class ModelMoleculesTests(ModelTest):

    def test_can_get_molecules(self):
        model = Model(self.chains[0], self.ligands[0], self.waters[0])
        self.assertEqual(model.molecules(), {self.chains[0], self.ligands[0], self.waters[0]})



class ModelResiduesTests(ModelTest):

    def test_can_get_residues(self):
        model = Model(*self.chains)
        residues = [Mock() for _ in range(6)]
        self.chains[0].residues.return_value = residues[:3]
        self.chains[1].residues.return_value = residues[3:]
        self.assertEqual(model.residues(), set(residues))



class ModelAtomsTests(ModelTest):

    @patch("atomium.structures.Model.molecules")
    def test_can_get_atoms(self, mock_mol):
        molecules = [Mock() for _ in range(3)]
        molecules[0]._atoms = {1: 2, 3: 4}
        molecules[1]._atoms = {5: 6}
        molecules[2]._atoms.side_effect = Exception
        molecules[2]._residues = {7: Mock(_atoms={8: 9}), 10: Mock(_atoms={10: 11})}
        mock_mol.return_value = molecules
        model = Model()
        self.assertEqual(model.atoms(), {2, 4, 6, 9, 11})



class ModelDehydrationTests(ModelTest):

    def test_can_dehydrate_model(self):
        model = Model(*(self.chains + self.ligands + self.waters))
        self.assertEqual(len(model._waters), 2)
        model.dehydrate()
        self.assertEqual(len(model._waters), 0)
        self.assertEqual(len(model._ligands), 2)
        self.assertEqual(len(model._chains), 2)



class ModelAddingTests(ModelTest):

    def test_can_add_molecule(self):
        model = Model()
        model.add(self.chains[0])
        self.assertEqual(len(model._chains), 1)
        model.add(self.ligands[0])
        self.assertEqual(len(model._ligands), 1)
        model.add(self.waters[0])
        self.assertEqual(len(model._waters), 1)



class ModelRemovingTests(ModelTest):

    def test_can_remove_molecule(self):
        model = Model(*(self.chains + self.ligands + self.waters))
        model.remove(self.chains[0])
        self.assertEqual(len(model._chains), 1)
        model.remove(self.ligands[0])
        self.assertEqual(len(model._ligands), 1)
        model.remove(self.waters[0])
        self.assertEqual(len(model._waters), 1)
