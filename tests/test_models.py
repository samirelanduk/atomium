from unittest import TestCase
import unittest.mock
from molecupy.structures import Model, AtomicStructure, SmallMolecule, Chain

class ModelTest(TestCase):

    def setUp(self):
        self.small_molecule1 = unittest.mock.Mock(spec=SmallMolecule)
        self.small_molecule1._model = None
        self.small_molecule1.molecule_id.return_value = "A100"
        self.small_molecule1.molecule_name.return_value = "MOL"
        self.small_molecule2 = unittest.mock.Mock(spec=SmallMolecule)
        self.small_molecule2._model = None
        self.small_molecule2.molecule_id.return_value = "A101"
        self.small_molecule2.molecule_name.return_value = "HET"
        self.chain1 = unittest.mock.Mock(spec=Chain)
        self.chain1._model = None
        self.chain1.chain_id.return_value = "A"
        self.chain2 = unittest.mock.Mock(spec=Chain)
        self.chain2._model = None
        self.chain2.chain_id.return_value = "B"



class ModelCreationTest(ModelTest):

    def test_can_create_chain(self):
        model = Model()
        self.assertIsInstance(model, AtomicStructure)
        self.assertEqual(model._atoms, set())


    def test_model_repr(self):
        model = Model()
        self.assertEqual(str(model), "<Model (0 atoms)>")



class ModelSmallMoleculeTests(ModelTest):

    def test_can_add_small_molecules(self):
        model = Model()
        self.assertEqual(model.small_molecules(), set())
        model.add_small_molecule(self.small_molecule1)
        self.assertEqual(model.small_molecules(), set([self.small_molecule1]))
        model.add_small_molecule(self.small_molecule2)
        self.assertEqual(
         model.small_molecules(),
         set([self.small_molecule1, self.small_molecule2])
        )


    def test_must_use_method_to_add_small_molecule(self):
        model = Model()
        self.assertEqual(model.small_molecules(), set())
        model.small_molecules().add(self.small_molecule1)
        self.assertEqual(model.small_molecules(), set())


    def test_can_remove_small_molecules(self):
        model = Model()
        model.add_small_molecule(self.small_molecule1)
        self.assertEqual(model.small_molecules(), set([self.small_molecule1]))
        model.remove_small_molecule(self.small_molecule1)
        self.assertEqual(model.small_molecules(), set())


    def test_small_molecule_knows_about_model(self):
        model = Model()
        self.assertIs(self.small_molecule1._model, None)
        model.add_small_molecule(self.small_molecule1)
        self.assertIs(self.small_molecule1._model, model)
        model.remove_small_molecule(self.small_molecule1)
        self.assertIs(self.small_molecule1._model, None)


    def test_can_only_add_small_molecules(self):
        model = Model()
        with self.assertRaises(TypeError):
            model.add_small_molecule("molecule")


    def test_can_get_small_molecule_by_id(self):
        model = Model()
        model.add_small_molecule(self.small_molecule1)
        model.add_small_molecule(self.small_molecule2)
        self.assertIs(model.get_small_molecule_by_id("A100"), self.small_molecule1)
        self.assertIs(model.get_small_molecule_by_id("A101"), self.small_molecule2)
        self.assertIs(model.get_small_molecule_by_id("A102"), None)


    def test_can_only_get_small_molecule_with_str_id(self):
        model = Model()
        with self.assertRaises(TypeError):
            model.get_small_molecule_by_id(100)


    def test_can_get_small_molecule_by_name(self):
        model = Model()
        model.add_small_molecule(self.small_molecule1)
        model.add_small_molecule(self.small_molecule2)
        self.assertIs(model.get_small_molecule_by_name("MOL"), self.small_molecule1)
        self.assertIs(model.get_small_molecule_by_name("HET"), self.small_molecule2)
        self.assertIs(model.get_small_molecule_by_name("ABC"), None)


    def test_can_only_get_small_molecule_with_str_name(self):
        model = Model()
        with self.assertRaises(TypeError):
            model.get_small_molecule_by_name(100)


    def test_can_get_small_molecules_by_name(self):
        model = Model()
        model.add_small_molecule(self.small_molecule1)
        model.add_small_molecule(self.small_molecule2)
        self.assertEqual(
         model.get_small_molecules_by_name("MOL"),
         set([self.small_molecule1])
        )
        self.small_molecule2.molecule_name.return_value = "MOL"
        self.assertEqual(
         model.get_small_molecules_by_name("MOL"),
         set([self.small_molecule1, self.small_molecule2])
        )
        self.assertEqual(model.get_small_molecules_by_name("ABC"), set())


    def test_can_only_get_small_molecules_with_str_name(self):
        model = Model()
        with self.assertRaises(TypeError):
            model.get_small_molecules_by_name(100)



class ModelChainTests(ModelTest):

    def test_can_add_chains(self):
        model = Model()
        self.assertEqual(model.chains(), set())
        model.add_chain(self.chain1)
        self.assertEqual(model.chains(), set([self.chain1]))
        model.add_chain(self.chain2)
        self.assertEqual(
         model.chains(),
         set([self.chain1, self.chain2])
        )


    def test_must_use_method_to_add_chain(self):
        model = Model()
        self.assertEqual(model.chains(), set())
        model.chains().add(self.chain1)
        self.assertEqual(model.chains(), set())


    def test_can_remove_chains(self):
        model = Model()
        model.add_chain(self.chain1)
        self.assertEqual(model.chains(), set([self.chain1]))
        model.remove_chain(self.chain1)
        self.assertEqual(model.chains(), set())


    def test_chain_knows_about_model(self):
        model = Model()
        self.assertIs(self.chain1._model, None)
        model.add_chain(self.chain1)
        self.assertIs(self.chain1._model, model)
        model.remove_chain(self.chain1)
        self.assertIs(self.chain1._model, None)


    def test_can_only_add_chains(self):
        model = Model()
        with self.assertRaises(TypeError):
            model.add_chain("chain")


    def test_can_get_chain_by_id(self):
        model = Model()
        model.add_chain(self.chain1)
        model.add_chain(self.chain2)
        self.assertIs(model.get_chain_by_id("A"), self.chain1)
        self.assertIs(model.get_chain_by_id("B"), self.chain2)
        self.assertIs(model.get_chain_by_id("C"), None)


    def test_can_only_get_chain_with_str_id(self):
        model = Model()
        with self.assertRaises(TypeError):
            model.get_chain_by_id(100)
