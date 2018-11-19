from unittest import TestCase
from unittest.mock import Mock, patch, PropertyMock, MagicMock
from atomium.structures import Het

class HetTest(TestCase):

    def setUp(self):
        class Structure(Het): pass
        self.het = Structure()
        self.patch = patch("atomium.structures.StructureSet")
        self.mock_set = self.patch.start()
        self.mock_set.return_value = Mock(structures=[])


    def tearDown(self):
        self.patch.stop()



class HetInitTests(HetTest):

    def test_can_initialise_het_atoms(self):
        atoms = [Mock(_id=20), Mock(_id=4)]
        Het.__init__(self.het, *atoms)
        self.assertEqual(self.het._atoms, self.mock_set.return_value)
        self.mock_set.assert_called_with(*atoms)
        for atom in atoms: self.assertIs(atom._structure, self.het)



class HetContainerTests(HetTest):

    def test_hets_are_containers_of_atoms(self):
        self.het._atoms = Mock(structures="AB")
        self.assertIn("A", self.het)
        self.assertIn("B", self.het)
        self.assertNotIn("C", self.het)



class HetChainTests(HetTest):

    def test_can_refer_to_chain(self):
        self.het._chain = "CHAIN"
        self.assertIs(self.het._chain, self.het.chain)



class HetModelTests(HetTest):

    def test_can_refer_to_model(self):
        self.het._chain = Mock(_model="MODEL")
        self.assertEqual(self.het.model, "MODEL")


    def test_can_refer_to_model_when_no_chain(self):
        self.het._chain = None
        self.assertIsNone(self.het.model)



class HetAtomAdditionTests(HetTest):

    def test_can_add_atoms(self):
        self.het._atoms = Mock()
        atom = Mock(id=5)
        self.het.add(atom)
        self.het._atoms.add.assert_called_with(atom)
        self.assertIs(atom._structure, self.het)



class HetAtomRemovalTests(HetTest):

    def test_can_remove_atoms(self):
        atom = Mock(id=5)
        self.het._atoms = Mock()
        self.het.remove(atom)
        self.het._atoms.remove.assert_called_with(atom)
        self.assertIsNone(atom._structure)



class HetNearbyStructuresTests(HetTest):

    def test_can_get_nearby_atoms(self):
        atoms = [Mock(), Mock(), Mock()]
        atoms[0].nearby_structures.return_value = {1, 2}
        atoms[1].nearby_structures.return_value = {3, 1, 4}
        atoms[2].nearby_structures.return_value = {1, 2, 3, 4, 9}
        self.het._atoms = Mock(structures={atoms[0], atoms[1], atoms[2]})
        self.assertEqual(self.het.nearby_structures(1, a=2), {1, 2, 3, 4, 9})
        for a in atoms:
            a.nearby_structures.assert_called_with(1, a=2)
