from unittest import TestCase
from unittest.mock import Mock, patch, PropertyMock, MagicMock
from atomium.structures import Chain, Residue, AtomStructure, Molecule

class ChainTest(TestCase):

    def setUp(self):
        self.residues = [Mock(Residue, _id=3), Mock(Residue, _id=4)]



class ChainCreationTests(ChainTest):

    def test_can_create_empty_chain(self):
        chain = Chain()
        self.assertEqual(chain._residues, {})
        self.assertIsInstance(chain, AtomStructure)
        self.assertIsInstance(chain, Molecule)
        self.assertIsNone(chain._model)
        self.assertIsNone(chain._id)
        self.assertIsNone(chain._internal_id)
        self.assertEqual(chain._sequence, "")


    def test_can_create_populated_chain(self):
        chain = Chain(*self.residues, id="A", internal_id="C", sequence="FG")
        self.assertEqual(chain._residues, {3: self.residues[0], 4: self.residues[1]})
        self.assertEqual(chain._id, "A")
        self.assertEqual(chain._internal_id, "C")
        self.assertEqual(chain._sequence, "FG")
        for res in self.residues: self.assertIs(res._chain, chain)



class ChainReprTests(ChainTest):

    def test_can_get_chain_repr(self):
        chain = Chain(*self.residues, id="A", internal_id="C", sequence="FG")
        self.assertEqual(repr(chain), "<Chain A (2 residues)>")



class ChainLenTests(ChainTest):

    def test_can_get_chain_len(self):
        chain = Chain(*self.residues, id="A", internal_id="C", sequence="FG")
        self.assertEqual(len(chain), 2)



class ChainIterTests(ChainTest):

    def test_can_get_chain_iter(self):
        chain = Chain(*self.residues, id="A", internal_id="C", sequence="FG")
        for res in chain: self.assertIn(res, self.residues)



class ChainIndexingTests(ChainTest):

    @patch("atomium.structures.Chain.residues")
    def test_can_get_chain_index(self, mock_res):
        mock_res.return_value = self.residues
        chain = Chain(*self.residues, id="A", internal_id="C", sequence="FG")
        self.assertIs(chain[0], self.residues[0])
        self.assertIs(chain[1], self.residues[1])



class ChainContainerTests(ChainTest):

    def test_chain_is_container_of_residues(self):
        chain = Chain(*self.residues, id="A", internal_id="C", sequence="FG")
        self.assertIn(self.residues[0], chain)


    @patch("atomium.structures.Chain.atoms")
    def test_chain_is_container_of_atoms(self, mock_atoms):
        mock_atoms.return_value = [1, 2]
        chain = Chain(*self.residues, id="A", internal_id="C", sequence="FG")
        self.assertIn(1, chain)
        self.assertIn(2, chain)
        self.assertNotIn(3, chain)



class ChainLengthTests(ChainTest):

    @patch("atomium.structures.Chain.__len__")
    def test_chain_length(self, mock_len):
        mock_len.return_value = 100
        chain = Chain(*self.residues, id="A", internal_id="C", sequence="FG")
        self.assertEqual(chain.length, 100)



class ChainResiduesProperty(ChainTest):

    def test_chain_residues(self):
        chain = Chain(*self.residues, id="A", internal_id="C", sequence="FG")
        self.assertEqual(chain.residues(), tuple(self.residues))



class ChainLigandsProperty(ChainTest):

    def test_chain_ligands(self):
        chain = Chain(*self.residues, id="A", internal_id="C", sequence="FG")
        ligands = {
         1: Mock(_chain=chain), 2: Mock(_chain=chain), 3: Mock(_chain=2)
        }
        model = Mock(_ligands=ligands)
        chain._model = model
        self.assertEqual(chain.ligands(), {ligands[1], ligands[2]})


    def test_chain_ligands_when_no_model(self):
        chain = Chain(*self.residues, id="A", internal_id="C", sequence="FG")
        self.assertEqual(chain.ligands(), set())



class ChainAtomsProperty(ChainTest):

    def test_chain_atoms(self):
        chain = Chain(*self.residues, id="A", internal_id="C", sequence="FG")
        self.residues[0]._atoms = {1: 2, 3: 4}
        self.residues[1]._atoms = {5: 6, 7: 8}
        self.assertEqual(chain.atoms(), {2, 4, 6, 8})



class ChainSequenceProperty(ChainTest):

    def test_can_get_chain(self):
        chain = Chain(*self.residues, id="A", internal_id="C", sequence="FG")
        self.assertIs(chain.sequence, chain._sequence)


    def test_can_set_chain(self):
        chain = Chain(*self.residues, id="A", internal_id="C", sequence="FG")
        chain.sequence = "123"
        self.assertEqual(chain._sequence, "123")



class ChainCopyingTests(ChainTest):

    @patch("atomium.structures.Chain.residues")
    def test_can_copy_chain(self, mock_res):
        mock_res.return_value = self.residues
        chain = Chain(*self.residues, id="A", internal_id="C", sequence="FG")
        p = patch("atomium.structures.Chain")
        mock_chain = p.start()
        try:
            copy = chain.copy()
            mock_chain.assert_called_with(
             self.residues[0].copy.return_value, self.residues[1].copy.return_value,
             id="A", internal_id="C", sequence="FG"
            )
        finally: p.stop()
