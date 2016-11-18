from unittest import TestCase
from unittest.mock import Mock
from molecupy.structures import Chain, Complex, ResiduicStructure, Residue, Atom
from molecupy.exceptions import DuplicateChainsError

class ComplexTest(TestCase):

    def setUp(self):
        self.chains = [Mock(Chain, _complex=None) for _ in range(3)]
        self.residues = [Mock(ResiduicStructure) for _ in range(6)]
        self.atoms = [Mock(Atom) for _ in range(18)]
        for index, chain in enumerate(self.chains):
            chain.residues.return_value = set(
             self.residues[index * 2: (index * 2) + 2]
            )
            chain.chain_id.return_value = chr(index + 65)
        for index, residue in enumerate(self.residues):
            residue.atoms.return_value = set(
             self.atoms[index * 3: (index * 3) + 3]
            )


class ComplexCreationTests(ComplexTest):

    def test_can_create_complex(self):
        complex_ = Complex("1", "A Complex", *self.chains)
        self.assertIsInstance(complex_, ResiduicStructure)
        self.assertEqual(complex_.residues(), set(self.residues))
        self.assertEqual(complex_.atoms(), set(self.atoms))
        self.assertEqual(complex_._complex_id, "1")
        self.assertEqual(complex_._complex_name, "A Complex")
        self.assertEqual(complex_._chains, set(self.chains))
        self.assertEqual(complex_._model, None)


    def test_complex_id_must_be_str(self):
        with self.assertRaises(TypeError):
            Complex(1, "A Complex", *self.chains)


    def test_complex_name_must_be_str(self):
        with self.assertRaises(TypeError):
            Complex("1", 100, *self.chains)


    def test_complex_chains_must_be_chains(self):
        with self.assertRaises(TypeError):
            Complex("1", "A Complex", "Chain 1", "Chain 2")


    def test_complex_updates_chains(self):
        complex_ = Complex("1", "A Complex", *self.chains)
        for chain in self.chains:
            self.assertIs(chain._complex, complex_)


    def test_repr(self):
        complex_ = Complex("1", "A Complex", *self.chains)
        self.assertEqual(str(complex_), "<Complex 'A Complex' (3 chains)>")



class ComplexPropertyTests(ComplexTest):

    def test_can_get_complex_properties(self):
        complex_ = Complex("1", "A Complex", *self.chains)
        self.assertEqual(complex_.complex_id(), "1")
        self.assertEqual(complex_.complex_name(), "A Complex")
        self.assertEqual(complex_.chains(), set(self.chains))
        self.assertEqual(complex_.model(), None)


    def test_can_change_complex_name(self):
        complex_ = Complex("1", "A Complex", *self.chains)
        self.assertEqual(complex_.complex_name(), "A Complex")
        complex_.complex_name("Complex2")
        self.assertEqual(complex_.complex_name(), "Complex2")


    def test_complex_name_can_only_be_changed_to_str(self):
        complex_ = Complex("1", "A Complex", *self.chains)
        with self.assertRaises(TypeError):
            complex_.complex_name(100)


    def test_complex_chains_is_read_only(self):
        chain4 = Mock(Chain)
        complex_ = Complex("1", "A Complex", *self.chains)
        self.assertEqual(len(complex_.chains()), 3)
        complex_.chains().add(chain4)
        self.assertEqual(len(complex_.chains()), 3)


    def test_can_add_chain(self):
        chain4 = Mock(Chain)
        complex_ = Complex("1", "A Complex", *self.chains)
        complex_.add_chain(chain4)
        self.assertEqual(len(complex_.chains()), 4)
        self.assertIn(chain4, complex_.chains())


    def test_can_only_add_chains(self):
        complex_ = Complex("1", "A Complex", *self.chains)
        with self.assertRaises(TypeError):
            complex_.add_chain("chain4")


    def test_can_only_add_chains_if_id_is_unique(self):
        complex_ = Complex("1", "A Complex", *self.chains[:-1])
        self.chains[-1].chain_id.return_value = "B"
        with self.assertRaises(DuplicateChainsError):
            complex_.add_chain(self.chains[-1])


    def test_can_remove_chains(self):
        complex_ = Complex("1", "A Complex", *self.chains)
        complex_.remove_chain(self.chains[1])
        self.assertEqual(len(complex_.chains()), 2)
        self.assertNotIn(self.chains[1], complex_.chains())
