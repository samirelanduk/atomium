from unittest import TestCase
from unittest.mock import Mock
from molecupy.structures import Chain, Complex, ResiduicStructure, Residue, Atom

class ComplexTest(TestCase):

    def setUp(self):
        self.chains = [Mock(Chain, _complex=None) for _ in range(3)]
        self.residues = [Mock(ResiduicStructure) for _ in range(6)]
        self.atoms = [Mock(Atom) for _ in range(18)]
        for index, chain in enumerate(self.chains):
            chain.residues.return_value = set(
             self.residues[index * 2: (index * 2) + 2]
            )
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
