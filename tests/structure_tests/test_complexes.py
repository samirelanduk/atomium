from unittest import TestCase
from unittest.mock import Mock
from molecupy.structures import Chain, Complex

class ComplexTest(TestCase):

    def setUp(self):
        self.chains = [Mock(Chain) for _ in range(3)]


class ComplexCreationTests(ComplexTest):

    def test_can_create_complex(self):
        complex_ = Complex("1", "A Complex", *self.chains)
        self.assertEqual(complex_._complex_id, "1")
        self.assertEqual(complex_._complex_name, "A Complex")
        self.assertEqual(complex_._chains, set(self.chains))
