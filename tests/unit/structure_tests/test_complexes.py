from unittest import TestCase
from unittest.mock import patch, Mock
from atomium.structures.models import Complex
from atomium.structures.molecules import AtomicStructure
from atomium.structures.atoms import Atom

class ComplexTest(TestCase):

    def setUp(self):
        self.atom1, self.atom2, self.atom3 = Mock(Atom), Mock(Atom), Mock(Atom)
        self.atoms = [self.atom1, self.atom2, self.atom3]
        def mock_init(obj, *args, **kwargs):
            obj._atoms = set(args)
        self.patch1 = patch("atomium.structures.molecules.AtomicStructure.__init__")
        self.mock_init = self.patch1.start()
        self.mock_init.side_effect = mock_init



class ComplexCreationTests(ComplexTest):

    @patch("atomium.structures.molecules.AtomicStructure.__init__")
    def test_model_is_atomic_structure(self, mock_init):
        mock_init.side_effect = self.mock_init
        cmplx = Complex(*self.atoms)
        self.assertIsInstance(cmplx, AtomicStructure)
        self.assertTrue(mock_init.called)
        self.assertIsNone(cmplx._id)
        self.assertIsNone(cmplx._name)


    def test_can_create_complex_with_id(self):
        complx = Complex(self.atom1, self.atom2, self.atom3, id="1")
        self.assertEqual(complx._id, "1")


    def test_complex_id_must_be_str(self):
        with self.assertRaises(TypeError):
            Complex(self.atom1, self.atom2, self.atom3, id=1000)


    def test_can_create_complex_with_name(self):
        complx = Complex(self.atom1, self.atom2, self.atom3, name="HEAVY")
        self.assertEqual(complx._name, "HEAVY")


    def test_complex_name_must_be_str(self):
        with self.assertRaises(TypeError):
            Complex(self.atom1, self.atom2, self.atom3, name=1000)


    def test_atoms_are_linked_to_complex(self):
        cmplx = Complex(*self.atoms)
        self.assertIs(self.atom1._complex, cmplx)
        self.assertIs(self.atom2._complex, cmplx)
        self.assertIs(self.atom3._complex, cmplx)



class ComplexReprTests(ComplexTest):

    def test_complex_repr(self):
        cmplx = Complex(*self.atoms)
        self.assertEqual(str(cmplx), "<Complex (3 atoms)>")



class ComplexIdTests(ComplexTest):

    def test_complex_id_property(self):
        complex = Complex(self.atom1, self.atom2, self.atom3, id="B10C")
        self.assertIs(complex._id, complex.id)



class ComplexNameTests(ComplexTest):

    def test_complex_name_property(self):
        complex = Complex(self.atom1, self.atom2, self.atom3, name="VAL")
        self.assertIs(complex._name, complex.name)


    def test_can_update_complex_name(self):
        complex = Complex(self.atom1, self.atom2, self.atom3, name="VAL")
        complex.name = "HIS"
        self.assertEqual(complex._name, "HIS")


    def test_complex_name_must_be_str(self):
        complex = Complex(self.atom1, self.atom2, self.atom3, name="VAL")
        with self.assertRaises(TypeError):
            complex.name = 10
