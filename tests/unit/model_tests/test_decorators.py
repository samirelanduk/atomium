from unittest import TestCase
from unittest.mock import patch, Mock, PropertyMock
from atomium.models.atoms import atom_query
from atomium.models.molecules import lower, upper

class AtomQueryTests(TestCase):

    def setUp(self):
        self.atoms = [
         Mock(_element="C", _name="CA", _id=15, mass=12),
         Mock(_element="C", _name="CA", _id=16, mass=12),
         Mock(_element="P", _name="PB", _id=17, mass=36),
         Mock(_element="H", _name="H1", _id=18, mass=1),
         Mock(_element="ZN", _name="ZN", _id=19, mass=58),
        ]
        def func(a, b, c=20):
            """"""
            return self.atoms
        self.func = atom_query(func)


    def test_atoms_decorator_can_do_nothing(self):
        self.assertEqual(self.func(1, 2, c=3), set(self.atoms))


    def test_can_get_by_attribute(self):
        self.assertEqual(self.func(1, 2, element="C", c=3), set(self.atoms[:2]))
        self.assertEqual(self.func(1, 2, name="CA", c=3), set(self.atoms[:2]))
        self.assertEqual(self.func(1, 2, id=18, c=3), {self.atoms[3]})


    def test_can_get_by_property(self):
        self.assertEqual(self.func(1, 2, mass=12, c=3), set(self.atoms[:2]))


    def test_can_get_by_regex(self):
        self.assertEqual(self.func(1, 2, element_regex=r"C|P", c=3), set(self.atoms[:3]))
        self.assertEqual(self.func(1, 2, name_regex=r"^.\d$", c=3), {self.atoms[3]})


    def test_can_get_by_numerical_measure(self):
        self.assertEqual(self.func(1, 2, id__gt=17, c=3), set(self.atoms[3:]))
        self.assertEqual(self.func(1, 2, mass__ge=36, c=3), {self.atoms[2], self.atoms[4]})
        self.assertEqual(self.func(1, 2, mass__ne=12, c=3), set(self.atoms[2:]))


    def test_can_combine_queries(self):
        self.assertEqual(
         self.func(1, 2, id__ne=15, element_regex=r"^.$", c=3), set(self.atoms[1:4])
        )



class LowerTests(TestCase):

    def test_can_create_function_for_getting_lower_structures(self):
        obj = Mock()
        obj._get.return_value = [1, 2, 3]
        func1, func2 = lower("AAA")
        self.assertEqual(func1(obj, "E", f="g"), [1, 2, 3])
        obj._get.assert_called_with("AAA", "E", f="g")
        self.assertEqual(func2(obj, "E", f="g"), 1)
        obj._get.assert_called_with("AAA", "E", f="g")



class UpperTests(TestCase):

    def test_can_create_function_for_getting_upper_structure(self):
        obj = Mock()
        obj._get.return_value = {1}
        func1 = upper("AAA")
        self.assertEqual(func1(obj), 1)
        obj._get.assert_called_with("AAA")
        obj._get.return_value = {1, 2}
        self.assertEqual(func1(obj), None)
        obj._get.assert_called_with("AAA")
