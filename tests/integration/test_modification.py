from unittest import TestCase
import atomium

class TrimTests(TestCase):

    def test_atom_trimming(self):
        model = atomium.open("tests/integration/files/1cbn.cif").model
        atom = model.atom(1)
        atom.trim(3)
        self.assertEqual(list(atom.location), [16.864, 14.059, 3.442])
        atom.trim(2)
        self.assertEqual(list(atom.location), [16.86, 14.06, 3.44])
        atom.trim(1)
        self.assertEqual(list(atom.location), [16.9, 14.1, 3.4])
        atom.trim(0)
        self.assertEqual(list(atom.location), [17, 14, 3])
    

    def test_structure_trimming(self):
        model = atomium.open("tests/integration/files/1cbn.cif").model
        mol = model.non_polymer()
        for places in [3, 2, 1]:
            mol.trim(places)
            for atom in mol.atoms():
                for pos in atom.location:
                    self.assertLessEqual(len(str(float(pos)).split(".")[-1]), places)