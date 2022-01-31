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



class RemovalTests(TestCase):

    def test_residue_can_remove_contents(self):
        model = atomium.open("tests/integration/files/1cbn.cif").model
        residue = model.residue("A.5")
        atom = residue.atom()
        self.assertEqual(len(residue.atoms()), 14)
        self.assertIn(atom, residue)
        residue.remove(atom)
        self.assertNotIn(atom, residue)
        self.assertEqual(len(residue.atoms()), 13)
    

    def test_non_polymer_can_remove_contents(self):
        model = atomium.open("tests/integration/files/1cbn.cif").model
        molecule = model.non_polymer("B")
        atom = molecule.atom()
        self.assertEqual(len(molecule.atoms()), 3)
        self.assertIn(atom, molecule)
        molecule.remove(atom)
        self.assertNotIn(atom, molecule)
        self.assertEqual(len(molecule.atoms()), 2)
    

    def test_water_can_remove_contents(self):
        model = atomium.open("tests/integration/files/1lol.cif").model
        molecule = model.water("G")
        atom = molecule.atom()
        self.assertEqual(len(molecule.atoms()), 96)
        self.assertIn(atom, molecule)
        molecule.remove(atom)
        self.assertNotIn(atom, molecule)
        self.assertEqual(len(molecule.atoms()), 95)
    

    def test_polymer_can_remove_contents(self):
        model = atomium.open("tests/integration/files/1cbn.cif").model
        polymer = model.polymer("A")

        # Remove residue
        residue = polymer.residue("A.20")
        self.assertEqual(len(polymer.residues()), 46)
        self.assertIn(residue, polymer)
        polymer.remove(residue)
        self.assertNotIn(residue, polymer)
        self.assertEqual(len(polymer.residues()), 45)

        # Remove atom
        atom = polymer.atom()
        self.assertEqual(len(polymer.atoms()), 630)
        self.assertIn(atom, polymer)
        polymer.remove(atom)
        self.assertNotIn(atom, polymer)
        self.assertEqual(len(polymer.atoms()), 629)
    

    def test_branched_polymer_can_remove_contents(self):
        model = atomium.open("tests/integration/files/6xlu.cif").model
        branched_polymer = model.branched_polymer("D")

        # Remove residue
        residue = branched_polymer.residue("D.1")
        self.assertEqual(len(branched_polymer.residues()), 2)
        self.assertIn(residue, branched_polymer)
        branched_polymer.remove(residue)
        self.assertNotIn(residue, branched_polymer)
        self.assertEqual(len(branched_polymer.residues()), 1)

        # Remove atom
        atom = branched_polymer.atom()
        self.assertEqual(len(branched_polymer.atoms()), 14)
        self.assertIn(atom, branched_polymer)
        branched_polymer.remove(atom)
        self.assertNotIn(atom, branched_polymer)
        self.assertEqual(len(branched_polymer.atoms()), 13)
        
    

    def test_model_can_remove_contents(self):
        # Remove polymer
        model = atomium.open("tests/integration/files/6xlu.cif").model
        polymer = model.polymer("B")
        self.assertEqual(len(model.polymers()), 3)
        self.assertIn(polymer, model)
        model.remove(polymer)
        self.assertNotIn(polymer, model)
        self.assertEqual(len(model.polymers()), 2)

        # Remove branched-polymer
        branched_polymer = model.branched_polymer("F")
        self.assertEqual(len(model.branched_polymers()), 15)
        self.assertIn(branched_polymer, model)
        model.remove(branched_polymer)
        self.assertNotIn(branched_polymer, model)
        self.assertEqual(len(model.branched_polymers()), 14)

        # Remove water
        water = model.water("YA")
        self.assertEqual(len(model.waters()), 3)
        self.assertIn(water, model)
        model.remove(water)
        self.assertNotIn(water, model)
        self.assertEqual(len(model.waters()), 2)

        # Remove non-polymer
        non_polymer = model.non_polymer("V")
        self.assertEqual(len(model.non_polymers()), 32)
        self.assertIn(non_polymer, model)
        model.remove(non_polymer)
        self.assertNotIn(non_polymer, model)
        self.assertEqual(len(model.non_polymers()), 31)

        # Remove residue
        residue = model.residue("A.20")
        self.assertEqual(len(model.residues()), 2145)
        self.assertIn(residue, model)
        model.remove(residue)
        self.assertNotIn(residue, model)
        self.assertEqual(len(model.residues()), 2144)

        # Remove atom
        atom = model.atom(1)
        self.assertEqual(len(model.atoms()), 17528)
        self.assertIn(atom, model)
        model.remove(atom)
        self.assertNotIn(atom, model)
        self.assertEqual(len(model.atoms()), 17527)

        

    '''def test_can_remove_atom(self):
        # Delete atom from polymer
        model = atomium.open("tests/integration/files/1cbn.cif").model
        atom = model.atom(1)
        self.assertEqual(len(model.atoms()), 640)
        atom.remove()
        self.assertEqual(len(model.atoms()), 639)

        # Delete atom fron non-polymer
        mol = model.non_polymer()
        atom = mol.atom()
        self.assertEqual(len(mol.atoms()), 3)
        atom.remove()
        self.assertEqual(len(mol.atoms()), 2)'''