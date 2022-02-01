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
        self.assertIsNone(atom.residue)
    

    def test_non_polymer_can_remove_contents(self):
        model = atomium.open("tests/integration/files/1cbn.cif").model
        molecule = model.non_polymer("B")
        atom = molecule.atom()
        self.assertEqual(len(molecule.atoms()), 3)
        self.assertIn(atom, molecule)
        molecule.remove(atom)
        self.assertNotIn(atom, molecule)
        self.assertEqual(len(molecule.atoms()), 2)
        self.assertIsNone(atom.non_polymer)
    

    def test_water_can_remove_contents(self):
        model = atomium.open("tests/integration/files/1lol.cif").model
        molecule = model.water("G")
        atom = molecule.atom()
        self.assertEqual(len(molecule.atoms()), 96)
        self.assertIn(atom, molecule)
        molecule.remove(atom)
        self.assertNotIn(atom, molecule)
        self.assertEqual(len(molecule.atoms()), 95)
        self.assertIsNone(atom.water)
    

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
        self.assertIsNone(residue.polymer)

        # Remove atom
        atom = polymer.atom()
        self.assertEqual(len(polymer.atoms()), 630)
        self.assertIn(atom, polymer)
        polymer.remove(atom)
        self.assertNotIn(atom, polymer)
        self.assertEqual(len(polymer.atoms()), 629)
        self.assertIsNone(atom.polymer)
    

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
        self.assertIsNone(residue.branched_polymer)

        # Remove atom
        atom = branched_polymer.atom()
        self.assertEqual(len(branched_polymer.atoms()), 14)
        self.assertIn(atom, branched_polymer)
        branched_polymer.remove(atom)
        self.assertNotIn(atom, branched_polymer)
        self.assertEqual(len(branched_polymer.atoms()), 13)
        self.assertIsNone(atom.branched_polymer)
        
    

    def test_model_can_remove_contents(self):
        # Remove polymer
        model = atomium.open("tests/integration/files/6xlu.cif").model
        polymer = model.polymer("B")
        self.assertEqual(len(model.polymers()), 3)
        self.assertIn(polymer, model)
        model.remove(polymer)
        self.assertNotIn(polymer, model)
        self.assertEqual(len(model.polymers()), 2)
        self.assertIsNone(polymer.model)

        # Remove branched-polymer
        branched_polymer = model.branched_polymer("F")
        self.assertEqual(len(model.branched_polymers()), 15)
        self.assertIn(branched_polymer, model)
        model.remove(branched_polymer)
        self.assertNotIn(branched_polymer, model)
        self.assertEqual(len(model.branched_polymers()), 14)
        self.assertIsNone(branched_polymer.model)

        # Remove water
        water = model.water("YA")
        self.assertEqual(len(model.waters()), 3)
        self.assertIn(water, model)
        model.remove(water)
        self.assertNotIn(water, model)
        self.assertEqual(len(model.waters()), 2)
        model.dehydrate()
        self.assertEqual(len(model.waters()), 0)
        self.assertIsNone(water.model)

        # Remove non-polymer
        non_polymer = model.non_polymer("V")
        self.assertEqual(len(model.non_polymers()), 32)
        self.assertIn(non_polymer, model)
        model.remove(non_polymer)
        self.assertNotIn(non_polymer, model)
        self.assertEqual(len(model.non_polymers()), 31)
        self.assertIsNone(non_polymer.model)

        # Remove residue
        residue = model.residue("A.20")
        self.assertEqual(len(model.residues()), 2145)
        self.assertIn(residue, model)
        model.remove(residue)
        self.assertNotIn(residue, model)
        self.assertEqual(len(model.residues()), 2144)
        self.assertIsNone(residue.model)

        # Remove atom
        atom = model.atom(1)
        self.assertEqual(len(model.atoms()), 17361)
        self.assertIn(atom, model)
        model.remove(atom)
        self.assertNotIn(atom, model)
        self.assertEqual(len(model.atoms()), 17360)
        self.assertIsNone(atom.model)



class CopyingTests(TestCase):

    def test_atom_copying(self):
        model = atomium.open("tests/integration/files/1cbn.cif").model
        atom = model.atom(1)
        residue = atom.residue
        copy = atom.copy()
        self.assertIsNone(copy.residue)
        self.assertNotIn(copy, residue)
        self.assertEqual(atom.id, copy.id)
        self.assertEqual(atom.element, copy.element)
        self.assertEqual(atom.name, copy.name)
        self.assertEqual(atom.bvalue, copy.bvalue)
        self.assertEqual(atom.anisotropy, copy.anisotropy)
        self.assertEqual(list(atom.location), list(copy.location))
    

    def test_residue_copying(self):
        model = atomium.open("tests/integration/files/1cbn.cif").model
        residue = model.residue("A.20")
        polymer = residue.polymer
        copy = residue.copy()
        self.assertIsNone(copy.polymer)
        self.assertNotIn(copy, polymer)
        self.assertEqual(residue.id, copy.id)
        self.assertEqual(residue.name, copy.name)
        self.assertEqual(residue.number, copy.number)
        self.assertEqual(residue.formula, copy.formula)
        self.assertEqual(
            len(residue.atoms() | copy.atoms()),
            len(residue.atoms()) + len(copy.atoms())
        )
    

    def test_non_polymer_copying(self):
        model = atomium.open("tests/integration/files/1lol.cif").model
        molecule = model.non_polymer("C")
        model = molecule.model
        copy = molecule.copy()
        self.assertIsNone(copy.model)
        self.assertNotIn(copy, model)
        self.assertEqual(molecule.id, copy.id)
        self.assertEqual(molecule.name, copy.name)
        self.assertEqual(molecule.entity_name, copy.entity_name)
        self.assertEqual(molecule.auth_id, copy.auth_id)
        self.assertEqual(molecule.formula, copy.formula)
        self.assertEqual(
            len(molecule.atoms() | copy.atoms()),
            len(molecule.atoms()) + len(copy.atoms())
        )
    

    def test_non_polymer_copying(self):
        model = atomium.open("tests/integration/files/1lol.cif").model
        molecule = model.water("G")
        model = molecule.model
        copy = molecule.copy()
        self.assertIsNone(copy.model)
        self.assertNotIn(copy, model)
        self.assertEqual(molecule.id, copy.id)
        self.assertEqual(molecule.name, copy.name)
        self.assertEqual(molecule.entity_name, copy.entity_name)
        self.assertEqual(molecule.auth_id, copy.auth_id)
        self.assertEqual(molecule.formula, copy.formula)
        self.assertEqual(
            len(molecule.atoms() | copy.atoms()),
            len(molecule.atoms()) + len(copy.atoms())
        )
    

    def test_polymer_copying(self):
        model = atomium.open("tests/integration/files/1cbn.cif").model
        molecule = model.polymer("A")
        model = molecule.model
        copy = molecule.copy()
        self.assertIsNone(copy.model)
        self.assertNotIn(copy, model)
        self.assertEqual(molecule.id, copy.id)
        self.assertEqual(molecule._helices, copy._helices)
        self.assertEqual(molecule._strands, copy._strands)
        self.assertEqual(molecule.entity_name, copy.entity_name)
        self.assertEqual(molecule.sequence, copy.sequence)
        self.assertEqual(molecule.present_sequence, copy.present_sequence)
        self.assertEqual(molecule.auth_id, copy.auth_id)
        self.assertEqual(molecule.formula, copy.formula)
        self.assertEqual(
            len(molecule.atoms() | copy.atoms()),
            len(molecule.atoms()) + len(copy.atoms())
        )
    

    def test_branched_polymer_copying(self):
        model = atomium.open("tests/integration/files/6xlu.cif").model
        molecule = model.branched_polymer("D")
        model = molecule.model
        copy = molecule.copy()
        self.assertIsNone(copy.model)
        self.assertNotIn(copy, model)
        self.assertEqual(molecule.id, copy.id)
        self.assertEqual(molecule.entity_name, copy.entity_name)
        self.assertEqual(molecule.auth_id, copy.auth_id)
        self.assertEqual(molecule.formula, copy.formula)
        self.assertEqual(
            len(molecule.atoms() | copy.atoms()),
            len(molecule.atoms()) + len(copy.atoms())
        )
        self.assertEqual(len(copy.residues()), len(molecule.residues()))
    

    def test_model_copying(self):
        model = atomium.open("tests/integration/files/6xlu.cif").model
        copy = model.copy()
        self.assertIs(model.file, copy.file)
        self.assertEqual(model.formula, copy.formula)
        self.assertEqual(
            len(model.atoms() | copy.atoms()),
            len(model.atoms()) + len(copy.atoms())
        )
        self.assertEqual(len(copy.molecules()), len(model.molecules()))



class TranslationTests(TestCase):

    def test_atom_translation(self):
        model = atomium.open("tests/integration/files/1cbn.cif").model
        atom = model.atom(1)
        atom.translate([2, -3, 4])
        self.assertEqual(list(atom.location), [18.864, 11.059, 7.442])
    

    def test_structure_translation(self):
        model = atomium.open("tests/integration/files/1lol.cif").model
        molecule = model.non_polymer("C")
        atom = molecule.atom(3192)
        molecule.translate([-10, 1, -2])
        self.assertEqual(list(atom.location), [-7.354, 46.112, 46.995])
        self.assertEqual(
            [round(v, 3) for v in molecule.center_of_mass],
            [-8.701, 45.529, 47.839]
        )



class RotationTests(TestCase):

    def test_atom_rotation(self):
        model = atomium.open("tests/integration/files/1cbn.cif").model
        atom = model.atom(1)
        atom.rotate([[0.9842, 0.1243, 0.1259], [0.0343, 0.5638, -0.8252], [-0.1737, 0.8164, 0.5506]])
        atom.trim(3)
        self.assertEqual(list(atom.location), [18.778, 5.665, 10.444])
    

    def test_structure_rotation(self):
        model = atomium.open("tests/integration/files/1lol.cif").model
        molecule = model.non_polymer("C")
        atom = molecule.atom(3192)
        molecule.rotate([[0.3535, 0.8838, 0.3061], [0.3535, 0.1767, -0.9186], [-0.8661, 0.433, -0.25]])
        molecule.trim(3)
        self.assertEqual(list(atom.location), [55.803, -36.1, 4.993])
        self.assertEqual(
            [round(v, 3) for v in molecule.center_of_mass],
            [55.069, -37.455, 5.696]
        )