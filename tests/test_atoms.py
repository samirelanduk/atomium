from unittest import TestCase
from molecupy import exceptions
from molecupy.atomic import Atom, Molecule

class AtomTest(TestCase):

    def check_valid_atom(self, atom, check_atom_id=False, check_atom_name=False):
        self.assertIsInstance(atom, Atom)
        self.assertIsInstance(atom.x, float)
        self.assertIsInstance(atom.y, float)
        self.assertIsInstance(atom.z, float)
        self.assertIsInstance(atom.element, str)
        self.assertIsInstance(atom.covalent_bonds, set)
        if check_atom_id:
            self.assertIsInstance(atom.atom_id, int)
        if check_atom_name:
            self.assertIsInstance(atom.atom_name, str)
        if atom.molecule:
            self.assertIsInstance(atom.molecule, Molecule)
        self.assertRegex(str(atom), r"<Atom \([a-zA-Z]{1,2}\)>")



class AtomCreationTests(AtomTest):

    def test_can_create_atom(self):
        atom = Atom(1.0, 2.0, 3.0, "C")
        self.check_valid_atom(atom)


    def test_coordinates_must_be_floats(self):
        with self.assertRaises(TypeError):
            atom = Atom("1", 2.0, 3.0, "C")
        with self.assertRaises(TypeError):
            atom = Atom(1.0, "2", 3.0, "C")
        with self.assertRaises(TypeError):
            atom = Atom(1.0, 2.0, "3", "C")


    def test_element_must_be_element(self):
        with self.assertRaises(TypeError):
            atom = Atom(1.0, 2.0, 3.0, None)
        with self.assertRaises(exceptions.InvalidElementError):
            atom = Atom(1.0, 2.0, 3.0, "")
        with self.assertRaises(exceptions.InvalidElementError):
            atom = Atom(1.0, 2.0, 3.0, "XXX")


    def test_can_create_atom_with_id(self):
        atom = Atom(1.0, 2.0, 3.0, "C", atom_id=1001)
        self.check_valid_atom(atom, check_atom_id=True)


    def test_atom_id_must_be_int(self):
        with self.assertRaises(TypeError):
            atom = Atom(1.0, 2.0, 3.0, "C", atom_id=1.1)
        with self.assertRaises(TypeError):
            atom = Atom(1.0, 2.0, 3.0, "C", atom_id="1001")


    def test_can_create_atom_with_name(self):
        atom = Atom(1.0, 2.0, 3.0, "C", atom_name="CA")
        self.check_valid_atom(atom, check_atom_name=True)


    def test_atom_name_must_be_str(self):
        with self.assertRaises(TypeError):
            atom = Atom(1.0, 2.0, 3.0, "C", atom_name=1001)



class AtomBehaviorTests(AtomTest):

    def test_can_get_atom_mass(self):
        lithium = Atom(1.0, 1.0, 1.0, "Li")
        sodium = Atom(1.0, 1.0, 1.0, "Na")
        iron = Atom(1.0, 1.0, 1.0, "Fe")
        uranium = Atom(1.0, 1.0, 1.0, "U")
        self.assertAlmostEqual(lithium.get_mass(), 7, delta=0.5)
        self.assertAlmostEqual(sodium.get_mass(), 23, delta=0.5)
        self.assertAlmostEqual(iron.get_mass(), 56, delta=0.5)
        self.assertAlmostEqual(uranium.get_mass(), 238, delta=0.5)


    def test_strange_elements_have_zero_mass(self):
        mysterium = Atom(1.0, 1.0, 1.0, "My")
        self.assertEqual(mysterium.get_mass(), 0)



class AtomInteractionTests(AtomTest):

    def test_can_determine_distance_between_atoms(self):
        atom1 = Atom(-0.791, 64.789, 30.59, "O") # Atom 2621 in 1LOL
        atom2 = Atom(5.132, 63.307, 56.785, "C") # Atom 1011 in 1LOL
        pymol_calculated_distance = 26.9
        self.assertAlmostEqual(
         atom1.distance_to(atom2),
         pymol_calculated_distance,
         delta=0.01
        )



class AtomConnectionTests(AtomTest):

    def test_can_bond_atoms_together(self):
        atom1 = Atom(1.0, 1.0, 1.0, "H")
        atom2 = Atom(1.0, 1.0, 2.0, "C")
        atom3 = Atom(1.0, 1.0, 3.0, "O")
        atom1.covalent_bond_to(atom2)
        atom3.covalent_bond_to(atom2)

        self.assertIn(atom2, atom1.get_covalent_bonded_atoms())
        self.assertNotIn(atom3, atom1.get_covalent_bonded_atoms())
        self.assertIn(atom1, atom2.get_covalent_bonded_atoms())
        self.assertIn(atom3, atom2.get_covalent_bonded_atoms())
        self.assertIn(atom2, atom3.get_covalent_bonded_atoms())
        self.assertNotIn(atom1, atom3.get_covalent_bonded_atoms())
        self.assertEqual(len(atom1.covalent_bonds), 1)
        self.assertEqual(len(atom2.covalent_bonds), 2)
        self.assertEqual(len(atom3.covalent_bonds), 1)
        self.assertIn(list(atom1.covalent_bonds)[0], atom2.covalent_bonds)
        self.assertIn(list(atom3.covalent_bonds)[0], atom2.covalent_bonds)


    def test_can_only_covalent_bond_atom_to_another_atom(self):
        atom1 = Atom(1.0, 1.0, 1.0, "H")
        atom2 = "Carbon"
        with self.assertRaises(TypeError):
            atom1.covalent_bond_to(atom2)


    def test_can_remove_covalent_bonds(self):
        atom1 = Atom(1.0, 1.0, 1.0, "H")
        atom2 = Atom(1.0, 1.0, 2.0, "C")
        atom1.covalent_bond_to(atom2)
        self.assertIn(atom2, atom1.get_covalent_bonded_atoms())
        self.assertIn(atom1, atom2.get_covalent_bonded_atoms())
        self.assertEqual(len(atom1.covalent_bonds), 1)
        self.assertEqual(len(atom2.covalent_bonds), 1)
        atom2.break_covalent_bond(list(atom2.covalent_bonds)[0])
        self.assertNotIn(atom2, atom1.get_covalent_bonded_atoms())
        self.assertNotIn(atom1, atom2.get_covalent_bonded_atoms())
        self.assertEqual(len(atom1.covalent_bonds), 0)
        self.assertEqual(len(atom2.covalent_bonds), 0)


    def test_can_break_bond_with_other_atom(self):
        atom1 = Atom(1.0, 1.0, 1.0, "H")
        atom2 = Atom(1.0, 1.0, 2.0, "C")
        atom1.covalent_bond_to(atom2)
        self.assertIn(atom2, atom1.get_covalent_bonded_atoms())
        self.assertIn(atom1, atom2.get_covalent_bonded_atoms())
        self.assertEqual(len(atom1.covalent_bonds), 1)
        self.assertEqual(len(atom2.covalent_bonds), 1)
        atom1.break_covalent_bond_with(atom2)
        self.assertNotIn(atom2, atom1.get_covalent_bonded_atoms())
        self.assertNotIn(atom1, atom2.get_covalent_bonded_atoms())
        self.assertEqual(len(atom1.covalent_bonds), 0)
        self.assertEqual(len(atom2.covalent_bonds), 0)


    def test_atom_can_get_all_atoms_covalently_accessible(self):
        atoms = [Atom(x / 10, x / 10, x / 10, "C") for x in range(11)]
        atoms[0].covalent_bond_to(atoms[1])
        atoms[1].covalent_bond_to(atoms[2])
        atoms[2].covalent_bond_to(atoms[3])
        atoms[3].covalent_bond_to(atoms[4])
        atoms[4].covalent_bond_to(atoms[5])
        atoms[5].covalent_bond_to(atoms[0])
        atoms[6].covalent_bond_to(atoms[3])
        atoms[7].covalent_bond_to(atoms[6])
        atoms[7].covalent_bond_to(atoms[8])
        atoms[7].covalent_bond_to(atoms[9])
        atoms[10].covalent_bond_to(atoms[8])
        atoms[10].covalent_bond_to(atoms[9])
        for atom in atoms:
            self.assertEqual(
             atom.get_covalent_accessible_atoms(),
             set([a for a in atoms if a is not atom])
            )
