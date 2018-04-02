from collections import Counter
from math import pi
from unittest import TestCase
from unittest.mock import Mock, patch
from atomium.structures.atoms import Atom
from atomium.structures.molecules import AtomicStructure, Residue

class AtomicStructureTest(TestCase):

    def setUp(self):
        self.atom1 = Mock(Atom)
        self.atom2 = Mock(Atom)
        self.atom3 = Mock(Atom)
        self.atom1.mass.return_value = 8
        self.atom2.mass.return_value = 5.1
        self.atom3.mass.return_value = 4
        self.atom1.element.return_value = "A"
        self.atom2.element.return_value = "B"
        self.atom3.element.return_value = "B"
        self.atom1.atom_id.return_value = 500
        self.atom2.atom_id.return_value = None
        self.atom3.atom_id.return_value = 700
        self.atom1.name.return_value = "CA"
        self.atom2.name.return_value = "CA"
        self.atom3.name.return_value = "NY"
        self.atoms = [self.atom1, self.atom2, self.atom3]



class AtomicStructureCreationTests(AtomicStructureTest):

    def test_can_create_atomic_structure_with_id_atoms(self):
        self.atom2.atom_id.return_value = 600
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(structure._id_atoms, {
         500: self.atom1, 600: self.atom2, 700: self.atom3
        })
        self.assertEqual(structure._atoms, set())


    def test_can_create_atomic_structure_with_no_id_atoms(self):
        self.atom1.atom_id.return_value = None
        self.atom3.atom_id.return_value = None
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(structure._id_atoms, {})
        self.assertEqual(structure._atoms, set(self.atoms))


    def test_can_create_atomic_structure_with_mixed_id_atoms(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(structure._atoms, set([self.atoms[1]]))
        self.assertEqual(structure._id_atoms, {500: self.atom1, 700: self.atom3})


    def test_atomic_structure_needs_atoms(self):
        with self.assertRaises(TypeError):
            AtomicStructure("self.atom1", self.atom2, self.atom3)


    def test_atomic_structure_will_accept_atomic_structures(self):
        structure = Mock(AtomicStructure)
        structure.atoms.return_value = set(self.atoms[1:])
        structure2 = AtomicStructure(self.atom1, structure)
        self.assertEqual(structure2._atoms, set([self.atoms[1]]))
        self.assertEqual(structure2._id_atoms, {500: self.atom1, 700: self.atom3})



class AtomicStructureReprTests(AtomicStructureTest):

    def test_atomic_structure_repr(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(str(structure), "<AtomicStructure (3 atoms)>")



class AtomicStructureContainerTests(AtomicStructureTest):

    def test_atomic_structure_is_container_of_its_atoms(self):
        structure = AtomicStructure(self.atom1, self.atom2)
        self.assertIn(self.atom1, structure)
        self.assertIn(self.atom2, structure)
        self.assertNotIn(self.atom3, structure)



class AtomicStructureAtomsTests(AtomicStructureTest):

    def test_can_get_all_atoms(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(structure.atoms(), set(self.atoms))



class AtomicStructureAtomTest(AtomicStructureTest):

    @patch("atomium.structures.molecules.AtomicStructure.atoms")
    def test_atom_calls_atoms(self, mock_atoms):
        mock_atoms.return_value = set([self.atoms[0]])
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        atom = structure.atom(fufu=500, element="A")
        mock_atoms.assert_called_with(fufu=500, element="A")
        self.assertIs(atom, self.atom1)


    @patch("atomium.structures.molecules.AtomicStructure.atoms")
    def test_atom_can_return_none(self, mock_atoms):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        mock_atoms.return_value = set()
        self.assertIs(structure.atom(element="C"), None)


    @patch("atomium.structures.molecules.AtomicStructure.atoms")
    def test_atom_can_get_atom_by_id(self, mock_atoms):
        mock_atoms.return_value = set([self.atoms[1]])
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        atom = structure.atom(atom_id=500)
        self.assertFalse(mock_atoms.called)
        self.assertIs(atom, self.atom1)
        atom = structure.atom(500)
        self.assertFalse(mock_atoms.called)
        self.assertIs(atom, self.atom1)



class AtomicStructureAtomAdditionTests(AtomicStructureTest):

    def test_can_add_atoms_to_structure(self):
        structure = AtomicStructure(self.atom1)
        self.assertEqual(structure._atoms, set())
        self.assertEqual(structure._id_atoms, {500: self.atom1})
        structure.add_atom(self.atom2)
        self.assertEqual(structure._atoms, set([self.atom2]))
        self.assertEqual(structure._id_atoms, {500: self.atom1})
        structure.add_atom(self.atom3)
        self.assertEqual(structure._atoms, set([self.atom2]))
        self.assertEqual(structure._id_atoms, {500: self.atom1, 700: self.atom3})


    def test_can_only_add_atoms(self):
        structure = AtomicStructure(self.atom1, self.atom2)
        with self.assertRaises(TypeError):
            structure.add_atom("atom")



class AtomicStructureAtomRemovalTests(AtomicStructureTest):

    def test_can_remove_atoms_from_structure(self):
        structure = AtomicStructure(self.atom1, self.atom2)
        self.assertEqual(structure._atoms, set([self.atom2]))
        self.assertEqual(structure._id_atoms, {500: self.atom1})
        structure.remove_atom(self.atom1)
        self.assertEqual(structure._atoms, set([self.atom2]))
        self.assertEqual(structure._id_atoms, {})
        structure.remove_atom(self.atom2)
        self.assertEqual(structure._atoms, set())
        self.assertEqual(structure._id_atoms, {})



class AtomicStructurePairwiseAtomTests(AtomicStructureTest):

    @patch("atomium.structures.molecules.AtomicStructure.atoms")
    def test_can_get_pairwise_atoms(self, mock_atoms):
        mock_atoms.return_value = set(self.atoms)
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        collected_atoms = []
        for pair in structure.pairwise_atoms(a=1, b=2):
            self.assertIsNot(pair[0], pair[1])
            collected_atoms.append(pair[0])
            collected_atoms.append(pair[1])
        for atom in self.atoms:
            self.assertEqual(collected_atoms.count(atom), 2)
        mock_atoms.assert_called_with(a=1, b=2)


    @patch("atomium.structures.molecules.AtomicStructure.atoms")
    def test_can_get_pairwise_atoms_one(self, mock_atoms):
        mock_atoms.return_value = set(self.atoms[:1])
        structure = AtomicStructure(self.atom1)
        self.assertEqual(len(list(structure.pairwise_atoms())), 0)



class AtomicStructureGridTests(AtomicStructureTest):

    def setUp(self):
        AtomicStructureTest.setUp(self)
        self.patch1 = patch("atomium.structures.molecules.AtomicStructure.atoms")
        self.mock_atoms = self.patch1.start()
        self.mock_atoms.return_value = set([self.atom1, self.atom2, self.atom3])
        self.atom1.location.return_value = (1, 1.1, 3)
        self.atom2.location.return_value = (-1, -2, -3)
        self.atom3.location.return_value = (1.5, -2.4, 1)


    def tearDown(self):
        self.patch1.stop()


    def test_can_get_grid(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        grid = list(structure.grid())
        self.mock_atoms.assert_called_with()
        self.atom1.location.assert_called_with()
        self.atom2.location.assert_called_with()
        self.atom3.location.assert_called_with()
        self.assertEqual(grid, [(x, y, z) for x in range(-1, 3)
         for y in range(-3, 3) for z in range(-3, 4)])


    def test_can_vary_grid_size(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        grid = list(structure.grid(size=0.5))
        self.mock_atoms.assert_called_with()
        self.atom1.location.assert_called_with()
        self.atom2.location.assert_called_with()
        self.atom3.location.assert_called_with()
        self.assertEqual(grid, [(x, y, z)
         for x in [-1, -0.5, 0, 0.5, 1, 1.5]
         for y in [-2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5]
         for z in [-3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3]])


    def test_can_get_grid_with_margin(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        grid = list(structure.grid(margin=5))
        self.mock_atoms.assert_called_with()
        self.atom1.location.assert_called_with()
        self.atom2.location.assert_called_with()
        self.atom3.location.assert_called_with()
        self.assertEqual(grid, [(x, y, z) for x in range(-6, 8)
         for y in range(-8, 8) for z in range(-8, 9)])



class AtomicStructureMassTests(AtomicStructureTest):

    @patch("atomium.structures.molecules.AtomicStructure.atoms")
    def test_structure_mass_is_sum_of_atom_masses(self, mock_atoms):
        mock_atoms.return_value = set(self.atoms)
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(structure.mass(), 17.1)



class AtomicStructureChargeTests(AtomicStructureTest):

    @patch("atomium.structures.molecules.AtomicStructure.atoms")
    def test_structure_charge_is_sum_of_atom_charges(self, mock_atoms):
        mock_atoms.return_value = set(self.atoms)
        self.atom1.charge.return_value = 0.2
        self.atom2.charge.return_value = -1.4
        self.atom3.charge.return_value = 0.6
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertAlmostEqual(structure.charge(), -0.6, delta=0.000005)



class AtomicStructureFormulaTests(AtomicStructureTest):

    @patch("atomium.structures.molecules.AtomicStructure.atoms")
    def test_can_get_formula(self, mock_atoms):
        mock_atoms.return_value = set(self.atoms)
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(structure.formula(), Counter({"A":1, "B":2}))



class AtomicStructureRoundingTests(AtomicStructureTest):

    @patch("atomium.structures.molecules.AtomicStructure.atoms")
    def test_can_round_atom_locations(self, mock_atoms):
        mock_atoms.return_value = set(self.atoms)
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        structure.round(10)
        self.atom1.round.assert_called_with(10)
        self.atom2.round.assert_called_with(10)
        self.atom3.round.assert_called_with(10)



class AtomSphereTests(AtomicStructureTest):

    def setUp(self):
        AtomicStructureTest.setUp(self)
        self.atom1.distance_to.return_value = 5
        self.atom2.distance_to.return_value = 10
        self.atom3.distance_to.return_value = 15
        self.patch1 = patch("atomium.structures.molecules.AtomicStructure.atoms")
        self.mock_atoms = self.patch1.start()
        self.mock_atoms.return_value = set([self.atom1, self.atom2, self.atom3])


    def tearDown(self):
        self.patch1.stop()


    def test_coordinates_must_be_numbers(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        with self.assertRaises(TypeError):
            structure.atoms_in_sphere("0", 0, 0, 10)
        with self.assertRaises(TypeError):
            structure.atoms_in_sphere(0, "0", 0, 10)
        with self.assertRaises(TypeError):
            structure.atoms_in_sphere(0, 0, "0", 10)


    def test_radius_must_be_positive_number(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        with self.assertRaises(TypeError):
            structure.atoms_in_sphere(0, 0, 0, "10")
        with self.assertRaises(ValueError):
            structure.atoms_in_sphere(0, 0, 0, -10)


    def test_can_get_atoms_in_sphere(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        atoms = structure.atoms_in_sphere(1, 2, 3, 10)
        self.assertEqual(atoms, set([self.atom1, self.atom2]))
        self.mock_atoms.assert_called_with()
        self.atom1.distance_to.assert_called_with((1, 2, 3))
        self.atom2.distance_to.assert_called_with((1, 2, 3))
        self.atom3.distance_to.assert_called_with((1, 2, 3))



class AtomicStructurePairingTests(AtomicStructureTest):

    def setUp(self):
        AtomicStructureTest.setUp(self)
        self.atoms = [Mock(Atom) for _ in range(10)]
        self.other_atoms = [Mock(Atom) for _ in range(10)]
        self.structure1 = AtomicStructure(*self.atoms)
        self.patch1 = patch("atomium.structures.molecules.AtomicStructure.atoms")
        self.mock_atoms = self.patch1.start()
        self.mock_atoms.return_value = set(self.atoms)
        self.structure2 = Mock(AtomicStructure)
        self.structure2.atoms.return_value = set(self.other_atoms)
        for i, atom1, atom2 in zip(range(10), self.atoms, self.other_atoms):
            atom1.element.return_value = atom2.element.return_value = chr(i + 65)
            atom1.bonds.return_value = atom2.bonds.return_value = []


    def tearDown(self):
        self.patch1.stop()


    def test_pairing_needs_structure(self):
        with self.assertRaises(TypeError):
            self.structure1.pairing_with("string")


    def test_paired_structure_must_be_equal_length(self):
        self.structure2.atoms.return_value = set(self.other_atoms[:-1])
        with self.assertRaises(ValueError):
            self.structure1.pairing_with(self.structure2)


    def test_can_pair_by_element(self):
        self.assertEqual(self.structure1.pairing_with(self.structure2), {
         self.atoms[0]: self.other_atoms[0], self.atoms[1]: self.other_atoms[1],
         self.atoms[2]: self.other_atoms[2], self.atoms[3]: self.other_atoms[3],
         self.atoms[4]: self.other_atoms[4], self.atoms[5]: self.other_atoms[5],
         self.atoms[6]: self.other_atoms[6], self.atoms[7]: self.other_atoms[7],
         self.atoms[8]: self.other_atoms[8], self.atoms[9]: self.other_atoms[9]
        })


    def test_can_pair_by_name(self):
        elements = ["A", "A", "B", "B", "C", "C", "D", "D", "E", "E"]
        names = ["1", "2", "1", "2", "1", "2", "1", "2", "1", "2"]
        for i, atom1, atom2 in zip(range(10), self.atoms, self.other_atoms):
            atom1.element.return_value = atom2.element.return_value = elements[i]
            atom1.name.return_value = atom2.name.return_value = names[i]
        self.assertEqual(self.structure1.pairing_with(self.structure2), {
         self.atoms[0]: self.other_atoms[0], self.atoms[1]: self.other_atoms[1],
         self.atoms[2]: self.other_atoms[2], self.atoms[3]: self.other_atoms[3],
         self.atoms[4]: self.other_atoms[4], self.atoms[5]: self.other_atoms[5],
         self.atoms[6]: self.other_atoms[6], self.atoms[7]: self.other_atoms[7],
         self.atoms[8]: self.other_atoms[8], self.atoms[9]: self.other_atoms[9]
        })


    def test_can_pair_by_bond_count(self):
        elements = ["A", "A", "A", "A", "A", "B", "B", "B", "B", "B"]
        names = ["A1", "A1", "A1", "A2", "A2", "B1", "B1", "B2", "B2", "B2"]
        bond_counts = [1, 2, 3, 1, 2, 1, 2, 1, 2, 3]
        for i, atom1, atom2 in zip(range(10), self.atoms, self.other_atoms):
            atom1.element.return_value = atom2.element.return_value = elements[i]
            atom1.name.return_value = atom2.name.return_value = names[i]
            atom1.bonds.return_value = atom2.bonds.return_value = [
             "bond" for _ in range(bond_counts[i])
            ]
        self.assertEqual(self.structure1.pairing_with(self.structure2), {
         self.atoms[0]: self.other_atoms[0], self.atoms[1]: self.other_atoms[1],
         self.atoms[2]: self.other_atoms[2], self.atoms[3]: self.other_atoms[3],
         self.atoms[4]: self.other_atoms[4], self.atoms[5]: self.other_atoms[5],
         self.atoms[6]: self.other_atoms[6], self.atoms[7]: self.other_atoms[7],
         self.atoms[8]: self.other_atoms[8], self.atoms[9]: self.other_atoms[9]
        })


    def test_can_pair_by_id(self):
        elements = ["A", "A", "A", "A", "A", "B", "B", "B", "B", "B"]
        names = ["A1", "A1", "A1", "A2", "A2", "B1", "B1", "B2", "B2", "B2"]
        bond_counts = [1, 1, 3, 1, 2, 1, 2, 1, 2, 2]
        ids = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        for i, atom1, atom2 in zip(range(10), self.atoms, self.other_atoms):
            atom1.element.return_value = atom2.element.return_value = elements[i]
            atom1.name.return_value = atom2.name.return_value = names[i]
            atom1.bonds.return_value = atom2.bonds.return_value = [
             "bond" for _ in range(bond_counts[i])
            ]
            atom1.atom_id.return_value = atom2.atom_id.return_value = ids[i]
        self.assertEqual(self.structure1.pairing_with(self.structure2), {
         self.atoms[0]: self.other_atoms[0], self.atoms[1]: self.other_atoms[1],
         self.atoms[2]: self.other_atoms[2], self.atoms[3]: self.other_atoms[3],
         self.atoms[4]: self.other_atoms[4], self.atoms[5]: self.other_atoms[5],
         self.atoms[6]: self.other_atoms[6], self.atoms[7]: self.other_atoms[7],
         self.atoms[8]: self.other_atoms[8], self.atoms[9]: self.other_atoms[9]
        })


    def test_can_pair_by_memory_address(self):
        elements = ["A"] * 10
        names = ["A"] * 10
        bond_counts = [0] * 10
        ids = [1] * 10
        for i, atom1, atom2 in zip(range(10), self.atoms, self.other_atoms):
            atom1.element.return_value = atom2.element.return_value = elements[i]
            atom1.name.return_value = atom2.name.return_value = names[i]
            atom1.bonds.return_value = atom2.bonds.return_value = [
             "bond" for _ in range(bond_counts[i])
            ]
            atom1.atom_id.return_value = atom2.atom_id.return_value = ids[i]
        self.assertEqual(
         self.structure1.pairing_with(self.structure2), {a1: a2 for a1, a2 in zip(
          sorted(self.atoms, key=lambda a: id(a)),
          sorted(self.other_atoms, key=lambda a: id(a))
        )})



class AtomicStructureTestRmsdTests(AtomicStructureTest):

    @patch("atomium.structures.molecules.AtomicStructure.pairing_with")
    def test_can_get_rmsd(self, mock_pair):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        other = Mock(AtomicStructure)
        other_atoms = [Mock(), Mock(), Mock()]
        self.atom1.distance_to.return_value = 5
        self.atom2.distance_to.return_value = 1
        self.atom3.distance_to.return_value = -1
        mock_pair.return_value = {
         self.atom1: other_atoms[0],
         self.atom2: other_atoms[1],
         self.atom3: other_atoms[2]
        }
        rmsd = structure.rmsd_with(other)
        mock_pair.assert_called_with(other)
        self.atom1.distance_to.assert_called_with(other_atoms[0])
        self.atom2.distance_to.assert_called_with(other_atoms[1])
        self.atom3.distance_to.assert_called_with(other_atoms[2])
        self.assertEqual(rmsd, 3)



class AtomicStructureOrientationTests(AtomicStructureTest):

    def setUp(self):
        AtomicStructureTest.setUp(self)
        self.patch1 = patch("atomium.structures.molecules.AtomicStructure.translate")
        self.patch2 = patch("atomium.structures.molecules.AtomicStructure.rotate")
        self.mock_translate = self.patch1.start()
        self.mock_rotate = self.patch2.start()
        self.atom1.location.return_value = (1, 2, 3)
        self.atom2.location.return_value = (1, 1, 1)
        self.atom3.location.return_value = (7, 8, 9)


    def tearDown(self):
        self.patch1.stop()
        self.patch2.stop()


    def test_can_orient_to_origin(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        structure.orient(self.atom1)
        self.mock_translate.assert_called_with(-1, -2, -3)


    def test_atom1_must_be_atom(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        with self.assertRaises(TypeError):
            structure.orient("atom")


    def test_can_orient_onto_x_axis(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        structure.orient(self.atom1, atom2=self.atom2, axis="x")
        self.mock_translate.assert_called_with(-1, -2, -3)
        first_rotation_call = self.mock_rotate.call_args_list[0][0]
        self.assertAlmostEqual(first_rotation_call[0], pi / 4, delta=0.00005)
        self.assertEqual(first_rotation_call[1], "y")
        second_rotation_call = self.mock_rotate.call_args_list[1][0]
        self.assertAlmostEqual(second_rotation_call[0], -0.9553166, delta=0.005)
        self.assertEqual(second_rotation_call[1], "z")


    def test_can_orient_onto_y_axis(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        structure.orient(self.atom1, atom2=self.atom2, axis="y")
        self.mock_translate.assert_called_with(-1, -2, -3)
        first_rotation_call = self.mock_rotate.call_args_list[0][0]
        self.assertAlmostEqual(first_rotation_call[0], pi / 4, delta=0.00005)
        self.assertEqual(first_rotation_call[1], "z")
        second_rotation_call = self.mock_rotate.call_args_list[1][0]
        self.assertAlmostEqual(second_rotation_call[0], -0.9553166, delta=0.005)
        self.assertEqual(second_rotation_call[1], "x")


    def test_can_orient_onto_z_axis(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        structure.orient(self.atom1, atom2=self.atom2, axis="z")
        self.mock_translate.assert_called_with(-1, -2, -3)
        first_rotation_call = self.mock_rotate.call_args_list[0][0]
        self.assertAlmostEqual(first_rotation_call[0], pi / 4, delta=0.00005)
        self.assertEqual(first_rotation_call[1], "x")
        second_rotation_call = self.mock_rotate.call_args_list[1][0]
        self.assertAlmostEqual(second_rotation_call[0], -0.9553166, delta=0.005)
        self.assertEqual(second_rotation_call[1], "y")


    def test_atom2_must_be_atom(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        with self.assertRaises(TypeError):
            structure.orient(self.atom2, atom2="atom")


    def test_axis_must_be_valid(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        with self.assertRaises(ValueError) as e:
            structure.orient(self.atom1, atom2=self.atom2, axis="7")
        self.assertIn("axis", str(e.exception))


    def test_can_rotate_onto_plane(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        structure.orient(
         self.atom1, atom2=self.atom2, axis="z", atom3=self.atom3, plane="xz"
        )
        self.mock_translate.assert_called_with(-1, -2, -3)
        first_rotation_call = self.mock_rotate.call_args_list[0][0]
        self.assertAlmostEqual(first_rotation_call[0], pi / 4, delta=0.00005)
        self.assertEqual(first_rotation_call[1], "x")
        second_rotation_call = self.mock_rotate.call_args_list[1][0]
        self.assertAlmostEqual(second_rotation_call[0], -0.9553166, delta=0.005)
        self.assertEqual(second_rotation_call[1], "y")
        third_rotation_call = self.mock_rotate.call_args_list[2][0]
        self.assertAlmostEqual(third_rotation_call[0], -0.8519663, delta=0.005)
        self.assertEqual(third_rotation_call[1], "z")


    def test_atom3_must_be_atom(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        with self.assertRaises(TypeError):
            structure.orient(self.atom1, atom2=self.atom2, atom3="atom3")


    def test_plane_must_be_valid(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        with self.assertRaises(ValueError) as e:
            structure.orient(self.atom1, atom2=self.atom2, atom3=self.atom3, plane="dtg")
        self.assertIn("plane", str(e.exception))
        with self.assertRaises(ValueError) as e:
            structure.orient(self.atom1, atom2=self.atom2, atom3=self.atom3, plane="xyz")
        self.assertIn("plane", str(e.exception))
        with self.assertRaises(ValueError) as e:
            structure.orient(self.atom1, atom2=self.atom2, atom3=self.atom3, plane="yz")
        self.assertIn("x-axis", str(e.exception))




class AtomicStructureTranslationTests(AtomicStructureTest):

    @patch("atomium.structures.molecules.AtomicStructure.atoms")
    def test_structure_translation(self, mock_atoms):
        mock_atoms.return_value = set(self.atoms)
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        structure.translate(5, 4, -2)
        for atom in self.atoms:
            atom.translate.assert_called_with(5, 4, -2)



class AtomicStructureRotationTests(AtomicStructureTest):

    @patch("atomium.structures.molecules.AtomicStructure.atoms")
    def test_structure_rotation(self, mock_atoms):
        mock_atoms.return_value = set(self.atoms)
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        structure.rotate(0.1, "z")
        for atom in self.atoms:
            atom.rotate.assert_called_with(0.1, "z")



class AtomicStructureCenterOfMassTests(AtomicStructureTest):

    @patch("atomium.structures.molecules.AtomicStructure.atoms")
    @patch("atomium.structures.molecules.AtomicStructure.mass")
    def test_can_get_center_of_mass_when_equal_mass(self, mock_mass, mock_atoms):
        mock_atoms.return_value = set(self.atoms[:2])
        mock_mass.return_value = 20
        self.atom1._x, self.atom1._y, self.atom1._z = 0, 0, 0
        self.atom2._x, self.atom2._y, self.atom2._z = 1, 1, 1
        self.atom1.mass.return_value = 10
        self.atom2.mass.return_value = 10
        structure = AtomicStructure(self.atom1, self.atom2)
        self.assertEqual(structure.center_of_mass(), (0.5, 0.5, 0.5))


    @patch("atomium.structures.molecules.AtomicStructure.atoms")
    @patch("atomium.structures.molecules.AtomicStructure.mass")
    def test_can_get_center_of_mass_when_unequal_mass(self, mock_mass, mock_atoms):
        mock_atoms.return_value = set(self.atoms[:2])
        mock_mass.return_value = 40
        self.atom1._x, self.atom1._y, self.atom1._z = 0, 0, 0
        self.atom2._x, self.atom2._y, self.atom2._z = 1, 1, 1
        self.atom1.mass.return_value = 10
        self.atom2.mass.return_value = 30
        structure = AtomicStructure(self.atom1, self.atom2)
        self.assertEqual(structure.center_of_mass(), (0.75, 0.75, 0.75))



class AtomicStructureRadiusOfGyrationTests(AtomicStructureTest):

    @patch("atomium.structures.molecules.AtomicStructure.atoms")
    @patch("atomium.structures.molecules.AtomicStructure.center_of_mass")
    def test_can_get_radius_of_gyration(self, mock_center, mock_atoms):
        mock_atoms.return_value = set(self.atoms[:2])
        mock_center.return_value = (5, 0, 0)
        self.atom1.distance_to.return_value = 5
        self.atom2.distance_to.return_value = 5
        self.atom1._x, self.atom1._y, self.atom1._z = 0, 0, 0
        self.atom2._x, self.atom2._y, self.atom2._z = 10, 0, 0
        structure = AtomicStructure(self.atom1, self.atom2)
        self.assertEqual(structure.radius_of_gyration(), 5)
        self.atom1.distance_to.assert_called_with((5, 0, 0))
        self.atom2.distance_to.assert_called_with((5, 0, 0))



'''class AtomicStructurePatternMatchingTests(AtomicStructureTest):

    def test_pattern_must_be_structure(self):
        structure = AtomicStructure(self.atom1, self.atom2)
        with self.assertRaises(TypeError):
            structure.find_pattern("pattern")'''



class AtomicStructureCopyTests(AtomicStructureTest):

    def test_can_create_copy_of_atomic_structure(self):
        new_atoms = [Mock(Atom), Mock(Atom)]
        self.atom1.copy.return_value, self.atom2.copy.return_value = new_atoms
        new_atoms[0].atom_id.return_value = None
        new_atoms[1].atom_id.return_value = None
        structure = AtomicStructure(self.atom1, self.atom2)
        copy = structure.copy()
        self.assertEqual(copy._atoms, set(new_atoms))



class AtomicStructureToStringTests(AtomicStructureTest):

    @patch("atomium.files.pdb2pdbdict.structure_to_pdb_dict")
    @patch("atomium.files.pdbdict2pdbstring.pdb_dict_to_pdb_string")
    def test_can_save_as_pdb_string_with_description(self, mock_string, mock_dict):
        structure = AtomicStructure(*self.atoms)
        pdb_dict = {}
        mock_string.return_value = "filecontents"
        mock_dict.return_value = pdb_dict
        s = structure.to_file_string("pdb", description="TTT")
        mock_dict.assert_called_with(structure)
        mock_string.assert_called_with({"title": "TTT"})
        self.assertEqual(s, "filecontents")


    @patch("atomium.files.xyz2xyzdict.structure_to_xyz_dict")
    @patch("atomium.files.xyzdict2xyzstring.xyz_dict_to_xyz_string")
    def test_can_save_as_xyz_string_with_description(self, mock_string, mock_dict):
        structure = AtomicStructure(*self.atoms)
        xyz_dict = {}
        mock_string.return_value = "filecontents"
        mock_dict.return_value = xyz_dict
        s = structure.to_file_string("xyz", description="DDD")
        mock_dict.assert_called_with(structure)
        mock_string.assert_called_with({"title": "DDD"})
        self.assertEqual(s, "filecontents")


    def test_invalid_file_format_is_error(self):
        structure = AtomicStructure(*self.atoms)
        with self.assertRaises(ValueError):
            structure.to_file_string("nosuchfile")



class AtomicStructureSavingTests(AtomicStructureTest):

    @patch("atomium.structures.molecules.AtomicStructure.to_file_string")
    @patch("atomium.files.utilities.string_to_file")
    def test_saving_uses_correct_functions(self, mock_save, mock_string):
        mock_string.return_value = "filestring"
        structure = AtomicStructure(*self.atoms)
        structure.save("file.xyz", "a description")
        mock_string.assert_called_with("xyz", "a description")
        mock_save.assert_called_with("filestring", "file.xyz")


    @patch("atomium.structures.molecules.AtomicStructure.to_file_string")
    @patch("atomium.files.utilities.string_to_file")
    def test_file_format_extracted(self, mock_save, mock_string):
        mock_string.return_value = "filestring"
        structure = AtomicStructure(*self.atoms)
        structure.save("path/to/file.dfghdfg", "a description")
        mock_string.assert_called_with("dfghdfg", "a description")
        mock_save.assert_called_with("filestring", "path/to/file.dfghdfg")
