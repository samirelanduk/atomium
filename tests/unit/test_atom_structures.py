from unittest import TestCase
from unittest.mock import Mock, patch, PropertyMock, MagicMock
from atomium.structures import AtomStructure

class AtomStructureTest(TestCase):

    def setUp(self):
        self.atoms = [
         Mock(_id=10, mass=2, charge=-1, element="C"),
         Mock(_id=20, mass=3, charge=-1, element="N"),
         Mock(_id=30, mass=0.5, charge=4, element="C"),
        ]
        class Structure(AtomStructure):
            atoms = lambda s, *args, **kwargs: self.atoms
        self.structure = Structure()
        self.structure._id, self.structure._name = 1, 2



class AtomStructureIdTests(AtomStructureTest):

    def test_can_get_id(self):
        self.assertIs(self.structure.id, self.structure._id)


    def test_no_id_is_graceful(self):
        del self.structure._id
        with self.assertRaises(AttributeError) as e:
            self.structure.id
        self.assertIn("Structures don't", str(e.exception))



class AtomStructureNameTests(AtomStructureTest):

    def test_can_get_name(self):
        self.assertIs(self.structure.name, self.structure._name)


    def test_no_name_is_graceful(self):
        del self.structure._name
        with self.assertRaises(AttributeError) as e:
            self.structure.name
        self.assertIn("Structures don't", str(e.exception))


    def test_can_set_name(self):
        self.structure.name = "AQ"
        self.assertEqual(self.structure._name, "AQ")



class AtomStructureMassTests(AtomStructureTest):

    def test_can_get_mass(self):
        self.assertEqual(self.structure.mass, 5.5)



class AtomStructureChargeTests(AtomStructureTest):

    def test_can_get_charge(self):
        self.assertEqual(self.structure.charge, 2)



class AtomStructureFormulaTests(AtomStructureTest):

    def test_can_get_formula(self):
        self.assertEqual(self.structure.formula, {"C": 2, "N": 1})



class AtomicStructureCenterOfMassTests(AtomStructureTest):

    @patch("atomium.structures.AtomStructure.mass", new_callable=PropertyMock)
    def test_can_get_center_of_mass_when_equal_mass(self, mock_mass):
        mock_mass.return_value = 20
        self.atoms[0].x, self.atoms[0].y, self.atoms[0].z = 0, 0, 0
        self.atoms[1].x, self.atoms[1].y, self.atoms[1].z = 1, 1, 1
        self.atoms[0].mass = 10
        self.atoms[1].mass = 10
        self.structure.atoms = lambda: set(self.atoms[:2])
        self.assertEqual(self.structure.center_of_mass, (0.5, 0.5, 0.5))


    @patch("atomium.structures.AtomStructure.mass", new_callable=PropertyMock)
    def test_can_get_center_of_mass_when_unequal_mass(self, mock_mass):
        mock_mass.return_value = 40
        self.atoms[0].x, self.atoms[0].y, self.atoms[0].z = 0, 0, 0
        self.atoms[1].x, self.atoms[1].y, self.atoms[1].z = 1, 1, 1
        self.atoms[0].mass = 10
        self.atoms[1].mass = 30
        self.structure.atoms = lambda: set(self.atoms[:2])
        self.assertEqual(self.structure.center_of_mass, (0.75, 0.75, 0.75))



class AtomStructureRadiusOfGyrationTests(AtomStructureTest):

    @patch("atomium.structures.AtomStructure.center_of_mass", new_callable=PropertyMock)
    def test_can_get_radius_of_gyration(self, mock_center):
        mock_center.return_value = (5, 0, 0)
        self.atoms[0].distance_to.return_value = 5
        self.atoms[1].distance_to.return_value = 5
        self.atoms[0].x, self.atoms[0].y, self.atoms[0].z = 0, 0, 0
        self.atoms[1].x, self.atoms[1].y, self.atoms[1].z = 10, 0, 0
        self.structure.atoms = lambda: set(self.atoms[:2])
        self.assertEqual(self.structure.radius_of_gyration, 5)
        self.atoms[0].distance_to.assert_called_with((5, 0, 0))
        self.atoms[1].distance_to.assert_called_with((5, 0, 0))



class AtomStructurePairwiseAtomTests(AtomStructureTest):

    def test_can_get_pairwise_atoms(self):
        collected_atoms = []
        for pair in self.structure.pairwise_atoms(a=1, b=2):
            collected_atoms += list(pair)
        for atom in self.atoms:
            self.assertEqual(collected_atoms.count(atom), 2)


    def test_can_get_pairwise_atoms_one(self):
        self.structure.atoms = lambda: set(self.atoms[:1])
        self.assertEqual(len(list(self.structure.pairwise_atoms())), 0)



class AtomStructureTranslationTests(AtomStructureTest):

    @patch("atomium.structures.Atom.translate_atoms")
    @patch("atomium.structures.AtomStructure.trim")
    def test_structure_translation(self, mock_trim, mock_trans):
        self.structure.translate(5, 4, -2, trim=9)
        mock_trans.assert_called_with((5, 4, -2), *self.atoms)
        mock_trim.assert_called_with(9)


    @patch("atomium.structures.Atom.translate_atoms")
    @patch("atomium.structures.AtomStructure.trim")
    def test_structure_translation_with_vector(self, mock_trim, mock_trans):
        self.structure.translate([5, 4, -2], trim=9)
        mock_trans.assert_called_with([5, 4, -2], *self.atoms)
        mock_trim.assert_called_with(9)



class AtomStructureTransformationTests(AtomStructureTest):

    @patch("atomium.structures.Atom.transform_atoms")
    @patch("atomium.structures.AtomStructure.trim")
    def test_structure_transformation(self, mock_trim, mock_trans):
        self.structure.transform([1, 2], trim=9)
        mock_trans.assert_called_with([1, 2], *self.atoms)
        mock_trim.assert_called_with(9)



class AtomStructureRotationTests(AtomStructureTest):

    @patch("atomium.structures.Atom.rotate_atoms")
    @patch("atomium.structures.AtomStructure.trim")
    def test_structure_transformation(self, mock_trim, mock_rot):
        self.structure.rotate(1, 2, trim=9)
        mock_rot.assert_called_with(1, 2, *self.atoms)
        mock_trim.assert_called_with(9)



class AtomStructureTrimmingTests(AtomStructureTest):

    def test_can_trim_structure(self):
        self.structure.trim(108)
        self.atoms[0].trim.assert_called_with(108)
        self.atoms[1].trim.assert_called_with(108)
        self.atoms[2].trim.assert_called_with(108)



class AtomStructurePairingTests(AtomStructureTest):

    def setUp(self):
        AtomStructureTest.setUp(self)
        self.atoms = [Mock() for _ in range(10)]
        self.other_atoms = [Mock() for _ in range(10)]
        self.structure1 = self.structure
        self.structure.atoms = lambda: set(self.atoms)
        self.structure2 = Mock(AtomStructure)
        self.structure2.atoms = MagicMock()
        self.structure2.atoms.return_value = set(self.other_atoms)
        for i, atom1, atom2 in zip(range(10), self.atoms, self.other_atoms):
            atom1._element = atom2._element = chr(i + 65)


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
            atom1._element = atom2._element = elements[i]
            atom1._name = atom2._name = names[i]
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
        for i, atom1, atom2 in zip(range(10), self.atoms, self.other_atoms):
            atom1._element = atom2._element = elements[i]
            atom1._name = atom2._name = names[i]
        self.assertEqual(
         self.structure1.pairing_with(self.structure2), {a1: a2 for a1, a2 in zip(
          sorted(self.atoms, key=lambda a: id(a)),
          sorted(self.other_atoms, key=lambda a: id(a))
        )})


    def test_common_ids_matched_first(self):
        self.atoms[0]._id, self.other_atoms[1]._id = 1, 1
        self.assertEqual(self.structure1.pairing_with(self.structure2), {
         self.atoms[0]: self.other_atoms[1], self.atoms[1]: self.other_atoms[0],
         self.atoms[2]: self.other_atoms[2], self.atoms[3]: self.other_atoms[3],
         self.atoms[4]: self.other_atoms[4], self.atoms[5]: self.other_atoms[5],
         self.atoms[6]: self.other_atoms[6], self.atoms[7]: self.other_atoms[7],
         self.atoms[8]: self.other_atoms[8], self.atoms[9]: self.other_atoms[9]
        })



class AtomStructureTestRmsdTests(AtomStructureTest):

    @patch("atomium.structures.AtomStructure.pairing_with")
    @patch("atomium.structures.AtomStructure.center_of_mass", new_callable=PropertyMock)
    def test_can_get_rmsd(self, mock_center, mock_pair):
        other = Mock(AtomStructure)
        other.center_of_mass = [100, 110, 120]
        mock_center.return_value = [1000, 1010, 1020]
        self.atoms[0].location = [10, 20, 30]
        self.atoms[1].location = [40, 50, 60]
        self.atoms[2].location = [70, 80, 90]
        other_atoms = [Mock(location=[1,2,3]), Mock(location=[4,5,6]), Mock(location=[7,8,9])]
        mock_pair.return_value = {
         self.atoms[0]: other_atoms[0],
         self.atoms[1]: other_atoms[1],
         self.atoms[2]: other_atoms[2]
        }
        rmsd = self.structure.rmsd_with(other)
        mock_pair.assert_called_with(other)
        self.assertAlmostEqual(rmsd, 1480.95, delta=0.01)



class AtomStructureGridTests(AtomStructureTest):

    def setUp(self):
        AtomStructureTest.setUp(self)
        self.atoms[0].location = (1, 1.1, 3)
        self.atoms[1].location = (-1, -2, -3)
        self.atoms[2].location = (1.5, -2.4, 1)


    def test_can_get_grid(self):
        grid = list(self.structure.grid())
        self.assertEqual(grid, [(x, y, z) for x in range(-1, 3)
         for y in range(-3, 3) for z in range(-3, 4)])


    def test_can_vary_grid_size(self):
        grid = list(self.structure.grid(size=0.5))
        self.assertEqual(grid, [(x, y, z)
         for x in [-1, -0.5, 0, 0.5, 1, 1.5]
         for y in [-2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5]
         for z in [-3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3]])


    def test_can_get_grid_with_margin(self):
        grid = list(self.structure.grid(margin=5))
        self.assertEqual(grid, [(x, y, z) for x in range(-6, 8)
         for y in range(-8, 8) for z in range(-8, 9)])



class AtomsInSphereTests(AtomStructureTest):

    def test_can_get_atoms_in_sphere(self):
        self.atoms[0].distance_to.return_value = 5
        self.atoms[1].distance_to.return_value = 10
        self.atoms[2].distance_to.return_value = 15
        atoms = self.structure.atoms_in_sphere((1, 2, 3), 10)
        self.assertEqual(atoms, {self.atoms[0], self.atoms[1]})
        self.atoms[0].distance_to.assert_called_with((1, 2, 3))
        self.atoms[1].distance_to.assert_called_with((1, 2, 3))
        self.atoms[2].distance_to.assert_called_with((1, 2, 3))



class AtomStructureNearbyAtomsTests(AtomStructureTest):

    def test_can_get_nearby_atoms(self):
        self.structure.atoms = lambda: set(self.atoms)
        self.atoms[0].nearby_atoms.return_value = {1, 2, self.atoms[1]}
        self.atoms[1].nearby_atoms.return_value = {3, 1, 4}
        self.atoms[2].nearby_atoms.return_value = {1, 9, self.atoms[0]}
        self.assertEqual(self.structure.nearby_atoms(1, a=2), {1, 2, 3, 4, 9})
        for a in self.atoms:
            a.nearby_atoms.assert_called_with(1, a=2)



class AtomStructureEquivalenceTests(AtomStructureTest):

    @patch("atomium.structures.AtomStructure.pairing_with")
    def test_structures_not_equivalent_if_cant_pair(self, mock_pair):
        other = Mock()
        mock_pair.side_effect = Exception
        self.assertFalse(self.structure.equivalent_to(other))


    @patch("atomium.structures.AtomStructure.pairing_with")
    def test_structures_not_equivalent_if_atoms_not_equivalent(self, mock_pair):
        other = Mock()
        atoms = [Mock(), Mock(), Mock(), Mock()]
        atoms[0].equivalent_to.return_value = True
        atoms[2].equivalent_to.return_value = False
        mock_pair.return_value = {
         atoms[0]: atoms[1], atoms[2]: atoms[3]
        }
        self.assertFalse(self.structure.equivalent_to(other))
        atoms[2].equivalent_to.return_value = True
        self.assertTrue(self.structure.equivalent_to(other))
