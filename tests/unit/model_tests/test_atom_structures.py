from unittest import TestCase
from unittest.mock import patch, Mock, PropertyMock
from atomium.models.structures import AtomStructure
from atomium.models.atoms import Atom

class AtomStructureTest(TestCase):

    def setUp(self):
        self.atom1 = Mock(Atom, mass=8, element="A", name="CA")
        self.atom2 = Mock(Atom, mass=5.1, element="B", name="CA")
        self.atom3 = Mock(Atom, mass=4, element="C", name="NY")
        self.atom1.id, self.atom2.id, self.atom3.id = 500, 500, 700
        self.atoms = [self.atom1, self.atom2, self.atom3]



class AtomStructureCreationTests(AtomStructureTest):

    def test_can_create_empty_atom_structure(self):
        structure = AtomStructure()
        self.assertEqual(structure._atoms, set())
        self.assertEqual(structure._id_atoms, {})
        self.assertEqual(structure._id, None)
        self.assertEqual(structure._name, None)


    def test_can_create_atom_structure_with_atoms(self):
        structure = AtomStructure(*self.atoms)
        self.assertEqual(structure._atoms, set(self.atoms))
        self.assertEqual(structure._id_atoms, {
         500: {self.atom1, self.atom2}, 700: {self.atom3}
        })


    def test_atom_structure_will_accept_atom_structures(self):
        structure = Mock(AtomStructure)
        structure._atoms = set(self.atoms[1:])
        structure2 = AtomStructure(self.atom1, structure)
        self.assertEqual(structure2._id_atoms, {
         500: {self.atom1, self.atom2}, 700: {self.atom3}
        })
        self.assertEqual(structure2._atoms, set(self.atoms))


    def test_atoms_are_updated(self):
        try:
            AtomStructure.__name__ = "AtomStructure"
            structure = AtomStructure(self.atom1)
            self.assertFalse(hasattr(self.atom1, "_atomstructure"))
            AtomStructure.__name__ = "Model"
            structure = AtomStructure(self.atom1)
            self.assertIs(self.atom1._model, structure)
            AtomStructure.__name__ = "Ligand"
            structure = AtomStructure(self.atom1)
            self.assertIs(self.atom1._ligand, structure)
            AtomStructure.__name__ = "Residue"
            structure = AtomStructure(self.atom1)
            self.assertIs(self.atom1._residue, structure)
            AtomStructure.__name__ = "Chain"
            structure = AtomStructure(self.atom1)
            self.assertIs(self.atom1._chain, structure)
        finally: AtomStructure.__name__ = "AtomStructure"


    def test_can_create_structure_with_id(self):
        structure = AtomStructure(id=100)
        self.assertEqual(structure._id, "100")


    def test_can_create_structure_with_name(self):
        structure = AtomStructure(name=100)
        self.assertEqual(structure._name, "100")



class AtomStructureReprTests(AtomStructureTest):

    def test_atomic_structure_repr(self):
        structure = AtomStructure(self.atoms[0])
        AtomStructure.__name__ = "OBJ"
        try:
            self.assertEqual(repr(structure), "<OBJ (1 atom)>")
        finally: AtomStructure.__name__ = "AtomStructure"


    def test_atomic_structure_repr_full(self):
        structure = AtomStructure(*self.atoms, name="N", id="T")
        AtomStructure.__name__ = "OBJ"
        try:
            self.assertEqual(repr(structure), "<OBJ N (T, 3 atoms)>")
        finally: AtomStructure.__name__ = "AtomStructure"



class AtomStructureContainerTests(AtomStructureTest):

    def test_atom_structure_is_container_of_its_atoms(self):
        structure = AtomStructure(self.atom1, self.atom2)
        self.assertIn(self.atom1, structure)
        self.assertIn(self.atom2, structure)
        self.assertNotIn(self.atom3, structure)


    def test_structure_is_container_of_other_structures(self):
        structure = AtomStructure(self.atom1, self.atom2)
        other = Mock()
        other._atoms = {self.atom2}
        self.assertIn(other, structure)
        other._atoms.add(self.atom3)
        self.assertNotIn(other, structure)



class AtomStructureAtomsTests(AtomStructureTest):

    def test_can_get_all_atoms(self):
        structure = AtomStructure(*self.atoms)
        self.assertEqual(structure.atoms(), set(self.atoms))



class AtomStructureAtomTest(AtomStructureTest):

    @patch("atomium.models.molecules.AtomStructure.atoms")
    def test_atom_calls_atoms(self, mock_atoms):
        mock_atoms.return_value = {self.atoms[0]}
        structure = AtomStructure(self.atom1, self.atom2, self.atom3)
        atom = structure.atom(fufu=500, element="A")
        mock_atoms.assert_called_with(fufu=500, element="A")
        self.assertIs(atom, self.atom1)


    @patch("atomium.models.molecules.AtomStructure.atoms")
    def test_atom_can_return_none(self, mock_atoms):
        structure = AtomStructure(self.atom1, self.atom2, self.atom3)
        mock_atoms.return_value = set()
        self.assertIs(structure.atom(element="C"), None)


    @patch("atomium.models.molecules.AtomStructure.atoms")
    def test_atom_can_get_atom_by_id(self, mock_atoms):
        mock_atoms.return_value = {self.atoms[1]}
        structure = AtomStructure(self.atom1, self.atom2, self.atom3)
        atom = structure.atom(id=700)
        self.assertFalse(mock_atoms.called)
        self.assertIs(atom, self.atom3)
        atom = structure.atom(500)
        self.assertFalse(mock_atoms.called)
        self.assertIn(atom, (self.atom1, self.atom2))
        self.assertIsNone(structure.atom(id=100))



class AtomStructurePairwiseAtomTests(AtomStructureTest):

    @patch("atomium.models.molecules.AtomStructure.atoms")
    def test_can_get_pairwise_atoms(self, mock_atoms):
        mock_atoms.return_value = set(self.atoms)
        structure = AtomStructure(self.atom1, self.atom2, self.atom3)
        collected_atoms = []
        for pair in structure.pairwise_atoms(a=1, b=2):
            collected_atoms += list(pair)
        for atom in self.atoms:
            self.assertEqual(collected_atoms.count(atom), 2)
        mock_atoms.assert_called_with(a=1, b=2)


    @patch("atomium.models.molecules.AtomStructure.atoms")
    def test_can_get_pairwise_atoms_one(self, mock_atoms):
        mock_atoms.return_value = set(self.atoms[:1])
        structure = AtomStructure(self.atom1)
        self.assertEqual(len(list(structure.pairwise_atoms())), 0)



class AtomStructureAdditionTests(AtomStructureTest):

    def test_can_add_atom_to_structure(self):
        structure = AtomStructure(self.atom1)
        self.assertEqual(structure._id_atoms, {500: {self.atom1}})
        self.assertEqual(structure._atoms, {self.atom1})
        structure.add(self.atom2)
        self.assertEqual(structure._id_atoms, {500: {self.atom1, self.atom2}})
        self.assertEqual(structure._atoms, {self.atom1, self.atom2})
        structure.add(self.atom3)
        self.assertEqual(structure._id_atoms, {
         500: {self.atom1, self.atom2}, 700: {self.atom3}
        })
        self.assertEqual(structure._atoms, set(self.atoms))


    def test_can_update_atom_awareness(self):
        try:
            AtomStructure.__name__ = "AtomStructure"
            structure = AtomStructure()
            structure.add(self.atom1)
            self.assertFalse(hasattr(self.atom1, "_atomstructure"))
            AtomStructure.__name__ = "Model"
            structure.add(self.atom1)
            self.assertIs(self.atom1._model, structure)
            AtomStructure.__name__ = "Ligand"
            structure.add(self.atom1)
            self.assertIs(self.atom1._ligand, structure)
            AtomStructure.__name__ = "Residue"
            structure.add(self.atom1)
            self.assertIs(self.atom1._residue, structure)
            AtomStructure.__name__ = "Chain"
            structure.add(self.atom1)
            self.assertIs(self.atom1._chain, structure)
        finally: AtomStructure.__name__ = "AtomStructure"


    def test_can_add_structure(self):
        structure = AtomStructure(self.atom1, self.atom2)
        other = Mock(AtomStructure)
        other._atoms = {self.atom3}
        structure.add(other)
        self.assertEqual(structure._id_atoms, {
         500: {self.atom1, self.atom2}, 700: {self.atom3}
        })
        self.assertEqual(structure._atoms, set(self.atoms))



class AtomStructureRemovalTests(AtomStructureTest):

    def test_can_remove_atom_from_structure(self):
        structure = AtomStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(structure._id_atoms, {
         500: {self.atom1, self.atom2}, 700: {self.atom3}
        })
        self.assertEqual(structure._atoms, set(self.atoms))
        structure.remove(self.atom1)
        self.assertEqual(structure._id_atoms, {
         500: {self.atom2}, 700: {self.atom3}
        })
        self.assertEqual(structure._atoms, set(self.atoms[1:]))
        mock_atom = Mock(Atom)
        mock_atom.id = 300
        structure.remove(mock_atom)
        structure.remove(self.atom2)
        self.assertEqual(structure._id_atoms, {700: {self.atom3}})
        self.assertEqual(structure._atoms, {self.atom3})
        structure.remove(self.atom3)
        self.assertEqual(structure._id_atoms, {})
        self.assertEqual(structure._atoms, set())


    def test_can_update_atom_awareness(self):
        self.atom1._model, self.atom2._model, self.atom3._model = 5, 5, 5
        try:
            structure = AtomStructure(self.atom1, self.atom2, self.atom3)
            structure.remove(self.atom2)
            AtomStructure.__name__ = "Model"
            structure.remove(self.atom1)
            self.assertIs(self.atom1._model, None)
        finally: AtomStructure.__name__ = "AtomStructure"


    def test_can_remove_structure(self):
        structure = AtomStructure(self.atom1, self.atom2, self.atom3)
        other = Mock(AtomStructure)
        other._atoms = {self.atom3}
        structure.remove(other)
        self.assertEqual(structure._id_atoms, {500: {self.atom1, self.atom2}})
        self.assertEqual(structure._atoms, {self.atom1, self.atom2})



class AtomStructureGettingTests(AtomStructureTest):

    def test_can_get_objects(self):
        structure = AtomStructure(self.atom1, self.atom2, self.atom3)
        self.atom1._aaa = 10; self.atom2._aaa = 11; self.atom3._aaa = None
        self.assertEqual(structure._get("aaa"), {10, 11})
        self.atom3._aaa = 10
        self.assertEqual(structure._get("aaa"), {10, 11})


    def test_can_get_objects_by_id(self):
        objects = [Mock(_id="AAA"), Mock(_id="BBB"), Mock(_id="CCC")]
        self.atom1._aaa, self.atom2._aaa, self.atom3._aaa = objects
        structure = AtomStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(structure._get("aaa"), set(objects))
        self.assertEqual(structure._get("aaa", id="AAA"), {objects[0]})
        self.assertEqual(structure._get("aaa", id_regex="(AAA)|(BBB)"), set(objects[:2]))


    def test_can_get_objects_by_name(self):
        objects = [Mock(_name="AAA"), Mock(_name="BBB"), Mock(_name="CCC")]
        self.atom1._aaa, self.atom2._aaa, self.atom3._aaa = objects
        structure = AtomStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(structure._get("aaa"), set(objects))
        self.assertEqual(structure._get("aaa", name="AAA"), {objects[0]})
        self.assertEqual(structure._get("aaa", name_regex="(AAA)|(BBB)"), set(objects[:2]))


    def test_can_get_objects_by_water_status(self):
        objects = [Mock(_name="HOH"), Mock(_name="BBB"), Mock(_name="WAT")]
        self.atom1._ligand, self.atom2._ligand, self.atom3._ligand = objects
        structure = AtomStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(structure._get("ligand"), set(objects))
        self.assertEqual(structure._get("ligand", water=False), {objects[1]})
        self.atom1._aaa, self.atom2._aaa, self.atom3._aaa = objects
        self.assertEqual(structure._get("aaa", water=False), set(objects))



class AtomStructureIdTests(AtomStructureTest):

    def test_structure_id_property(self):
        structure = AtomStructure(self.atom1, self.atom2, self.atom3, id="B10C")
        self.assertIs(structure._id, structure.id)



class AtomStructureNameTests(AtomStructureTest):

    def test_structure_name_property(self):
        structure = AtomStructure(self.atom1, self.atom2, self.atom3, name="VAL")
        self.assertIs(structure._name, structure.name)


    def test_can_update_structure_name(self):
        structure = AtomStructure(self.atom1, self.atom2, self.atom3, name="VAL")
        structure.name = "HIS"
        self.assertEqual(structure._name, "HIS")



class AtomStructureTrimmingTests(AtomStructureTest):

    def test_can_trim_structure(self):
        structure = AtomStructure(self.atom1, self.atom2, self.atom3)
        structure.trim(108)
        self.atom1.trim.assert_called_with(108)
        self.atom2.trim.assert_called_with(108)
        self.atom3.trim.assert_called_with(108)



class AtomStructureTranslationTests(AtomStructureTest):

    def test_structure_translation(self):
        structure = AtomStructure(self.atom1, self.atom2, self.atom3)
        structure.translate(5, 4, -2, a=1, b=2)
        for atom in self.atoms:
            atom.translate.assert_called_with(5, 4, -2, a=1, b=2)



class AtomStructureTransformationTests(AtomStructureTest):

    def test_structure_transformation(self):
        structure = AtomStructure(self.atom1, self.atom2, self.atom3)
        structure.transform(5, 4, -2, a=1, b=2)
        for atom in self.atoms:
            atom.transform.assert_called_with(5, 4, -2, a=1, b=2)



class AtomStructureRotationTests(AtomStructureTest):

    def test_structure_rotation(self):
        structure = AtomStructure(self.atom1, self.atom2, self.atom3)
        structure.rotate(5, 4, -2, a=1, b=2)
        for atom in self.atoms:
            atom.rotate.assert_called_with(5, 4, -2, a=1, b=2)



class AtomStructureMassTests(AtomStructureTest):

    def test_structure_mass_is_sum_of_atom_masses(self):
        self.atom1.mass = 0.2
        self.atom2.mass = 16
        self.atom3.mass = 0.9
        structure = AtomStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(structure.mass, 17.1)



class AtomStructureChargeTests(AtomStructureTest):

    def test_structure_charge_is_sum_of_atom_charges(self):
        self.atom1.charge = 0.2
        self.atom2.charge = -1.4
        self.atom3.charge = 0.6
        structure = AtomStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(structure.charge, -0.6)



class AtomStructureFormulaTests(AtomStructureTest):

    def test_can_get_formula(self):
        self.atom1.element = "A"
        self.atom2.element = "B"
        self.atom3.element = "B"
        structure = AtomStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(structure.formula, {"A": 1, "B": 2})



class AtomStructureCenterOfMassTests(AtomStructureTest):

    @patch("atomium.models.structures.AtomStructure.mass", new_callable=PropertyMock)
    def test_can_get_center_of_mass_when_equal_mass(self, mock_mass):
        mock_mass.return_value = 20
        self.atom1.x, self.atom1.y, self.atom1.z = 0, 0, 0
        self.atom2.x, self.atom2.y, self.atom2.z = 1, 1, 1
        self.atom1.mass = 10
        self.atom2.mass = 10
        structure = AtomStructure(self.atom1, self.atom2)
        self.assertEqual(structure.center_of_mass, (0.5, 0.5, 0.5))


    @patch("atomium.models.structures.AtomStructure.mass", new_callable=PropertyMock)
    def test_can_get_center_of_mass_when_unequal_mass(self, mock_mass):
        mock_mass.return_value = 40
        self.atom1.x, self.atom1.y, self.atom1.z = 0, 0, 0
        self.atom2.x, self.atom2.y, self.atom2.z = 1, 1, 1
        self.atom1.mass = 10
        self.atom2.mass = 30
        structure = AtomStructure(self.atom1, self.atom2)
        self.assertEqual(structure.center_of_mass, (0.75, 0.75, 0.75))



class AtomStructureRadiusOfGyrationTests(AtomStructureTest):

    @patch("atomium.models.structures.AtomStructure.center_of_mass", new_callable=PropertyMock)
    def test_can_get_radius_of_gyration(self, mock_center):
        mock_center.return_value = (5, 0, 0)
        self.atom1.distance_to.return_value = 5
        self.atom2.distance_to.return_value = 5
        self.atom1.x, self.atom1.y, self.atom1.z = 0, 0, 0
        self.atom2.x, self.atom2.y, self.atom2.z = 10, 0, 0
        structure = AtomStructure(self.atom1, self.atom2)
        self.assertEqual(structure.radius_of_gyration, 5)
        self.atom1.distance_to.assert_called_with((5, 0, 0))
        self.atom2.distance_to.assert_called_with((5, 0, 0))



class AtomStructurePairingTests(AtomStructureTest):

    def setUp(self):
        AtomStructureTest.setUp(self)
        self.atoms = [Mock(Atom) for _ in range(10)]
        self.other_atoms = [Mock(Atom) for _ in range(10)]
        self.structure1 = AtomStructure(*self.atoms)
        self.structure2 = Mock(AtomStructure)
        self.structure2._atoms = set(self.other_atoms)
        for i, atom1, atom2 in zip(range(10), self.atoms, self.other_atoms):
            atom1.element = atom2.element = chr(i + 65)
            atom1.bonded_atoms = atom2.bonded_atoms = []


    def test_paired_structure_must_be_equal_length(self):
        self.structure2._atoms = set(self.other_atoms[:-1])
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
            atom1.element = atom2.element = elements[i]
            atom1.name = atom2.name = names[i]
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
            atom1.element = atom2.element = elements[i]
            atom1.name = atom2.name = names[i]
            atom1.bonded_atoms = atom2.bonded_atoms = [
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
            atom1.element = atom2.element = elements[i]
            atom1.name = atom2.name = names[i]
            atom1.bonds = atom2.bonds = [
             "bond" for _ in range(bond_counts[i])
            ]
            atom1.id = atom2.id = ids[i]
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
            atom1.element = atom2.element = elements[i]
            atom1.name = atom2.name = names[i]
            atom1.bonds = atom2.bonds = [
             "bond" for _ in range(bond_counts[i])
            ]
            atom1.id = atom2.id = ids[i]
        self.assertEqual(
         self.structure1.pairing_with(self.structure2), {a1: a2 for a1, a2 in zip(
          sorted(self.atoms, key=lambda a: id(a)),
          sorted(self.other_atoms, key=lambda a: id(a))
        )})



class AtomStructureSuperimposingTests(AtomStructureTest):

    @patch("atomium.models.structures.AtomStructure.pairing_with")
    @patch("atomium.models.structures.AtomStructure.center_of_mass", new_callable=PropertyMock)
    @patch("atomium.models.structures.AtomStructure.translate")
    def test_can_superimpose(self, mock_tran, mock_cntr, mock_pair):
        m1, m2, m3 = Mock(), Mock(), Mock()
        mock_pair.return_value = {self.atom1: m1, self.atom2: m2, self.atom3: m3}
        other = Mock()
        other.center_of_mass = (10, 20, 30)
        mock_cntr.return_value = (1, 2, 3)
        self.atom1.location, m1.x, m1.y, m1.z = (0.1, 0.2, 0.3), 0.6, 0.7, 0.8
        self.atom2.location, m2.x, m2.y, m2.z = (0.2, 0.3, 0.4), 0.7, 0.8, 0.9
        self.atom3.location, m3.x, m3.y, m3.z = (0.3, 0.4, 0.5), 0.8, 0.9, 1
        structure = AtomStructure(self.atom1, self.atom2, self.atom3)
        structure.superimpose_onto(other)
        mock_tran.assert_any_call(-1, -2, -3)
        mock_pair.assert_called_with(other)
        self.assertAlmostEqual(self.atom1.move_to.call_args_list[0][0][0], -0.135, delta=0.005)
        self.assertAlmostEqual(self.atom1.move_to.call_args_list[0][0][1], -0.208, delta=0.005)
        self.assertAlmostEqual(self.atom1.move_to.call_args_list[0][0][2], -0.280, delta=0.005)
        self.assertAlmostEqual(self.atom2.move_to.call_args_list[0][0][0], -0.138, delta=0.005)
        self.assertAlmostEqual(self.atom2.move_to.call_args_list[0][0][1], -0.286, delta=0.005)
        self.assertAlmostEqual(self.atom2.move_to.call_args_list[0][0][2], -0.434, delta=0.005)
        self.assertAlmostEqual(self.atom3.move_to.call_args_list[0][0][0], -0.142, delta=0.005)
        self.assertAlmostEqual(self.atom3.move_to.call_args_list[0][0][1], -0.365, delta=0.005)
        self.assertAlmostEqual(self.atom3.move_to.call_args_list[0][0][2], -0.589, delta=0.005)
        mock_tran.assert_any_call(10, 20, 30)



class AtomStructureTestRmsdTests(AtomStructureTest):

    @patch("atomium.models.structures.AtomStructure.pairing_with")
    def test_can_get_rmsd(self, mock_pair):
        structure = AtomStructure(self.atom1, self.atom2, self.atom3)
        other = Mock(AtomStructure)
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


    @patch("atomium.models.structures.AtomStructure.pairing_with")
    @patch("atomium.models.structures.AtomStructure.superimpose_onto")
    def test_can_get_rmsd_after_superposition(self, mock_onto, mock_pair):
        structure = AtomStructure(self.atom1, self.atom2, self.atom3)
        other = Mock(AtomStructure)
        other_atoms = [Mock(), Mock(), Mock()]
        self.atom1.distance_to.return_value = 5
        self.atom2.distance_to.return_value = 1
        self.atom3.distance_to.return_value = -1
        self.atom1.location = (1, 2, 3)
        self.atom2.location = (4, 5, 6)
        self.atom3.location = (7, 8, 9)
        mock_pair.return_value = {
         self.atom1: other_atoms[0],
         self.atom2: other_atoms[1],
         self.atom3: other_atoms[2]
        }
        rmsd = structure.rmsd_with(other, superimpose=True)
        mock_onto.assert_called_with(other)
        self.atom1.move_to.assert_called_with(1, 2, 3)
        self.atom2.move_to.assert_called_with(4, 5, 6)
        self.atom3.move_to.assert_called_with(7, 8, 9)



class AtomStructureGridTests(AtomStructureTest):

    def setUp(self):
        AtomStructureTest.setUp(self)
        self.atom1.location = (1, 1.1, 3)
        self.atom2.location = (-1, -2, -3)
        self.atom3.location = (1.5, -2.4, 1)


    def test_can_get_grid(self):
        structure = AtomStructure(self.atom1, self.atom2, self.atom3)
        grid = list(structure.grid())
        self.assertEqual(grid, [(x, y, z) for x in range(-1, 3)
         for y in range(-3, 3) for z in range(-3, 4)])


    def test_can_vary_grid_size(self):
        structure = AtomStructure(self.atom1, self.atom2, self.atom3)
        grid = list(structure.grid(size=0.5))
        self.assertEqual(grid, [(x, y, z)
         for x in [-1, -0.5, 0, 0.5, 1, 1.5]
         for y in [-2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5]
         for z in [-3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3]])


    def test_can_get_grid_with_margin(self):
        structure = AtomStructure(self.atom1, self.atom2, self.atom3)
        grid = list(structure.grid(margin=5))
        self.assertEqual(grid, [(x, y, z) for x in range(-6, 8)
         for y in range(-8, 8) for z in range(-8, 9)])



class AtomSphereTests(AtomStructureTest):

    def setUp(self):
        AtomStructureTest.setUp(self)
        self.atom1.distance_to.return_value = 5
        self.atom2.distance_to.return_value = 10
        self.atom3.distance_to.return_value = 15


    def test_can_get_atoms_in_sphere(self):
        structure = AtomStructure(self.atom1, self.atom2, self.atom3)
        atoms = structure.atoms_in_sphere(1, 2, 3, 10)
        self.assertEqual(atoms, {self.atom1, self.atom2})
        self.atom1.distance_to.assert_called_with((1, 2, 3))
        self.atom2.distance_to.assert_called_with((1, 2, 3))
        self.atom3.distance_to.assert_called_with((1, 2, 3))



class NearbyAtomTests(AtomStructureTest):

    def test_can_get_nearby_atoms(self):
        structure = AtomStructure(self.atom1, self.atom2, self.atom3)
        self.atom1.nearby_atoms.return_value = {1, 2, 3, self.atom2}
        self.atom2.nearby_atoms.return_value = {1, 3, 7, self.atom1, self.atom3}
        self.atom3.nearby_atoms.return_value = {1, 3, 7, 9, self.atom2}
        self.assertEqual(structure.nearby_atoms(2, a=4), {1, 2, 3, 7, 9})
        self.atom1.nearby_atoms.assert_called_with(2, a=4)
        self.atom2.nearby_atoms.assert_called_with(2, a=4)
        self.atom3.nearby_atoms.assert_called_with(2, a=4)



class NearbyResidueTests(AtomStructureTest):

    @patch("atomium.structures.AtomStructure._get")
    def test_can_get_nearby_residues(self, mock_get):
        mock_get.side_effect = [{1, 2}, {3}]
        structure = AtomStructure(self.atom1, self.atom2, self.atom3)
        self.atom1.nearby_residues.return_value = {1, 2, 3, 6}
        self.atom2.nearby_residues.return_value = {1, 3, 7}
        self.atom3.nearby_residues.return_value = {1, 3, 7, 9}
        self.assertEqual(structure.nearby_residues(2, a=4), {6, 7, 9})
        self.atom1.nearby_residues.assert_called_with(2, a=4)
        self.atom2.nearby_residues.assert_called_with(2, a=4)
        self.atom3.nearby_residues.assert_called_with(2, a=4)
        mock_get.assert_any_call("ligand")
        mock_get.assert_any_call("residue")



class AtomStructureCopyTests(AtomStructureTest):

    def test_can_create_copy_of_atomic_structure(self):
        new_atoms = [Mock(Atom), Mock(Atom)]
        self.atom1.copy.return_value, self.atom2.copy.return_value = new_atoms
        new_atoms[0].id = 0
        new_atoms[1].id = 0
        structure = AtomStructure(self.atom1, self.atom2)
        structure._id, structure._name = 10, 20
        copy = structure.copy()
        self.assertEqual(copy._id_atoms, {0: set(new_atoms)})
        self.assertEqual(copy._atoms, set(new_atoms))
        self.assertEqual(copy._id, 10)
        self.assertEqual(copy._name, 20)



class AtomStructureToStringTests(AtomStructureTest):

    @patch("atomium.files.pdb2pdbdict.structure_to_pdb_dict")
    @patch("atomium.files.pdbdict2pdbstring.pdb_dict_to_pdb_string")
    def test_can_save_as_pdb_string_with_description(self, mock_string, mock_dict):
        structure = AtomStructure(*self.atoms)
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
        structure = AtomStructure(*self.atoms)
        xyz_dict = {}
        mock_string.return_value = "filecontents"
        mock_dict.return_value = xyz_dict
        s = structure.to_file_string("xyz", description="DDD")
        mock_dict.assert_called_with(structure)
        mock_string.assert_called_with({"title": "DDD"})
        self.assertEqual(s, "filecontents")


    def test_invalid_file_format_is_error(self):
        structure = AtomStructure(*self.atoms)
        with self.assertRaises(ValueError):
            structure.to_file_string("nosuchfile")



class AtomStructureSavingTests(AtomStructureTest):

    @patch("atomium.models.structures.AtomStructure.to_file_string")
    @patch("atomium.files.utilities.string_to_file")
    def test_saving_uses_correct_functions(self, mock_save, mock_string):
        mock_string.return_value = "filestring"
        structure = AtomStructure(*self.atoms)
        structure.save("file.xyz", "a description")
        mock_string.assert_called_with("xyz", "a description")
        mock_save.assert_called_with("filestring", "file.xyz")


    @patch("atomium.models.structures.AtomStructure.to_file_string")
    @patch("atomium.files.utilities.string_to_file")
    def test_file_format_extracted(self, mock_save, mock_string):
        mock_string.return_value = "filestring"
        structure = AtomStructure(*self.atoms)
        structure.save("path/to/file.dfghdfg", "a description")
        mock_string.assert_called_with("dfghdfg", "a description")
        mock_save.assert_called_with("filestring", "path/to/file.dfghdfg")
