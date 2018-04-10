import math
from unittest import TestCase
from unittest.mock import Mock, patch, PropertyMock
from atomium.structures.atoms import Atom
from atomium.structures.molecules import AtomicStructure, Molecule, Residue
from atomium.structures.chains import Chain

class AtomicStructureTest(TestCase):

    def setUp(self):
        self.atom1, self.atom2, self.atom3 = Mock(Atom), Mock(Atom), Mock(Atom)
        self.atom1.mass, self.atom2.mass, self.atom3.mass = 8, 5.1, 4
        self.atom1.element, self.atom2.element, self.atom3.element = "A", "B", "B"
        self.atom1.id, self.atom2.id, self.atom3.id = 500, 500, 700
        self.atom1.name, self.atom2.name, self.atom3.name = "CA", "CA", "NY"
        self.atoms = [self.atom1, self.atom2, self.atom3]



class AtomicStructureCreationTests(AtomicStructureTest):

    def test_can_create_atomic_structure_with_id_atoms(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(structure._atoms, set(self.atoms))
        self.assertEqual(structure._id_atoms, {
         500: {self.atom1, self.atom2}, 700: {self.atom3}
        })


    def test_atomic_structure_needs_atoms(self):
        with self.assertRaises(TypeError):
            AtomicStructure("self.atom1", self.atom2, self.atom3)


    def test_atomic_structure_will_accept_atomic_structures(self):
        structure = Mock(AtomicStructure)
        structure._atoms = set(self.atoms[1:])
        structure2 = AtomicStructure(self.atom1, structure)
        self.assertEqual(structure2._id_atoms, {
         500: {self.atom1, self.atom2}, 700: {self.atom3}
        })
        self.assertEqual(structure2._atoms, set(self.atoms))



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
        mock_atoms.return_value = {self.atoms[0]}
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
        mock_atoms.return_value = {self.atoms[1]}
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        atom = structure.atom(id=700)
        self.assertFalse(mock_atoms.called)
        self.assertIs(atom, self.atom3)
        atom = structure.atom(500)
        self.assertFalse(mock_atoms.called)
        self.assertIn(atom, (self.atom1, self.atom2))
        self.assertIsNone(structure.atom(id=100))



class AtomicStructureAtomAdditionTests(AtomicStructureTest):

    def test_can_add_atoms_to_structure(self):
        structure = AtomicStructure(self.atom1)
        self.assertEqual(structure._id_atoms, {500: {self.atom1}})
        self.assertEqual(structure._atoms, {self.atom1})
        structure.add_atom(self.atom2)
        self.assertEqual(structure._id_atoms, {500: {self.atom1, self.atom2}})
        self.assertEqual(structure._atoms, {self.atom1, self.atom2})
        structure.add_atom(self.atom3)
        self.assertEqual(structure._id_atoms, {
         500: {self.atom1, self.atom2}, 700: {self.atom3}
        })
        self.assertEqual(structure._atoms, set(self.atoms))


    def test_can_only_add_atoms(self):
        structure = AtomicStructure(self.atom1, self.atom2)
        with self.assertRaises(TypeError):
            structure.add_atom("atom")


    def test_can_update_atom_awareness(self):
        try:
            structure = AtomicStructure()
            structure.add_atom(self.atom1)
            self.assertFalse(hasattr(structure, "_atomicstructure"))
            AtomicStructure.__name__ = "Model"
            structure.add_atom(self.atom1)
            self.assertIs(self.atom1._model, structure)
            AtomicStructure.__name__ = "Molecule"
            structure.add_atom(self.atom1)
            self.assertIs(self.atom1._molecule, structure)
            AtomicStructure.__name__ = "Chain"
            structure.add_atom(self.atom1)
            self.assertIs(self.atom1._chain, structure)
            AtomicStructure.__name__ = "Model"
            structure.add_atom(self.atom1)
            self.assertIs(self.atom1._model, structure)
        finally: AtomicStructure.__name__ = "AtomicStructure"



class AtomicStructureAtomRemovalTests(AtomicStructureTest):

    def test_can_remove_atoms_from_structure(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(structure._id_atoms, {
         500: {self.atom1, self.atom2}, 700: {self.atom3}
        })
        self.assertEqual(structure._atoms, set(self.atoms))
        structure.remove_atom(self.atom1)
        self.assertEqual(structure._id_atoms, {
         500: {self.atom2}, 700: {self.atom3}
        })
        self.assertEqual(structure._atoms, set(self.atoms[1:]))
        mock_atom = Mock(Atom)
        mock_atom.id = 300
        structure.remove_atom(mock_atom)
        structure.remove_atom(self.atom2)
        self.assertEqual(structure._id_atoms, {700: {self.atom3}})
        self.assertEqual(structure._atoms, {self.atom3})
        structure.remove_atom(self.atom3)
        self.assertEqual(structure._id_atoms, {})
        self.assertEqual(structure._atoms, set())


    def test_can_update_atom_awareness(self):
        self.atom1._model, self.atom2._model, self.atom3._model = 5, 5, 5
        try:
            structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
            structure.remove_atom(self.atom2)
            AtomicStructure.__name__ = "Model"
            structure.remove_atom(self.atom1)
            self.assertIs(self.atom1._model, None)
        finally: AtomicStructure.__name__ = "AtomicStructure"



class AtomicStructurePairwiseAtomTests(AtomicStructureTest):

    @patch("atomium.structures.molecules.AtomicStructure.atoms")
    def test_can_get_pairwise_atoms(self, mock_atoms):
        mock_atoms.return_value = set(self.atoms)
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        collected_atoms = []
        for pair in structure.pairwise_atoms(a=1, b=2):
            collected_atoms += list(pair)
        for atom in self.atoms:
            self.assertEqual(collected_atoms.count(atom), 2)
        mock_atoms.assert_called_with(a=1, b=2)


    @patch("atomium.structures.molecules.AtomicStructure.atoms")
    def test_can_get_pairwise_atoms_one(self, mock_atoms):
        mock_atoms.return_value = set(self.atoms[:1])
        structure = AtomicStructure(self.atom1)
        self.assertEqual(len(list(structure.pairwise_atoms())), 0)



class StructureAdditionTests(AtomicStructureTest):

    @patch("atomium.structures.molecules.AtomicStructure.add_atom")
    def test_can_add_structure(self, mock_add):
        structure = AtomicStructure(self.atom1, self.atom2)
        other = Mock(AtomicStructure)
        other._atoms = {self.atom3}
        structure.add(other)
        mock_add.assert_called_with(self.atom3)


    def test_can_only_add_structures(self):
        structure = AtomicStructure(self.atom1, self.atom2)
        with self.assertRaises(TypeError):
            structure.add("self.chain4")



class StructureRemovalTests(AtomicStructureTest):

    @patch("atomium.structures.molecules.AtomicStructure.remove_atom")
    def test_can_remove_structure(self, mock_remove):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        other = Mock(AtomicStructure)
        other._atoms = {self.atom3}
        structure.remove(other)
        mock_remove.assert_called_with(self.atom3)


    def test_can_only_remove_structures(self):
        structure = AtomicStructure(self.atom1, self.atom2)
        with self.assertRaises(TypeError):
            structure.remove("self.chain4")



class StructureResiduesTests(AtomicStructureTest):

    def setUp(self):
        AtomicStructureTest.setUp(self)
        self.structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.res1, self.res2, self.res3 = Mock(), Mock(), Mock()
        self.atom1.residue, self.atom2.residue, self.atom3.residue = (
         self.res1, self.res2, self.res3
        )
        self.res1.id, self.res2.id, self.res3.id = "A1", "A2", "A3"
        self.res1.name, self.res2.name, self.res3.name = "VAL", "TYR", "TYR"


    def test_can_get_residues(self):
        self.assertEqual(
         self.structure.residues(),
         set([self.res1, self.res2, self.res3])
        )


    def test_can_filter_none_from_residues(self):
        self.atom3.residue = None
        self.assertEqual(self.structure.residues(), {self.res1, self.res2})


    def test_can_get_residues_by_id(self):
        self.assertEqual(
         self.structure.residues(id="A1"), {self.res1}
        )
        self.assertEqual(
         self.structure.residues(id="A2"), {self.res2}
        )
        self.assertEqual(self.structure.residues(id="A5"), set())


    def test_can_get_residues_by_name(self):
        self.assertEqual(
         self.structure.residues(name="VAL"), set([self.res1])
        )
        self.assertEqual(
         self.structure.residues(name="TYR"), set([self.res2, self.res3])
        )
        self.assertEqual(self.structure.residues(name="GLY"), set())



class StructureResidueTest(AtomicStructureTest):

    def setUp(self):
        AtomicStructureTest.setUp(self)
        self.structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.res1, self.res2, self.res3 = Mock(), Mock(), Mock()
        self.res1.id, self.res2.id, self.res3.id = "A1", "A2", "A3"
        self.res1.name, self.res2.name, self.res3.name = "VAL", "TYR", "TYR"


    @patch("atomium.structures.molecules.AtomicStructure.residues")
    def test_residue_calls_residues(self, mock_residues):
        mock_residues.return_value = set([self.res3])
        residue = self.structure.residue(name="A")
        mock_residues.assert_called_with(name="A")
        self.assertIs(residue, self.res3)


    @patch("atomium.structures.molecules.AtomicStructure.residues")
    def test_residue_can_return_none(self, mock_residues):
        mock_residues.return_value = set()
        self.assertIs(self.structure.residue(name="C"), None)



class StructureMoleculesTests(AtomicStructureTest):

    def setUp(self):
        AtomicStructureTest.setUp(self)
        self.structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.mol1, self.mol2, self.mol3 = Mock(), Mock(Residue), Mock(Chain)
        self.atom1.molecule, self.atom2.molecule, self.atom3.molecule = (
         self.mol1, self.mol2, self.mol3
        )
        self.mol1.id, self.mol2.id, self.mol3.id = "A1", "A2", "A3"
        self.mol1.name, self.mol2.name, self.mol3.name = "VAL", "TYR", "TYR"


    def test_can_get_molecules(self):
        self.assertEqual(
         self.structure.molecules(),
         set([self.mol1, self.mol2, self.mol3])
        )


    def test_can_filter_none_from_molecules(self):
        self.atom3.molecule = None
        self.assertEqual(self.structure.molecules(), {self.mol1, self.mol2})


    def test_can_get_molecules_by_id(self):
        self.assertEqual(
         self.structure.molecules(id="A1"), {self.mol1}
        )
        self.assertEqual(
         self.structure.molecules(id="A2"), {self.mol2}
        )
        self.assertEqual(self.structure.molecules(id="A5"), set())


    def test_can_get_molecules_by_name(self):
        self.assertEqual(
         self.structure.molecules(name="VAL"), set([self.mol1])
        )
        self.assertEqual(
         self.structure.molecules(name="TYR"), set([self.mol2, self.mol3])
        )
        self.assertEqual(self.structure.molecules(name="GLY"), set())


    def test_can_filter_out_specific_molecule_types(self):
        self.assertEqual(
         self.structure.molecules(generic=True), set([self.mol1])
        )


    def test_can_filter_out_water(self):
        self.mol2.name = "WAT"
        self.assertEqual(
         self.structure.molecules(water=False), {self.mol1, self.mol3}
        )
        self.mol2.name = "HOH"
        self.assertEqual(
         self.structure.molecules(water=False), {self.mol1, self.mol3}
        )



class StructureMoleculeTest(AtomicStructureTest):

    def setUp(self):
        AtomicStructureTest.setUp(self)
        self.structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.res1, self.res2, self.res3 = Mock(), Mock(), Mock()
        self.res1.id, self.res2.id, self.res3.id = "A1", "A2", "A3"
        self.res1.name, self.res2.name, self.res3.name = "VAL", "TYR", "TYR"


    @patch("atomium.structures.molecules.AtomicStructure.molecules")
    def test_molecule_calls_molecules(self, mock_molecules):
        mock_molecules.return_value = set([self.res3])
        molecule = self.structure.molecule(name="A")
        mock_molecules.assert_called_with(name="A")
        self.assertIs(molecule, self.res3)


    @patch("atomium.structures.molecules.AtomicStructure.molecules")
    def test_molecule_can_return_none(self, mock_molecules):
        mock_molecules.return_value = set()
        self.assertIs(self.structure.molecule(name="C"), None)



class StructureChainsTests(AtomicStructureTest):

    def setUp(self):
        AtomicStructureTest.setUp(self)
        self.structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.chain1, self.chain2, self.chain3 = Mock(), Mock(), Mock()
        self.atom1.chain, self.atom2.chain, self.atom3.chain = (
         self.chain1, self.chain2, self.chain3
        )
        self.chain1.id, self.chain2.id, self.chain3.id = "A1", "A2", "A3"
        self.chain1.name, self.chain2.name, self.chain3.name = "VAL", "TYR", "TYR"


    def test_can_get_chains(self):
        self.assertEqual(
         self.structure.chains(),
         set([self.chain1, self.chain2, self.chain3])
        )


    def test_can_filter_none_from_chains(self):
        self.atom3.chain = None
        self.assertEqual(self.structure.chains(), {self.chain1, self.chain2})


    def test_can_get_chains_by_id(self):
        self.assertEqual(
         self.structure.chains(id="A1"), {self.chain1}
        )
        self.assertEqual(
         self.structure.chains(id="A2"), {self.chain2}
        )
        self.assertEqual(self.structure.chains(id="A5"), set())


    def test_can_get_chains_by_name(self):
        self.assertEqual(
         self.structure.chains(name="VAL"), set([self.chain1])
        )
        self.assertEqual(
         self.structure.chains(name="TYR"), set([self.chain2, self.chain3])
        )
        self.assertEqual(self.structure.chains(name="GLY"), set())



class StructureChainTest(AtomicStructureTest):

    def setUp(self):
        AtomicStructureTest.setUp(self)
        self.structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.chain1, self.chain2, self.chain3 = Mock(), Mock(), Mock()
        self.chain1.id, self.chain2.id, self.chain3.id = "A1", "A2", "A3"
        self.chain1.name, self.chain2.name, self.chain3.name = "VAL", "TYR", "TYR"


    @patch("atomium.structures.molecules.AtomicStructure.chains")
    def test_chain_calls_chains(self, mock_chains):
        mock_chains.return_value = set([self.chain3])
        chain = self.structure.chain(name="A")
        mock_chains.assert_called_with(name="A")
        self.assertIs(chain, self.chain3)


    @patch("atomium.structures.molecules.AtomicStructure.chains")
    def test_chain_can_return_none(self, mock_chains):
        mock_chains.return_value = set()
        self.assertIs(self.structure.chain(name="C"), None)



class StructureComplexesTests(AtomicStructureTest):

    def setUp(self):
        AtomicStructureTest.setUp(self)
        self.structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.complex1, self.complex2, self.complex3 = Mock(), Mock(), Mock()
        self.atom1.complex, self.atom2.complex, self.atom3.complex = (
         self.complex1, self.complex2, self.complex3
        )
        self.complex1.id, self.complex2.id, self.complex3.id = "1", "2", "3"
        self.complex1.name, self.complex2.name, self.complex3.name = "AA", "CC", "CC"


    def test_can_get_complexes(self):
        self.assertEqual(
         self.structure.complexes(),
         set([self.complex1, self.complex2, self.complex3])
        )


    def test_can_filter_none_from_complexes(self):
        self.atom3.complex = None
        self.assertEqual(self.structure.complexes(), {self.complex1, self.complex2})


    def test_can_get_complexes_by_id(self):
        self.assertEqual(
         self.structure.complexes(id="1"), {self.complex1}
        )
        self.assertEqual(
         self.structure.complexes(id="2"), {self.complex2}
        )
        self.assertEqual(self.structure.complexes(id="4"), set())


    def test_can_get_complexes_by_name(self):
        self.assertEqual(
         self.structure.complexes(name="AA"), set([self.complex1])
        )
        self.assertEqual(
         self.structure.complexes(name="CC"), set([self.complex2, self.complex3])
        )
        self.assertEqual(self.structure.complexes(name="DD"), set())



class StructureComplexTest(AtomicStructureTest):

    def setUp(self):
        AtomicStructureTest.setUp(self)
        self.structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.complex1, self.complex2, self.complex3 = Mock(), Mock(), Mock()
        self.complex1.id, self.complex2.id, self.complex3.id = "1", "2", "3"
        self.complex1.name, self.complex2.name, self.complex3.name = "AA", "BB", "CC"


    @patch("atomium.structures.molecules.AtomicStructure.complexes")
    def test_complex_calls_complexes(self, mock_complexes):
        mock_complexes.return_value = set([self.complex3])
        complex = self.structure.complex(name="1")
        mock_complexes.assert_called_with(name="1")
        self.assertIs(complex, self.complex3)


    @patch("atomium.structures.molecules.AtomicStructure.complexes")
    def test_complex_can_return_none(self, mock_complexes):
        mock_complexes.return_value = set()
        self.assertIs(self.structure.complex(name="AA"), None)



class AtomicStructureTrimmingTests(AtomicStructureTest):

    @patch("atomium.structures.molecules.AtomicStructure.atoms")
    def test_can_trim_structure(self, mock_atoms):
        mock_atoms.return_value = set(self.atoms)
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        structure.trim(108)
        mock_atoms.assert_called_with()
        self.atom1.trim.assert_called_with(108)
        self.atom2.trim.assert_called_with(108)
        self.atom3.trim.assert_called_with(108)



class AtomicStructureTranslationTests(AtomicStructureTest):

    def test_structure_translation(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        structure.translate(5, 4, -2, a=1, b=2)
        for atom in self.atoms:
            atom.translate.assert_called_with(5, 4, -2, a=1, b=2)



class AtomicStructureRotationTests(AtomicStructureTest):

    @patch("atomium.structures.molecules.Atom.generate_rotation_matrix")
    @patch("atomium.structures.molecules.AtomicStructure.trim")
    def test_structure_rotation(self, mock_trim, mock_matrix):
        matrix = Mock()
        mock_matrix.return_value = matrix
        matrix.dot.side_effect = ([10, 20, 30], [40, 50, 60], [70, 80, 90])
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        structure._atoms = self.atoms
        structure.rotate(0.1, "z")
        mock_matrix.assert_called_with(None, 0.1, "z")
        for i, atom in enumerate(self.atoms):
            for j, coord in enumerate((atom._x, atom._y, atom._z)):
                self.assertEqual(coord, 10 + (30 * i) + (10 * j))
        mock_trim.assert_called_with(12)


    @patch("atomium.structures.molecules.Atom.generate_rotation_matrix")
    @patch("atomium.structures.molecules.AtomicStructure.trim")
    def test_structure_rotation_varied_trim(self, mock_trim, mock_matrix):
        matrix = Mock()
        mock_matrix.return_value = matrix
        matrix.dot.side_effect = ([10, 20, 30], [40, 50, 60], [70, 80, 90])
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        structure._atoms = self.atoms
        structure.rotate(0.1, "z", trim=9)
        mock_trim.assert_called_with(9)


    @patch("atomium.structures.molecules.Atom.generate_rotation_matrix")
    @patch("atomium.structures.molecules.AtomicStructure.trim")
    def test_structure_rotation_in_degrees(self, mock_trim, mock_matrix):
        matrix = Mock()
        mock_matrix.return_value = matrix
        matrix.dot.side_effect = ([10, 20, 30], [40, 50, 60], [70, 80, 90])
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        structure._atoms = self.atoms
        structure.rotate(0.1, "z", degrees=True)
        mock_matrix.assert_called_with(None, math.radians(0.1), "z")


    def test_rotation_axis_must_be_valid(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        with self.assertRaises(ValueError):
            structure.rotate(0.1, "s")



class AtomicStructureMassTests(AtomicStructureTest):

    def test_structure_mass_is_sum_of_atom_masses(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(structure.mass, 17.1)



class AtomicStructureChargeTests(AtomicStructureTest):

    def test_structure_charge_is_sum_of_atom_charges(self):
        self.atom1.charge = 0.2
        self.atom2.charge = -1.4
        self.atom3.charge = 0.6
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(structure.charge, -0.6)



class AtomicStructureFormulaTests(AtomicStructureTest):

    def test_can_get_formula(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        self.assertEqual(structure.formula, {"A":1, "B":2})



class AtomicStructureCenterOfMassTests(AtomicStructureTest):

    @patch("atomium.structures.molecules.AtomicStructure.mass", new_callable=PropertyMock)
    def test_can_get_center_of_mass_when_equal_mass(self, mock_mass):
        mock_mass.return_value = 20
        self.atom1.x, self.atom1.y, self.atom1.z = 0, 0, 0
        self.atom2.x, self.atom2.y, self.atom2.z = 1, 1, 1
        self.atom1.mass = 10
        self.atom2.mass = 10
        structure = AtomicStructure(self.atom1, self.atom2)
        self.assertEqual(structure.center_of_mass, (0.5, 0.5, 0.5))


    @patch("atomium.structures.molecules.AtomicStructure.mass", new_callable=PropertyMock)
    def test_can_get_center_of_mass_when_unequal_mass(self, mock_mass):
        mock_mass.return_value = 40
        self.atom1.x, self.atom1.y, self.atom1.z = 0, 0, 0
        self.atom2.x, self.atom2.y, self.atom2.z = 1, 1, 1
        self.atom1.mass = 10
        self.atom2.mass = 30
        structure = AtomicStructure(self.atom1, self.atom2)
        self.assertEqual(structure.center_of_mass, (0.75, 0.75, 0.75))



class AtomicStructureRadiusOfGyrationTests(AtomicStructureTest):

    @patch("atomium.structures.molecules.AtomicStructure.center_of_mass", new_callable=PropertyMock)
    def test_can_get_radius_of_gyration(self, mock_center):
        mock_center.return_value = (5, 0, 0)
        self.atom1.distance_to.return_value = 5
        self.atom2.distance_to.return_value = 5
        self.atom1.x, self.atom1.y, self.atom1.z = 0, 0, 0
        self.atom2.x, self.atom2.y, self.atom2.z = 10, 0, 0
        structure = AtomicStructure(self.atom1, self.atom2)
        self.assertEqual(structure.radius_of_gyration, 5)
        self.atom1.distance_to.assert_called_with((5, 0, 0))
        self.atom2.distance_to.assert_called_with((5, 0, 0))



class AtomicStructurePairingTests(AtomicStructureTest):

    def setUp(self):
        AtomicStructureTest.setUp(self)
        self.atoms = [Mock(Atom) for _ in range(10)]
        self.other_atoms = [Mock(Atom) for _ in range(10)]
        self.structure1 = AtomicStructure(*self.atoms)
        self.structure2 = Mock(AtomicStructure)
        self.structure2._atoms = set(self.other_atoms)
        for i, atom1, atom2 in zip(range(10), self.atoms, self.other_atoms):
            atom1.element = atom2.element = chr(i + 65)
            atom1.bonds = atom2.bonds = []


    def test_pairing_needs_structure(self):
        with self.assertRaises(TypeError):
            self.structure1.pairing_with("string")


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
            atom1.bonds = atom2.bonds = [
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



class AtomicStructureSuperimposingTests(AtomicStructureTest):

    @patch("atomium.structures.molecules.AtomicStructure.pairing_with")
    @patch("atomium.structures.molecules.AtomicStructure.center_of_mass", new_callable=PropertyMock)
    @patch("atomium.structures.molecules.AtomicStructure.translate")
    def test_can_superimpose(self, mock_tran, mock_cntr, mock_pair):
        m1, m2, m3 = Mock(), Mock(), Mock()
        mock_pair.return_value = {self.atom1: m1, self.atom2: m2, self.atom3: m3}
        other = Mock()
        other.center_of_mass = (10, 20, 30)
        mock_cntr.return_value = (1, 2, 3)
        self.atom1.location, m1.x, m1.y, m1.z = (0.1, 0.2, 0.3), 0.6, 0.7, 0.8
        self.atom2.location, m2.x, m2.y, m2.z = (0.2, 0.3, 0.4), 0.7, 0.8, 0.9
        self.atom3.location, m3.x, m3.y, m3.z = (0.3, 0.4, 0.5), 0.8, 0.9, 1
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
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


    @patch("atomium.structures.molecules.AtomicStructure.pairing_with")
    @patch("atomium.structures.molecules.AtomicStructure.superimpose_onto")
    def test_can_get_rmsd_after_superposition(self, mock_onto, mock_pair):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        other = Mock(AtomicStructure)
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



class AtomicStructureCopyTests(AtomicStructureTest):

    def test_can_create_copy_of_atomic_structure(self):
        new_atoms = [Mock(Atom), Mock(Atom)]
        self.atom1.copy.return_value, self.atom2.copy.return_value = new_atoms
        new_atoms[0].id = 0
        new_atoms[1].id = 0
        structure = AtomicStructure(self.atom1, self.atom2)
        structure._id, structure._name = 10, 20
        copy = structure.copy()
        self.assertEqual(copy._id_atoms, {0: set(new_atoms)})
        self.assertEqual(copy._atoms, set(new_atoms))
        self.assertEqual(copy._id, 10)
        self.assertEqual(copy._name, 20)



class AtomicStructureGridTests(AtomicStructureTest):

    def setUp(self):
        AtomicStructureTest.setUp(self)
        self.atom1.location = (1, 1.1, 3)
        self.atom2.location = (-1, -2, -3)
        self.atom3.location = (1.5, -2.4, 1)


    def test_can_get_grid(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        grid = list(structure.grid())
        self.assertEqual(grid, [(x, y, z) for x in range(-1, 3)
         for y in range(-3, 3) for z in range(-3, 4)])


    def test_can_vary_grid_size(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        grid = list(structure.grid(size=0.5))
        self.assertEqual(grid, [(x, y, z)
         for x in [-1, -0.5, 0, 0.5, 1, 1.5]
         for y in [-2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5]
         for z in [-3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3]])


    def test_can_get_grid_with_margin(self):
        structure = AtomicStructure(self.atom1, self.atom2, self.atom3)
        grid = list(structure.grid(margin=5))
        self.assertEqual(grid, [(x, y, z) for x in range(-6, 8)
         for y in range(-8, 8) for z in range(-8, 9)])



class AtomSphereTests(AtomicStructureTest):

    def setUp(self):
        AtomicStructureTest.setUp(self)
        self.atom1.distance_to.return_value = 5
        self.atom2.distance_to.return_value = 10
        self.atom3.distance_to.return_value = 15


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
        self.assertEqual(atoms, {self.atom1, self.atom2})
        self.atom1.distance_to.assert_called_with((1, 2, 3))
        self.atom2.distance_to.assert_called_with((1, 2, 3))
        self.atom3.distance_to.assert_called_with((1, 2, 3))



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
