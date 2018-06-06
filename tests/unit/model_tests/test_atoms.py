import math
from unittest import TestCase
from unittest.mock import patch, Mock, PropertyMock
from atomium.models.atoms import Atom

class AtomCreationTests(TestCase):

    def test_can_create_atom(self):
        atom = Atom("C")
        self.assertEqual(atom._element, "C")
        self.assertEqual((atom._x, atom._y, atom._z), (0, 0, 0))
        self.assertEqual(atom._id, 0)
        self.assertEqual(atom._name, None)
        self.assertEqual(atom._charge, 0)
        self.assertEqual(atom._bfactor, 0)
        self.assertEqual(atom._anisotropy, [0, 0, 0, 0, 0, 0])
        self.assertEqual(atom._bonded_atoms, set())
        self.assertEqual(atom._residue, None)
        self.assertEqual(atom._ligand, None)
        self.assertEqual(atom._chain, None)
        self.assertEqual(atom._complex, None)
        self.assertEqual(atom._model, None)


    def test_can_create_atom_with_location(self):
        atom = Atom("C", x=2.3, y=5, z=-9)
        self.assertEqual((atom._x, atom._y, atom._z), (2.3, 5, -9))


    def test_can_create_atom_with_id(self):
        atom = Atom("C", id=2.0)
        self.assertEqual(atom._id, 2)


    def test_can_create_atom_with_name(self):
        atom = Atom("C", name="N2")
        self.assertEqual(atom._name, "N2")


    def test_can_create_atom_with_charge(self):
        atom = Atom("C", charge=-1)
        self.assertEqual(atom._charge, -1)


    def test_can_create_atom_with_bfactor(self):
        atom = Atom("C", bfactor=1.1)
        self.assertEqual(atom._bfactor, 1.1)


    def test_can_create_atom_with_anisotropy(self):
        atom = Atom("C", anisotropy=[2, 5, 3, 2, 5, 6])
        self.assertEqual(atom._anisotropy, [2, 5, 3, 2, 5, 6])



class AtomReprTests(TestCase):

    def test_atom_repr(self):
        atom = Atom("C", 2, 3, 5, id=15, name="CA")
        self.assertEqual(repr(atom), "<C(CA) Atom 15 at (2, 3, 5)>")


    def test_basic_atom_repr(self):
        atom = Atom("C")
        self.assertEqual(repr(atom), "<C Atom at (0, 0, 0)>")



class AtomStrTests(TestCase):

    def test_atom_str(self):
        atom = Atom("C", 2, 3, 5, id=15, name="CA")
        self.assertEqual(str(atom), "<Atom 15 (CA)>")


    def test_basic_atom_str(self):
        atom = Atom("C")
        self.assertEqual(str(atom), "<Atom (C)>")



class AtomElementTests(TestCase):

    def test_element_property(self):
        atom = Atom("C", 2, 3, 5)
        self.assertIs(atom._element, atom.element)


    def test_can_update_element(self):
        atom = Atom("C", 2, 3, 5)
        atom.element = "N"
        self.assertEqual(atom._element, "N")



class AtomXTests(TestCase):

    def test_x_property(self):
        atom = Atom("C", 2, 3, 5)
        self.assertIs(atom._x, atom.x)


    def test_can_update_x(self):
        atom = Atom("C", 2, 3, 5)
        atom.x = 6
        self.assertEqual(atom._x, 6)



class AtomYTests(TestCase):

    def test_y_property(self):
        atom = Atom("C", 2, 3, 5)
        self.assertIs(atom._y, atom.y)


    def test_can_update_y(self):
        atom = Atom("C", 2, 3, 5)
        atom.y = 6
        self.assertEqual(atom._y, 6)



class AtomZTests(TestCase):

    def test_x_property(self):
        atom = Atom("C", 2, 3, 5)
        self.assertIs(atom._z, atom.z)


    def test_can_update_x(self):
        atom = Atom("C", 2, 3, 5)
        atom.z = 6
        self.assertEqual(atom._z, 6)



class AtomIdTests(TestCase):

    def test_id_property(self):
        atom = Atom("C", id=100)
        self.assertIs(atom._id, atom.id)



class AtomNameTests(TestCase):

    def test_name_property(self):
        atom = Atom("C", name="CA")
        self.assertIs(atom._name, atom.name)


    def test_can_update_name(self):
        atom = Atom("C", name="CA")
        atom.name = "CB"
        self.assertEqual(atom._name, "CB")



class AtomChargeTests(TestCase):

    def test_charge_property(self):
        atom = Atom("C", charge=0.5)
        self.assertIs(atom._charge, atom.charge)


    def test_can_update_charge(self):
        atom = Atom("C")
        atom.charge = 6
        self.assertEqual(atom._charge, 6)



class AtomBfactorTests(TestCase):

    def test_bfactor_property(self):
        atom = Atom("C", bfactor=0.5)
        self.assertIs(atom._bfactor, atom.bfactor)


    def test_can_update_bfactor(self):
        atom = Atom("C")
        atom.bfactor = 6
        self.assertEqual(atom._bfactor, 6)



class AtomAnisotropyTests(TestCase):

    def test_anisotropy_property(self):
        atom = Atom("C", anisotropy=(1, 2, 3, 4, 5, 6))
        self.assertIs(atom._anisotropy, atom.anisotropy)



class AtomLocationTests(TestCase):

    def test_atom_location(self):
        atom = Atom("C", 20, 30, 50)
        self.assertEqual(atom.location, (20, 30, 50))



class AtomTrimmingTests(TestCase):

    def test_can_not_round_atom_location(self):
        atom = Atom("C", 0.1111, 0.4444, 0.8888)
        atom.trim(None)
        self.assertEqual((atom._x, atom._y, atom._z), (0.1111, 0.4444, 0.8888))


    def test_can_round_atom_location(self):
        atom = Atom("C", 0.1111, 0.4444, 0.8888)
        atom.trim(3)
        self.assertEqual((atom._x, atom._y, atom._z), (0.111, 0.444, 0.889))
        atom.trim(2)
        self.assertEqual((atom._x, atom._y, atom._z), (0.11, 0.44, 0.89))
        atom.trim(1)
        self.assertEqual((atom._x, atom._y, atom._z), (0.1, 0.4, 0.9))
        atom.trim(0)
        self.assertEqual((atom._x, atom._y, atom._z), (0, 0, 1))



class AtomTranslationTests(TestCase):

    @patch("atomium.models.atoms.Atom.trim")
    def test_can_translate_atoms(self, mock_trim):
        atom = Atom("C", 20, 30, 50)
        atom.translate(5, -4, 12)
        self.assertEqual(atom._x, 25)
        self.assertEqual(atom._y, 26)
        self.assertEqual(atom._z, 62)
        mock_trim.assert_called_with(12)


    @patch("atomium.models.atoms.Atom.trim")
    def test_can_translate_atoms_with_iterable(self, mock_trim):
        atom = Atom("C", 20, 30, 50)
        atom.translate((5, -4, 12))
        self.assertEqual(atom._x, 25)
        self.assertEqual(atom._y, 26)
        self.assertEqual(atom._z, 62)
        mock_trim.assert_called_with(12)


    @patch("atomium.models.atoms.Atom.trim")
    def test_can_translate_atoms_and_specify_trim(self, mock_trim):
        atom = Atom("C", 20, 30, 50)
        atom.translate(5, -4, 12, trim=18)
        mock_trim.assert_called_with(18)
        atom.translate((5, -4, 12), trim=None)
        mock_trim.assert_called_with(None)



class AtomMovingTests(TestCase):

    def test_can_move_atom(self):
        atom = Atom("C", 20, 30, 50)
        atom.move_to(1, 2, 3)
        self.assertEqual(atom._x, 1)
        self.assertEqual(atom._y, 2)
        self.assertEqual(atom._z, 3)



class AtomRotationTests(TestCase):

    @patch("atomium.models.atoms.Atom.trim")
    def test_can_rotate_atoms(self, mock_trim):
        atom = Atom("C", 20, 30, 50)
        atom.rotate([[1, 0, 0], [0, 0.7071, -0.7071], [0, 0.7071, 0.7071]], trim=12)
        self.assertEqual(atom._x, 20)
        self.assertAlmostEqual(atom._y, -14.142, delta=3)
        self.assertAlmostEqual(atom._z, 56.569, delta=3)
        mock_trim.assert_called_with(12)


    @patch("atomium.models.atoms.Atom.trim")
    def test_can_rotate_atoms_and_specify_trim(self, mock_trim):
        atom = Atom("C", 20, 30, 50)
        atom.rotate([[1, 0, 0], [0, 0.7, -0.7], [0, 0.7, 0.7]], trim=1)
        mock_trim.assert_called_with(1)



class AtomMassTests(TestCase):

    def test_known_element_mass(self):
        atom = Atom("C", 2, 3, 5)
        self.assertAlmostEqual(atom.mass, 12, delta=0.1)
        atom._element = "H"
        self.assertAlmostEqual(atom.mass, 1, delta=0.1)


    def test_atom_mass_case_insensitive(self):
        atom = Atom("he", 2, 3, 5)
        self.assertAlmostEqual(atom.mass, 4, delta=0.1)
        atom = Atom("He", 2, 3, 5)
        self.assertAlmostEqual(atom.mass, 4, delta=0.1)
        atom = Atom("hE", 2, 3, 5)
        self.assertAlmostEqual(atom.mass, 4, delta=0.1)
        atom = Atom("HE", 2, 3, 5)
        self.assertAlmostEqual(atom.mass, 4, delta=0.1)


    def test_unknown_atom_mass(self):
        atom = Atom("XX", 2, 3, 5)
        self.assertEqual(atom.mass, 0)



class AtomMetalTests(TestCase):

    def test_some_atoms_are_metals(self):
        self.assertFalse(Atom("C", 2, 3, 5).is_metal)
        self.assertFalse(Atom("c", 2, 3, 5).is_metal)
        self.assertTrue(Atom("FE", 2, 3, 5).is_metal)
        self.assertTrue(Atom("zn", 2, 3, 5).is_metal)



class AtomDistanceToTests(TestCase):

    def test_can_get_distance_between_atoms(self):
        atom1 = Atom("C", 4, 8, 3)
        atom2 = Mock(Atom)
        atom2.location = (2, 3, 5)
        self.assertAlmostEqual(atom1.distance_to(atom2), 5.744, delta=0.001)


    def test_atom_distance_can_be_zero(self):
        atom1 = Atom("C", 4, 8, 3)
        atom2 = Mock(Atom)
        atom2.location = (4, 8, 3)
        self.assertEqual(atom1.distance_to(atom2), 0)


    def test_can_get_distance_to_xyz_tuple(self):
        atom1 = Atom("C", 4, 8, 3)
        atom2 = (2, 3, 5)
        self.assertAlmostEqual(atom1.distance_to(atom2), 5.744, delta=0.001)



class AtomLigandTests(TestCase):

    def test_ligand_property(self):
        ligand = Mock()
        atom = Atom("C", 2, 3, 5)
        atom._ligand = ligand
        self.assertIs(atom.ligand, ligand)



class AtomResidueTests(TestCase):

    def test_residue_property(self):
        residue = Mock()
        atom = Atom("C", 2, 3, 5)
        atom._residue = residue
        self.assertIs(atom.residue, residue)



class AtomChainTests(TestCase):

    def test_chain_property(self):
        chain = Mock()
        atom = Atom("C", 2, 3, 5)
        atom._chain = chain
        self.assertIs(atom.chain, chain)



class AtomComplexTests(TestCase):

    def test_complex_property(self):
        complex = Mock()
        atom = Atom("C", 2, 3, 5)
        atom._complex = complex
        self.assertIs(atom.complex, complex)



class AtomModelTests(TestCase):

    def test_model_property(self):
        model = Mock()
        atom = Atom("C", 2, 3, 5)
        atom._model = model
        self.assertIs(atom.model, model)



class AtomBondedAtomTests(TestCase):

    def test_bonded_atoms_property(self):
        atom = Atom("C", 2, 3, 5)
        atom._bonded_atoms = [1, 2, 3]
        self.assertEqual(atom.bonded_atoms, {1, 2, 3})



class AtomBondingTests(TestCase):

    def test_can_bond_to_atom(self):
        atom = Atom("C")
        atom2 = Mock(_bonded_atoms=set())
        atom.bond_to(atom2)
        self.assertEqual(atom._bonded_atoms, {atom2})
        self.assertEqual(atom2._bonded_atoms, {atom})



class AtomUnbondingTests(TestCase):

    def test_can_unbond_from_atom(self):
        atom = Atom("C")
        atom2 = Mock(_bonded_atoms={atom})
        atom._bonded_atoms = {atom2}
        atom.unbond_from(atom2)
        self.assertEqual(atom._bonded_atoms, set())
        self.assertEqual(atom2._bonded_atoms, set())


    def test_can_handle_unbonding_from_unconnected_atom(self):
        atom = Atom("C")
        atom2 = Mock(_bonded_atoms=set())
        atom.unbond_from(atom2)
        self.assertEqual(atom._bonded_atoms, set())
        self.assertEqual(atom2._bonded_atoms, set())



class NearbyAtomTests(TestCase):

    def test_atom_with_no_model_has_no_nearby_atoms(self):
        atom = Atom("C", 4, 8, 3)
        self.assertEqual(atom.nearby_atoms(cutoff=1), set())
        self.assertEqual(atom.nearby_atoms(cutoff=100), set())


    @patch("atomium.models.atoms.Atom.location", new_callable=PropertyMock)
    def test_can_get_nearby_atoms(self, mock_loc):
        model = Mock()
        mock_loc.return_value = (1, 2, 3)
        atom = Atom("C", 4, 8, 3)
        atom._model = model
        model.atoms_in_sphere.return_value = [1, 2, atom, 4]
        atoms = atom.nearby_atoms(4, a=1, b=2)
        model.atoms_in_sphere.assert_called_with(1, 2, 3, 4, a=1, b=2)
        self.assertEqual(atoms, [1, 2, 4])



class NearbyResidueests(TestCase):

    def setUp(self):
        self.atoms = [Mock(), Mock(), Mock(), Mock(), Mock()]
        self.residues = [Mock(), Mock()]
        self.atoms[0].residue = self.residues[0]
        self.atoms[1].residue = self.residues[0]
        self.atoms[2].residue = self.residues[1]
        self.atoms[3].residue = None
        self.atoms[4].residue = None
        for atom in self.atoms: atom.ligand = None
        self.atoms[4].ligand = "LIG"


    @patch("atomium.models.atoms.Atom.nearby_atoms")
    def test_can_get_nearby_residues(self, mock_atoms):
        mock_atoms.return_value = self.atoms
        atom = Atom("C", 4, 8, 3)
        self.assertEqual(atom.nearby_residues(3, a=1), set(self.residues))
        mock_atoms.assert_called_with(3, a=1)


    @patch("atomium.models.atoms.Atom.nearby_atoms")
    def test_can_get_nearby_residues_and_ligands(self, mock_atoms):
        mock_atoms.return_value = self.atoms
        atom = Atom("C", 4, 8, 3)
        self.assertEqual(
         atom.nearby_residues(3, ligands=True, a=1), set(self.residues + ["LIG"])
        )
        mock_atoms.assert_called_with(3, a=1)



class AtomCopyingTests(TestCase):

    def test_can_make_atom_copy(self):
        atom = Atom("he", 2, 3, 5, 10, "H1", 0.5, 0.4, (1, 1, 1, 1, 1, 1))
        patcher = patch("atomium.models.atoms.Atom")
        mock_atom = patcher.start()
        try:
            copy = atom.copy()
            self.assertIs(copy, mock_atom.return_value)
            mock_atom.assert_called_with(element="he", x=2, y=3, z=5, id=10,
             name="H1", charge=0.5, bfactor=0.4, anisotropy=[1, 1, 1, 1, 1, 1])
        finally:
            patcher.stop()
