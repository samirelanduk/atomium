from unittest import TestCase
from unittest.mock import Mock, patch
from atomium.structures.molecules import Molecule, AtomicStructure
from atomium.structures.atoms import Atom

class MoleculeTest(TestCase):

    def setUp(self):
        self.atom1, self.atom2, self.atom3 = Mock(Atom), Mock(Atom), Mock(Atom)
        self.atoms = [self.atom1, self.atom2, self.atom3]
        def mock_init(obj, *args, **kwargs):
            obj._atoms = set(args[:1])
            obj._id_atoms = {i: arg for i, arg in enumerate(args[1:])}
        self.mock_init = mock_init



class MoleculeCreationTests(MoleculeTest):

    @patch("atomium.structures.molecules.AtomicStructure.__init__")
    def test_can_create_molecule(self, mock_init):
        mock_init.side_effect = self.mock_init
        mol = Molecule(self.atom1, self.atom2, self.atom3)
        self.assertIsInstance(mol, AtomicStructure)
        mock_init.assert_called_with(mol, self.atom1, self.atom2, self.atom3)
        self.assertEqual(mol._id, None)
        self.assertEqual(mol._name, None)


    def test_can_create_molecule_with_id(self):
        mol = Molecule(self.atom1, self.atom2, self.atom3, molecule_id="A100")
        self.assertEqual(mol._id, "A100")


    def test_molecule_id_must_be_str(self):
        with self.assertRaises(TypeError):
            Molecule(self.atom1, self.atom2, self.atom3, molecule_id=1000)


    def test_can_create_molecule_with_name(self):
        mol = Molecule(self.atom1, self.atom2, self.atom3, name="HIS")
        self.assertEqual(mol._name, "HIS")


    def test_molecule_name_must_be_str(self):
        with self.assertRaises(TypeError):
            Molecule(self.atom1, self.atom2, self.atom3, name=1000)


    def test_atoms_are_linked_to_molecule(self):
        self.atom1.atom_id.return_value = 100
        mol = Molecule(self.atom1, self.atom2, self.atom3)
        self.assertIs(self.atom1._molecule, mol)
        self.assertIs(self.atom2._molecule, mol)
        self.assertIs(self.atom3._molecule, mol)



class MoleculeReprTests(MoleculeTest):

    def test_molecule_repr_no_id_or_name(self):
        mol = Molecule(self.atom1, self.atom2, self.atom3)
        self.assertEqual(str(mol), "<Molecule (3 atoms)>")


    def test_molecule_repr_id_no_name(self):
        mol = Molecule(self.atom1, self.atom2, self.atom3, molecule_id="B10A")
        self.assertEqual(str(mol), "<Molecule B10A (3 atoms)>")


    def test_molecule_repr_name_no_id(self):
        mol = Molecule(self.atom1, self.atom2, self.atom3, name="GLY")
        self.assertEqual(str(mol), "<Molecule (GLY, 3 atoms)>")


    def test_molecule_repr_id_and_name(self):
        mol = Molecule(self.atom1, self.atom2, molecule_id="B10B", name="GLY")
        self.assertEqual(str(mol), "<Molecule B10B (GLY, 2 atoms)>")



class MoleculeIdTests(MoleculeTest):

    def test_molecule_id_property(self):
        mol = Molecule(self.atom1, self.atom2, self.atom3, molecule_id="B10C")
        self.assertIs(mol._id, mol.molecule_id())



class MoleculeNameTests(MoleculeTest):

    def test_molecule_name_property(self):
        mol = Molecule(self.atom1, self.atom2, self.atom3, name="VAL")
        self.assertIs(mol._name, mol.name())


    def test_can_update_molecule_name(self):
        mol = Molecule(self.atom1, self.atom2, self.atom3, name="VAL")
        mol.name("HIS")
        self.assertEqual(mol._name, "HIS")


    def test_molecule_name_must_be_str(self):
        mol = Molecule(self.atom1, self.atom2, self.atom3, name="VAL")
        with self.assertRaises(TypeError):
            mol.name(10)



class MoleculeAtomAdditionTests(MoleculeTest):

    def test_adding_atom_updates_atom(self):
        mol = Molecule(self.atom1, self.atom2)
        mol.add_atom(self.atom3)
        self.assertIs(self.atom3._molecule, mol)



class MoleculeAtomRemovalTests(MoleculeTest):

    def test_removing_atom_updates_atom(self):
        mol = Molecule(self.atom1, self.atom2, self.atom3)
        mol.remove_atom(self.atom3)
        self.assertIs(self.atom3._molecule, None)



class MoleculeModelTests(MoleculeTest):

    @patch("atomium.structures.molecules.Molecule.atoms")
    def test_can_get_model(self, mock_atoms):
        mock_atoms.return_value = set(self.atoms)
        model = Mock()
        self.atom1.model.return_value = model
        self.atom2.model.return_value = model
        self.atom3.model.return_value = model
        mol = Molecule(self.atom1, self.atom2, self.atom3)
        self.assertIs(mol.model(), model)


    @patch("atomium.structures.molecules.Molecule.atoms")
    def test_can_get_no_model(self, mock_atoms):
        mock_atoms.return_value = set(self.atoms)
        self.atom1.model.return_value = None
        self.atom2.model.return_value = None
        self.atom3.model.return_value = None
        mol = Molecule(self.atom1, self.atom2, self.atom3)
        self.assertIs(mol.model(), None)



class MoleculeSiteTests(MoleculeTest):

    @patch("atomium.structures.Molecule.atoms")
    @patch("atomium.structures.chains.Site")
    def test_can_get_site(self, mock_site, mock_atoms):
        residues = [Mock(), Mock(), Mock(), Mock(), Mock(), Mock(), Mock()]
        atoms = []
        for residue in residues:
            atoms += [Mock(), Mock(), Mock()]
            atoms[-3].residue.return_value = residue
            atoms[-2].residue.return_value = residue
            atoms[-1].residue.return_value = residue
        mock_atoms.return_value = set(self.atoms)
        self.atom1.nearby.return_value = set([self.atom2, atoms[1], atoms[3]])
        self.atom2.nearby.return_value = set([self.atom1, atoms[3], atoms[6]])
        self.atom3.nearby.return_value = set([self.atom3, atoms[6], atoms[-1]])
        self.atom1.residue.return_value = None
        self.atom2.residue.return_value = None
        self.atom3.residue.return_value = None
        site = Mock()
        mock_site.return_value = site
        mol = Molecule(self.atom1, self.atom2, self.atom3)
        returned_site = mol.site()
        self.assertIs(site, returned_site)
        mock_atoms.assert_called_with(exclude="H")
        self.atom1.nearby.assert_called_with(4, exclude="H")
        self.atom2.nearby.assert_called_with(4, exclude="H")
        self.atom3.nearby.assert_called_with(4, exclude="H")
        site_args, site_kwargs = mock_site.call_args_list[0]
        self.assertEqual(set(site_args), set([
         residues[0], residues[1], residues[2], residues[6]
        ]))
        self.assertEqual(site_kwargs, {"ligand": mol})
