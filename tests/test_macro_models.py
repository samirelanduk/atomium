import copy
from unittest import TestCase
from molecupy import exceptions
from molecupy.molecules import Atom, Molecule, Model
from molecupy.macromolecules import Residue, ResiduicStructure, Chain, MacroModel, Complex, Site

class MacroModelTest(TestCase):

    def setUp(self):
        self.atom1 = Atom(1.0, 1.0, 1.0, "H", atom_id=1, atom_name="H1")
        self.atom2 = Atom(1.0, 1.0, 2.0, "C", atom_id=2, atom_name="CA")
        self.atom3 = Atom(1.0, 1.0, 3.0, "O", atom_id=3, atom_name="OX1")
        self.atom2.covalent_bond_to(self.atom1)
        self.atom2.covalent_bond_to(self.atom3)
        self.residue1 = Residue(1, "MON1", self.atom1, self.atom2, self.atom3)
        self.atom4 = Atom(1.0, 1.0, 4.0, "H", atom_id=4, atom_name="H1")
        self.atom5 = Atom(1.0, 1.0, 5.0, "C", atom_id=5, atom_name="CA")
        self.atom6 = Atom(1.0, 1.0, 6.0, "O", atom_id=6, atom_name="OX1")
        self.atom5.covalent_bond_to(self.atom4)
        self.atom5.covalent_bond_to(self.atom6)
        self.residue2 = Residue(2, "MON2", self.atom4, self.atom5, self.atom6)
        self.atom7 = Atom(1.0, 1.0, 7.0, "H", atom_id=7, atom_name="H1")
        self.atom8 = Atom(1.0, 1.0, 8.0, "C", atom_id=8, atom_name="CA")
        self.atom9 = Atom(1.0, 1.0, 9.0, "O", atom_id=9, atom_name="OX1")
        self.atom8.covalent_bond_to(self.atom7)
        self.atom8.covalent_bond_to(self.atom9)
        self.residue3 = Residue(3, "MON3", self.atom7, self.atom8, self.atom9)
        self.atom10 = Atom(1.0, 1.0, 10.0, "H", atom_id=10, atom_name="H1")
        self.atom11 = Atom(1.0, 1.0, 11.0, "C", atom_id=11, atom_name="CA")
        self.atom12 = Atom(1.0, 1.0, 12.0, "O", atom_id=12, atom_name="OX1")
        self.atom11.covalent_bond_to(self.atom10)
        self.atom11.covalent_bond_to(self.atom12)
        self.residue4 = Residue(4, "MON4", self.atom10, self.atom11, self.atom12)
        self.residue1.connect_to(self.residue2, self.atom3, self.atom4)
        self.residue3.connect_to(self.residue4, self.atom9, self.atom10)


    def check_valid_macro_model(self, macro_model):
        self.assertIsInstance(macro_model, MacroModel)
        self.assertIsInstance(macro_model, Model)
        self.assertIsInstance(macro_model.atoms, set)
        self.assertIsInstance(macro_model._molecules, set)
        self.assertIsInstance(macro_model._small_molecules, set)
        self.assertIsInstance(macro_model._chains, set)
        self.assertIsInstance(macro_model._complexes, set)
        self.assertIsInstance(macro_model._sites, set)
        self.assertRegex(str(macro_model), r"<MacroModel \((\d+) atoms\)>")



class MacroModelCreationTests(MacroModelTest):

    def test_can_make_macro_model(self):
        macro_model = MacroModel()
        self.check_valid_macro_model(macro_model)



class SmallMoleculeAdditionTests(MacroModelTest):

    def setUp(self):
        MacroModelTest.setUp(self)
        self.molecule1 = Molecule(self.atom1, self.atom2, self.atom3)
        self.molecule2 = Molecule(self.atom4, self.atom5, self.atom6)


    def test_can_add_small_molecules(self):
        macro_model = MacroModel()
        macro_model.add_small_molecule(self.molecule1)
        self.assertEqual(len(macro_model.get_small_molecules()), 1)
        self.assertEqual(len(macro_model.get_molecules()), 1)
        self.assertEqual(self.molecule1.model, macro_model)
        self.assertEqual(len(macro_model.atoms), 3)
        self.assertIn(self.molecule1, macro_model)
        macro_model.add_small_molecule(self.molecule2)
        self.assertEqual(len(macro_model.get_small_molecules()), 2)
        self.assertEqual(len(macro_model.get_molecules()), 2)
        self.assertEqual(self.molecule2.model, macro_model)
        self.assertEqual(len(macro_model.atoms), 6)
        self.assertIn(self.molecule2, macro_model)


    def test_molecules_must_be_molecules(self):
        macro_model = MacroModel()
        with self.assertRaises(TypeError):
            macro_model.add_chain("mol")


    def test_molecules_must_be_unique(self):
        macro_model = MacroModel()
        self.molecule1.molecule_id = 23
        self.molecule2.molecule_id = 23
        macro_model.add_small_molecule(self.molecule1)
        with self.assertRaises(exceptions.DuplicateMoleculeError):
            macro_model.add_small_molecule(self.molecule2)


    def test_can_get_small_molecules_by_id(self):
        macro_model = MacroModel()
        self.molecule1.molecule_id = 1
        self.molecule2.molecule_id = 2
        macro_model.add_small_molecule(self.molecule1)
        macro_model.add_small_molecule(self.molecule2)
        self.assertEqual(
         macro_model.get_small_molecule_by_id(1),
         self.molecule1
        )
        self.assertEqual(
         macro_model.get_small_molecule_by_id(2),
         self.molecule2
        )
        self.assertEqual(
         macro_model.get_small_molecule_by_id(3),
         None
        )


    def test_can_only_search_by_numeric_id(self):
        macro_model = MacroModel()
        with self.assertRaises(TypeError):
            macro_model.get_small_molecule_by_id(None)


    def test_can_get_small_molecules_by_name(self):
        macro_model = MacroModel()
        self.molecule1.molecule_name = "STK"
        self.molecule2.molecule_name = "LAN"
        macro_model.add_small_molecule(self.molecule1)
        macro_model.add_small_molecule(self.molecule2)
        self.assertEqual(
         macro_model.get_small_molecule_by_name("LAN"),
         self.molecule2
        )
        self.assertEqual(
         macro_model.get_small_molecule_by_name("GJY"),
         None
        )


    def test_can_get_multiple_small_molecules_by_name(self):
        macro_model = MacroModel()
        self.molecule1.molecule_name = "STK"
        self.molecule2.molecule_name = "STK"
        macro_model.add_small_molecule(self.molecule1)
        macro_model.add_small_molecule(self.molecule2)
        self.assertEqual(
         macro_model.get_small_molecules_by_name("STK"),
         set([self.molecule1, self.molecule2])
        )
        self.assertEqual(macro_model.get_small_molecules_by_name("GJY"), set())


    def test_can_only_search_by_string_name(self):
        macro_model = MacroModel()
        with self.assertRaises(TypeError):
            macro_model.get_small_molecule_by_name(None)
        with self.assertRaises(TypeError):
            macro_model.get_small_molecules_by_name(None)



class ChainAdditionTests(MacroModelTest):

    def setUp(self):
        MacroModelTest.setUp(self)
        self.chain1 = Chain(self.residue1, self.residue2)
        self.chain2 = Chain(self.residue3, self.residue4)


    def test_can_add_chains(self):
        macro_model = MacroModel()
        macro_model.add_chain(self.chain1)
        self.assertEqual(len(macro_model.get_chains()), 1)
        self.assertEqual(len(macro_model.get_molecules()), 1)
        self.assertEqual(self.chain1.model, macro_model)
        self.assertEqual(len(macro_model.atoms), 6)
        self.assertIn(self.chain1, macro_model)
        macro_model.add_chain(self.chain2)
        self.assertEqual(len(macro_model.get_chains()), 2)
        self.assertEqual(len(macro_model.get_molecules()), 2)
        self.assertEqual(self.chain2.model, macro_model)
        self.assertEqual(len(macro_model.atoms), 12)
        self.assertIn(self.chain2, macro_model)


    def test_chains_must_be_chains(self):
        macro_model = MacroModel()
        with self.assertRaises(TypeError):
            macro_model.add_chain("chain")


    def test_chains_must_be_unique(self):
        macro_model = MacroModel()
        self.chain1.chain_id = "A"
        self.chain2.chain_id = "A"
        macro_model.add_chain(self.chain1)
        with self.assertRaises(exceptions.DuplicateChainError):
            macro_model.add_chain(self.chain2)


    def test_can_get_chains_by_id(self):
        macro_model = MacroModel()
        self.chain1.chain_id = "A"
        self.chain2.chain_id = "B"
        macro_model.add_chain(self.chain1)
        macro_model.add_chain(self.chain2)
        self.assertEqual(
         macro_model.get_chain_by_id("A"),
         self.chain1
        )
        self.assertEqual(
         macro_model.get_chain_by_id("B"),
         self.chain2
        )
        self.assertEqual(
         macro_model.get_chain_by_id("C"),
         None
        )


    def test_can_only_search_by_str_id(self):
        macro_model = MacroModel()
        with self.assertRaises(TypeError):
            macro_model.get_chain_by_id(None)



class ComplexAdditionTests(MacroModelTest):

    def setUp(self):
        MacroModelTest.setUp(self)
        self.chain1 = Chain(self.residue1, self.residue2)
        self.chain2 = Chain(self.residue3, self.residue4)
        self.complex = Complex(self.chain1, self.chain2)


    def test_can_add_complexes(self):
        macro_model = MacroModel()
        macro_model.add_complex(self.complex)
        self.assertEqual(len(macro_model.get_complexes()), 1)
        self.assertEqual(len(macro_model.get_chains()), 2)
        self.assertEqual(len(macro_model.get_molecules()), 2)
        self.assertEqual(len(macro_model.get_small_molecules()), 0)
        self.assertEqual(self.complex.model, macro_model)
        self.assertEqual(self.chain1.model, macro_model)
        self.assertEqual(self.chain2.model, macro_model)
        self.assertEqual(len(macro_model.atoms), 12)
        self.assertIn(self.complex, macro_model)


    def test_complexes_must_be_complexes(self):
        macro_model = MacroModel()
        with self.assertRaises(TypeError):
            macro_model.add_complex("complex")


    def test_complexes_must_be_unique(self):
        self.complex.complex_id = 1
        complex2 = copy.deepcopy(self.complex)
        macro_model = MacroModel()
        macro_model.add_complex(self.complex)
        with self.assertRaises(exceptions.DuplicateComplexError):
            macro_model.add_complex(complex2)



class SiteAdditionTests(MacroModelTest):

    def setUp(self):
        MacroModelTest.setUp(self)
        self.site1 = Site(self.residue1, self.residue2)
        self.site2 = Site(self.residue3, self.residue4)


    def test_can_add_sites(self):
        macro_model = MacroModel()
        macro_model.add_site(self.site1)
        self.assertEqual(len(macro_model.get_sites()), 1)
        self.assertEqual(len(macro_model.get_molecules()), 0)
        self.assertEqual(self.site1.model, macro_model)
        self.assertEqual(len(macro_model.atoms), 6)
        self.assertIn(self.site1, macro_model)
        macro_model.add_site(self.site2)
        self.assertEqual(len(macro_model.get_sites()), 2)
        self.assertEqual(len(macro_model.get_molecules()), 0)
        self.assertEqual(self.site2.model, macro_model)
        self.assertEqual(len(macro_model.atoms), 12)
        self.assertIn(self.site2, macro_model)


    def test_sites_must_be_sites(self):
        macro_model = MacroModel()
        with self.assertRaises(TypeError):
            macro_model.add_chain("site")


    def test_sites_must_be_unique(self):
        macro_model = MacroModel()
        self.site1.site_id = 1
        self.site2.site_id = 1
        macro_model.add_site(self.site1)
        with self.assertRaises(exceptions.DuplicateSiteError):
            macro_model.add_site(self.site2)


    def test_can_get_sites_by_id(self):
        macro_model = MacroModel()
        self.site1.site_id = 1
        self.site2.site_id = 2
        macro_model.add_site(self.site1)
        macro_model.add_site(self.site2)
        self.assertEqual(
         macro_model.get_site_by_id(1),
         self.site1
        )
        self.assertEqual(
         macro_model.get_site_by_id(2),
         self.site2
        )
        self.assertEqual(
         macro_model.get_site_by_id(3),
         None
        )


    def test_can_only_search_by_int_id(self):
        macro_model = MacroModel()
        with self.assertRaises(TypeError):
            macro_model.get_site_by_id(None)


    def test_can_sites_by_name(self):
        macro_model = MacroModel()
        self.site1.site_name = "AB1"
        self.site2.site_name = "AB2"
        macro_model.add_site(self.site1)
        macro_model.add_site(self.site2)
        self.assertEqual(
         macro_model.get_site_by_name("AB1"),
         self.site1
        )
        self.assertEqual(
         macro_model.get_site_by_name("AB3"),
         None
        )


    def test_can_get_multiple_sites_by_name(self):
        macro_model = MacroModel()
        self.site1.site_name = "AB1"
        self.site2.site_name = "AB1"
        macro_model.add_site(self.site1)
        macro_model.add_site(self.site2)
        self.assertEqual(
         macro_model.get_sites_by_name("AB1"),
         set([self.site1, self.site2])
        )
        self.assertEqual(macro_model.get_sites_by_name("AB3"), set())


    def test_can_only_search_by_string_name(self):
        macro_model = MacroModel()
        with self.assertRaises(TypeError):
            macro_model.get_site_by_name(None)
        with self.assertRaises(TypeError):
            macro_model.get_sites_by_name(None)



class BindSiteDetectionTests(MacroModelTest):

    def test_can_get_bind_site(self):
        self.residue2.connect_to(self.residue3, self.atom6, self.atom7)
        chain = Chain(self.residue1, self.residue2, self.residue3, self.residue4)
        atom1 = Atom(1.0, 3.5, 2.5, "C", atom_id=13, atom_name="C")
        atom2 = Atom(1.0, 3.5, 4.0, "C", atom_id=14, atom_name="C")
        atom3 = Atom(1.0, 3.5, 5.3, "C", atom_id=15, atom_name="C")
        atom4 = Atom(1.0, 4.5, 6.5, "C", atom_id=16, atom_name="C")
        atom5 = Atom(1.0, 4.5, 8.0, "C", atom_id=17, atom_name="C")
        atom6 = Atom(1.0, 4.5, 9.0, "C", atom_id=18, atom_name="C")
        atom7 = Atom(1.0, 4.5, 10.2, "C", atom_id=19, atom_name="C")
        atom8 = Atom(1.0, 4.5, 11.5, "C", atom_id=20, atom_name="C")
        atom9 = Atom(1.0, 4.0, 13.0, "C", atom_id=21, atom_name="C")
        atom10 = Atom(1.0, 2.5, 12.8, "C", atom_id=22, atom_name="C")
        atom1.covalent_bond_to(atom2)
        atom2.covalent_bond_to(atom3)
        atom3.covalent_bond_to(atom4)
        atom4.covalent_bond_to(atom5)
        atom5.covalent_bond_to(atom6)
        atom6.covalent_bond_to(atom7)
        atom7.covalent_bond_to(atom8)
        atom8.covalent_bond_to(atom9)
        atom9.covalent_bond_to(atom10)
        molecule = Molecule(atom1, atom2, atom3, atom4, atom5, atom6, atom7, atom8, atom9, atom10)
        macro_model = MacroModel()
        macro_model.add_chain(chain)
        macro_model.add_small_molecule(molecule)
        self.assertEqual(
         macro_model.get_adjacent_residues(molecule),
         set([self.residue1, self.residue4, self.residue2])
        )
