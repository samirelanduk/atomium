from unittest import TestCase
from molecupy import exceptions
from molecupy.structures import *

class ModelTest(TestCase):

    def check_valid_model(self, model):
        self.assertIsInstance(model, PdbModel)
        self.assertIsInstance(model, AtomicStructure)
        self.assertIsInstance(model.atoms, set)
        self.assertIsInstance(model.small_molecules, set)
        self.assertIsInstance(model.chains, set)
        self.assertIsInstance(model.sites, set)
        self.assertRegex(str(model), r"<Model \((\d+) atoms\)>")



class ModelCreationTests(ModelTest):

    def test_can_make_model(self):
        model = PdbModel()
        self.check_valid_model(model)



class SmallMoleculeAdditionTests(ModelTest):

    def setUp(self):
        self.atom1 = PdbAtom(1.0, 1.0, 1.0, "H", 1, "H1")
        self.atom2 = PdbAtom(1.0, 1.0, 2.0, "C", 2, "CA")
        self.molecule1 = PdbSmallMolecule("A1", "MOL1", self.atom1, self.atom2)
        self.atom3 = PdbAtom(1.0, 1.0, 3.0, "H", 3, "H1")
        self.atom4 = PdbAtom(1.0, 1.0, 4.0, "C", 4, "CA")
        self.molecule2 = PdbSmallMolecule("A2", "MOL2", self.atom1, self.atom2)
        self.model = PdbModel()


    def test_can_add_small_molecules(self):
        self.model.add_small_molecule(self.molecule1)
        self.assertIn(self.molecule1, self.model.small_molecules)


    def test_cannot_add_two_small_molecules_with_same_id(self):
        self.model.add_small_molecule(self.molecule1)
        self.molecule2.molecule_id = "A1"
        with self.assertRaises(exceptions.DuplicateSmallMoleculeError):
            self.model.add_small_molecule(self.molecule2)
        self.molecule2.molecule_id = "A2"
        self.model.add_small_molecule(self.molecule2)
        self.assertIn(self.molecule1, self.model.small_molecules)
        self.assertIn(self.molecule2, self.model.small_molecules)


    def test_add_small_molecule_atoms_to_model(self):
        self.model.add_small_molecule(self.molecule1)
        self.assertEqual(self.model.atoms, self.molecule1.atoms)
        self.model.add_small_molecule(self.molecule2)
        self.assertEqual(
         self.model.atoms,
         set(list(self.molecule1.atoms) + list(self.molecule2.atoms))
        )


    def test_small_molecules_must_be_small_molecules(self):
        with self.assertRaises(TypeError):
            self.model.add_small_molecule("mol")


    def test_can_get_small_molecules_by_id(self):
        self.model.add_small_molecule(self.molecule1)
        self.model.add_small_molecule(self.molecule2)
        self.assertEqual(
         self.model.get_small_molecule_by_id("A1"),
         self.molecule1
        )
        self.assertEqual(
         self.model.get_small_molecule_by_id("A2"),
         self.molecule2
        )
        self.assertEqual(
         self.model.get_small_molecule_by_id("A3"),
         None
        )


    def test_can_only_search_by_str_id(self):
        with self.assertRaises(TypeError):
            self.model.get_small_molecule_by_id(None)


    def test_can_get_small_molecules_by_name(self):
        self.model.add_small_molecule(self.molecule1)
        self.model.add_small_molecule(self.molecule2)
        self.assertEqual(
         self.model.get_small_molecule_by_name("MOL1"),
         self.molecule1
        )
        self.assertEqual(
         self.model.get_small_molecule_by_name("MOLX"),
         None
        )


    def test_can_get_multiple_small_molecules_by_name(self):
        self.model.add_small_molecule(self.molecule1)
        self.model.add_small_molecule(self.molecule2)
        self.molecule2.molecule_name = "MOL1"
        self.assertEqual(
         self.model.get_small_molecules_by_name("MOL1"),
         set([self.molecule1, self.molecule2])
        )
        self.assertEqual(self.model.get_small_molecules_by_name("XXX"), set())


    def test_can_only_search_by_string_name(self):
        with self.assertRaises(TypeError):
            self.model.get_small_molecule_by_name(None)
        with self.assertRaises(TypeError):
            self.model.get_small_molecules_by_name(None)



class ChainAdditionTests(ModelTest):

    def setUp(self):
        self.atom1 = PdbAtom(1.0, 1.0, 1.0, "H", 1, "H1")
        self.atom2 = PdbAtom(1.0, 1.0, 2.0, "C", 2, "CA")
        self.atom3 = PdbAtom(1.0, 1.0, 3.0, "O", 3, "OX1")
        self.residue1 = PdbResidue("A1", "ARG", self.atom1, self.atom2, self.atom3)
        self.atom4 = PdbAtom(1.0, 1.0, 4.0, "H", 4, "H1")
        self.atom5 = PdbAtom(1.0, 1.0, 5.0, "C", 5, "CA")
        self.atom6 = PdbAtom(1.0, 1.0, 6.0, "O", 6, "OX1")
        self.residue2 = PdbResidue("A2", "HST", self.atom4, self.atom5, self.atom6)
        self.atom7 = PdbAtom(1.0, 1.0, 7.0, "H", 7, "H1")
        self.atom8 = PdbAtom(1.0, 1.0, 8.0, "C", 8, "CA")
        self.atom9 = PdbAtom(1.0, 1.0, 9.0, "O", 9, "OX1")
        self.residue3 = PdbResidue("A3", "TRP", self.atom7, self.atom8, self.atom9)
        self.atom10 = PdbAtom(1.0, 1.0, 10.0, "H", 10, "H1")
        self.atom11 = PdbAtom(1.0, 1.0, 11.0, "C", 11, "CA")
        self.atom12 = PdbAtom(1.0, 1.0, 12.0, "O", 12, "OX1")
        self.residue4 = PdbResidue("A4", "TRP", self.atom10, self.atom11, self.atom12)
        self.chain1 = PdbChain("A", self.residue1, self.residue2)
        self.chain2 = PdbChain("B", self.residue3, self.residue4)
        self.model = PdbModel()


    def test_can_add_chains(self):
        self.model.add_chain(self.chain1)
        self.assertIn(self.chain1, self.model.chains)


    def test_cannot_add_two_chains_with_same_id(self):
        self.model.add_chain(self.chain1)
        self.chain2.chain_id = "A"
        with self.assertRaises(exceptions.DuplicateChainError):
            self.model.add_chain(self.chain2)
        self.chain2.chain_id = "B"
        self.model.add_chain(self.chain2)
        self.assertIn(self.chain1, self.model.chains)
        self.assertIn(self.chain2, self.model.chains)


    def test_add_chain_atoms_to_model(self):
        self.model.add_chain(self.chain1)
        self.assertEqual(self.model.atoms, self.chain1.atoms)
        self.model.add_chain(self.chain2)
        self.assertEqual(
         self.model.atoms,
         set(list(self.chain1.atoms) + list(self.chain2.atoms))
        )


    def test_chains_must_be_chains(self):
        with self.assertRaises(TypeError):
            self.model.add_chain(self.residue1)


    def test_can_get_chains_by_id(self):
        self.model.add_chain(self.chain1)
        self.model.add_chain(self.chain2)
        self.assertEqual(
         self.model.get_chain_by_id("A"),
         self.chain1
        )
        self.assertEqual(
         self.model.get_chain_by_id("B"),
         self.chain2
        )
        self.assertEqual(
         self.model.get_chain_by_id("C"),
         None
        )


    def test_can_only_search_by_str_id(self):
        with self.assertRaises(TypeError):
            self.model.get_chain_by_id(None)



class SiteAdditionTests(ModelTest):

    def setUp(self):
        self.atom1 = PdbAtom(1.0, 1.0, 1.0, "H", 1, "H1")
        self.atom2 = PdbAtom(1.0, 1.0, 2.0, "C", 2, "CA")
        self.atom3 = PdbAtom(1.0, 1.0, 3.0, "O", 3, "OX1")
        self.residue1 = PdbResidue("A1", "ARG", self.atom1, self.atom2, self.atom3)
        self.atom4 = PdbAtom(1.0, 1.0, 4.0, "H", 4, "H1")
        self.atom5 = PdbAtom(1.0, 1.0, 5.0, "C", 5, "CA")
        self.atom6 = PdbAtom(1.0, 1.0, 6.0, "O", 6, "OX1")
        self.residue2 = PdbResidue("A2", "HST", self.atom4, self.atom5, self.atom6)
        self.atom7 = PdbAtom(1.0, 1.0, 7.0, "H", 7, "H1")
        self.atom8 = PdbAtom(1.0, 1.0, 8.0, "C", 8, "CA")
        self.atom9 = PdbAtom(1.0, 1.0, 9.0, "O", 9, "OX1")
        self.residue3 = PdbResidue("A3", "TRP", self.atom7, self.atom8, self.atom9)
        self.atom10 = PdbAtom(1.0, 1.0, 10.0, "H", 10, "H1")
        self.atom11 = PdbAtom(1.0, 1.0, 11.0, "C", 11, "CA")
        self.atom12 = PdbAtom(1.0, 1.0, 12.0, "O", 12, "OX1")
        self.residue4 = PdbResidue("A4", "TRP", self.atom10, self.atom11, self.atom12)
        self.site1 = PdbSite("AB1", self.residue1, self.residue2)
        self.site2 = PdbSite("AB2", self.residue3, self.residue4)
        self.model = PdbModel()


    def test_can_add_sites(self):
        self.model.add_site(self.site1)
        self.assertIn(self.site1, self.model.sites)


    def test_cannot_add_two_sites_with_same_id(self):
        self.model.add_site(self.site1)
        self.site2.site_id = "AB1"
        with self.assertRaises(exceptions.DuplicateSiteError):
            self.model.add_site(self.site2)
        self.site2.site_id = "AB2"
        self.model.add_site(self.site2)
        self.assertIn(self.site1, self.model.sites)
        self.assertIn(self.site2, self.model.sites)


    def test_add_site_atoms_to_model(self):
        self.model.add_site(self.site1)
        self.assertEqual(self.model.atoms, self.site1.atoms)
        self.model.add_site(self.site2)
        self.assertEqual(
         self.model.atoms,
         set(list(self.site1.atoms) + list(self.site2.atoms))
        )


    def test_sites_must_be_sites(self):
        with self.assertRaises(TypeError):
            self.model.add_site(self.residue1)


    def test_can_get_sites_by_id(self):
        self.model.add_site(self.site1)
        self.model.add_site(self.site2)
        self.assertEqual(
         self.model.get_site_by_id("AB1"),
         self.site1
        )
        self.assertEqual(
         self.model.get_site_by_id("AB2"),
         self.site2
        )
        self.assertEqual(
         self.model.get_site_by_id("AB3"),
         None
        )


    def test_can_only_search_by_str_id(self):
        with self.assertRaises(TypeError):
            self.model.get_site_by_id(None)
