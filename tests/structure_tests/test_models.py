from unittest import TestCase
import unittest.mock
from molecupy.structures import Model, AtomicStructure, SmallMolecule, Chain
from molecupy.structures import BindSite, Atom, GhostAtom, Complex, Residue
from molecupy.exceptions import DuplicateSmallMoleculesError, DuplicateChainsError
from molecupy.exceptions import DuplicateBindSitesError, DuplicateComplexesError
from molecupy.pdb.pdbdatafile import PdbDataFile

#TODO: Make sure user can't add SmallMolecule, Chain etc. if there are Atom ID clashes

class ModelTest(TestCase):

    def setUp(self):
        self.small_molecule1 = unittest.mock.Mock(spec=SmallMolecule)
        self.small_molecule1._model = None
        self.small_molecule1.molecule_id.return_value = "A100"
        self.small_molecule1.molecule_name.return_value = "MOL"
        self.small_molecule1.atoms.return_value = set()
        self.small_molecule2 = unittest.mock.Mock(spec=SmallMolecule)
        self.small_molecule2._model = None
        self.small_molecule2.molecule_id.return_value = "A101"
        self.small_molecule2.molecule_name.return_value = "HET"
        self.small_molecule2.atoms.return_value = set()
        self.atoms = [GhostAtom("A", i + 1, "ATM") for i in range(10)]
        self.residues = [
         Residue("A%i" % (i + 1), "RES", *self.atoms[i:i+1]) for i in range(10)
        ]
        self.chain1 = unittest.mock.Mock(spec=Chain)
        self.chain1._model = None
        self.chain1.chain_id.return_value = "A"
        self.chain1.residues.return_value = self.residues
        self.chain1.atoms.return_value = set(self.atoms)
        self.chain2 = unittest.mock.Mock(spec=Chain)
        self.chain2._model = None
        self.chain2.chain_id.return_value = "B"
        self.chain2.atoms.return_value = set()
        self.site1 = unittest.mock.Mock(spec=BindSite)
        self.site1._model = None
        self.site1.site_id.return_value = "AA"
        self.site1.atoms.return_value = set()
        self.site2 = unittest.mock.Mock(spec=BindSite)
        self.site2._model = None
        self.site2.site_id.return_value = "BB"
        self.site2.atoms.return_value = set()
        self.complex1 = unittest.mock.Mock(Complex)
        self.complex1.complex_id.return_value = "1"
        self.complex1.complex_name.return_value = "COM1"
        self.complex1._model = None
        self.complex2 = unittest.mock.Mock(Complex)
        self.complex2.complex_id.return_value = "2"
        self.complex2.complex_name.return_value = "COM2"
        self.complex2._model = None
        self.small_molecule1.atoms.return_value = set(self.atoms)
        self.complex1.chains.return_value = set([self.chain1, self.chain2])
        self.complex2.chains.return_value = set([self.chain1, self.chain2])



class ModelCreationTest(ModelTest):

    def test_can_create_model(self):
        model = Model()
        self.assertIsInstance(model, AtomicStructure)
        self.assertEqual(model._atoms, set())


    def test_model_repr(self):
        model = Model()
        self.assertEqual(str(model), "<Model (0 atoms)>")



class ModelSmallMoleculeTests(ModelTest):

    def test_can_add_small_molecules(self):
        model = Model()
        self.assertEqual(model.small_molecules(), set())
        model.add_small_molecule(self.small_molecule1)
        self.assertEqual(model.small_molecules(), set([self.small_molecule1]))
        model.add_small_molecule(self.small_molecule2)
        self.assertEqual(
         model.small_molecules(),
         set([self.small_molecule1, self.small_molecule2])
        )


    def test_must_use_method_to_add_small_molecule(self):
        model = Model()
        self.assertEqual(model.small_molecules(), set())
        model.small_molecules().add(self.small_molecule1)
        self.assertEqual(model.small_molecules(), set())


    def test_cannot_have_duplicate_small_molecules(self):
        model = Model()
        model.add_small_molecule(self.small_molecule1)
        self.small_molecule2.molecule_id.return_value = "A100"
        with self.assertRaises(DuplicateSmallMoleculesError):
            model.add_small_molecule(self.small_molecule2)


    def test_can_remove_small_molecules(self):
        model = Model()
        model.add_small_molecule(self.small_molecule1)
        self.assertEqual(model.small_molecules(), set([self.small_molecule1]))
        model.remove_small_molecule(self.small_molecule1)
        self.assertEqual(model.small_molecules(), set())


    def test_small_molecule_knows_about_model(self):
        model = Model()
        self.assertIs(self.small_molecule1._model, None)
        model.add_small_molecule(self.small_molecule1)
        self.assertIs(self.small_molecule1._model, model)
        model.remove_small_molecule(self.small_molecule1)
        self.assertIs(self.small_molecule1._model, None)


    def test_can_only_add_small_molecules(self):
        model = Model()
        with self.assertRaises(TypeError):
            model.add_small_molecule("molecule")


    def test_can_get_small_molecule_by_id(self):
        model = Model()
        model.add_small_molecule(self.small_molecule1)
        model.add_small_molecule(self.small_molecule2)
        self.assertIs(model.get_small_molecule_by_id("A100"), self.small_molecule1)
        self.assertIs(model.get_small_molecule_by_id("A101"), self.small_molecule2)
        self.assertIs(model.get_small_molecule_by_id("A102"), None)


    def test_can_only_get_small_molecule_with_str_id(self):
        model = Model()
        with self.assertRaises(TypeError):
            model.get_small_molecule_by_id(100)


    def test_can_get_small_molecule_by_name(self):
        model = Model()
        model.add_small_molecule(self.small_molecule1)
        model.add_small_molecule(self.small_molecule2)
        self.assertIs(model.get_small_molecule_by_name("MOL"), self.small_molecule1)
        self.assertIs(model.get_small_molecule_by_name("HET"), self.small_molecule2)
        self.assertIs(model.get_small_molecule_by_name("ABC"), None)


    def test_can_only_get_small_molecule_with_str_name(self):
        model = Model()
        with self.assertRaises(TypeError):
            model.get_small_molecule_by_name(100)


    def test_can_get_small_molecules_by_name(self):
        model = Model()
        model.add_small_molecule(self.small_molecule1)
        model.add_small_molecule(self.small_molecule2)
        self.assertEqual(
         model.get_small_molecules_by_name("MOL"),
         set([self.small_molecule1])
        )
        self.small_molecule2.molecule_name.return_value = "MOL"
        self.assertEqual(
         model.get_small_molecules_by_name("MOL"),
         set([self.small_molecule1, self.small_molecule2])
        )
        self.assertEqual(model.get_small_molecules_by_name("ABC"), set())


    def test_can_only_get_small_molecules_with_str_name(self):
        model = Model()
        with self.assertRaises(TypeError):
            model.get_small_molecules_by_name(100)


    def test_can_duplicate_small_molecules(self):
        model = Model()
        model.add_small_molecule(self.small_molecule1)
        self.assertEqual(model.small_molecules(), set([self.small_molecule1]))
        model.duplicate_small_molecule(self.small_molecule1)
        self.assertEqual(len(model.small_molecules()), 2)


    def test_can_only_duplicate_present_small_molecule(self):
        model = Model()
        model.add_chain(self.chain1)
        with self.assertRaises(TypeError):
            model.duplicate_small_molecule(self.chain1)
        with self.assertRaises(ValueError):
            model.duplicate_small_molecule(self.small_molecule1)


    def test_duplicate_small_molecules_have_unique_id(self):
        model = Model()
        model.add_small_molecule(self.small_molecule1)
        model.add_small_molecule(self.small_molecule2)
        new_molecule = model.duplicate_small_molecule(self.small_molecule1)
        self.assertEqual(new_molecule.molecule_id(), "A102")
        new_molecule = model.duplicate_small_molecule(new_molecule)
        self.assertEqual(new_molecule.molecule_id(), "A103")


    def test_duplicate_small_molecules_can_be_given_id(self):
        model = Model()
        model.add_small_molecule(self.small_molecule1)
        new_molecule = model.duplicate_small_molecule(
         self.small_molecule1, molecule_id="A123"
        )
        self.assertEqual(new_molecule.molecule_id(), "A123")


    def test_dupicate_small_molecule_id_must_be_str(self):
        model = Model()
        model.add_small_molecule(self.small_molecule1)
        with self.assertRaises(TypeError):
            new_molecule = model.duplicate_small_molecule(
             self.small_molecule1, molecule_id=123
            )


    def test_duplicate_small_molecule_ids_must_be_valid(self):
        model = Model()
        model.add_small_molecule(self.small_molecule1)
        model.add_small_molecule(self.small_molecule2)
        with self.assertRaises(ValueError):
            new_molecule = model.duplicate_small_molecule(
             self.small_molecule1, molecule_id="A101"
            )
        with self.assertRaises(ValueError):
            new_molecule = model.duplicate_small_molecule(
             self.small_molecule1, molecule_id="105"
            )
        with self.assertRaises(ValueError):
            new_molecule = model.duplicate_small_molecule(
             self.small_molecule1, molecule_id="A_105"
            )


    def test_duplicate_small_molecules_have_same_name(self):
        model = Model()
        model.add_small_molecule(self.small_molecule1)
        new_molecule = model.duplicate_small_molecule(self.small_molecule1)
        self.assertEqual(new_molecule.molecule_name(), "MOL")


    def test_duplicate_small_molecules_have_distinct_atoms(self):
        model = Model()
        model.add_small_molecule(self.small_molecule1)
        new_molecule = model.duplicate_small_molecule(self.small_molecule1)
        self.assertNotEqual(new_molecule.atoms(atom_type="all"), self.small_molecule1.atoms())
        for atom in new_molecule.atoms(atom_type="all"):
            self.assertNotIn(atom, self.small_molecule1.atoms())


    def test_duplicte_small_molecules_have_atoms_with_ok_ids(self):
        model = Model()
        model.add_small_molecule(self.small_molecule1)
        existing_model_ids = [atom.atom_id() for atom in model.atoms(atom_type="all")]
        new_molecule = model.duplicate_small_molecule(self.small_molecule1)
        for atom in new_molecule.atoms(atom_type="all"):
            self.assertNotIn(atom.atom_id(), existing_model_ids)


    def test_duplicate_small_molecule_has_ghost_and_localised_atoms(self):
        real_atoms = [Atom(1.0, 1.0, 1.0, "A", i + 1, "ATM") for i in range(10, 20)]
        self.small_molecule1.atoms.return_value = self.atoms + real_atoms
        model = Model()
        model.add_small_molecule(self.small_molecule1)
        self.assertEqual(len(model.atoms()), 10)
        self.assertEqual(len(model.atoms(atom_type="all")), 20)
        new_molecule = model.duplicate_small_molecule(self.small_molecule1)
        self.assertEqual(len(new_molecule.atoms()), 10)
        self.assertEqual(len(new_molecule.atoms(atom_type="all")), 20)



class ModelChainTests(ModelTest):

    def test_can_add_chains(self):
        model = Model()
        self.assertEqual(model.chains(), set())
        model.add_chain(self.chain1)
        self.assertEqual(model.chains(), set([self.chain1]))
        model.add_chain(self.chain2)
        self.assertEqual(
         model.chains(),
         set([self.chain1, self.chain2])
        )


    def test_must_use_method_to_add_chain(self):
        model = Model()
        self.assertEqual(model.chains(), set())
        model.chains().add(self.chain1)
        self.assertEqual(model.chains(), set())


    def test_cannot_have_duplicate_chains(self):
        model = Model()
        model.add_chain(self.chain1)
        self.chain2.chain_id.return_value = "A"
        with self.assertRaises(DuplicateChainsError):
            model.add_chain(self.chain2)


    def test_can_remove_chains(self):
        model = Model()
        model.add_chain(self.chain1)
        self.assertEqual(model.chains(), set([self.chain1]))
        model.remove_chain(self.chain1)
        self.assertEqual(model.chains(), set())


    def test_chain_knows_about_model(self):
        model = Model()
        self.assertIs(self.chain1._model, None)
        model.add_chain(self.chain1)
        self.assertIs(self.chain1._model, model)
        model.remove_chain(self.chain1)
        self.assertIs(self.chain1._model, None)


    def test_can_only_add_chains(self):
        model = Model()
        with self.assertRaises(TypeError):
            model.add_chain("chain")


    def test_can_get_chain_by_id(self):
        model = Model()
        model.add_chain(self.chain1)
        model.add_chain(self.chain2)
        self.assertIs(model.get_chain_by_id("A"), self.chain1)
        self.assertIs(model.get_chain_by_id("B"), self.chain2)
        self.assertIs(model.get_chain_by_id("C"), None)


    def test_can_only_get_chain_with_str_id(self):
        model = Model()
        with self.assertRaises(TypeError):
            model.get_chain_by_id(100)


    def test_can_duplicate_chains(self):
        model = Model()
        model.add_chain(self.chain1)
        self.assertEqual(model.chains(), set([self.chain1]))
        model.duplicate_chain(self.chain1)
        self.assertEqual(len(model.chains()), 2)


    def test_can_only_duplicate_chains(self):
        model = Model()
        model.add_small_molecule(self.small_molecule1)
        with self.assertRaises(TypeError):
            model.duplicate_chain(self.small_molecule1)
        with self.assertRaises(ValueError):
            model.duplicate_chain(self.chain1)


    def test_duplicated_chains_have_unique_ids(self):
        model = Model()
        model.add_chain(self.chain1)
        model.add_chain(self.chain2)
        new_chain = model.duplicate_chain(self.chain1)
        self.assertEqual(new_chain.chain_id(), "C")
        new_chain = model.duplicate_chain(new_chain)
        self.assertEqual(new_chain.chain_id(), "D")


    def test_duplicate_chains_can_be_given_id(self):
        model = Model()
        model.add_chain(self.chain1)
        new_chain = model.duplicate_chain(self.chain1, chain_id="Q")
        self.assertEqual(new_chain.chain_id(), "Q")


    def test_duplicate_chain_id_must_be_str(self):
        model = Model()
        model.add_chain(self.chain1)
        with self.assertRaises(TypeError):
            new_chain = model.duplicate_chain(self.chain1, chain_id=12)


    def test_duplicate_chain_must_be_valid(self):
        model = Model()
        model.add_chain(self.chain1)
        model.add_chain(self.chain2)
        with self.assertRaises(ValueError):
            new_chain = model.duplicate_chain(self.chain1, chain_id="B")
        with self.assertRaises(ValueError):
            new_chain = model.duplicate_chain(self.chain1, chain_id="CD")
        with self.assertRaises(ValueError):
            new_chain = model.duplicate_chain(self.chain1, chain_id="1")


    def test_duplicate_chains_have_distinct_residues(self):
        model = Model()
        model.add_chain(self.chain1)
        new_chain = model.duplicate_chain(self.chain1)
        self.assertNotEqual(self.chain1.residues(), new_chain.residues())
        for residue in new_chain.residues():
            self.assertNotIn(residue, self.chain1.residues())


    def test_duplicate_chains_have_missing_and_present_residues(self):
        model = Model()
        model.add_chain(self.chain1)
        for index, residue in enumerate(self.residues[::2]):
            residue.add_atom(Atom(1.0, 1.0, 1.0, "Z", index + 100, "Z"))
        new_chain = model.duplicate_chain(self.chain1)
        self.assertEqual(len(new_chain.residues()), 10)
        self.assertEqual(len(new_chain.residues(include_missing=False)), 5)


    def test_duplicate_chains_have_distinct_atoms(self):
        model = Model()
        model.add_chain(self.chain1)
        new_chain = model.duplicate_chain(self.chain1)
        self.assertNotEqual(new_chain.atoms(atom_type="all"), self.chain1.atoms())
        for atom in new_chain.atoms(atom_type="all"):
            self.assertNotIn(atom, self.chain1.atoms())


    def test_duplicate_chains_have_atoms_with_ok_ids(self):
        model = Model()
        model.add_chain(self.chain1)
        existing_model_ids = [atom.atom_id() for atom in model.atoms(atom_type="all")]
        new_chain = model.duplicate_chain(self.chain1)
        for atom in new_chain.atoms(atom_type="all"):
            self.assertNotIn(atom.atom_id(), existing_model_ids)


    def test_duplicate_chain_has_ghost_and_localised_atoms(self):
        real_atoms = [Atom(1.0, 1.0, 1.0, "A", i + 1, "ATM") for i in range(10, 20)]
        for index, residue in enumerate(self.residues):
            residue.add_atom(real_atoms[index])
        self.chain1.atoms.return_value = real_atoms + self.atoms
        model = Model()
        model.add_chain(self.chain1)
        self.assertEqual(len(model.atoms()), 10)
        self.assertEqual(len(model.atoms(atom_type="all")), 20)
        new_chain = model.duplicate_chain(self.chain1)
        self.assertEqual(len(new_chain.atoms(atom_type="all")), 20)
        self.assertEqual(len(new_chain.atoms()), 10)



class ModelBindSiteTests(ModelTest):

    def test_can_add_bind_sites(self):
        model = Model()
        self.assertEqual(model.bind_sites(), set())
        model.add_bind_site(self.site1)
        self.assertEqual(model.bind_sites(), set([self.site1]))
        model.add_bind_site(self.site2)
        self.assertEqual(
         model.bind_sites(),
         set([self.site1, self.site2])
        )


    def test_must_use_method_to_add_bind_site(self):
        model = Model()
        self.assertEqual(model.bind_sites(), set())
        model.bind_sites().add(self.site1)
        self.assertEqual(model.bind_sites(), set())


    def test_cannot_have_duplicate_sites(self):
        model = Model()
        model.add_bind_site(self.site1)
        self.site2.site_id.return_value = "AA"
        with self.assertRaises(DuplicateBindSitesError):
            model.add_bind_site(self.site2)


    def test_can_remove_sites(self):
        model = Model()
        model.add_bind_site(self.site1)
        self.assertEqual(model.bind_sites(), set([self.site1]))
        model.remove_bind_site(self.site1)
        self.assertEqual(model.bind_sites(), set())


    def test_site_knows_about_model(self):
        model = Model()
        self.assertIs(self.site1._model, None)
        model.add_bind_site(self.site1)
        self.assertIs(self.site1._model, model)
        model.remove_bind_site(self.site1)
        self.assertIs(self.site1._model, None)


    def test_can_only_add_bind_sites(self):
        model = Model()
        with self.assertRaises(TypeError):
            model.add_bind_site("site")


    def test_can_get_site_by_id(self):
        model = Model()
        model.add_bind_site(self.site1)
        model.add_bind_site(self.site2)
        self.assertIs(model.get_bind_site_by_id("AA"), self.site1)
        self.assertIs(model.get_bind_site_by_id("BB"), self.site2)
        self.assertIs(model.get_bind_site_by_id("CC"), None)


    def test_can_only_get_site_with_str_id(self):
        model = Model()
        with self.assertRaises(TypeError):
            model.get_bind_site_by_id(100)



class ModelComplexTests(ModelTest):

    def test_can_add_complexes(self):
        model = Model()
        self.assertEqual(model.complexes(), set())
        model.add_complex(self.complex1)
        self.assertEqual(model.complexes(), set([self.complex1]))
        model.add_complex(self.complex2)
        self.assertEqual(model.complexes(), set([self.complex1, self.complex2]))


    def test_must_use_method_to_get_complexes(self):
        model = Model()
        self.assertEqual(model.complexes(), set())
        model.complexes().add(self.complex1)
        self.assertEqual(model.complexes(), set())


    def test_cannot_have_duplicate_complexes(self):
        model = Model()
        model.add_complex(self.complex1)
        self.complex2.complex_id.return_value = "1"
        with self.assertRaises(DuplicateComplexesError):
            model.add_complex(self.complex2)


    def test_can_only_add_complexes(self):
        model = Model()
        with self.assertRaises(TypeError):
            model.add_complex("site")


    def test_can_remove_complexes(self):
        model = Model()
        model.add_complex(self.complex1)
        self.assertEqual(model.complexes(), set([self.complex1]))
        model.remove_complex(self.complex1)
        self.assertEqual(model.complexes(), set())


    def test_complex_knows_about_model(self):
        model = Model()
        self.assertIs(self.complex1._model, None)
        model.add_complex(self.complex1)
        self.assertIs(self.complex1._model, model)
        model.remove_complex(self.complex1)
        self.assertIs(self.complex1._model, None)


    def test_can_get_complex_by_id(self):
        model = Model()
        model.add_complex(self.complex1)
        model.add_complex(self.complex2)
        self.assertIs(model.get_complex_by_id("1"), self.complex1)
        self.assertIs(model.get_complex_by_id("2"), self.complex2)
        self.assertIs(model.get_complex_by_id("3"), None)


    def test_can_only_get_complex_with_str_id(self):
        model = Model()
        with self.assertRaises(TypeError):
            model.get_complex_by_id(100)


    def test_can_get_complex_by_name(self):
        model = Model()
        model.add_complex(self.complex1)
        model.add_complex(self.complex2)
        self.assertIs(model.get_complex_by_name("COM1"), self.complex1)
        self.assertIs(model.get_complex_by_name("COM2"), self.complex2)
        self.assertIs(model.get_complex_by_name("COM3"), None)


    def test_can_only_get_complex_str_name(self):
        model = Model()
        with self.assertRaises(TypeError):
            model.get_complex_by_name(100)


    def test_can_get_complexes_by_name(self):
        model = Model()
        model.add_complex(self.complex1)
        model.add_complex(self.complex2)
        self.assertEqual(
         model.get_complexes_by_name("COM1"),
         set([self.complex1])
        )
        self.complex2.complex_name.return_value = "COM1"
        self.assertEqual(
         model.get_complexes_by_name("COM1"),
         set([self.complex1, self.complex2])
        )
        self.assertEqual(model.get_complexes_by_name("ABC"), set())


    def test_can_only_get_complexes_with_str_name(self):
        model = Model()
        with self.assertRaises(TypeError):
            model.get_complexes_by_name(100)


    def test_can_duplicate_complexes(self):
        model = Model()
        model.add_complex(self.complex1)
        self.assertEqual(model.complexes(), set([self.complex1]))
        model.duplicate_complex(self.complex1)
        self.assertEqual(len(model.complexes()), 2)


    def test_can_only_duplicate_complexes(self):
        model = Model()
        model.add_small_molecule(self.small_molecule1)
        with self.assertRaises(TypeError):
            model.duplicate_complex(self.small_molecule1)
        with self.assertRaises(ValueError):
            model.duplicate_complex(self.complex1)


    def test_duplicated_complexes_have_unique_ids(self):
        model = Model()
        model.add_complex(self.complex1)
        model.add_complex(self.complex2)
        new_complex = model.duplicate_complex(self.complex1)
        self.assertEqual(new_complex.complex_id(), "1_1")
        new_complex = model.duplicate_complex(new_complex)
        self.assertEqual(new_complex.complex_id(), "1_2")
        new_complex = model.duplicate_complex(self.complex1)
        self.assertEqual(new_complex.complex_id(), "1_3")
        new_complex = model.duplicate_complex(self.complex2)
        self.assertEqual(new_complex.complex_id(), "2_1")
        new_complex = model.duplicate_complex(new_complex)
        self.assertEqual(new_complex.complex_id(), "2_2")


    def test_duplicate_complexes_can_be_given_id(self):
        model = Model()
        model.add_complex(self.complex1)
        new_complex = model.duplicate_complex(self.complex1, complex_id=".+.")
        self.assertEqual(new_complex.complex_id(), ".+.")


    def test_dupicate_complex_id_must_be_str(self):
        model = Model()
        model.add_complex(self.complex1)
        with self.assertRaises(TypeError):
            new_complex = model.duplicate_complex(self.complex1, complex_id=100)


    def test_duplicate_complex_id_must_be_valid(self):
        model = Model()
        model.add_complex(self.complex1)
        model.add_complex(self.complex2)
        with self.assertRaises(ValueError):
            new_complex = model.duplicate_complex(self.complex1, complex_id="2")


    def test_duplicate_complex_name_stays_the_same(self):
        model = Model()
        model.add_complex(self.complex1)
        new_complex = model.duplicate_complex(self.complex1)
        self.assertEqual(new_complex.complex_name(), "COM1")



class ModelAtomsTests(ModelTest):

    def test_model_has_atoms(self):
        model = Model()
        self.small_molecule1.atoms.return_value = set(["1", "2"])
        self.small_molecule2.atoms.return_value = set(["3", "4"])
        self.chain1.atoms.return_value = set(["5", "6"])
        self.chain2.atoms.return_value = set(["7", "8"])
        self.site1.atoms.return_value = set(["9", "10"])
        self.site2.atoms.return_value = set(["11", "12"])
        model.add_small_molecule(self.small_molecule1)
        model.add_small_molecule(self.small_molecule2)
        model.add_chain(self.chain1)
        model.add_chain(self.chain2)
        model.add_bind_site(self.site1)
        model.add_bind_site(self.site2)
        self.assertEqual(model._atoms, set([str(i) for i in range(1, 13)]))
        self.assertEqual(
         model.atoms(atom_type="all"),
         set([str(i) for i in range(1, 13)])
        )



class ConversionToPdbDataFileTests(ModelTest):

    def test_model_can_produce_pdbdatafile(self):
        model = Model()
        data_file = model.pdb_data_file()
        self.assertIsInstance(data_file, PdbDataFile)


    def test_can_add_atoms_to_pdb_data_file(self):
        model = Model()
        atom1 = unittest.mock.Mock(Atom)
        atom1.x.return_value, atom1.y.return_value, atom1.z.return_value = (
         1.2, 2.3, 3.4
        )
        atom1.element.return_value, atom1.atom_id.return_value, atom1.atom_name.return_value = (
         "G", 23, "GX"
        )
        atom2 = unittest.mock.Mock(Atom)
        atom2.x.return_value, atom2.y.return_value, atom2.z.return_value = (
         11.2, 11.3, 34.4
        )
        atom2.element.return_value, atom2.atom_id.return_value, atom2.atom_name.return_value = (
         "Y", 38, "YT"
        )
        model._atoms = set((atom1, atom2))
        data_file = model.pdb_data_file()
        self.assertEqual(
         data_file.atoms(),
         [
          {
           "atom_id": 23,
           "atom_name": "GX",
           "alt_loc": None,
           "residue_name": None,
           "chain_id": None,
           "residue_id": None,
           "insert_code": None,
           "x": 1.2,
           "y": 2.3,
           "z": 3.4,
           "occupancy": None,
           "temperature_factor": None,
           "element": "G",
           "charge": None,
           "model_id": 1
          }, {
           "atom_id": 38,
           "atom_name": "YT",
           "alt_loc": None,
           "residue_name": None,
           "chain_id": None,
           "residue_id": None,
           "insert_code": None,
           "x": 11.2,
           "y": 11.3,
           "z": 34.4,
           "occupancy": None,
           "temperature_factor": None,
           "element": "Y",
           "charge": None,
           "model_id": 1
          }
         ]
        )
