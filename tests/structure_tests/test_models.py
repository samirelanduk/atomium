from unittest import TestCase
import unittest.mock
from molecupy.structures import Model, AtomicStructure, SmallMolecule, Chain
from molecupy.structures import BindSite, Atom, Complex
from molecupy.exceptions import DuplicateSmallMoleculesError, DuplicateChainsError
from molecupy.exceptions import DuplicateBindSitesError, DuplicateComplexesError
from molecupy.pdb.pdbdatafile import PdbDataFile

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
        self.chain1 = unittest.mock.Mock(spec=Chain)
        self.chain1._model = None
        self.chain1.chain_id.return_value = "A"
        self.chain1.atoms.return_value = set()
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
