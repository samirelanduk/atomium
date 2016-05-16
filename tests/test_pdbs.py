import datetime
from unittest import TestCase
from molecupy.pdbfile import PdbFile
from molecupy.pdbdatafile import PdbDataFile
from molecupy.pdb import Pdb
from molecupy.structures import PdbModel, PdbSmallMolecule
from molecupy.exceptions import *
from molecupy import pdb_from_string, get_pdb_remotely, get_pdb_from_file

class PdbTest(TestCase):

    def setUp(self):
        self.empty = Pdb(PdbDataFile(PdbFile("")))



class PdbPropertiesTests(PdbTest):

    def test_has_pdb_file(self):
        self.assertIsInstance(self.empty.data_file, PdbDataFile)


    def test_repr(self):
        self.assertRegex(
         str(self.empty),
         r"<Pdb \(([^\s]{4})\)>"
        )


    def test_other_properties(self):
        self.assertEqual(
         self.empty.classification,
         self.empty.data_file.classification
        )
        self.assertEqual(
         self.empty.deposition_date,
         self.empty.data_file.deposition_date
        )
        self.assertEqual(
         self.empty.pdb_code,
         self.empty.data_file.pdb_code
        )
        self.assertEqual(
         self.empty.is_obsolete,
         self.empty.data_file.is_obsolete
        )
        self.assertEqual(
         self.empty.obsolete_date,
         self.empty.data_file.obsolete_date
        )
        self.assertEqual(
         self.empty.replacement_code,
         self.empty.data_file.replacement_code
        )
        self.assertEqual(
         self.empty.title,
         self.empty.data_file.title
        )
        self.assertEqual(
         self.empty.split_codes,
         self.empty.data_file.split_codes
        )
        self.assertEqual(
         self.empty.caveat,
         self.empty.data_file.caveat
        )
        self.assertEqual(
         self.empty.keywords,
         self.empty.data_file.keywords
        )
        self.assertEqual(
         self.empty.experimental_techniques,
         self.empty.data_file.experimental_techniques
        )
        self.assertEqual(
         self.empty.model_count,
         self.empty.data_file.model_count
        )
        self.assertEqual(
         self.empty.model_annotations,
         self.empty.data_file.model_annotations
        )
        self.assertEqual(
         self.empty.revisions,
         self.empty.data_file.revisions
        )
        self.assertEqual(
         self.empty.supercedes,
         self.empty.data_file.supercedes
        )
        self.assertEqual(
         self.empty.supercede_date,
         self.empty.data_file.supercede_date
        )
        self.assertEqual(
         self.empty.journal,
         self.empty.data_file.journal
        )



class ModelTests(PdbTest):

    def setUp(self):
        PdbTest.setUp(self)
        self.multi_model = Pdb(PdbDataFile(PdbFile(
         "MODEL        1\n"
         "ATOM    107  N   GLY A  13      12.681  37.302 -25.211 1.000 15.56           N\n"
         "ENDMDL\n"
         "MODEL        2\n"
         "ATOM    107  N   GLY A  13      12.681  37.302 -25.211 1.000 15.56           N\n"
         "ENDMDL"
        )))
        self.single_model = Pdb(PdbDataFile(PdbFile(
         "ATOM    107  N   GLY A  13      12.681  37.302 -25.211 1.000 15.56           N"
        )))


    def test_single_model_case(self):
        self.assertEqual(len(self.single_model.models), 1)
        self.assertEqual(len(self.empty.models), 1)
        self.assertIsInstance(self.single_model.models[0], PdbModel)
        self.assertIsInstance(self.single_model.models[0], PdbModel)


    def test_multi_models(self):
        self.assertEqual(len(self.multi_model.models), 2)
        self.assertIsInstance(self.multi_model.models[0], PdbModel)
        self.assertIsInstance(self.multi_model.models[1], PdbModel)
        self.assertIsNot(self.multi_model.models[0], self.multi_model.models[1])


    def test_model(self):
        self.assertIs(self.empty.model, self.empty.models[0])
        self.assertIs(self.single_model.model, self.single_model.models[0])
        self.assertIs(self.multi_model.model, self.multi_model.models[0])



class SmallMoleculeTests(PdbTest):
    def setUp(self):
        PdbTest.setUp(self)
        self.one_het = Pdb(PdbDataFile(PdbFile(
         "HETATM 3194  C1  BU2 A5001       2.646  45.112  48.995  1.00 43.24           C\n"
         "HETATM 3195  O1  BU2 A5001       1.781  45.484  47.929  1.00 42.82           O\n"
         "HETATM 3196  C2  BU2 A5001       1.922  45.088  50.288  1.00 44.82           C\n"
         "HETATM 3197  C3  BU2 A5001       0.706  44.197  50.309  1.00 43.92           C\n"
         "HETATM 3198  O3  BU2 A5001       1.101  42.889  50.701  1.00 45.94           O\n"
         "HETATM 3199  C4  BU2 A5001      -0.456  44.629  51.162  1.00 42.35           C"
        )))
        self.two_hets = Pdb(PdbDataFile(PdbFile(
         "HETATM 3194  C1  BU2 A5001       2.646  45.112  48.995  1.00 43.24           C\n"
         "HETATM 3195  O1  BU2 A5001       1.781  45.484  47.929  1.00 42.82           O\n"
         "HETATM 3196  C2  BU2 A5001       1.922  45.088  50.288  1.00 44.82           C\n"
         "HETATM 3197  C3  BU2 A5001       0.706  44.197  50.309  1.00 43.92           C\n"
         "HETATM 3198  O3  BU2 A5001       1.101  42.889  50.701  1.00 45.94           O\n"
         "HETATM 3199  C4  BU2 A5001      -0.456  44.629  51.162  1.00 42.35           C\n"
         "HETATM 3224  C1  BU2 B5002     -14.563  61.208  49.005  1.00 45.50           C\n"
         "HETATM 3225  O1  BU2 B5002     -15.048  61.106  50.333  1.00 45.86           O\n"
         "HETATM 3226  C2  BU2 B5002     -15.717  61.232  48.004  1.00 45.69           C\n"
         "HETATM 3227  C3  BU2 B5002     -16.272  59.866  47.666  1.00 45.71           C\n"
         "HETATM 3228  O3  BU2 B5002     -17.648  59.990  47.338  1.00 48.42           O\n"
         "HETATM 3229  C4  BU2 B5002     -15.600  59.073  46.589  1.00 44.02           C"
        )))
        self.two_models = Pdb(PdbDataFile(PdbFile(
         "MODEL        1\n"
         "HETATM 3194  C1  BU2 A5001       2.646  45.112  48.995  1.00 43.24           C\n"
         "HETATM 3195  O1  BU2 A5001       1.781  45.484  47.929  1.00 42.82           O\n"
         "HETATM 3196  C2  BU2 A5001       1.922  45.088  50.288  1.00 44.82           C\n"
         "HETATM 3197  C3  BU2 A5001       0.706  44.197  50.309  1.00 43.92           C\n"
         "HETATM 3198  O3  BU2 A5001       1.101  42.889  50.701  1.00 45.94           O\n"
         "HETATM 3199  C4  BU2 A5001      -0.456  44.629  51.162  1.00 42.35           C\n"
         "ENDMDL\n"
         "MODEL        2\n"
         "HETATM 3194  C1  BU2 A5001       2.646  45.112  48.995  1.00 43.24           C\n"
         "HETATM 3195  O1  BU2 A5001       1.781  45.484  47.929  1.00 42.82           O\n"
         "HETATM 3196  C2  BU2 A5001       1.922  45.088  50.288  1.00 44.82           C\n"
         "HETATM 3197  C3  BU2 A5001       0.706  44.197  50.309  1.00 43.92           C\n"
         "HETATM 3198  O3  BU2 A5001       1.101  42.889  50.701  1.00 45.94           O\n"
         "HETATM 3199  C4  BU2 A5001      -0.456  44.629  51.162  1.00 42.35           C\n"
         "ENDMDL"
        )))


    def test_het_in_model(self):
        self.assertEqual(len(self.one_het.model.small_molecules), 1)
        self.assertIsInstance(
         list(self.one_het.model.small_molecules)[0],
         PdbSmallMolecule
        )
        self.assertEqual(
         list(self.one_het.model.small_molecules)[0].molecule_id,
         5001
        )
        self.assertEqual(
         list(self.one_het.model.small_molecules)[0].molecule_name,
         "BU2"
        )
        self.assertEqual(
         len(list(self.one_het.model.small_molecules)[0].atoms),
         6
        )


    def test_multiple_hets(self):
        self.assertEqual(len(self.two_hets.model.small_molecules), 2)
        for small_molecule in self.two_hets.model.small_molecules:
            self.assertIsInstance(small_molecule, PdbSmallMolecule)


    def test_multi_model_hets(self):
        self.assertEqual(len(self.two_models.models[0].small_molecules), 1)
        self.assertEqual(len(self.two_models.models[1].small_molecules), 1)
        self.assertIsNot(
         list(self.two_models.models[0].small_molecules)[0],
         list(self.two_models.models[1].small_molecules)[0]
        )
        self.assertEqual(
         len(list(self.two_models.model.small_molecules)[0].atoms),
         6
        )
        self.assertEqual(
         len(list(self.two_models.model.small_molecules)[0].atoms),
         6
        )





class PdbFromStringTests(TestCase):

    def test_can_get_pdb_from_string(self):
        pdb = pdb_from_string("TITLE     CRYSTAL")
        self.assertIsInstance(pdb, Pdb)
        self.assertEqual(pdb.title, "CRYSTAL")



class PdbFromFileTests(TestCase):

    def test_can_get_pdb_from_file(self):
        pdb = get_pdb_from_file("tests/1AOI.pdb")
        self.assertIsInstance(pdb, Pdb)
        self.assertEqual(pdb.classification, "DNA BINDING PROTEIN/DNA")



class PdbFromRemotelyTests(TestCase):

    def test_can_get_pdb_remotely(self):
        pdb = get_pdb_remotely("1NVQ")
        self.assertIsInstance(pdb, Pdb)
        self.assertEqual(pdb.pdb_code, "1NVQ")


    def test_invalid_code_raises_error(self):
        with self.assertRaises(InvalidPdbCodeError):
            pdb = get_pdb_remotely("XXXX")
