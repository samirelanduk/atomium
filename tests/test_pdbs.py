import datetime
from unittest import TestCase
from molecupy.pdbfile import PdbFile
from molecupy.pdbdatafile import PdbDataFile
from molecupy.pdb import Pdb
from molecupy.structures import PdbModel, PdbSmallMolecule, PdbChain
from molecupy.exceptions import *

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
         "A5001"
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



class ChainTests(PdbTest):

    def setUp(self):
        PdbTest.setUp(self)
        self.one_chain = Pdb(PdbDataFile(PdbFile(
         "ATOM      1  N   VAL A  11       3.696  33.898  63.219  1.00 21.50           N\n"
         "ATOM      2  CA  VAL A  11       3.198  33.218  61.983  1.00 19.76           C\n"
         "ATOM      3  C   VAL A  11       3.914  31.863  61.818  1.00 19.29           C\n"
         "ATOM      4  O   VAL A  11       5.132  31.792  61.932  1.00 19.78           O\n"
         "ATOM      5  CB  VAL A  11       3.431  34.149  60.743  1.00 22.70           C\n"
         "ATOM      6  CG1 VAL A  11       3.512  33.359  59.474  1.00 20.55           C\n"
         "ATOM      7  CG2 VAL A  11       2.283  35.168  60.648  1.00 21.37           C\n"
         "ATOM      8  N   MET A  12       3.155  30.797  61.557  1.00 17.03           N\n"
         "ATOM      9  CA  MET A  12       3.728  29.464  61.400  1.00 17.91           C\n"
         "ATOM     10  C   MET A  12       4.757  29.459  60.275  1.00 17.01           C\n"
         "ATOM     11  O   MET A  12       4.454  29.842  59.143  1.00 16.20           O\n"
         "ATOM     12  CB  MET A  12       2.609  28.448  61.115  1.00 17.66           C\n"
         "ATOM     13  CG  MET A  12       3.046  26.992  61.089  1.00 19.46           C\n"
         "ATOM     14  SD  MET A  12       1.652  25.909  60.639  1.00 21.70           S\n"
         "ATOM     15  CE  MET A  12       2.419  24.308  60.655  1.00 21.05           C"
        )))
        self.two_chains = Pdb(PdbDataFile(PdbFile(
         "ATOM      1  N   VAL A  11       3.696  33.898  63.219  1.00 21.50           N\n"
         "ATOM      2  CA  VAL A  11       3.198  33.218  61.983  1.00 19.76           C\n"
         "ATOM      3  C   VAL A  11       3.914  31.863  61.818  1.00 19.29           C\n"
         "ATOM      4  O   VAL A  11       5.132  31.792  61.932  1.00 19.78           O\n"
         "ATOM      5  CB  VAL A  11       3.431  34.149  60.743  1.00 22.70           C\n"
         "ATOM      6  CG1 VAL A  11       3.512  33.359  59.474  1.00 20.55           C\n"
         "ATOM      7  CG2 VAL A  11       2.283  35.168  60.648  1.00 21.37           C\n"
         "ATOM      8  N   MET A  12       3.155  30.797  61.557  1.00 17.03           N\n"
         "ATOM      9  CA  MET A  12       3.728  29.464  61.400  1.00 17.91           C\n"
         "ATOM     10  C   MET A  12       4.757  29.459  60.275  1.00 17.01           C\n"
         "ATOM     11  O   MET A  12       4.454  29.842  59.143  1.00 16.20           O\n"
         "ATOM     12  CB  MET A  12       2.609  28.448  61.115  1.00 17.66           C\n"
         "ATOM     13  CG  MET A  12       3.046  26.992  61.089  1.00 19.46           C\n"
         "ATOM     14  SD  MET A  12       1.652  25.909  60.639  1.00 21.70           S\n"
         "ATOM     15  CE  MET A  12       2.419  24.308  60.655  1.00 21.05           C\n"
         "ATOM   1559  N   VAL B1011     -26.384  61.433  36.898  1.00 39.30           N\n"
         "ATOM   1560  CA  VAL B1011     -26.779  61.969  35.563  1.00 40.04           C\n"
         "ATOM   1561  C   VAL B1011     -28.230  62.451  35.541  1.00 39.30           C\n"
         "ATOM   1562  O   VAL B1011     -28.472  63.639  35.306  1.00 39.01           O\n"
         "ATOM   1563  CB  VAL B1011     -26.576  60.922  34.442  1.00 39.96           C\n"
         "ATOM   1564  CG1 VAL B1011     -27.078  61.464  33.120  1.00 41.05           C\n"
         "ATOM   1565  CG2 VAL B1011     -25.101  60.574  34.320  1.00 42.51           C\n"
         "ATOM   1566  N   MET B1012     -29.202  61.566  35.777  1.00 37.52           N\n"
         "ATOM   1567  CA  MET B1012     -30.582  62.053  35.752  1.00 35.97           C\n"
         "ATOM   1568  C   MET B1012     -30.760  63.085  36.844  1.00 33.78           C\n"
         "ATOM   1569  O   MET B1012     -30.490  62.813  38.016  1.00 32.48           O\n"
         "ATOM   1570  CB  MET B1012     -31.622  60.946  35.948  1.00 37.52           C\n"
         "ATOM   1571  CG  MET B1012     -33.059  61.521  35.956  1.00 39.27           C\n"
         "ATOM   1572  SD  MET B1012     -34.423  60.340  35.684  1.00 42.85           S\n"
         "ATOM   1573  CE  MET B1012     -34.740  59.810  37.386  1.00 39.87           C"
        )))
        self.two_models = Pdb(PdbDataFile(PdbFile(
         "MODEL        1\n"
         "ATOM      1  N   VAL A  11       3.696  33.898  63.219  1.00 21.50           N\n"
         "ATOM      2  CA  VAL A  11       3.198  33.218  61.983  1.00 19.76           C\n"
         "ATOM      3  C   VAL A  11       3.914  31.863  61.818  1.00 19.29           C\n"
         "ATOM      4  O   VAL A  11       5.132  31.792  61.932  1.00 19.78           O\n"
         "ATOM      5  CB  VAL A  11       3.431  34.149  60.743  1.00 22.70           C\n"
         "ATOM      6  CG1 VAL A  11       3.512  33.359  59.474  1.00 20.55           C\n"
         "ATOM      7  CG2 VAL A  11       2.283  35.168  60.648  1.00 21.37           C\n"
         "ATOM      8  N   MET A  12       3.155  30.797  61.557  1.00 17.03           N\n"
         "ATOM      9  CA  MET A  12       3.728  29.464  61.400  1.00 17.91           C\n"
         "ATOM     10  C   MET A  12       4.757  29.459  60.275  1.00 17.01           C\n"
         "ATOM     11  O   MET A  12       4.454  29.842  59.143  1.00 16.20           O\n"
         "ATOM     12  CB  MET A  12       2.609  28.448  61.115  1.00 17.66           C\n"
         "ATOM     13  CG  MET A  12       3.046  26.992  61.089  1.00 19.46           C\n"
         "ATOM     14  SD  MET A  12       1.652  25.909  60.639  1.00 21.70           S\n"
         "ATOM     15  CE  MET A  12       2.419  24.308  60.655  1.00 21.05           C\n"
         "ENDMDL\n"
         "MODEL        2\n"
         "ATOM      1  N   VAL A  11       3.696  33.898  63.219  1.00 21.50           N\n"
         "ATOM      2  CA  VAL A  11       3.198  33.218  61.983  1.00 19.76           C\n"
         "ATOM      3  C   VAL A  11       3.914  31.863  61.818  1.00 19.29           C\n"
         "ATOM      4  O   VAL A  11       5.132  31.792  61.932  1.00 19.78           O\n"
         "ATOM      5  CB  VAL A  11       3.431  34.149  60.743  1.00 22.70           C\n"
         "ATOM      6  CG1 VAL A  11       3.512  33.359  59.474  1.00 20.55           C\n"
         "ATOM      7  CG2 VAL A  11       2.283  35.168  60.648  1.00 21.37           C\n"
         "ATOM      8  N   MET A  12       3.155  30.797  61.557  1.00 17.03           N\n"
         "ATOM      9  CA  MET A  12       3.728  29.464  61.400  1.00 17.91           C\n"
         "ATOM     10  C   MET A  12       4.757  29.459  60.275  1.00 17.01           C\n"
         "ATOM     11  O   MET A  12       4.454  29.842  59.143  1.00 16.20           O\n"
         "ATOM     12  CB  MET A  12       2.609  28.448  61.115  1.00 17.66           C\n"
         "ATOM     13  CG  MET A  12       3.046  26.992  61.089  1.00 19.46           C\n"
         "ATOM     14  SD  MET A  12       1.652  25.909  60.639  1.00 21.70           S\n"
         "ATOM     15  CE  MET A  12       2.419  24.308  60.655  1.00 21.05           C\n"
         "ENDMDL"
        )))


    def test_chain_in_model(self):
        self.assertEqual(len(self.one_chain.model.chains), 1)
        self.assertIsInstance(
         list(self.one_chain.model.chains)[0],
         PdbChain
        )
        self.assertEqual(
         list(self.one_chain.model.chains)[0].chain_id,
         "A"
        )
        self.assertEqual(
         len(list(self.one_chain.model.chains)[0].residues),
         2
        )
        self.assertEqual(
         len(list(self.one_chain.model.chains)[0].atoms),
         15
        )


    def test_multiple_chains(self):
        self.assertEqual(len(self.two_chains.model.chains), 2)
        for chain in self.two_chains.model.chains:
            self.assertIsInstance(chain, PdbChain)


    def test_multi_model_chains(self):
        self.assertEqual(len(self.two_models.models[0].chains), 1)
        self.assertEqual(len(self.two_models.models[1].chains), 1)
        self.assertIsNot(
         list(self.two_models.models[0].chains)[0],
         list(self.two_models.models[1].chains)[0]
        )
        self.assertEqual(
         len(list(self.two_models.model.chains)[0].atoms),
         15
        )
        self.assertEqual(
         len(list(self.two_models.model.chains)[0].atoms),
         15
        )
