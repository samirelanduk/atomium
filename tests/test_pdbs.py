import datetime
from unittest import TestCase
from molecupy.pdbfile import PdbFile
from molecupy.pdbdatafile import PdbDataFile
from molecupy.pdb import Pdb, _residue_id_is_greater_than_residue_id
from molecupy.structures import PdbModel, PdbSmallMolecule, PdbChain, PdbSite
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
         "ATOM      8  N   MET A  12A      3.155  30.797  61.557  1.00 17.03           N\n"
         "ATOM      9  CA  MET A  12A      3.728  29.464  61.400  1.00 17.91           C\n"
         "ATOM     10  C   MET A  12A      4.757  29.459  60.275  1.00 17.01           C\n"
         "ATOM     11  O   MET A  12A      4.454  29.842  59.143  1.00 16.20           O\n"
         "ATOM     12  CB  MET A  12A      2.609  28.448  61.115  1.00 17.66           C\n"
         "ATOM     13  CG  MET A  12A      3.046  26.992  61.089  1.00 19.46           C\n"
         "ATOM     14  SD  MET A  12A      1.652  25.909  60.639  1.00 21.70           S\n"
         "ATOM     15  CE  MET A  12A      2.419  24.308  60.655  1.00 21.05           C"
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


    def test_residue_ids(self):
        self.assertEqual(
         self.one_chain.model.get_chain_by_id("A").residues[0].residue_id,
         "A11"
        )
        self.assertEqual(
         self.one_chain.model.get_chain_by_id("A").residues[1].residue_id,
         "A12A"
        )



class ConnectionTests(PdbTest):

    def setUp(self):
        self.pdb = Pdb(PdbDataFile(PdbFile(
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
         "ATOM     16  N   ASN A  13       5.980  29.039  60.600  1.00 16.04           N\n"
         "ATOM     17  CA  ASN A  13       7.092  28.983  59.649  1.00 16.16           C\n"
         "ATOM     18  C   ASN A  13       7.504  30.310  58.987  1.00 16.90           C\n"
         "ATOM     19  O   ASN A  13       8.250  30.315  57.989  1.00 13.64           O\n"
         "ATOM     20  CB  ASN A  13       6.809  27.938  58.561  1.00 16.22           C\n"
         "ATOM     21  CG  ASN A  13       6.952  26.527  59.080  1.00 19.20           C\n"
         "HETATM 3194  C1  BU2 A5001       2.646  45.112  48.995  1.00 43.24           C\n"
         "HETATM 3195  O1  BU2 A5001       1.781  45.484  47.929  1.00 42.82           O\n"
         "HETATM 3196  C2  BU2 A5001       1.922  45.088  50.288  1.00 44.82           C\n"
         "HETATM 3197  C3  BU2 A5001       0.706  44.197  50.309  1.00 43.92           C\n"
         "HETATM 3198  O3  BU2 A5001       1.101  42.889  50.701  1.00 45.94           O\n"
         "HETATM 3199  C4  BU2 A5001      -0.456  44.629  51.162  1.00 42.35           C\n"
         "CONECT 3194 3195 3196\n"
         "CONECT 3195 3194\n"
         "CONECT 3196 3194 3197\n"
         "CONECT 3197 3196 3198 3199\n"
         "CONECT 3198 3197\n"
         "CONECT 3199 3197\n"
        )))
        self.missing_residues_pdb = Pdb(PdbDataFile(PdbFile(
         "REMARK 465\n"
         "REMARK 465 MISSING RESIDUES\n"
         "REMARK 465 THE FOLLOWING RESIDUES WERE NOT LOCATED IN THE\n"
         "REMARK 465 EXPERIMENT. (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN\n"
         "REMARK 465 IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)\n"
         "REMARK 465\n"
         "REMARK 465   M RES C SSSEQI\n"
         "REMARK 465     LEU A     1\n"
         "REMARK 465     ARG A     2\n"
         "REMARK 465     SER A     3\n"
         "REMARK 465     ARG A     4\n"
         "REMARK 465     ARG A     5\n"
         "REMARK 465     VAL A     6\n"
         "REMARK 465     ASP A     7\n"
         "REMARK 465     VAL A     8\n"
         "REMARK 465     MET A     9\n"
         "REMARK 465     ASP A    10\n"
         "REMARK 465     VAL A   14\n"
         "REMARK 465     GLY A   16\n"
         "REMARK 465     GLY A   17A\n"
         "REMARK 465     GLY B   1\n"
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
         "ATOM     16  N   ASN A  13       5.980  29.039  60.600  1.00 16.04           N\n"
         "ATOM     17  CA  ASN A  13       7.092  28.983  59.649  1.00 16.16           C\n"
         "ATOM     18  C   ASN A  13       7.504  30.310  58.987  1.00 16.90           C\n"
         "ATOM     19  O   ASN A  13       8.250  30.315  57.989  1.00 13.64           O\n"
         "ATOM     20  CB  ASN A  13       6.809  27.938  58.561  1.00 16.22           C\n"
         "ATOM     21  CG  ASN A  13       6.952  26.527  59.080  1.00 19.20           C\n"
         "ATOM     22  OD1 ASN A  13       7.957  26.203  59.717  1.00 19.60           O\n"
         "ATOM     23  ND2 ASN A  13       5.960  25.681  58.818  1.00 13.22           N\n"
         "ATOM     35  N   LEU A  15       5.991  32.147  57.112  1.00 12.76           N\n"
         "ATOM     36  CA  LEU A  15       5.548  32.181  55.716  1.00 13.64           C\n"
         "ATOM     37  C   LEU A  15       4.045  32.424  55.621  1.00 13.49           C\n"
         "ATOM     38  O   LEU A  15       3.270  31.589  56.067  1.00 14.06           O\n"
         "ATOM     39  CB  LEU A  15       5.858  30.841  55.045  1.00 12.60           C\n"
         "ATOM     40  CG  LEU A  15       5.461  30.647  53.577  1.00 12.53           C\n"
         "ATOM     41  CD1 LEU A  15       6.056  31.744  52.709  1.00 12.38           C\n"
         "ATOM     42  CD2 LEU A  15       5.953  29.270  53.111  1.00 14.79           C\n"
         "ATOM     51  N   LEU A  17       0.690  32.856  53.234  1.00 14.97           N\n"
         "ATOM     52  CA  LEU A  17       0.201  32.561  51.891  1.00 12.72           C\n"
         "ATOM     53  C   LEU A  17      -0.680  33.711  51.411  1.00 14.43           C\n"
         "ATOM     54  O   LEU A  17      -1.652  34.080  52.071  1.00 10.95           O\n"
         "ATOM     55  CB  LEU A  17      -0.634  31.277  51.879  1.00 12.78           C\n"
         "ATOM     56  CG  LEU A  17      -1.426  31.024  50.587  1.00 14.29           C\n"
         "ATOM     57  CD1 LEU A  17      -0.459  30.721  49.426  1.00 12.74           C\n"
         "ATOM     58  CD2 LEU A  17      -2.376  29.850  50.805  1.00 14.20           C\n"
        )))


    def test_can_bond_atoms_together_from_conect(self):
        self.assertTrue(set((
          self.pdb.model.get_atom_by_id(3195),
          self.pdb.model.get_atom_by_id(3196)
        )).issubset(
         self.pdb.model.get_atom_by_id(3194).get_covalent_bonded_atoms()
        ))
        self.assertTrue(set((
          self.pdb.model.get_atom_by_id(3194),
        )).issubset(
         self.pdb.model.get_atom_by_id(3195).get_covalent_bonded_atoms()
        ))
        self.assertTrue(set((
          self.pdb.model.get_atom_by_id(3194),
          self.pdb.model.get_atom_by_id(3197)
        )).issubset(
         self.pdb.model.get_atom_by_id(3196).get_covalent_bonded_atoms()
        ))
        self.assertTrue(set((
          self.pdb.model.get_atom_by_id(3196),
          self.pdb.model.get_atom_by_id(3198),
          self.pdb.model.get_atom_by_id(3199)
        )).issubset(
         self.pdb.model.get_atom_by_id(3197).get_covalent_bonded_atoms()
        ))
        self.assertTrue(set((
          self.pdb.model.get_atom_by_id(3197),
        )).issubset(
         self.pdb.model.get_atom_by_id(3198).get_covalent_bonded_atoms()
        ))
        self.assertTrue(set((
          self.pdb.model.get_atom_by_id(3197),
        )).issubset(
         self.pdb.model.get_atom_by_id(3199).get_covalent_bonded_atoms()
        ))


    def test_can_connect_residue_atoms(self):
        self.assertTrue(set((
          self.pdb.model.get_atom_by_id(9),
        )).issubset(
         self.pdb.model.get_atom_by_id(8).get_covalent_bonded_atoms()
        ))
        self.assertTrue(set((
          self.pdb.model.get_atom_by_id(8),
          self.pdb.model.get_atom_by_id(10),
          self.pdb.model.get_atom_by_id(12)
        )).issubset(
         self.pdb.model.get_atom_by_id(9).get_covalent_bonded_atoms()
        ))
        self.assertTrue(set((
          self.pdb.model.get_atom_by_id(9),
          self.pdb.model.get_atom_by_id(11)
        )).issubset(
         self.pdb.model.get_atom_by_id(10).get_covalent_bonded_atoms()
        ))
        self.assertTrue(set((
          self.pdb.model.get_atom_by_id(10),
        )).issubset(
         self.pdb.model.get_atom_by_id(11).get_covalent_bonded_atoms()
        ))
        self.assertTrue(set((
          self.pdb.model.get_atom_by_id(9),
          self.pdb.model.get_atom_by_id(13)
        )).issubset(
         self.pdb.model.get_atom_by_id(12).get_covalent_bonded_atoms()
        ))
        self.assertTrue(set((
          self.pdb.model.get_atom_by_id(12),
          self.pdb.model.get_atom_by_id(14)
        )).issubset(
         self.pdb.model.get_atom_by_id(13).get_covalent_bonded_atoms()
        ))
        self.assertTrue(set((
          self.pdb.model.get_atom_by_id(13),
          self.pdb.model.get_atom_by_id(15)
        )).issubset(
         self.pdb.model.get_atom_by_id(14).get_covalent_bonded_atoms()
        ))
        self.assertTrue(set((
          self.pdb.model.get_atom_by_id(14),
        )).issubset(
         self.pdb.model.get_atom_by_id(15).get_covalent_bonded_atoms()
        ))


    def test_can_connect_residues(self):
        self.assertIn(
         self.pdb.model.get_atom_by_id(8),
         self.pdb.model.get_atom_by_id(3).get_covalent_bonded_atoms()
        )
        self.assertIn(
         self.pdb.model.get_atom_by_id(3),
         self.pdb.model.get_atom_by_id(8).get_covalent_bonded_atoms()
        )
        self.assertIn(
         self.pdb.model.get_atom_by_id(16),
         self.pdb.model.get_atom_by_id(10).get_covalent_bonded_atoms()
        )
        self.assertIn(
         self.pdb.model.get_atom_by_id(10),
         self.pdb.model.get_atom_by_id(16).get_covalent_bonded_atoms()
        )
        self.assertIs(
         self.pdb.model.get_chain_by_id("A").residues[0].upstream_residue,
         None
        )
        self.assertIs(
         self.pdb.model.get_chain_by_id("A").residues[0].downstream_residue,
         self.pdb.model.get_chain_by_id("A").residues[1]
        )
        self.assertIs(
         self.pdb.model.get_chain_by_id("A").residues[1].upstream_residue,
         self.pdb.model.get_chain_by_id("A").residues[0]
        )
        self.assertIs(
         self.pdb.model.get_chain_by_id("A").residues[1].downstream_residue,
         self.pdb.model.get_chain_by_id("A").residues[2]
        )
        self.assertIs(
         self.pdb.model.get_chain_by_id("A").residues[2].upstream_residue,
         self.pdb.model.get_chain_by_id("A").residues[1]
        )
        self.assertIs(
         self.pdb.model.get_chain_by_id("A").residues[2].downstream_residue,
         None
        )


    def test_next_residue_id_prediction(self):
        self.assertFalse(_residue_id_is_greater_than_residue_id("A1", "A1"))
        self.assertFalse(_residue_id_is_greater_than_residue_id("A19", "A20"))
        self.assertTrue(_residue_id_is_greater_than_residue_id("A20", "A19"))
        self.assertTrue(_residue_id_is_greater_than_residue_id("A19A", "A19"))
        self.assertTrue(_residue_id_is_greater_than_residue_id("A19F", "A19E"))
        self.assertFalse(_residue_id_is_greater_than_residue_id("A19", "A19A"))
        self.assertFalse(_residue_id_is_greater_than_residue_id("A19A", "A20"))


    def test_chains_list_missing_residues(self):
        self.assertEqual(
         self.missing_residues_pdb.model.get_chain_by_id("A").missing_residues,
         ["A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10", "A14", "A16", "A17A"]
        )


    def test_can_detect_missing_residues(self):
        self.assertIs(
         self.missing_residues_pdb.model.get_chain_by_id("A").residues[2].downstream_residue,
         None
        )
        self.assertNotIn(
         self.missing_residues_pdb.model.get_atom_by_id(35),
         self.missing_residues_pdb.model.get_atom_by_id(18).get_covalent_bonded_atoms()
        )


    def test_can_form_disulphide_bonds(self):
        ssbond_pdb = Pdb(PdbDataFile(PdbFile(
         "SSBOND   1 CYS A  171    CYS B  876                          1555   1555  2.05\n"
         "ATOM   1286  SG  CYS A 171      24.158   4.518 -40.070  1.00 28.25           S\n"
         "ATOM   3510  SG  CYS B 876      25.429   4.012 -41.598  1.00 19.93           S"
        )))
        self.assertIn(
         ssbond_pdb.model.get_atom_by_id(1286),
         ssbond_pdb.model.get_atom_by_id(3510).get_covalent_bonded_atoms()
        )
        self.assertIn(
         ssbond_pdb.model.get_atom_by_id(3510),
         ssbond_pdb.model.get_atom_by_id(1286).get_covalent_bonded_atoms()
        )


    def test_can_make_bonds_from_link_records(self):
        link_pdb = Pdb(PdbDataFile(PdbFile(
         "LINK         OD2 ASP A  10                NA    NA A 489     1555   1555  2.69\n"
         "LINK         O   TYR A  15                NA    NA A 489     1555   1555  2.78\n"
         "ATOM     43  OD2 ASP A  10      13.704   1.167 -28.173  1.00 27.76           O\n"
         "ATOM     73  O   TYR A  15      12.612  -3.402 -26.932  1.00 26.47           O\n"
         "HETATM 3780 NA    NA A 489      12.993  -0.702 -26.370  1.00 22.88          NA"
        )))
        self.assertIn(
         link_pdb.model.get_atom_by_id(3780),
         link_pdb.model.get_atom_by_id(43).get_covalent_bonded_atoms()
        )
        self.assertIn(
         link_pdb.model.get_atom_by_id(3780),
         link_pdb.model.get_atom_by_id(73).get_covalent_bonded_atoms()
        )



class SiteTests(PdbTest):

    def setUp(self):
        PdbTest.setUp(self)
        self.one_site = Pdb(PdbDataFile(PdbFile(
         "SITE     1 AC1  2 VAL A  11  MET A  12\n"
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
         "HETATM 3194  C1  BU2 A5001       2.646  45.112  48.995  1.00 43.24           C\n"
         "HETATM 3195  O1  BU2 A5001       1.781  45.484  47.929  1.00 42.82           O"
        )))
        self.two_sites = Pdb(PdbDataFile(PdbFile(
         "REMARK 800\n"
         "REMARK 800 SITE\n"
         "REMARK 800 SITE_IDENTIFIER: AC1\n"
         "REMARK 800 EVIDENCE_CODE: SOFTWARE\n"
         "REMARK 800 SITE_DESCRIPTION: BINDING SITE FOR RESIDUE BU2 A 5001\n"
         "REMARK 800 SITE_IDENTIFIER: AC2\n"
         "REMARK 800 EVIDENCE_CODE: SOFTWARE\n"
         "REMARK 800 SITE_DESCRIPTION: BINDING SITE FOR RESIDUE BU2 B 5002\n"
         "SITE     1 AC1  2 VAL A  11  MET A  12\n"
         "SITE     2 AC2  2 VAL B1011  MET B1012\n"
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
         "ATOM   1573  CE  MET B1012     -34.740  59.810  37.386  1.00 39.87           C\n"
         "HETATM 3194  C1  BU2 A5001       2.646  45.112  48.995  1.00 43.24           C\n"
         "HETATM 3195  O1  BU2 A5001       1.781  45.484  47.929  1.00 42.82           O\n"
         "HETATM 3224  C1  BU2 B5002     -14.563  61.208  49.005  1.00 45.50           C\n"
         "HETATM 3225  O1  BU2 B5002     -15.048  61.106  50.333  1.00 45.86           O\n"
        )))
        self.two_models = Pdb(PdbDataFile(PdbFile(
         "SITE     1 AC1  2 VAL A  11  MET A  12\n"
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
         "HETATM 3194  C1  BU2 A5001       2.646  45.112  48.995  1.00 43.24           C\n"
         "HETATM 3195  O1  BU2 A5001       1.781  45.484  47.929  1.00 42.82           O\n"
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
         "HETATM 3194  C1  BU2 A5001       2.646  45.112  48.995  1.00 43.24           C\n"
         "HETATM 3195  O1  BU2 A5001       1.781  45.484  47.929  1.00 42.82           O\n"
         "ENDMDL"
        )))


    def test_site_in_model(self):
        self.assertEqual(len(self.one_site.model.sites), 1)
        self.assertIsInstance(
         list(self.one_site.model.sites)[0],
         PdbSite
        )
        self.assertEqual(
         list(self.one_site.model.sites)[0].site_id,
         "AC1"
        )
        self.assertEqual(
         len(list(self.one_site.model.sites)[0].residues),
         2
        )
        self.assertEqual(
         len(list(self.one_site.model.sites)[0].atoms),
         15
        )


    def test_multiple_sites(self):
        self.assertEqual(len(self.two_sites.model.sites), 2)
        for site in self.two_sites.model.sites:
            self.assertIsInstance(site, PdbSite)


    def test_multi_model_sites(self):
        self.assertEqual(len(self.two_models.models[0].sites), 1)
        self.assertEqual(len(self.two_models.models[1].sites), 1)
        self.assertIsNot(
         list(self.two_models.models[0].sites)[0],
         list(self.two_models.models[1].sites)[0]
        )
        self.assertEqual(
         len(list(self.two_models.model.sites)[0].atoms),
         15
        )
        self.assertEqual(
         len(list(self.two_models.model.sites)[0].atoms),
         15
        )


    def test_sites_have_ligands(self):
        self.assertEqual(
         self.two_sites.model.get_site_by_id("AC1").ligand,
         self.two_sites.model.get_small_molecule_by_id("A5001")
        )
        self.assertEqual(
         self.two_sites.model.get_site_by_id("AC2").ligand,
         self.two_sites.model.get_small_molecule_by_id("B5002")
        )


    def test_site_ligand_mapping_when_chain_conflated(self):
        one_site = Pdb(PdbDataFile(PdbFile(
         "REMARK 800\n"
         "REMARK 800 SITE\n"
         "REMARK 800 SITE_IDENTIFIER: AC1\n"
         "REMARK 800 EVIDENCE_CODE: SOFTWARE\n"
         "REMARK 800 SITE_DESCRIPTION: BINDING SITE FOR RESIDUE BU2 A5001\n"
         "SITE     1 AC1  2 VAL A  11  MET A  12\n"
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
         "HETATM 3194  C1  BU2 A5001       2.646  45.112  48.995  1.00 43.24           C\n"
         "HETATM 3195  O1  BU2 A5001       1.781  45.484  47.929  1.00 42.82           O"
        )))
        self.assertIs(
         one_site.model.get_small_molecule_by_id("A5001").get_binding_site(),
         one_site.model.get_site_by_id("AC1")
        )
