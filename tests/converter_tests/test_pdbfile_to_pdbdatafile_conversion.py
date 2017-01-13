import datetime
from unittest import TestCase
from unittest.mock import Mock
from molecupy.pdb.pdbfile import PdbFile, PdbRecord
from molecupy.pdb.pdbdatafile import PdbDataFile
from molecupy.converters.pdbfile2pdbdatafile import pdb_data_file_from_pdb_file
from molecupy.converters.pdbfile2pdbdatafile import date_from_string
from molecupy.converters.pdbfile2pdbdatafile import merge_records
from molecupy.converters.pdbfile2pdbdatafile import records_to_token_value_dicts

class PdbFile2PdbDataFileTest(TestCase):

    def setUp(self):
        self.empty = pdb_data_file_from_pdb_file(PdbFile(""))



class BasicPdbDataFileCreationTests(PdbFile2PdbDataFileTest):

    def test_can_create_pdb_data_file(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(""))
        self.assertIsInstance(data_file, PdbDataFile)


    def test_can_only_convert_pdb_files(self):
        with self.assertRaises(TypeError):
            pdb_data_file_from_pdb_file("PDB file")



class HeaderRecordProcessingTests(PdbFile2PdbDataFileTest):

    def test_missing_header_processing(self):
        self.assertEqual(self.empty._classification, None)
        self.assertEqual(self.empty._deposition_date, None)
        self.assertEqual(self.empty._pdb_code, None)


    def test_header_processing(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "HEADER    LYASE                                   06-MAY-02   1LOL"
        ))
        self.assertEqual(data_file._classification, "LYASE")
        self.assertEqual(
         data_file._deposition_date,
         datetime.datetime(2002, 5, 6).date()
        )
        self.assertEqual(data_file._pdb_code, "1LOL")



class ObslteRecordProcessingTests(PdbFile2PdbDataFileTest):

    def test_missing_header_processing(self):
        self.assertFalse(self.empty._is_obsolete)
        self.assertEqual(self.empty._obsolete_date, None)
        self.assertEqual(self.empty._replacement_code, None)


    def test_obslte_processing(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "OBSLTE     30-SEP-93 1LOL      1SAM"
        ))
        self.assertTrue(data_file._is_obsolete)
        self.assertEqual(
         data_file._obsolete_date,
         datetime.datetime(1993, 9, 30).date()
        )
        self.assertEqual(
         data_file._replacement_code,
         "1SAM"
        )



class TitleRecordProcessingTests(PdbFile2PdbDataFileTest):

    def test_missing_title_processing(self):
        self.assertEqual(self.empty._title, None)


    def test_title_processing(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "TITLE     CRYSTAL STRUCTURE OF OROTIDINE MONOPHOSPHATE DECARBOXYLASE\n"
         "TITLE    2 COMPLEX WITH XMP"
        ))
        self.assertEqual(
         data_file._title,
         "CRYSTAL STRUCTURE OF OROTIDINE MONOPHOSPHATE DECARBOXYLASE COMPLEX WITH XMP"
        )



class SplitRecordProcessingTests(PdbFile2PdbDataFileTest):

    def test_missing_split_processing(self):
        self.assertEqual(self.empty._split_codes, [])


    def test_split_processing(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "SPLIT      1VOQ 1VOR 1VOS 1VOU 1VOV 1VOW 1VOX 1VOY 1VP0 1VOZ 1VOY 1VP0 1VOZ 1VOZ\n"
         "SPLIT      1VOQ 1VOR 1VOS 1VOU 1VOV 1VOW 1VOX 1VOY 1VP0 1VOZ"
        ))
        self.assertEqual(
         data_file._split_codes,
         [
          "1VOQ", "1VOR", "1VOS", "1VOU", "1VOV", "1VOW",
          "1VOX", "1VOY", "1VP0", "1VOZ", "1VOY", "1VP0",
          "1VOZ", "1VOZ", "1VOQ", "1VOR", "1VOS", "1VOU",
          "1VOV", "1VOW", "1VOX", "1VOY", "1VP0", "1VOZ"
         ]
        )



class CaveatRecordProcessingTests(PdbFile2PdbDataFileTest):

    def test_missing_caveat_processing(self):
        self.assertEqual(self.empty._caveat, None)


    def test_caveat_processing(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "CAVEAT     1SAM    THE CRYSTAL TRANSFORMATION IS IN ERROR BUT IS\n"
         "CAVEAT   2 1SAM    UNCORRECTABLE AT THIS TIME"
        ))
        self.assertEqual(
         data_file._caveat,
         "THE CRYSTAL TRANSFORMATION IS IN ERROR BUT IS UNCORRECTABLE AT THIS TIME"
        )



class CompndRecordProcessingTests(PdbFile2PdbDataFileTest):

    def test_missing_compnd_processing(self):
        self.assertEqual(self.empty._compounds, [])


    def test_compnd_processing_single_compound(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "COMPND    MOL_ID: 1;\n"
         "COMPND   2 MOLECULE: OROTIDINE 5'-MONOPHOSPHATE DECARBOXYLASE;\n"
         "COMPND   3 CHAIN: A, B;\n"
         "COMPND   4 SYNONYM: OMP DECARBOXYLASE, OMPDCASE, OMPDECASE;\n"
         "COMPND   5 EC: 4.1.1.23;\n"
         "COMPND   6 ENGINEERED: YES;"
        ))
        self.assertEqual(
         data_file._compounds,
         [{
          "MOL_ID": 1,
          "MOLECULE": "OROTIDINE 5'-MONOPHOSPHATE DECARBOXYLASE",
          "CHAIN": ["A", "B"],
          "SYNONYM": [
           "OMP DECARBOXYLASE",
           "OMPDCASE",
           "OMPDECASE"
          ],
          "EC": "4.1.1.23",
          "ENGINEERED": True
         }]
        )


    def test_compnd_processing_multiple_compound(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "COMPND    MOL_ID: 1;\n"
         "COMPND   2 MOLECULE: OROTIDINE 5'-MONOPHOSPHATE DECARBOXYLASE;\n"
         "COMPND   3 CHAIN: A, B;\n"
         "COMPND   4 SYNONYM: OMP DECARBOXYLASE, OMPDCASE, OMPDECASE;\n"
         "COMPND   5 EC: 4.1.1.23;\n"
         "COMPND   6 ENGINEERED: YES;\n"
         "COMPND   7 MOL_ID: 2;\n"
         "COMPND   8 MOLECULE: OROTIDINE 5'-MONOPHOSPHATE DECARBOXYLASE\n"
         "COMPND   9 PLUS;"
        ))
        self.assertEqual(
         data_file._compounds,
         [{
          "MOL_ID": 1,
          "MOLECULE": "OROTIDINE 5'-MONOPHOSPHATE DECARBOXYLASE",
          "CHAIN": ["A", "B"],
          "SYNONYM": [
           "OMP DECARBOXYLASE",
           "OMPDCASE",
           "OMPDECASE"
          ],
          "EC": "4.1.1.23",
          "ENGINEERED": True
         }, {
          "MOL_ID": 2,
          "MOLECULE": "OROTIDINE 5'-MONOPHOSPHATE DECARBOXYLASE PLUS"
         }]
        )



class SourceRecordProcessingTests(PdbFile2PdbDataFileTest):

    def test_missing_source_processing(self):
        self.assertEqual(self.empty._sources, [])


    def test_source_processing_single_source(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
          "SOURCE    MOL_ID: 1;\n"
          "SOURCE   2 ORGANISM_SCIENTIFIC: METHANOTHERMOBACTER\n"
          "SOURCE   3 THERMAUTOTROPHICUS STR. DELTA H;\n"
          "SOURCE   4 ORGANISM_TAXID: 187420;\n"
          "SOURCE   5 STRAIN: DELTA H;\n"
          "SOURCE   6 EXPRESSION_SYSTEM: ESCHERICHIA COLI;\n"
          "SOURCE   7 EXPRESSION_SYSTEM_TAXID: 562;\n"
          "SOURCE   8 EXPRESSION_SYSTEM_PLASMID: PET15B\n"
        ))
        self.assertEqual(
         data_file._sources,
         [{
          "MOL_ID": 1,
          "ORGANISM_SCIENTIFIC": "METHANOTHERMOBACTER THERMAUTOTROPHICUS STR. DELTA H",
          "ORGANISM_TAXID": 187420,
          "STRAIN": "DELTA H",
          "EXPRESSION_SYSTEM": "ESCHERICHIA COLI",
          "EXPRESSION_SYSTEM_TAXID": 562,
          "EXPRESSION_SYSTEM_PLASMID": "PET15B"
         }]
        )


    def test_source_processing_multiple_sources(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
          "SOURCE    MOL_ID: 1;\n"
          "SOURCE   2 ORGANISM_SCIENTIFIC: METHANOTHERMOBACTER\n"
          "SOURCE   3 THERMAUTOTROPHICUS STR. DELTA H;\n"
          "SOURCE   4 ORGANISM_TAXID: 187420;\n"
          "SOURCE   5 STRAIN: DELTA H;\n"
          "SOURCE   6 EXPRESSION_SYSTEM: ESCHERICHIA COLI;\n"
          "SOURCE   7 EXPRESSION_SYSTEM_TAXID: 562;\n"
          "SOURCE   8 MOL_ID: 2;\n"
          "SOURCE   9 ORGANISM_SCIENTIFIC: METHANOTHERMOBACTER;\n"
        ))
        self.assertEqual(
         data_file._sources,
         [{
          "MOL_ID": 1,
          "ORGANISM_SCIENTIFIC": "METHANOTHERMOBACTER THERMAUTOTROPHICUS STR. DELTA H",
          "ORGANISM_TAXID": 187420,
          "STRAIN": "DELTA H",
          "EXPRESSION_SYSTEM": "ESCHERICHIA COLI",
          "EXPRESSION_SYSTEM_TAXID": 562
         }, {
          "MOL_ID": 2,
          "ORGANISM_SCIENTIFIC": "METHANOTHERMOBACTER"
         }]
        )



class DateFromStringTests(PdbFile2PdbDataFileTest):

    def test_can_get_date_from_string(self):
        self.assertEqual(
         date_from_string("01-JAN-00"),
         datetime.datetime(2000, 1, 1).date()
        )
        self.assertEqual(
         date_from_string("28-SEP-99"),
         datetime.datetime(1999, 9, 28).date()
        )


    def test_date_conversion_will_return_none_if_given_nothing(self):
        self.assertEqual(date_from_string(""), None)
        self.assertEqual(date_from_string(None), None)



class RecordMergingTests(TestCase):

    def setUp(self):
        self.records = [PdbRecord(l) for l in [
         "0123456789",
         "abcdefghij",
         "0123456789"
        ]]

        self.punc_records = [PdbRecord(l) for l in [
         "0123, 456789",
         "abcd  efghij",
         "0123; 456789"
        ]]


    def test_can_merge_records(self):
        self.assertEqual(
         merge_records(self.records, 5),
         "56789 fghij 56789"
        )
        self.assertEqual(
         merge_records(self.records, 8),
         "89 ij 89"
        )


    def test_can_vary_join(self):
        self.assertEqual(
         merge_records(self.records, 5, join=""),
         "56789fghij56789"
        )
        self.assertEqual(
         merge_records(self.records, 8, join="."),
         "89.ij.89"
        )


    def test_can_condense(self):
        self.assertEqual(
         merge_records(self.punc_records, 2),
         "23,456789 cd efghij 23;456789"
        )


    def test_can_ignore_condensors(self):
        self.assertEqual(
         merge_records(self.punc_records, 2, dont_condense=","),
         "23, 456789 cd efghij 23;456789"
        )
        self.assertEqual(
         merge_records(self.punc_records, 2, dont_condense=";"),
         "23,456789 cd efghij 23; 456789"
        )
        self.assertEqual(
         merge_records(self.punc_records, 2, dont_condense=";,"),
         "23, 456789 cd efghij 23; 456789"
        )
        self.assertEqual(
         merge_records(self.punc_records, 2, dont_condense=";, "),
         "23, 456789 cd  efghij 23; 456789"
        )



class RecordsToDictTests(TestCase):

    def test_can_make_dicts(self):
        records = [PdbRecord(l) for l in [
         "COMPND    MOL_ID: A;",
         "COMPND   2 MOLECULE: MOLNAME;",
         "COMPND   3 CHAIN_: CHAINS;",
         "COMPND   4 MOL_ID: B;",
         "COMPND   5 MOLECULE: MOLNAME2;",
         "COMPND   6 CHAIN_: CHAINS2;"
        ]]
        self.assertEqual(
         records_to_token_value_dicts(records),
         [
          {"MOL_ID": "A", "MOLECULE": "MOLNAME", "CHAIN_": "CHAINS"},
          {"MOL_ID": "B", "MOLECULE": "MOLNAME2", "CHAIN_": "CHAINS2"}
         ]
        )


    def test_can_detect_numeric_fields(self):
        records = [PdbRecord(l) for l in [
         "COMPND    MOL_ID: 1;",
         "COMPND   2 MOLECULE: MOLNAME;",
         "COMPND   3 CHAIN_: CHAINS;"
        ]]
        self.assertEqual(
         records_to_token_value_dicts(records),
         [
          {"MOL_ID": 1, "MOLECULE": "MOLNAME", "CHAIN_": "CHAINS"}
         ]
        )


    def test_can_detect_boolean_fields(self):
        records = [PdbRecord(l) for l in [
         "COMPND    MOL_ID: 1;",
         "COMPND   2 MOLECULE: YES;",
         "COMPND   3 CHAIN_: NO;"
        ]]
        self.assertEqual(
         records_to_token_value_dicts(records),
         [
          {"MOL_ID": 1, "MOLECULE": True, "CHAIN_": False}
         ]
        )


    def test_can_split_chains_and_synonyms(self):
        records = [PdbRecord(l) for l in [
         "COMPND    MOL_ID: 1;",
         "COMPND   2 CHAIN: A,B;",
         "COMPND   2 SYNONYM: BEEP, BOOP;"
        ]]
        self.assertEqual(
         records_to_token_value_dicts(records),
         [
          {"MOL_ID": 1, "CHAIN": ["A", "B"], "SYNONYM": ["BEEP", "BOOP"]}
         ]
        )


    def test_can_account_for_people_not_knowing_how_to_format_source_records_semi_colons(self):
        records = [PdbRecord(l) for l in [
         "COMPND    MOL_ID: 1;",
         "COMPND   2 FIELD: VALUE;",
         "COMPND   2 FIELD2: VALUE2; EXTRA;",
         "COMPND   2 FIELD3: VALUE3;"
        ]]
        self.assertEqual(
         records_to_token_value_dicts(records),
         [
          {"MOL_ID": 1, "FIELD": "VALUE", "FIELD2": "VALUE2; EXTRA", "FIELD3": "VALUE3"}
         ]
        )
