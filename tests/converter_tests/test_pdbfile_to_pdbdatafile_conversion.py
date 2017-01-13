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



class KeywdsRecordProcessingTests(PdbFile2PdbDataFileTest):

    def test_missing_keywds_processing(self):
        self.assertEqual(self.empty._keywords, [])


    def test_keywds_processing(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "KEYWDS    TIM BARREL, LYASE"
        ))
        self.assertEqual(
         data_file._keywords,
         ["TIM BARREL", "LYASE"]
        )



class ExpdtaRecordProcessingTests(PdbFile2PdbDataFileTest):

    def test_missing_expdta_processing(self):
        self.assertEqual(self.empty._experimental_techniques, [])


    def test_expdta_processing(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "EXPDTA    NEUTRON DIFFRACTION; X-RAY DIFFRACTION"
        ))
        self.assertEqual(
         data_file._experimental_techniques,
         ["NEUTRON DIFFRACTION", "X-RAY DIFFRACTION"]
        )



class NummdlRecordProcessingTests(PdbFile2PdbDataFileTest):

    def test_missing_nummdl_processing(self):
        self.assertEqual(self.empty._model_count, 1)


    def test_nummdl_processing(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "NUMMDL    2"
        ))
        self.assertEqual(data_file._model_count, 2)



class MdltypRecordProcessingTests(PdbFile2PdbDataFileTest):

    def test_missing_mdltyp_processing(self):
        self.assertEqual(self.empty._model_annotations, [])


    def test_mdltyp_processing(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "MDLTYP    CA ATOMS ONLY, CHAIN A, B, C, D, E, F, G, H, I, J, K ; P ATOMS ONLY,\n"
         "MDLTYP   2 CHAIN X, Y, Z"
        ))
        self.assertEqual(
         data_file._model_annotations,
         [
          "CA ATOMS ONLY, CHAIN A, B, C, D, E, F, G, H, I, J, K",
          "P ATOMS ONLY, CHAIN X, Y, Z"
         ]
        )



class AuthorRecordProcessingTests(PdbFile2PdbDataFileTest):

    def test_missing_author_processing(self):
        self.assertEqual(self.empty._authors, [])


    def test_author_processing(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "AUTHOR    M.B.BERRY,B.MEADOR,T.BILDERBACK,P.LIANG,M.GLASER,\n"
         "AUTHOR   2 G.N.PHILLIPS JR.,T.L.ST. STEVENS"
        ))
        self.assertEqual(
         data_file._authors,
         [
          "M.B.BERRY", "B.MEADOR", "T.BILDERBACK", "P.LIANG", "M.GLASER",
          "G.N.PHILLIPS JR.", "T.L.ST. STEVENS"
         ]
        )



class RevdatRecordProcessingTests(PdbFile2PdbDataFileTest):

    def test_missing_revdat_processing(self):
        self.assertEqual(self.empty._revisions, [])


    def test_revdat_processing(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "REVDAT   4 1 24-FEB-09 1LOL    1       VERSN  COMPND EXPDTA CAVEAT\n"
         "REVDAT   4 2                   1       SOURCE JRNL\n"
         "REVDAT   3   01-APR-03 1LOL    1       JRNL\n"
         "REVDAT   2   14-AUG-02 1LOL    1       DBREF\n"
         "REVDAT   1   07-AUG-02 1LOL    0"
        ))
        self.assertEqual(
         data_file._revisions,
         [{
          "number": 1, "date": datetime.datetime(2002, 8, 7).date(),
          "type": 0, "records": []
         }, {
          "number": 2, "date": datetime.datetime(2002, 8, 14).date(),
          "type": 1, "records": ["DBREF"]
         }, {
          "number": 3, "date": datetime.datetime(2003, 4, 1).date(),
          "type": 1, "records": ["JRNL"]
         }, {
          "number": 4, "date": datetime.datetime(2009, 2, 24).date(),
          "type": 1, "records": ["VERSN", "COMPND", "EXPDTA", "CAVEAT", "SOURCE", "JRNL"]
         }]
        )



class SprsdeRecordProcessingTests(PdbFile2PdbDataFileTest):

    def test_missing_sprsde_processing(self):
        self.assertEqual(self.empty._supercedes, [])
        self.assertEqual(self.empty._supercede_date, None)


    def test_sprsde_processing(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "SPRSDE     27-FEB-95 1GDJ      1LH4 2LH4"
        ))
        self.assertEqual(
         data_file._supercedes,
         ["1LH4", "2LH4"]
        )
        self.assertEqual(
         data_file._supercede_date,
         datetime.datetime(1995, 2, 27).date()
        )



class JrnlRecordProcessingTests(PdbFile2PdbDataFileTest):

    def test_empty_jrnl_processing(self):
        self.assertEqual(self.empty._journal, None)


    def test_jrnl_authors_processing(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "JRNL        AUTH   N.WU,E.F.PAI"
        ))
        self.assertEqual(
         data_file._journal["authors"],
         ["N.WU", "E.F.PAI"]
        )


    def test_empty_jrnl_authors_processing(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile("JRNL"))
        self.assertEqual(data_file._journal["authors"], [])


    def test_jrnl_title_processing(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "JRNL        TITL   CRYSTAL STRUCTURES OF INHIBITOR COMPLEXES REVEAL\n"
         "JRNL        TITL 2 AN ALTERNATE BINDING MODE IN\n"
         "JRNL        TITL 3 OROTIDINE-5'-MONOPHOSPHATE DECARBOXYLASE."
        ))
        self.assertEqual(
         data_file._journal["title"],
         "CRYSTAL STRUCTURES OF INHIBITOR COMPLEXES REVEAL AN ALTERNA"
         "TE BINDING MODE IN OROTIDINE-5'-MONOPHOSPHATE DECARBOXYLASE."
        )


    def test_empty_jrnl_title_processing(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile("JRNL"))
        self.assertEqual(data_file._journal["title"], None)


    def test_jrnl_editors_processing(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "JRNL        EDIT   J.REN,C.NICHOLS,L.BIRD,P.CHAMBERLAIN,K.WEAVER,\n"
         "JRNL        EDIT 2 S.SHORT,D.I.STUART,D.K.STAMMERS"
        ))
        self.assertEqual(
         data_file._journal["editors"],
         [
          "J.REN", "C.NICHOLS", "L.BIRD", "P.CHAMBERLAIN", "K.WEAVER",
          "S.SHORT", "D.I.STUART", "D.K.STAMMERS"
         ]
        )


    def test_empty_jrnl_editors_processing(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile("JRNL"))
        self.assertEqual(data_file._journal["editors"], [])


    def test_jrnl_reference_processing(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "JRNL        REF    J.BIOL.CHEM.                  V. 277 28080 2002"
        ))
        self.assertEqual(
         data_file._journal["reference"],
         {
          "published": True, "publication": "J.BIOL.CHEM.",
          "volume": 277, "page": 28080, "year": 2002
         }
        )

    def test_jrnl_unpublished_reference_processing(self):
        data_file = pdb_data_file_from_pdb_file(
         PdbFile("JRNL        REF    TO BE PUBLISHED")
        )
        self.assertEqual(
         data_file._journal["reference"],
         {
          "published": False, "publication": None,
          "volume": None, "page": None, "year": None
         }
        )


    def test_empty_jrnl_reference_processing(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile("JRNL"))
        self.assertEqual(data_file._journal["reference"], None)


    def test_jrnl_publisher_processing(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "JRNL        PUBL   AMERICAN ASSOCIATION FOR THE ADVANCEMENT OF SCIENCE\n"
         "JRNL        PUBL 2 WASHINGTON, D.C."
        ))
        self.assertEqual(
         data_file._journal["publisher"],
         "AMERICAN ASSOCIATION FOR THE ADVANCEMENT OF SCIENCE WASHINGTON, D.C."
        )


    def test_empty_jrnl_publisher_processing(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile("JRNL"))
        self.assertEqual(data_file._journal["publisher"], None)


    def test_jrnl_referencenumber__processing(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "JRNL        REFN                   ISSN 0021-9258"
        ))
        self.assertEqual(
         data_file._journal["reference_number"],
         {"type": "ISSN", "value": "0021-9258"}
        )


    def test_empty_jrnl_reference_number_processing(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile("JRNL"))
        self.assertEqual(data_file._journal["reference_number"], None)


    def test_jrnl_pubmed_processing(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "JRNL        PMID   12011084"
        ))
        self.assertEqual(
         data_file._journal["pubmed"],
         "12011084"
        )


    def test_empty_jrnl_pubmed_processing(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile("JRNL"))
        self.assertEqual(data_file._journal["pubmed"], None)


    def test_jrnl_doi_processing(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "JRNL        DOI    10.1074/JBC.M202362200"
        ))
        self.assertEqual(
         data_file._journal["doi"],
         "10.1074/JBC.M202362200"
        )


    def test_empty_jrnl_doi_processing(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile("JRNL"))
        self.assertEqual(data_file._journal["doi"], None)


    def test_full_jrnl_processing(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "JRNL        AUTH   N.WU,E.F.PAI\n"
         "JRNL        TITL   CRYSTAL STRUCTURES OF INHIBITOR COMPLEXES REVEAL\n"
         "JRNL        TITL 2 AN ALTERNATE BINDING MODE IN\n"
         "JRNL        TITL 3 OROTIDINE-5'-MONOPHOSPHATE DECARBOXYLASE.\n"
         "JRNL        EDIT   J.REN,C.NICHOLS,L.BIRD,P.CHAMBERLAIN,K.WEAVER,\n"
         "JRNL        EDIT 2 S.SHORT,D.I.STUART,D.K.STAMMERS\n"
         "JRNL        REF    J.BIOL.CHEM.                  V. 277 28080 2002\n"
         "JRNL        PUBL   AMERICAN ASSOCIATION FOR THE ADVANCEMENT OF SCIENCE\n"
         "JRNL        PUBL 2 WASHINGTON, D.C.\n"
         "JRNL        REFN                   ISSN 0021-9258\n"
         "JRNL        PMID   12011084\n"
         "JRNL        DOI    10.1074/JBC.M202362200"
        ))
        self.assertEqual(
         data_file._journal,
         {
          "authors": ["N.WU", "E.F.PAI"],
          "title": "CRYSTAL STRUCTURES OF INHIBITOR COMPLEXES REVEAL AN ALTERNA"
          "TE BINDING MODE IN OROTIDINE-5'-MONOPHOSPHATE DECARBOXYLASE.",
          "editors": [
           "J.REN", "C.NICHOLS", "L.BIRD", "P.CHAMBERLAIN", "K.WEAVER",
           "S.SHORT", "D.I.STUART", "D.K.STAMMERS"
          ],
          "reference": {
           "published": True, "publication": "J.BIOL.CHEM.",
           "volume": 277, "page": 28080, "year": 2002
          },
          "publisher": "AMERICAN ASSOCIATION FOR THE ADVANCEMENT OF SCIENCE WASHINGTON, D.C.",
          "reference_number": {"type": "ISSN", "value": "0021-9258"},
          "pubmed": "12011084",
          "doi": "10.1074/JBC.M202362200"
         }
        )



class RemarkRecordTests(PdbFile2PdbDataFileTest):

    def test_missing_remark_processing(self):
        self.assertEqual(self.empty._remarks, [])


    def test_remark_processing(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "REMARK   2\n"
         "REMARK   2 RESOLUTION.    1.90 ANGSTROMS.\n"
         "REMARK 999\n"
         "REMARK 999  SEQUENCE\n"
         "REMARK 999 AUTHOR STATES THAT ALTHOUGH RESIDUES 1 AND 1001 ARE MET\n"
         "REMARK 999 AND RESIDUES 101 AND 1101 ARE ARG ACCORDING TO THE\n"
         "REMARK 999 SWISSPROT ENTRY, RESIDUES 1 AND 1001 WERE LEU AND RESIDUES\n"
         "REMARK 999 101 AND 1101 WERE PRO IN THE ORIGINAL CONSTRUCT CLONED\n"
         "REMARK 999 OF MT GENOMIC DNA."
        ))
        self.assertEqual(
         data_file._remarks,
         [{
          "number": 2,
          "content": "RESOLUTION.    1.90 ANGSTROMS."
         }, {
          "number": 999,
          "content": "SEQUENCE\n"
          "AUTHOR STATES THAT ALTHOUGH RESIDUES 1 AND 1001 ARE MET\n"
          "AND RESIDUES 101 AND 1101 ARE ARG ACCORDING TO THE\n"
          "SWISSPROT ENTRY, RESIDUES 1 AND 1001 WERE LEU AND RESIDUES\n"
          "101 AND 1101 WERE PRO IN THE ORIGINAL CONSTRUCT CLONED\n"
          "OF MT GENOMIC DNA."
         }]
        )



class DbrefRecordProcessingTests(PdbFile2PdbDataFileTest):

    def test_missing_dbref_processing(self):
        self.assertEqual(self.empty._dbreferences, [])


    def test_dbref_processing(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "DBREF  1LOL A    1   229  UNP    O26232   PYRF_METTH       1    228\n"
         "DBREF  1LOL B 1001  1229  UNP    O26232   PYRF_METTH       1    228"
        ))
        self.assertEqual(
         data_file._dbreferences,
         [{
          "chain_id": "A",
          "sequence_begin": 1,
          "insert_begin": "",
          "sequence_end": 229,
          "insert_end": "",
          "database": "UNP",
          "accession": "O26232",
          "db_id": "PYRF_METTH",
          "db_sequence_begin": 1,
          "db_insert_begin": None,
          "db_sequence_end": 228,
          "db_insert_end": None
         }, {
          "chain_id": "B",
          "sequence_begin": 1001,
          "insert_begin": "",
          "sequence_end": 1229,
          "insert_end": "",
          "database": "UNP",
          "accession": "O26232",
          "db_id": "PYRF_METTH",
          "db_sequence_begin": 1,
          "db_insert_begin": None,
          "db_sequence_end": 228,
          "db_insert_end": None
         }]
        )


    def test_long_dbref_processing(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "DBREF1 1LOL C   61   322  GB                   AE017221\n"
         "DBREF2 1LOL C     46197919                      1534489     1537377"
        ))
        self.assertEqual(
         data_file._dbreferences,
         [{
          "chain_id": "C",
          "sequence_begin": 61,
          "insert_begin": "",
          "sequence_end": 322,
          "insert_end": "",
          "database": "GB",
          "accession": "46197919",
          "db_id": "AE017221",
          "db_sequence_begin": 1534489,
          "db_insert_begin": None,
          "db_sequence_end": 1537377,
          "db_insert_end": None
         }]
        )


    def test_mixed_dbref_processing(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "DBREF  1LOL A    1   229  UNP    O26232   PYRF_METTH       1    228\n"
         "DBREF  1LOL B 1001  1229  UNP    O26232   PYRF_METTH       1    228\n"
         "DBREF1 1LOL C   61   322  GB                   AE017221\n"
         "DBREF2 1LOL C     46197919                      1534489     1537377"
        ))
        self.assertEqual(
         data_file._dbreferences,
         [{
          "chain_id": "A",
          "sequence_begin": 1,
          "insert_begin": "",
          "sequence_end": 229,
          "insert_end": "",
          "database": "UNP",
          "accession": "O26232",
          "db_id": "PYRF_METTH",
          "db_sequence_begin": 1,
          "db_insert_begin": None,
          "db_sequence_end": 228,
          "db_insert_end": None
         }, {
          "chain_id": "B",
          "sequence_begin": 1001,
          "insert_begin": "",
          "sequence_end": 1229,
          "insert_end": "",
          "database": "UNP",
          "accession": "O26232",
          "db_id": "PYRF_METTH",
          "db_sequence_begin": 1,
          "db_insert_begin": None,
          "db_sequence_end": 228,
          "db_insert_end": None
         }, {
          "chain_id": "C",
          "sequence_begin": 61,
          "insert_begin": "",
          "sequence_end": 322,
          "insert_end": "",
          "database": "GB",
          "accession": "46197919",
          "db_id": "AE017221",
          "db_sequence_begin": 1534489,
          "db_insert_begin": None,
          "db_sequence_end": 1537377,
          "db_insert_end": None
         }]
        )



class SeqadvRecordProcessingTests(PdbFile2PdbDataFileTest):

    def test_missing_seqadv_processing(self):
        self.assertEqual(self.empty._sequence_differences, [])


    def test_seqadv_processing(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "SEQADV 1LOL GLU A  229  UNP  O26232              INSERTION\n"
         "SEQADV 1LOL GLU B 1229  UNP  O26232              INSERTION"
        ))
        self.assertEqual(
         data_file._sequence_differences,
         [{
          "residue_name": "GLU",
          "chain_id": "A",
          "residue_id": 229,
          "insert_code": "",
          "database": "UNP",
          "accession": "O26232",
          "db_residue_name": None,
          "db_residue_id": None,
          "conflict": "INSERTION"
         }, {
          "residue_name": "GLU",
          "chain_id": "B",
          "residue_id": 1229,
          "insert_code": "",
          "database": "UNP",
          "accession": "O26232",
          "db_residue_name": None,
          "db_residue_id": None,
          "conflict": "INSERTION"
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
