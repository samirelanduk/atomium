import datetime
from unittest import TestCase
from molecupy.parsers.pdb.pdb_file import PdbFile
from molecupy.parsers.pdb.pdb_data_file import PdbDataFile, date_from_string

class PdbDataFileTest(TestCase):

    def setUp(self):
         with open("tests/parsers/pdb/test.pdb") as f:
            pdb_file = PdbFile(f.read())
            self.data_file = PdbDataFile(pdb_file)
            self.empty_data_file = PdbDataFile(PdbFile(""))


    def check_pdb_data_file(self, data_file):
        self.assertIsInstance(data_file, PdbDataFile)
        self.assertIsInstance(data_file.pdb_file, PdbFile)
        self.assertRegex(
         str(data_file),
         r"<PdbDataFile \(([^\s]{4})\)>"
        )


    def test_can_create_pdb_data_file(self):
        self.check_pdb_data_file(self.data_file)
        self.check_pdb_data_file(self.empty_data_file)


    def test_can_get_date_from_string(self):
        self.assertEqual(
         date_from_string("01-JAN-00"),
         datetime.datetime(2000, 1, 1).date()
        )
        self.assertEqual(
         date_from_string("28-SEP-99"),
         datetime.datetime(1999, 9, 28).date()
        )



class TitleSectionTest(PdbDataFileTest):

    def test_header(self):
        self.assertEqual(self.data_file.classification, "LYASE")
        self.assertEqual(self.empty_data_file.classification, None)
        self.assertEqual(
         self.data_file.deposition_date, datetime.datetime(2002, 5, 6).date()
        )
        self.assertEqual(self.empty_data_file.deposition_date, None)
        self.assertEqual(self.data_file.pdb_code, "1LOL")
        self.assertEqual(self.empty_data_file.pdb_code, None)


    def test_obslte(self):
        self.assertTrue(self.data_file.is_obsolete)
        self.assertFalse(self.empty_data_file.is_obsolete)
        self.assertEqual(
         self.data_file.obsolete_date,
         datetime.datetime(1993, 9, 30).date()
        )
        self.assertEqual(self.empty_data_file.obsolete_date, None)
        self.assertEqual(
         self.data_file.replacement_code,
         "1SAM"
        )
        self.assertEqual(self.empty_data_file.replacement_code, None)


    def test_title(self):
        self.assertEqual(
         self.data_file.title,
         "CRYSTAL STRUCTURE OF OROTIDINE MONOPHOSPHATE DECARBOXYLASE COMPLEX WITH XMP"
        )
        self.assertEqual(self.empty_data_file.title, None)


    def test_split(self):
        self.assertEqual(
         self.data_file.split_codes,
         [
          "1VOQ", "1VOR", "1VOS", "1VOU", "1VOV", "1VOW",
          "1VOX", "1VOY", "1VP0", "1VOZ", "1VOY", "1VP0",
          "1VOZ", "1VOZ", "1VOQ", "1VOR", "1VOS", "1VOU",
          "1VOV", "1VOW", "1VOX", "1VOY", "1VP0", "1VOZ"
         ]
        )
        self.assertEqual(self.empty_data_file.split_codes, [])


    def test_caveat(self):
        self.assertEqual(
         self.data_file.caveat,
         "THE CRYSTAL TRANSFORMATION IS IN ERROR BUT IS UNCORRECTABLE AT THIS TIME"
        )
        self.assertEqual(self.empty_data_file.caveat, None)


    def test_compnd(self):
        self.assertEqual(
         self.data_file.compounds,
         [
          {
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
          }
         ]
        )
        self.assertEqual(self.empty_data_file.compounds, [])


    def test_source(self):
        self.assertEqual(
         self.data_file.sources,
         [
          {
           "MOL_ID": 1,
           "ORGANISM_SCIENTIFIC": "METHANOTHERMOBACTER THERMAUTOTROPHICUS STR. DELTA H",
           "ORGANISM_TAXID": 187420,
           "STRAIN": "DELTA H",
           "EXPRESSION_SYSTEM": "ESCHERICHIA COLI",
           "EXPRESSION_SYSTEM_TAXID": 562,
           "EXPRESSION_SYSTEM_PLASMID": "PET15B"
          }
         ]
        )
        self.assertEqual(self.empty_data_file.sources, [])


    def test_keywds(self):
        self.assertEqual(
         self.data_file.keywords,
         ["TIM BARREL", "LYASE"]
        )
        self.assertEqual(self.empty_data_file.keywords, [])


    def test_expdata(self):
        self.assertEqual(
         self.data_file.experimental_techniques,
         ["NEUTRON DIFFRACTION", "X-RAY DIFFRACTION"]
        )
        self.assertEqual(self.empty_data_file.experimental_techniques, [])


    def test_nummdl(self):
        self.assertEqual(self.data_file.model_num, 2)
        self.assertEqual(self.empty_data_file.model_num, 1)


    def test_mdltyp(self):
        self.assertEqual(
         self.data_file.model_annotations,
         [
          "CA ATOMS ONLY, CHAIN A, B, C, D, E, F, G, H, I, J, K",
          "P ATOMS ONLY, CHAIN X, Y, Z"
         ]
        )
        self.assertEqual(self.empty_data_file.model_annotations, [])


    def test_author(self):
        self.assertEqual(
         self.data_file.authors,
         [
          "N.WU", "E.F.PAI"
         ]
        )
        self.assertEqual(self.empty_data_file.authors, [])


    def test_revdat(self):
        self.assertEqual(
         self.data_file.revisions,
         [
          {
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
          }
         ]
        )
        self.assertEqual(self.empty_data_file.revisions, [])


    def test_sprsde(self):
        self.assertEqual(
         self.data_file.supercedes,
         ["1LH4", "2LH4"]
        )
        self.assertEqual(self.empty_data_file.supercedes, [])
        self.assertEqual(
         self.data_file.supercede_date,
         datetime.datetime(1995, 2, 27).date()
        )
        self.assertEqual(self.empty_data_file.supercede_date, None)


    def test_jrnl(self):
        self.assertEqual(
         self.data_file.journal,
         {
          "authors": ["N.WU", "E.F.PAI"],
          "title": "CRYSTAL STRUCTURES OF INHIBITOR COMPLEXES REVEAL AN ALTERNA"
          "TE BINDING MODE IN OROTIDINE-5'-MONOPHOSPHATE DECARBOXYLASE.",
          "editors": [
           "J.REN", "C.NICHOLS", "L.BIRD", "P.CHAMBERLAIN", "K.WEAVER",
           "S.SHORT", "D.I.STUART", "D.K.STAMMERS"
          ],
          "reference": {"published": True, "publication": "J.BIOL.CHEM.", "volume": 277, "page": 28080, "year": 2002},
          "publisher": "AMERICAN ASSOCIATION FOR THE ADVANCEMENT OF SCIENCE WASHINGTON, D.C.",
          "reference_number": {"type": "ISSN", "value": "0021-9258"},
          "pubmed": "12011084",
          "doi": "10.1074/JBC.M202362200"
         }
        )
        self.assertEqual(self.empty_data_file.journal, None)


    def test_remark(self):
        self.assertIn(
         {
          "number": 2,
          "content": "RESOLUTION.    1.90 ANGSTROMS."
         }, self.data_file.remarks
        )
        self.assertIn(
         {
           "number": 999,
           "content": "SEQUENCE\n"
           "AUTHOR STATES THAT ALTHOUGH RESIDUES 1 AND 1001 ARE MET\n"
           "AND RESIDUES 101 AND 1101 ARE ARG ACCORDING TO THE\n"
           "SWISSPROT ENTRY, RESIDUES 1 AND 1001 WERE LEU AND RESIDUES\n"
           "101 AND 1101 WERE PRO IN THE ORIGINAL CONSTRUCT CLONED\n"
           "OF MT GENOMIC DNA."
          }, self.data_file.remarks
        )
        self.assertEqual(self.empty_data_file.remarks, [])



class PrimaryStructureSectionTest(PdbDataFileTest):

    def test_dbref(self):
        self.assertEqual(
         self.data_file.dbreferences,
         [
          {
           "chain_id": "A",
           "sequence_begin": 1,
           "insert_begin": None,
           "sequence_end": 229,
           "insert_end": None,
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
           "insert_begin": None,
           "sequence_end": 1229,
           "insert_end": None,
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
           "insert_begin": None,
           "sequence_end": 322,
           "insert_end": None,
           "database": "GB",
           "accession": "46197919",
           "db_id": "AE017221",
           "db_sequence_begin": 1534489,
           "db_insert_begin": None,
           "db_sequence_end": 1537377,
           "db_insert_end": None
          }
         ]
        )
        self.assertEqual(self.empty_data_file.dbreferences, [])


    def test_seqadv(self):
        self.assertEqual(
         self.data_file.sequence_differences,
         [
          {
           "residue_name": "GLU",
           "chain_id": "A",
           "residue_id": 229,
           "insert_code": None,
           "database": "UNP",
           "accession": "O26232",
           "db_residue_name": None,
           "db_residue_id": None,
           "conflict": "INSERTION"
          }, {
           "residue_name": "GLU",
           "chain_id": "B",
           "residue_id": 1229,
           "insert_code": None,
           "database": "UNP",
           "accession": "O26232",
           "db_residue_name": None,
           "db_residue_id": None,
           "conflict": "INSERTION"
          }
         ]
        )
        self.assertEqual(self.empty_data_file.sequence_differences, [])


    def test_seqres(self):
        self.assertEqual(
         self.data_file.residue_sequences,
         [
          {
           "chain_id": "A",
           "length": 229,
           "residues": [
            "LEU", "ARG", "SER", "ARG", "ARG", "VAL", "ASP", "VAL",
            "MET", "ASP", "VAL", "MET", "ASN", "ARG", "LEU", "ILE",
            "LEU", "ALA", "MET", "ASP", "LEU", "MET", "ASN", "ARG",
            "ASP", "ASP", "ALA", "LEU", "ARG", "VAL", "THR", "GLY",
            "GLU", "VAL", "ARG", "GLU", "TYR", "ILE", "ASP", "THR",
            "VAL", "LYS", "ILE", "GLY", "TYR", "PRO", "LEU", "VAL",
            "LEU", "SER", "GLU", "GLY", "MET", "ASP", "ILE", "ILE",
            "ALA", "GLU", "PHE", "ARG", "LYS", "ARG", "PHE", "GLY",
            "CYS", "ARG", "ILE", "ILE", "ALA", "ASP", "PHE", "LYS",
            "VAL", "ALA", "ASP", "ILE", "PRO", "GLU", "THR", "ASN",
            "GLU", "LYS", "ILE", "CYS", "ARG", "ALA", "THR", "PHE",
            "LYS", "ALA", "GLY", "ALA", "ASP", "ALA", "ILE", "ILE",
            "VAL", "HIS", "GLY", "PHE", "PRO", "GLY", "ALA", "ASP",
            "SER", "VAL", "ARG", "ALA", "CYS", "LEU", "ASN", "VAL",
            "ALA", "GLU", "GLU", "MET", "GLY", "ARG", "GLU", "VAL",
            "PHE", "LEU", "LEU", "THR", "GLU", "MET", "SER", "HIS",
            "PRO", "GLY", "ALA", "GLU", "MET", "PHE", "ILE", "GLN",
            "GLY", "ALA", "ALA", "ASP", "GLU", "ILE", "ALA", "ARG",
            "MET", "GLY", "VAL", "ASP", "LEU", "GLY", "VAL", "LYS",
            "ASN", "TYR", "VAL", "GLY", "PRO", "SER", "THR", "ARG",
            "PRO", "GLU", "ARG", "LEU", "SER", "ARG", "LEU", "ARG",
            "GLU", "ILE", "ILE", "GLY", "GLN", "ASP", "SER", "PHE",
            "LEU", "ILE", "SER", "PRO", "GLY", "VAL", "GLY", "ALA",
            "GLN", "GLY", "GLY", "ASP", "PRO", "GLY", "GLU", "THR",
            "LEU", "ARG", "PHE", "ALA", "ASP", "ALA", "ILE", "ILE",
            "VAL", "GLY", "ARG", "SER", "ILE", "TYR", "LEU", "ALA",
            "ASP", "ASN", "PRO", "ALA", "ALA", "ALA", "ALA", "ALA",
            "GLY", "ILE", "ILE", "GLU", "SER", "ILE", "LYS", "ASP",
            "LEU", "LEU", "ILE", "PRO", "GLU"
           ]
          }, {
           "chain_id": "B",
           "length": 229,
           "residues": [
            "LEU", "ARG", "SER", "ARG", "ARG", "VAL", "ASP", "VAL",
            "MET", "ASP", "VAL", "MET", "ASN", "ARG", "LEU", "ILE",
            "LEU", "ALA", "MET", "ASP", "LEU", "MET", "ASN", "ARG",
            "ASP", "ASP", "ALA", "LEU", "ARG", "VAL", "THR", "GLY",
            "GLU", "VAL", "ARG", "GLU", "TYR", "ILE", "ASP", "THR",
            "VAL", "LYS", "ILE", "GLY", "TYR", "PRO", "LEU", "VAL",
            "LEU", "SER", "GLU", "GLY", "MET", "ASP", "ILE", "ILE",
            "ALA", "GLU", "PHE", "ARG", "LYS", "ARG", "PHE", "GLY",
            "CYS", "ARG", "ILE", "ILE", "ALA", "ASP", "PHE", "LYS",
            "VAL", "ALA", "ASP", "ILE", "PRO", "GLU", "THR", "ASN",
            "GLU", "LYS", "ILE", "CYS", "ARG", "ALA", "THR", "PHE",
            "LYS", "ALA", "GLY", "ALA", "ASP", "ALA", "ILE", "ILE",
            "VAL", "HIS", "GLY", "PHE", "PRO", "GLY", "ALA", "ASP",
            "SER", "VAL", "ARG", "ALA", "CYS", "LEU", "ASN", "VAL",
            "ALA", "GLU", "GLU", "MET", "GLY", "ARG", "GLU", "VAL",
            "PHE", "LEU", "LEU", "THR", "GLU", "MET", "SER", "HIS",
            "PRO", "GLY", "ALA", "GLU", "MET", "PHE", "ILE", "GLN",
            "GLY", "ALA", "ALA", "ASP", "GLU", "ILE", "ALA", "ARG",
            "MET", "GLY", "VAL", "ASP", "LEU", "GLY", "VAL", "LYS",
            "ASN", "TYR", "VAL", "GLY", "PRO", "SER", "THR", "ARG",
            "PRO", "GLU", "ARG", "LEU", "SER", "ARG", "LEU", "ARG",
            "GLU", "ILE", "ILE", "GLY", "GLN", "ASP", "SER", "PHE",
            "LEU", "ILE", "SER", "PRO", "GLY", "VAL", "GLY", "ALA",
            "GLN", "GLY", "GLY", "ASP", "PRO", "GLY", "GLU", "THR",
            "LEU", "ARG", "PHE", "ALA", "ASP", "ALA", "ILE", "ILE",
            "VAL", "GLY", "ARG", "SER", "ILE", "TYR", "LEU", "ALA",
            "ASP", "ASN", "PRO", "ALA", "ALA", "ALA", "ALA", "ALA",
            "GLY", "ILE", "ILE", "GLU", "SER", "ILE", "LYS", "ASP",
            "LEU", "LEU", "ILE", "PRO", "GLU"
           ]
          }
         ]
        )
        self.assertEqual(self.empty_data_file.residue_sequences, [])


    def test_modres(self):
        self.assertEqual(
         self.data_file.modifies_residues,
         [
          {
           "residue_name": "ASP",
           "chain_id": "A",
           "residue_id": 10,
           "insert_code": None,
           "standard_resisdue_name": 'ASP',
           "comment": "GLYCOSYLATION SITE"
          }
         ]
        )
        self.assertEqual(self.empty_data_file.modifies_residues, [])



class HeterogenSectionTest(PdbDataFileTest):

    def test_het(self):
        self.assertEqual(
         self.data_file.hets,
         [
          {
           "het_name": "BU2",
           "chain_id": "A",
           "het_id": 5001,
           "insert_code": None,
           "atom_num": 6,
           "description": None
          }, {
           "het_name": "BU2",
           "chain_id": "B",
           "het_id": 5002,
           "insert_code": None,
           "atom_num": 6,
           "description": None
          }, {
           "het_name": "XMP",
           "chain_id": "A",
           "het_id": 2001,
           "insert_code": None,
           "atom_num": 24,
           "description": None
          }, {
           "het_name": "XMP",
           "chain_id": "B",
           "het_id": 2002,
           "insert_code": None,
           "atom_num": 24,
           "description": None
          }
         ]
        )
        self.assertEqual(self.empty_data_file.hets, [])


    def test_hetnam(self):
        self.assertEqual(
         self.data_file.het_names,
         {
          "BU2": "1,3-BUTANEDIOL",
          "XMP": "XANTHOSINE-5'-MONOPHOSPHATE"
         }
        )
        self.assertEqual(self.empty_data_file.het_names, {})


    def test_hetsyn(self):
        self.assertEqual(
         self.data_file.het_synonyms,
         {
          "BU2": ["BOOM BOOM BOMB", "WYRDSTUFF"],
          "XMP": ["5-MONOPHOSPHATE-9-BETA-D-RIBOFURANOSYL XANTHINE"]
         }
        )
        self.assertEqual(self.empty_data_file.het_synonyms, {})


    def test_formul(self):
        self.assertEqual(
         self.data_file.het_formulae,
         {
          "BU2": {"component_number": 3, "is_water": False, "formula": "2(C4 H10 O2)"},
          "XMP": {"component_number": 5, "is_water": False, "formula": "2(C10 H14 N4 O9 P 1+)"},
          "HOH": {"component_number": 7, "is_water": True, "formula": "180(H2 O)"}
         }
        )
        self.assertEqual(self.empty_data_file.het_formulae, {})



class SecondaryStructureSectionTest(PdbDataFileTest):

    def test_helix(self):
        self.assertEqual(
         self.data_file.helices,
         [
          {
           "helix_id": 1,
           "helix_name": "1",
           "start_residue_name": "VAL",
           "start_residue_chain_id": "A",
           "start_residue_id": 11,
           "start_residue_insert": None,
           "end_residue_name": "ASN",
           "end_residue_chain_id": "A",
           "end_residue_id": 13,
           "end_residue_insert": None,
           "helix_class": 5,
           "comment": None,
           "length": 3
          }, {
           "helix_id": 2,
           "helix_name": "2",
           "start_residue_name": "ASN",
           "start_residue_chain_id": "A",
           "start_residue_id": 23,
           "start_residue_insert": None,
           "end_residue_name": "ARG",
           "end_residue_chain_id": "A",
           "end_residue_id": 35,
           "end_residue_insert": None,
           "helix_class": 1,
           "comment": None,
           "length": 13
          }
         ]
        )
        self.assertEqual(self.empty_data_file.helices, [])


    def test_sheet(self):
        self.assertEqual(
         self.data_file.sheets,
         [
          {
           "sheet_id": "A",
           "strand_count": 2,
           "strands": [{
            "strand_id": 1,
            "start_residue_name": "LEU",
            "start_residue_chain_id": "A",
            "start_residue_id": 15,
            "start_residue_insert": None,
            "end_residue_name": "MET",
            "end_residue_chain_id": "A",
            "end_residue_id": 19,
            "end_residue_insert": None,
            "sense": 0,
            "current_atom": None,
            "current_residue_name": None,
            "current_chain_id": None,
            "current_residue_id": None,
            "current_insert": None,
            "previous_atom": None,
            "previous_residue_name": None,
            "previous_chain_id": None,
            "previous_residue_id": None,
            "previous_insert": None
           }, {
            "strand_id": 2,
            "start_residue_name": "THR",
            "start_residue_chain_id": "A",
            "start_residue_id": 40,
            "start_residue_insert": None,
            "end_residue_name": "GLY",
            "end_residue_chain_id": "A",
            "end_residue_id": 44,
            "end_residue_insert": None,
            "sense": 1,
            "current_atom": "O",
            "current_residue_name": "LYS",
            "current_chain_id": "A",
            "current_residue_id": 42,
            "current_insert": None,
            "previous_atom": "N",
            "previous_residue_name": "LEU",
            "previous_chain_id": "A",
            "previous_residue_id": 17,
            "previous_insert": None
           }]
          }
         ]
        )
        self.assertEqual(self.empty_data_file.sheets, [])



class ConnectivityAnnotationSectionTest(PdbDataFileTest):

    def test_ssbond(self):
        self.assertEqual(
         self.data_file.ss_bonds,
         [{
          "serial_num": 1,
          "residue_name_1": "CYS",
          "chain_id_1": "A",
          "residue_id_1": 123,
          "insert_code_1": None,
          "residue_name_2": "CYS",
          "chain_id_2": "A",
          "residue_id_2": 155,
          "insert_code_2": None,
          "symmetry_1": "1555",
          "symmetry_2": "1555",
          "length": 2.04
         }]
        )
        self.assertEqual(self.empty_data_file.ss_bonds, [])


    def test_link(self):
        self.assertEqual(
         self.data_file.links,
         [{
          "atom_1": "O",
          "alt_loc_1": None,
          "residue_name_1": "TYR",
          "chain_id_1": "A",
          "residue_id_1": 146,
          "insert_code_1": None,
          "atom_2": "K",
          "alt_loc_2": None,
          "residue_name_2": "K",
          "chain_id_2": "A",
          "residue_id_2": 501,
          "insert_code_2": None,
          "symmetry_1": "1555",
          "symmetry_2": "1555",
          "length": 2.75
         }, {
          "atom_1": "OG1",
          "alt_loc_1": None,
          "residue_name_1": "THR",
          "chain_id_1": "A",
          "residue_id_1": 143,
          "insert_code_1": None,
          "atom_2": "K",
          "alt_loc_2": None,
          "residue_name_2": "K",
          "chain_id_2": "A",
          "residue_id_2": 504,
          "insert_code_2": None,
          "symmetry_1": "1555",
          "symmetry_2": "1555",
          "length": None
         }]
        )
        self.assertEqual(self.empty_data_file.links, [])


    def test_cispep(self):
        self.assertEqual(
         self.data_file.cis_peptides,
         [{
          "serial_num": None,
          "residue_name_1": "ASP",
          "chain_id_1": "B",
          "residue_id_1": 1188,
          "insert_1": None,
          "residue_name_2": "PRO",
          "chain_id_2": "B",
          "residue_id_2": 1189,
          "insert_2": None,
          "model_number": 0,
          "angle": 0.35
         }]
        )



class MiscellaneousSectionTest(PdbDataFileTest):

    def test_site(self):
        self.assertEqual(
         self.data_file.sites,
         [
          {
           "site_id": "AC1",
           "residue_count": 6,
           "residues": [
            {"residue_name": "ASP", "chain": "A", "residue_id": 70, "insert_code": None},
            {"residue_name": "LYS", "chain": "A", "residue_id": 72, "insert_code": None},
            {"residue_name": "LEU", "chain": "A", "residue_id": 123, "insert_code": None},
            {"residue_name": "VAL", "chain": "A", "residue_id": 155, "insert_code": None},
            {"residue_name": "XMP", "chain": "A", "residue_id": 2001, "insert_code": None},
            {"residue_name": "HOH", "chain": "A", "residue_id": 3015, "insert_code": None}
           ]
          }, {
           "site_id": "AC3",
           "residue_count": 8,
           "residues": [
            {"residue_name": "ALA", "chain": "A", "residue_id": 18, "insert_code": None},
            {"residue_name": "ASP", "chain": "A", "residue_id": 20, "insert_code": None},
            {"residue_name": "LYS", "chain": "A", "residue_id": 42, "insert_code": None},
            {"residue_name": "ASP", "chain": "A", "residue_id": 70, "insert_code": None},
            {"residue_name": "MET", "chain": "A", "residue_id": 126, "insert_code": None},
            {"residue_name": "SER", "chain": "A", "residue_id": 127, "insert_code": None},
            {"residue_name": "SER", "chain": "A", "residue_id": 158, "insert_code": None},
            {"residue_name": "PRO", "chain": "A", "residue_id": 180, "insert_code": None}
           ]
          }
         ]
        )
        self.assertEqual(self.empty_data_file.sites, [])



class CrystalSectionTest(PdbDataFileTest):

    def test_crystal(self):
        self.assertEqual(self.data_file.crystal_a, 57.57)
        self.assertEqual(self.data_file.crystal_b, 55.482)
        self.assertEqual(self.data_file.crystal_c, 66.129)
        self.assertEqual(self.data_file.crystal_alpha, 90.0)
        self.assertEqual(self.data_file.crystal_beta, 94.28)
        self.assertEqual(self.data_file.crystal_gamma, 90.0)
        self.assertEqual(self.data_file.crystal_s_group, "P 1 21 1")
        self.assertEqual(self.data_file.crystal_z, 4)
        self.assertEqual(self.empty_data_file.crystal_a, None)
        self.assertEqual(self.empty_data_file.crystal_b, None)
        self.assertEqual(self.empty_data_file.crystal_c, None)
        self.assertEqual(self.empty_data_file.crystal_alpha, None)
        self.assertEqual(self.empty_data_file.crystal_beta, None)
        self.assertEqual(self.empty_data_file.crystal_gamma, None)
        self.assertEqual(self.empty_data_file.crystal_s_group, None)
        self.assertEqual(self.empty_data_file.crystal_z, None)


    def test_origx(self):
        self.assertEqual(self.data_file.crystal_o11, 0.963457)
        self.assertEqual(self.data_file.crystal_o12, 0.136613)
        self.assertEqual(self.data_file.crystal_o13, 0.230424)
        self.assertEqual(self.data_file.crystal_t1, 16.61)
        self.assertEqual(self.data_file.crystal_o21, -0.158977)
        self.assertEqual(self.data_file.crystal_o22, 0.983924)
        self.assertEqual(self.data_file.crystal_o23, 0.081383)
        self.assertEqual(self.data_file.crystal_t2, 13.72)
        self.assertEqual(self.data_file.crystal_o31, -0.215598)
        self.assertEqual(self.data_file.crystal_o32, -0.115048)
        self.assertEqual(self.data_file.crystal_o33, 0.969683)
        self.assertEqual(self.data_file.crystal_t3, 37.65)
        self.assertEqual(self.empty_data_file.crystal_o11, None)
        self.assertEqual(self.empty_data_file.crystal_o12, None)
        self.assertEqual(self.empty_data_file.crystal_o13, None)
        self.assertEqual(self.empty_data_file.crystal_t1, None)
        self.assertEqual(self.empty_data_file.crystal_o21, None)
        self.assertEqual(self.empty_data_file.crystal_o22, None)
        self.assertEqual(self.empty_data_file.crystal_o23, None)
        self.assertEqual(self.empty_data_file.crystal_t2, None)
        self.assertEqual(self.empty_data_file.crystal_o31, None)
        self.assertEqual(self.empty_data_file.crystal_o32, None)
        self.assertEqual(self.empty_data_file.crystal_o33, None)
        self.assertEqual(self.empty_data_file.crystal_t3, None)


    def test_scale(self):
        self.assertEqual(self.data_file.crystal_s11, 0.01737)
        self.assertEqual(self.data_file.crystal_s12, 0.0)
        self.assertEqual(self.data_file.crystal_s13, 0.001301)
        self.assertEqual(self.data_file.crystal_u1, 0.0)
        self.assertEqual(self.data_file.crystal_s21, 0.0)
        self.assertEqual(self.data_file.crystal_s22, 0.018024)
        self.assertEqual(self.data_file.crystal_s23, 0.0)
        self.assertEqual(self.data_file.crystal_u2, 0.0)
        self.assertEqual(self.data_file.crystal_s31, 0.0)
        self.assertEqual(self.data_file.crystal_s32, 0.0)
        self.assertEqual(self.data_file.crystal_s33, 0.015164)
        self.assertEqual(self.data_file.crystal_u3, 0.0)
        self.assertEqual(self.empty_data_file.crystal_s11, None)
        self.assertEqual(self.empty_data_file.crystal_s12, None)
        self.assertEqual(self.empty_data_file.crystal_s13, None)
        self.assertEqual(self.empty_data_file.crystal_u1, None)
        self.assertEqual(self.empty_data_file.crystal_s21, None)
        self.assertEqual(self.empty_data_file.crystal_s22, None)
        self.assertEqual(self.empty_data_file.crystal_s23, None)
        self.assertEqual(self.empty_data_file.crystal_u2, None)
        self.assertEqual(self.empty_data_file.crystal_s31, None)
        self.assertEqual(self.empty_data_file.crystal_s32, None)
        self.assertEqual(self.empty_data_file.crystal_s33, None)
        self.assertEqual(self.empty_data_file.crystal_u3, None)


    def test_mtrix(self):
        self.assertEqual(self.data_file.crystal_serial_1, 1)
        self.assertEqual(self.data_file.crystal_m11, -1.0)
        self.assertEqual(self.data_file.crystal_m12, 0.0)
        self.assertEqual(self.data_file.crystal_m13, 0.0)
        self.assertEqual(self.data_file.crystal_v1, 0.0)
        self.assertEqual(self.data_file.crystal_i_given_1, True)
        self.assertEqual(self.data_file.crystal_serial_2, 1)
        self.assertEqual(self.data_file.crystal_m21, 0.0)
        self.assertEqual(self.data_file.crystal_m22, 1.0)
        self.assertEqual(self.data_file.crystal_m23, 0.0)
        self.assertEqual(self.data_file.crystal_v2, 0.0)
        self.assertEqual(self.data_file.crystal_i_given_3, True)
        self.assertEqual(self.data_file.crystal_serial_3, 1)
        self.assertEqual(self.data_file.crystal_m31, 0.0)
        self.assertEqual(self.data_file.crystal_m32, 0.0)
        self.assertEqual(self.data_file.crystal_m33, -1.0)
        self.assertEqual(self.data_file.crystal_v3, 0.0)
        self.assertEqual(self.data_file.crystal_i_given_3, True)
        self.assertEqual(self.empty_data_file.crystal_serial_1, None)
        self.assertEqual(self.empty_data_file.crystal_m11, None)
        self.assertEqual(self.empty_data_file.crystal_m12, None)
        self.assertEqual(self.empty_data_file.crystal_m13, None)
        self.assertEqual(self.empty_data_file.crystal_v1, None)
        self.assertEqual(self.empty_data_file.crystal_i_given_3, False)
        self.assertEqual(self.empty_data_file.crystal_serial_1, None)
        self.assertEqual(self.empty_data_file.crystal_m21, None)
        self.assertEqual(self.empty_data_file.crystal_m22, None)
        self.assertEqual(self.empty_data_file.crystal_m23, None)
        self.assertEqual(self.empty_data_file.crystal_v2, None)
        self.assertEqual(self.empty_data_file.crystal_i_given_3, False)
        self.assertEqual(self.empty_data_file.crystal_serial_1, None)
        self.assertEqual(self.empty_data_file.crystal_m31, None)
        self.assertEqual(self.empty_data_file.crystal_m32, None)
        self.assertEqual(self.empty_data_file.crystal_m33, None)
        self.assertEqual(self.empty_data_file.crystal_v3, None)
        self.assertEqual(self.empty_data_file.crystal_i_given_3, False)



class CoordinateSectionTest(PdbDataFileTest):

    def test_model(self):
        self.assertEqual(
         self.data_file.models,
         [
          {
           "model_id": 1,
           "start_record": 462,
           "end_record": 469
          }, {
           "model_id": 2,
           "start_record": 470,
           "end_record": 477
          }
         ]
        )
        self.assertEqual(self.empty_data_file.models, [])


    def test_atom(self):
        self.assertEqual(
         self.data_file.atoms,
         [
          {
           "atom_id": 107,
           "atom_name": "N",
           "alt_loc": None,
           "residue_name": "GLY",
           "chain": "A",
           "residue_id": 13,
           "insert_code": None,
           "x": 12.681,
           "y": 37.302,
           "z": -25.211,
           "occupancy": 1.0,
           "temperature_factor": 15.56,
           "element": "N",
           "charge": None,
           "model_id": 1
          }, {
           "atom_id": 108,
           "atom_name": "CA",
           "alt_loc": None,
           "residue_name": "GLY",
           "chain": "A",
           "residue_id": 13,
           "insert_code": None,
           "x": 11.982,
           "y": 37.996,
           "z": -26.241,
           "occupancy": 1.0,
           "temperature_factor": 16.92,
           "element": "C",
           "charge": None,
           "model_id": 1
          }, {
           "atom_id": 107,
           "atom_name": "N",
           "alt_loc": None,
           "residue_name": "GLY",
           "chain": "A",
           "residue_id": 13,
           "insert_code": None,
           "x": 12.681,
           "y": 37.302,
           "z": -25.211,
           "occupancy": 1.0,
           "temperature_factor": 15.56,
           "element": "N",
           "charge": None,
           "model_id": 2
          }, {
           "atom_id": 108,
           "atom_name": "CA",
           "alt_loc": None,
           "residue_name": "GLY",
           "chain": "A",
           "residue_id": 13,
           "insert_code": None,
           "x": 11.982,
           "y": 37.996,
           "z": -26.241,
           "occupancy": 1.0,
           "temperature_factor": 16.92,
           "element": "C",
           "charge": None,
           "model_id": 2
          }
         ]
        )
        self.assertEqual(self.empty_data_file.atoms, [])


    def test_anisou(self):
        self.assertEqual(
         self.data_file.anisou,
         [
          {
           "atom_id": 107,
           "atom_name": "N",
           "alt_loc": None,
           "residue_name": "GLY",
           "chain": "A",
           "residue_id": 13,
           "insert_code": None,
           "u11": 2406,
           "u22": 1892,
           "u33": 1614,
           "u12": 198,
           "u13": 519,
           "u23": -328,
           "element": "N",
           "charge": None,
           "model_id": 1
          }, {
           "atom_id": 108,
           "atom_name": "CA",
           "alt_loc": None,
           "residue_name": "GLY",
           "chain": "A",
           "residue_id": 13,
           "insert_code": None,
           "u11": 2748,
           "u22": 2004,
           "u33": 1679,
           "u12": -21,
           "u13": 155,
           "u23": -419,
           "element": "C",
           "charge": None,
           "model_id": 1
          }, {
           "atom_id": 107,
           "atom_name": "N",
           "alt_loc": None,
           "residue_name": "GLY",
           "chain": "A",
           "residue_id": 13,
           "insert_code": None,
           "u11": 2406,
           "u22": 1892,
           "u33": 1614,
           "u12": 198,
           "u13": 519,
           "u23": -328,
           "element": "N",
           "charge": None,
           "model_id": 2
          }, {
           "atom_id": 108,
           "atom_name": "CA",
           "alt_loc": None,
           "residue_name": "GLY",
           "chain": "A",
           "residue_id": 13,
           "insert_code": None,
           "u11": 2748,
           "u22": 2004,
           "u33": 1679,
           "u12": -21,
           "u13": 155,
           "u23": -419,
           "element": "C",
           "charge": None,
           "model_id": 2
          }
         ]
        )
        self.assertEqual(self.empty_data_file.anisou, [])


    def test_ter(self):
        self.assertEqual(
         self.data_file.terminals,
         [
          {
           "atom_id": 109,
           "residue_name": "GLY",
           "chain": "A",
           "residue_id": 13,
           "insert_code": None,
           "model_id": 1
          }, {
           "atom_id": 109,
           "residue_name": "GLY",
           "chain": "A",
           "residue_id": 13,
           "insert_code": None,
           "model_id": 2
          }
         ]
        )
        self.assertEqual(self.empty_data_file.terminals, [])


    def test_hetatm(self):
        self.assertEqual(
         self.data_file.heteroatoms,
         [
          {
           "atom_id": 8237,
           "atom_name": "MG",
           "alt_loc": None,
           "residue_name": "MG",
           "chain": "A",
           "residue_id": 1001,
           "insert_code": None,
           "x": 13.872,
           "y": -2.555,
           "z": -29.045,
           "occupancy": 1.0,
           "temperature_factor": 27.36,
           "element": "MG",
           "charge": None,
           "model_id": 1
          }, {
           "atom_id": 8237,
           "atom_name": "MG",
           "alt_loc": None,
           "residue_name": "MG",
           "chain": "A",
           "residue_id": 1001,
           "insert_code": None,
           "x": 13.872,
           "y": -2.555,
           "z": -29.045,
           "occupancy": 1.0,
           "temperature_factor": 27.36,
           "element": "MG",
           "charge": None,
           "model_id": 2
          }
         ]
        )
        self.assertEqual(self.empty_data_file.heteroatoms, [])
