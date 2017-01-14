import datetime
from unittest import TestCase
from unittest.mock import patch
from molecupy.pdb.pdbfile import PdbFile
from molecupy.pdb.pdbdatafile import PdbDataFile
from molecupy.converters.pdbfile2pdbdatafile import pdb_data_file_from_pdb_file

class PdbDataFileTest(TestCase):

    def setUp(self):
        self.blank = PdbDataFile()



class PdbDataFileCreationTests(PdbDataFileTest):

    def test_can_create_empty_data_file(self):
        data_file = PdbDataFile()
        self.assertEqual(data_file._source, None)

        self.assertEqual(data_file._classification, None)
        self.assertEqual(data_file._deposition_date, None)
        self.assertEqual(data_file._pdb_code, None)
        self.assertFalse(data_file._is_obsolete)
        self.assertEqual(data_file._obsolete_date, None)
        self.assertEqual(data_file._replacement_code, None)
        self.assertEqual(data_file._title, None)
        self.assertEqual(data_file._split_codes, [])
        self.assertEqual(data_file._caveat, None)
        self.assertEqual(data_file._compounds, [])
        self.assertEqual(data_file._sources, [])
        self.assertEqual(data_file._keywords, [])
        self.assertEqual(data_file._experimental_techniques, [])
        self.assertEqual(data_file._model_count, 1)
        self.assertEqual(data_file._model_annotations, [])
        self.assertEqual(data_file._authors, [])
        self.assertEqual(data_file._revisions, [])
        self.assertEqual(data_file._supercedes, [])
        self.assertEqual(data_file._supercede_date, None)
        self.assertEqual(data_file._journal, None)
        self.assertEqual(data_file._remarks, [])

        self.assertEqual(data_file._dbreferences, [])
        self.assertEqual(data_file._sequence_differences, [])
        self.assertEqual(data_file._residue_sequences, [])
        self.assertEqual(data_file._modified_residues, [])

        self.assertEqual(data_file._hets, [])
        self.assertEqual(data_file._het_names, {})
        self.assertEqual(data_file._het_synonyms, {})
        self.assertEqual(data_file._formulae, {})

        self.assertEqual(data_file._helices, [])
        self.assertEqual(data_file._sheets, [])

        self.assertEqual(data_file._ss_bonds, [])
        self.assertEqual(data_file._links, [])
        self.assertEqual(data_file._cis_peptides, [])

        self.assertEqual(data_file._sites, [])

        self.assertEqual(data_file._crystal, None)
        self.assertEqual(data_file._origix, None)
        self.assertEqual(data_file._scale, None)
        self.assertEqual(data_file._matrix, None)

        self.assertEqual(
         data_file._models,
         [{"model_id": 1, "start_record": -1, "end_record": -1}]
        )
        self.assertEqual(data_file._atoms, [])
        self.assertEqual(data_file._anisou, [])
        self.assertEqual(data_file._termini, [])
        self.assertEqual(data_file._heteroatoms, [])

        self.assertEqual(data_file._connections, [])

        self.assertEqual(data_file._master, None)


    def test_repr(self):
        data_file = PdbDataFile()
        self.assertEqual(str(data_file), "<PdbDataFile (????)>")


    def test_repr_when_there_is_pdb_code(self):
        data_file = PdbDataFile()
        data_file._pdb_code = "ABCD"
        self.assertEqual(str(data_file), "<PdbDataFile (ABCD)>")



class GeneralPdbDataFilePropertyTests(PdbDataFileTest):

    def test_basic_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile())
        self.assertIs(data_file.source(), data_file._source)



class PdbDataFileConversionTests(PdbDataFileTest):

    @patch("molecupy.converters.pdbdatafile2pdbfile.pdb_file_from_pdb_data_file")
    def test_can_convert_to_pdb_file(self, mock_converter):
        value = "Return value"
        mock_converter.return_value = value
        self.assertIs(PdbDataFile().to_pdb_file(), value)



class TitleSectionPropertyTests(PdbDataFileTest):

    def test_header_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "HEADER    LYASE                                   06-MAY-02   1LOL"
        ))
        self.assertIs(data_file._classification, data_file.classification())
        self.assertIs(data_file._deposition_date, data_file.deposition_date())
        self.assertIs(data_file._pdb_code, data_file.pdb_code())


    def test_can_modify_header_properties(self):
        self.blank.classification("TEST CLASS")
        self.blank.deposition_date(datetime.datetime(2008, 1, 24).date())
        self.blank.pdb_code("1xxx")
        self.assertEqual(self.blank._classification, "TEST CLASS")
        self.assertEqual(
         self.blank._deposition_date,
         datetime.datetime(2008, 1, 24).date()
        )
        self.assertEqual(self.blank._pdb_code, "1xxx")


    def test_classification_must_be_str(self):
        with self.assertRaises(TypeError):
            self.blank.classification(1000)


    def test_classifcation_must_be_less_than_40_chars(self):
        with self.assertRaises(ValueError):
            self.blank.classification("-" * 41)


    def test_deposition_date_must_be_date(self):
        with self.assertRaises(TypeError):
            self.blank.deposition_date("1-1-91")


    def test_pdb_code_must_be_string(self):
        with self.assertRaises(TypeError):
            self.blank.pdb_code(1000)


    def test_pdb_code_must_be_4_chars(self):
        with self.assertRaises(ValueError):
            self.blank.pdb_code("1xx")
        with self.assertRaises(ValueError):
            self.blank.pdb_code("1xxxx")


    def test_obslte_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "OBSLTE     30-SEP-93 1LOL      1SAM"
        ))
        self.assertIs(data_file._is_obsolete, data_file.is_obsolete())
        self.assertIs(data_file._obsolete_date, data_file.obsolete_date())
        self.assertIs(data_file._replacement_code, data_file.replacement_code())


    def test_can_modify_obslte_properties(self):
        self.blank.is_obsolete(True)
        self.blank.obsolete_date(datetime.datetime(2008, 1, 24).date())
        self.blank.replacement_code("1xxx")
        self.assertTrue(self.blank._is_obsolete)
        self.assertEqual(
         self.blank._obsolete_date,
         datetime.datetime(2008, 1, 24).date()
        )
        self.assertEqual(self.blank._replacement_code, "1xxx")


    def test_is_obsolete_must_be_bool(self):
        with self.assertRaises(TypeError):
            self.blank.is_obsolete("yes")


    def test_obsolete_date_must_be_date(self):
        with self.assertRaises(TypeError):
            self.blank.obsolete_date("1-1-91")


    def test_replacement_code_must_be_string(self):
        with self.assertRaises(TypeError):
            self.blank.replacement_code(1000)


    def test_replacement_code_must_be_4_chars(self):
        with self.assertRaises(ValueError):
            self.blank.replacement_code("1xx")
        with self.assertRaises(ValueError):
            self.blank.replacement_code("1xxxx")


    def test_title_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "TITLE     CRYSTAL STRUCTURE OF OROTIDINE MONOPHOSPHATE DECARBOXYLASE\n"
         "TITLE    2 COMPLEX WITH XMP"
        ))
        self.assertIs(data_file._title, data_file.title())


    def test_can_modify_title_properties(self):
        self.blank.title("123" * 10)
        self.assertEqual(self.blank.title(), "123" * 10)


    def test_title_must_be_str(self):
        with self.assertRaises(TypeError):
            self.blank.title(100)


    def test_split_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "SPLIT      1VOQ 1VOR 1VOS 1VOU 1VOV 1VOW 1VOX 1VOY 1VP0 1VOZ 1VOY 1VP0 1VOZ 1VOZ\n"
         "SPLIT      1VOQ 1VOR 1VOS 1VOU 1VOV 1VOW 1VOX 1VOY 1VP0 1VOZ"
        ))
        self.assertIs(data_file._split_codes, data_file.split_codes())


    def test_caveat_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "TITLE     CRYSTAL STRUCTURE OF OROTIDINE MONOPHOSPHATE DECARBOXYLASE\n"
         "TITLE    2 COMPLEX WITH XMP"
        ))
        self.assertIs(data_file._caveat, data_file.caveat())


    def test_can_modify_caveat_properties(self):
        self.blank.caveat("123" * 10)
        self.assertEqual(self.blank.caveat(), "123" * 10)


    def test_caveat_must_be_str(self):
        with self.assertRaises(TypeError):
            self.blank.caveat(100)


    def test_compnd_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "COMPND    MOL_ID: 1;\n"
         "COMPND   2 MOLECULE: OROTIDINE 5'-MONOPHOSPHATE DECARBOXYLASE;\n"
        ))
        self.assertIs(data_file._compounds, data_file.compounds())


    def test_source_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "SOURCE    MOL_ID: 1;\n"
         "SOURCE   2 ORGANISM_SCIENTIFIC: METHANOTHERMOBACTER\n"
        ))
        self.assertIs(data_file._sources, data_file.sources())


    def test_keyword_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "SOURCE    MOL_ID: 1;\n"
         "SOURCE   2 ORGANISM_SCIENTIFIC: METHANOTHERMOBACTER\n"
        ))
        self.assertIs(data_file._keywords, data_file.keywords())


    def test_expdta_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "SOURCE    MOL_ID: 1;\n"
         "SOURCE   2 ORGANISM_SCIENTIFIC: METHANOTHERMOBACTER\n"
        ))
        self.assertIs(
         data_file._experimental_techniques,
         data_file.experimental_techniques()
        )


    def test_nummdl_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "NUMMDL    2"
        ))
        self.assertIs(data_file._model_count, data_file.model_count())


    def test_can_modify_nummdl_properties(self):
        self.blank.model_count(9)
        self.assertEqual(self.blank.model_count(), 9)


    def test_nummdl_must_be_int(self):
        with self.assertRaises(TypeError):
            self.blank.model_count("1")
        with self.assertRaises(TypeError):
            self.blank.model_count(9.8)


    def test_mdltyp_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "MDLTYP    CA ATOMS ONLY, CHAIN A, B, C, D, E, F, G, H, I, J, K ; P ATOMS ONLY,\n"
         "MDLTYP   2 CHAIN X, Y, Z"
        ))
        self.assertIs(
         data_file._model_annotations,
         data_file.model_annotations()
        )


    def test_author_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "AUTHOR    M.B.BERRY,B.MEADOR,T.BILDERBACK,P.LIANG,M.GLASER,\n"
         "AUTHOR   2 G.N.PHILLIPS JR.,T.L.ST. STEVENS"
        ))
        self.assertIs(data_file._authors, data_file.authors())


    def test_revdat_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "REVDAT   4 1 24-FEB-09 1LOL    1       VERSN  COMPND EXPDTA CAVEAT\n"
         "REVDAT   4 2                   1       SOURCE JRNL\n"
        ))
        self.assertIs(data_file._revisions, data_file.revisions())


    def test_sprsde_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "SPRSDE     27-FEB-95 1GDJ      1LH4 2LH4"
        ))
        self.assertIs(data_file._supercedes, data_file.supercedes())
        self.assertIs(data_file._supercede_date, data_file.supercede_date())


    def test_can_modify_sprsde_properties(self):
        self.blank.supercede_date(datetime.datetime(2008, 1, 24).date())
        self.assertEqual(
         self.blank.supercede_date(),
         datetime.datetime(2008, 1, 24).date()
        )


    def test_supercede_date_must_be_date(self):
        with self.assertRaises(TypeError):
            self.blank.supercede_date("1-1-91")


    def test_jrnl_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "JRNL        REF    J.BIOL.CHEM.                  V. 277 28080 2002"
        ))
        self.assertEqual(data_file._journal, data_file.journal())


    def test_can_modify_jrnl_properties(self):
        self.blank.journal({"authors": ["N.WU", "E.F.PAI"]})
        self.assertEqual(self.blank.journal(), {"authors": ["N.WU", "E.F.PAI"]})


    def test_journal_must_be_dict(self):
        with self.assertRaises(TypeError):
            self.blank.journal("aaa")


    def test_remark_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "REMARK   2\n"
         "REMARK   2 RESOLUTION.    1.90 ANGSTROMS."
        ))
        self.assertEqual(data_file._remarks, data_file.remarks())


    def test_can_get_remark_by_number(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "REMARK   2\n"
         "REMARK 999\n"
         "REMARK   2 RESOLUTION.    1.90 ANGSTROMS."
        ))
        self.assertEqual(
         data_file.get_remark_by_number(2),
         {
          "number": 2,
          "content": "RESOLUTION.    1.90 ANGSTROMS."
         }
        )
        self.assertEqual(data_file.get_remark_by_number(3), None)



class PrimaryStructureSectionPropertyTests(PdbDataFileTest):

    def test_dbref_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "DBREF1 1LOL C   61   322  GB                   AE017221\n"
         "DBREF2 1LOL C     46197919                      1534489     1537377"
        ))
        self.assertEqual(data_file._dbreferences, data_file.dbreferences())


    def test_seqadv_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "SEQADV 1LOL GLU A  229  UNP  O26232              INSERTION\n"
         "SEQADV 1LOL GLU B 1229  UNP  O26232              INSERTION"
        ))
        self.assertEqual(
         data_file._sequence_differences,
         data_file.sequence_differences()
        )


    def test_seqres_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "SEQRES   1 A    8  LEU ARG SER ARG ARG VAL ASP VAL MET ASP VAL MET ASN\n"
         "SEQRES   2 A    8  ARG LEU ILE\n"
        ))
        self.assertEqual(
         data_file._residue_sequences,
         data_file.residue_sequences()
        )



    def test_modres_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "MODRES 1LOL ASP A   10  ASP  GLYCOSYLATION SITE"
        ))
        self.assertEqual(
         data_file._modified_residues,
         data_file.modified_residues()
        )



class HeterogenSectionPropertyTests(PdbDataFileTest):

    def test_het_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "HET    BU2  A5001       6\n"
         "HET    BU2  B5002       6\n"
         "HET    XMP  A2001      24\n"
         "HET    XMP  B2002      24"
        ))
        self.assertEqual(data_file._hets, data_file.hets())


    def test_hetnam_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "HETNAM     BU2 1,3-BUTANEDIOL\n"
         "HETNAM     XMP XANTHOSINE-5'-MONOPHOSPHATE"
        ))
        self.assertEqual(data_file._het_names, data_file.het_names())


    def test_hetsyn_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "HETSYN     BU2 BOOM BOOM BOMB; WYRDSTUFF\n"
         "HETSYN     XMP 5-MONOPHOSPHATE-9-BETA-D-RIBOFURANOSYL XANTHINE"
        ))
        self.assertEqual(data_file._het_synonyms, data_file.het_synonyms())


    def test_formul_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "FORMUL   3  BU2    2(C4 H10 O2)\n"
        ))
        self.assertEqual(data_file._formulae, data_file.formulae())



class SecondaryStructurePropertyTests(PdbDataFileTest):

    def test_helix_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "HELIX    1   1 VAL A   11  ASN A   13  5                                   3\n"
         "HELIX    2   2 ASN A   23  ARG A   35  1                                  13"
        ))
        self.assertEqual(data_file._helices, data_file.helices())


    def test_sheet_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "SHEET    1   A 2 LEU A  15  MET A  19  0\n"
         "SHEET    2   A 2 THR A  40  GLY A  44  1  O  LYS A  42   N  LEU A  17"
        ))
        self.assertEqual(data_file._sheets, data_file.sheets())



class ConnectivityAnnotationPropertyTests(PdbDataFileTest):

    def test_ssbond_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "SSBOND   1 CYS A  123    CYS A  155                          1555   1555  2.04"
        ))
        self.assertEqual(data_file._ss_bonds, data_file.ss_bonds())


    def test_link_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "LINK         O   TYR A 146                 K     K A 501     1555   1555  2.75"
        ))
        self.assertEqual(data_file._links, data_file.links())


    def test_cispep_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "CISPEP     ASP B 1188    PRO B 1189          0         0.35"
        ))
        self.assertEqual(data_file._cis_peptides, data_file.cis_peptides())



class MiscellaneousSectionPropertyTests(PdbDataFileTest):

    def test_site_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "SITE     1 AC1  6 ASP A  70  LYS A  72  LEU A 123  VAL A 155"
        ))
        self.assertEqual(data_file._sites, data_file.sites())



class CrystallographySectionPropertyTests(PdbDataFileTest):

    def test_crystal_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "CRYST1   57.570   55.482   66.129  90.00  94.28  90.00 P 1 21 1      4"
        ))
        self.assertEqual(data_file._crystal, data_file.crystal())


    def test_can_modify_crystal_properties(self):
        self.blank.crystal({"a": 100.1})
        self.assertEqual(self.blank.crystal(), {"a": 100.1})


    def test_crystal_must_be_dict(self):
        with self.assertRaises(TypeError):
            self.blank.crystal("aaa")


    def test_origx_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "ORIGX1      0.963457  0.136613  0.230424       16.61000\n"
         "ORIGX2     -0.158977  0.983924  0.081383       13.72000\n"
         "ORIGX3     -0.215598 -0.115048  0.969683       37.65000"
        ))
        self.assertEqual(data_file._origx, data_file.origx())


    def test_can_modify_origix_properties(self):
        self.blank.origx({"o22": 100.1})
        self.assertEqual(self.blank.origx(), {"o22": 100.1})


    def test_origix_must_be_dict(self):
        with self.assertRaises(TypeError):
            self.blank.origx("aaa")


    def test_scale_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "SCALE1      0.963457  0.136613  0.230424       16.61000\n"
         "SCALE2     -0.158977  0.983924  0.081383       13.72000\n"
         "SCALE3     -0.215598 -0.115048  0.969683       37.65000"
        ))
        self.assertEqual(data_file._scale, data_file.scale())


    def test_can_modify_scale_properties(self):
        self.blank.scale({"s22": 100.1})
        self.assertEqual(self.blank.scale(), {"s22": 100.1})


    def test_scale_must_be_dict(self):
        with self.assertRaises(TypeError):
            self.blank.scale("aaa")


    def test_mtrix_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "MTRIX1   1 -1.000000  0.000000  0.000000        0.00000    1\n"
         "MTRIX2   1  0.000000  1.000000  0.000000        0.00000    1\n"
         "MTRIX3   1  0.000000  0.000000 -1.000000        0.00000    1"
        ))
        self.assertEqual(data_file._matrix, data_file.matrix())


    def test_can_modify_mtrx_properties(self):
        self.blank.matrix({"m22": 100.1})
        self.assertEqual(self.blank.matrix(), {"m22": 100.1})


    def test_scale_must_be_dict(self):
        with self.assertRaises(TypeError):
            self.blank.matrix("aaa")



class CoordinateSectionPropertyTests(PdbDataFileTest):

    def test_model_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "MODEL        1\n"
         "ATOM    107  N   GLY A  13      12.681  37.302 -25.211 1.000 15.56           N\n"
         "ENDMDL\n"
        ))
        self.assertEqual(data_file._models, data_file.models())


    def test_atom_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "ATOM    107  N   GLY A  13      12.681  37.302 -25.211 1.000 15.56           N\n"
         "ATOM    108  CA  GLY A  13      11.982  37.996 -26.241 1.000 16.92           C"
        ))
        self.assertEqual(data_file._atoms, data_file.atoms())


    def test_anisou_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "ANISOU  107  N   GLY A  13     2406   1892   1614    198    519   -328       N\n"
         "ANISOU  108  CA  GLY A  13     2748   2004   1679    -21    155   -419       C"
        ))
        self.assertEqual(data_file._anisou, data_file.anisou())


    def test_ter_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "TER     109      GLY A  13"
        ))
        self.assertEqual(data_file._termini, data_file.termini())


    def test_hetatm_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "HETATM 8237 MG    MG A1001      13.872  -2.555 -29.045  1.00 27.36          MG"
        ))
        self.assertEqual(data_file._heteroatoms, data_file.heteroatoms())



class ConnectionSectionPropertyTests(PdbDataFileTest):

    def test_conect_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "CONECT 1179  746 1184 1195 1203\n"
         "CONECT 1179 1211 1222"
        ))
        self.assertEqual(data_file._connections, data_file.connections())



class MasterSectionPropertyTests(PdbDataFileTest):

    def test_master_properties(self):
        data_file = pdb_data_file_from_pdb_file(PdbFile(
         "MASTER       40    0    0    0    0    0    0    6 2930    2    0   29"
        ))
        self.assertEqual(data_file._master, data_file.master())


    def test_can_modify_master_properties(self):
        self.blank.master({
         "remark_num": 40,
         "het_num": 0,
         "helix_num": 0,
         "sheet_num": 0,
         "site_num": 0,
         "crystal_num": 6,
         "coordinate_num": 2930,
         "ter_num": 2,
         "conect_num": 0,
         "seqres_num": 29
        })
        self.assertEqual(self.blank.master(), {
         "remark_num": 40,
         "het_num": 0,
         "helix_num": 0,
         "sheet_num": 0,
         "site_num": 0,
         "crystal_num": 6,
         "coordinate_num": 2930,
         "ter_num": 2,
         "conect_num": 0,
         "seqres_num": 29
        })


    def test_master_must_be_dict(self):
        with self.assertRaises(TypeError):
            self.blank.master("aaa")
