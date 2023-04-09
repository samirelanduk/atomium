import os
import shutil
from unittest import TestCase
import atomium

class FileToDictTests(TestCase):

    def test_1lol(self):
        d = atomium.open("tests/integration/files/1lol.pdb", dictionary=True)
        self.assertEqual(len(d.keys()), 40)

        # Metadata categories
        self.assertEqual(d["entry"], [{"id": "1LOL"}])
        self.assertEqual(d["pdbx_database_status"], [
            {"status_code": "REL", "entry_id": "1LOL", "recvd_initial_deposition_date": "2002-05-06"},
        ])
        self.assertEqual(d["struct_keywords"], [{"entry_id": "1LOL", "pdbx_keywords": "LYASE", "text": "TIM BARREL, LYASE"}])
        self.assertNotIn("pdbx_database_PDB_obs_spr", d)
        self.assertEqual(d["struct"], [{
            "entry_id": "1LOL",
            "title": "CRYSTAL STRUCTURE OF OROTIDINE MONOPHOSPHATE DECARBOXYLASE COMPLEX WITH XMP",
            "pdbx_descriptor": "?", "pdbx_model_details": "?",
            "pdbx_CASP_flag": "?", "pdbx_model_type_details": "?",
        }])
        self.assertNotIn("pdbx_database_related", d)
        self.assertNotIn("database_PDB_caveat", d)
        self.assertEqual(d["exptl"], [{"entry_id": "1LOL", "method": "X-RAY DIFFRACTION"}])
        self.assertNotIn("pdbx_coordinate_model", d)
        self.assertEqual(d["pdbx_audit_revision_history"], [{
            "ordinal": "1", "data_content_type": "Structure model",
            "major_revision": "1", "minor_revision": "0",
            "revision_date": "2002-08-07",
        }, {
            "ordinal": "2", "data_content_type": "Structure model",
            "major_revision": "2", "minor_revision": "0",
            "revision_date": "2002-08-14",
        }, {
            "ordinal": "3", "data_content_type": "Structure model",
            "major_revision": "3", "minor_revision": "0",
            "revision_date": "2003-04-01",
        }, {
            "ordinal": "4", "data_content_type": "Structure model",
            "major_revision": "4", "minor_revision": "0",
            "revision_date": "2009-02-24",
        }])

        # Entity categories
        self.assertEqual(d["entity"], [{
            "id": "1", "type": "polymer", "src_method": "man",
            "pdbx_description": "OROTIDINE 5'-MONOPHOSPHATE DECARBOXYLASE",
            "formula_weight": "?", "pdbx_number_of_molecules": "2",
            "pdbx_ec": "4.1.1.23", "pdbx_mutation": "?", "pdbx_fragment": "?",
            "details": "?",
        }, {
            "id": "2", "type": "non-polymer", "src_method": "syn",
            "pdbx_description": "1,3-BUTANEDIOL", "formula_weight": "?",
            "pdbx_number_of_molecules": "2", "pdbx_ec": "?", "pdbx_mutation": "?",
            "pdbx_fragment": "?", "details": "?",
        }, {
            "id": "3", "type": "non-polymer", "src_method": "syn",
            "pdbx_description": "XANTHOSINE-5'-MONOPHOSPHATE",
            "formula_weight": "?", "pdbx_number_of_molecules": "2", "pdbx_ec": "?",
            "pdbx_mutation": "?", "pdbx_fragment": "?", "details": "?",
        }, {
            "id": "4", "type": "water", "src_method": "nat",
            "pdbx_description": "water", "formula_weight": "?",
            "pdbx_number_of_molecules": "180", "pdbx_ec": "?", "pdbx_mutation": "?",
            "pdbx_fragment": "?", "details": "?",
        }])
        self.assertEqual(d["entity_name_com"], [{"entity_id": "1", "name": "OMP DECARBOXYLASE, OMPDCASE, OMPDECASE"}])
        self.assertEqual(d["entity_poly"], [{
            "entity_id": "1", "type": "polypeptide(L)", "nstd_linkage": "no",
            "nstd_monomer": "no",
            "pdbx_seq_one_letter_code": "LRSRRVDVMDVMNRLILAMDLMNRDDALRVTGEVREYIDTVKIGYPLVLSEGMDIIAEFRKRFGCRIIADFKVADIPETNEKICRATFKAGADAIIVHGFPGADSVRACLNVAEEMGREVFLLTEMSHPGAEMFIQGAADEIARMGVDLGVKNYVGPSTRPERLSRLREIIGQDSFLISPGVGAQGGDPGETLRFADAIIVGRSIYLADNPAAAAAGIIESIKDLLIPE",
            "pdbx_seq_one_letter_code_can": "LRSRRVDVMDVMNRLILAMDLMNRDDALRVTGEVREYIDTVKIGYPLVLSEGMDIIAEFRKRFGCRIIADFKVADIPETNEKICRATFKAGADAIIVHGFPGADSVRACLNVAEEMGREVFLLTEMSHPGAEMFIQGAADEIARMGVDLGVKNYVGPSTRPERLSRLREIIGQDSFLISPGVGAQGGDPGETLRFADAIIVGRSIYLADNPAAAAAGIIESIKDLLIPE",
            "pdbx_strand_id": "A,B", "pdbx_target_identifier": "?",
        }])
        self.assertEqual(len(d["entity_poly_seq"]), 229)
        self.assertEqual(d["entity_poly_seq"][0], {
            "entity_id": "1", "num": "1", "mon_id": "LEU", "hetero": "n",
        })
        self.assertEqual(d["entity_poly_seq"][-1], {
            "entity_id": "1", "num": "229", "mon_id": "GLU", "hetero": "n",
        })
        self.assertEqual(d["struct_ref"], [{
            "id": "1", "db_name": "UNP", "db_code": "PYRF_METTH",
            "pdbx_db_accession": "O26232", "pdbx_db_isoform": "?", "entity_id": "1",
            "pdbx_seq_one_letter_code": "?", "pdbx_align_begin": "?",
        }])
        self.assertEqual(d["struct_ref_seq"], [{
            "align_id": "1", "ref_id": "1", "pdbx_PDB_id_code": "1LOL",
            "pdbx_strand_id": "A", "seq_align_beg": "?",
            "pdbx_seq_align_beg_ins_code": "?", "seq_align_end": "?",
            "pdbx_seq_align_end_ins_code": "?", "pdbx_db_accession": "O26232",
            "db_align_beg": "1", "pdbx_db_align_beg_ins_code": "?",
            "db_align_end": "228", "pdbx_db_align_end_ins_code": "?",
            "pdbx_auth_seq_align_beg": "1", "pdbx_auth_seq_align_end": "229",
        }, {
            "align_id": "2", "ref_id": "1", "pdbx_PDB_id_code": "1LOL",
            "pdbx_strand_id": "B", "seq_align_beg": "?",
            "pdbx_seq_align_beg_ins_code": "?", "seq_align_end": "?",
            "pdbx_seq_align_end_ins_code": "?", "pdbx_db_accession": "O26232",
            "db_align_beg": "1", "pdbx_db_align_beg_ins_code": "?",
            "db_align_end": "228", "pdbx_db_align_end_ins_code": "?",
            "pdbx_auth_seq_align_beg": "1001", "pdbx_auth_seq_align_end": "1229",
        }])
        self.assertEqual(len(d["struct_ref_seq_dif"]), 8)
        self.assertEqual(d["struct_ref_seq_dif"][0], {
            "align_id": "1", "pdbx_pdb_id_code": "1LOL", "mon_id": "LEU",
            "pdbx_pdb_strand_id": "A", "seq_num": "?", "pdbx_pdb_ins_code": "?",
            "pdbx_seq_db_name": "UNP", "pdbx_seq_db_accession_code": "O26232",
            "db_mon_id": "MET", "pdbx_seq_db_seq_num": "1",
            "details": "SEE REMARK 999", "pdbx_auth_seq_num": "1",
            "pdbx_ordinal": "1",
        })
        self.assertEqual(d["struct_ref_seq_dif"][-1], {
            "align_id": "2", "pdbx_pdb_id_code": "1LOL", "mon_id": "GLU",
            "pdbx_pdb_strand_id": "B", "seq_num": "?", "pdbx_pdb_ins_code": "?",
            "pdbx_seq_db_name": "UNP", "pdbx_seq_db_accession_code": "O26232",
            "db_mon_id": "?", "pdbx_seq_db_seq_num": "?", "details": "INSERTION",
            "pdbx_auth_seq_num": "1229", "pdbx_ordinal": "8",
        })
        self.assertNotIn("pdbx_struct_mod_residue", d)
        self.assertEqual(len(d["chem_comp"]), 22)
        self.assertEqual(d["chem_comp"][0], {
            "id": "ALA", "type": "L-peptide linking", "mon_nstd_flag": "y",
            "name": "ALANINE", "pdbx_synonyms": "?", "formula": "C3 H7 N O2",
            "formula_weight": "89.093",
        })
        self.assertEqual(d["chem_comp"][4], {
            "id": "BU2", "type": "non-polymer", "mon_nstd_flag": ".",
            "name": "1,3-BUTANEDIOL", "pdbx_synonyms": "?", "formula": "C4 H10 O2",
            "formula_weight": "90.121",
        })
        self.assertEqual(d["chem_comp"][10], {
            "id": "HOH", "type": "non-polymer", "mon_nstd_flag": ".",
            "name": "WATER", "pdbx_synonyms": "?", "formula": "H2 O",
            "formula_weight": "18.015",
        })
        self.assertEqual(d["chem_comp"][-1], {
            "id": "XMP", "type": "non-polymer", "mon_nstd_flag": ".",
            "name": "XANTHOSINE-5'-MONOPHOSPHATE",
            "pdbx_synonyms": "5-MONOPHOSPHATE-9-BETA-D-RIBOFURANOSYL XANTHINE",
            "formula": "C10 H14 N4 O9 P 1+", "formula_weight": "365.213",
        })
        self.assertEqual(d["pdbx_entity_nonpoly"], [
            {"entity_id": "2", "name": "1,3-BUTANEDIOL", "comp_id": "BU2"},
            {"entity_id": "3", "name": "XANTHOSINE-5'-MONOPHOSPHATE", "comp_id": "XMP"},
            {"entity_id": "4", "name": "water", "comp_id": "HOH"},
        ])
        self.assertEqual(d["struct_asym"], [
            {"id": "A", "pdbx_blank_PDB_chainid_flag": "N", "pdbx_modified": "N", "entity_id": "1", "details": "?"},
            {"id": "B", "pdbx_blank_PDB_chainid_flag": "N", "pdbx_modified": "N", "entity_id": "1", "details": "?"},
            {"id": "C", "pdbx_blank_PDB_chainid_flag": "N", "pdbx_modified": "N", "entity_id": "2", "details": "?"},
            {"id": "D", "pdbx_blank_PDB_chainid_flag": "N", "pdbx_modified": "N", "entity_id": "3", "details": "?"},
            {"id": "E", "pdbx_blank_PDB_chainid_flag": "N", "pdbx_modified": "N", "entity_id": "2", "details": "?"},
            {"id": "F", "pdbx_blank_PDB_chainid_flag": "N", "pdbx_modified": "N", "entity_id": "3", "details": "?"},
            {"id": "G", "pdbx_blank_PDB_chainid_flag": "N", "pdbx_modified": "N", "entity_id": "4", "details": "?"},
            {"id": "H", "pdbx_blank_PDB_chainid_flag": "N", "pdbx_modified": "N", "entity_id": "4", "details": "?"},
        ])

        # Atom categories
        self.assertEqual(d["atom_type"], [
            {"symbol": "C"},
            {"symbol": "N"},
            {"symbol": "O"},
            {"symbol": "P"},
            {"symbol": "S"},
        ])
        self.assertEqual(len(d["atom_site"]), 3431)
        self.assertEqual(d["atom_site"][0], {
            "group_PDB": "ATOM", "id": "1", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "VAL",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "11",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "3.696", "Cartn_y": "33.898",
            "Cartn_z": "63.219", "occupancy": "1.00", "B_iso_or_equiv": "21.50",
            "pdbx_formal_charge": "?", "auth_seq_id": "11", "auth_comp_id": "VAL",
            "auth_asym_id": "A", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1557], {
            "group_PDB": "ATOM", "id": "1558", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "VAL",
            "label_asym_id": "B", "label_entity_id": "1", "label_seq_id": "11",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-26.384", "Cartn_y": "61.433",
            "Cartn_z": "36.898", "occupancy": "1.00", "B_iso_or_equiv": "39.30",
            "pdbx_formal_charge": "?", "auth_seq_id": "1011", "auth_comp_id": "VAL",
            "auth_asym_id": "B", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][3191], {
            "group_PDB": "HETATM", "id": "3192", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "BU2",
            "label_asym_id": "C", "label_entity_id": "2", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "2.646", "Cartn_y": "45.112",
            "Cartn_z": "48.995", "occupancy": "1.00", "B_iso_or_equiv": "43.24",
            "pdbx_formal_charge": "?", "auth_seq_id": "5001", "auth_comp_id": "BU2",
            "auth_asym_id": "A", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][3197], {
            "group_PDB": "HETATM", "id": "3198", "type_symbol": "P",
            "label_atom_id": "P", "label_alt_id": ".", "label_comp_id": "XMP",
            "label_asym_id": "D", "label_entity_id": "3", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "3.293", "Cartn_y": "36.948",
            "Cartn_z": "44.605", "occupancy": "1.00", "B_iso_or_equiv": "20.47",
            "pdbx_formal_charge": "?", "auth_seq_id": "2001", "auth_comp_id": "XMP",
            "auth_asym_id": "A", "auth_atom_id": "P", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][3221], {
            "group_PDB": "HETATM", "id": "3222", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "BU2",
            "label_asym_id": "E", "label_entity_id": "2", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-14.563", "Cartn_y": "61.208",
            "Cartn_z": "49.005", "occupancy": "1.00", "B_iso_or_equiv": "45.50",
            "pdbx_formal_charge": "?", "auth_seq_id": "5002", "auth_comp_id": "BU2",
            "auth_asym_id": "B", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][3227], {
            "group_PDB": "HETATM", "id": "3228", "type_symbol": "P",
            "label_atom_id": "P", "label_alt_id": ".", "label_comp_id": "XMP",
            "label_asym_id": "F", "label_entity_id": "3", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-23.846", "Cartn_y": "63.007",
            "Cartn_z": "53.020", "occupancy": "1.00", "B_iso_or_equiv": "33.88",
            "pdbx_formal_charge": "?", "auth_seq_id": "2002", "auth_comp_id": "XMP",
            "auth_asym_id": "B", "auth_atom_id": "P", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][3251], {
            "group_PDB": "HETATM", "id": "3252", "type_symbol": "O",
            "label_atom_id": "O", "label_alt_id": ".", "label_comp_id": "HOH",
            "label_asym_id": "G", "label_entity_id": "4", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-7.624", "Cartn_y": "43.710",
            "Cartn_z": "47.691", "occupancy": "1.00", "B_iso_or_equiv": "15.72",
            "pdbx_formal_charge": "?", "auth_seq_id": "3005", "auth_comp_id": "HOH",
            "auth_asym_id": "A", "auth_atom_id": "O", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][3347], {
            "group_PDB": "HETATM", "id": "3348", "type_symbol": "O",
            "label_atom_id": "O", "label_alt_id": ".", "label_comp_id": "HOH",
            "label_asym_id": "H", "label_entity_id": "4", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-9.220", "Cartn_y": "50.800",
            "Cartn_z": "49.092", "occupancy": "1.00", "B_iso_or_equiv": "17.69",
            "pdbx_formal_charge": "?", "auth_seq_id": "3001", "auth_comp_id": "HOH",
            "auth_asym_id": "B", "auth_atom_id": "O", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][-1], {
            "group_PDB": "HETATM", "id": "3431", "type_symbol": "O",
            "label_atom_id": "O", "label_alt_id": ".", "label_comp_id": "HOH",
            "label_asym_id": "H", "label_entity_id": "4", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-39.239", "Cartn_y": "51.357",
            "Cartn_z": "40.064", "occupancy": "1.00", "B_iso_or_equiv": "37.93",
            "pdbx_formal_charge": "?", "auth_seq_id": "3180", "auth_comp_id": "HOH",
            "auth_asym_id": "B", "auth_atom_id": "O", "pdbx_PDB_model_num": "1",
        })
        self.assertNotIn("atom_site_anisotrop", d)

        # Annotation categories
        self.assertEqual(len(d["struct_conf"]), 22)
        self.assertEqual(d["struct_conf"][0], {
            "conf_type_id": "HELX_P", "id": "HELX_P1", "pdbx_PDB_helix_id": "1",
            "beg_label_comp_id": "VAL", "beg_label_asym_id": "A",
            "beg_label_seq_id": "11", "pdbx_beg_PDB_ins_code": "?",
            "end_label_comp_id": "ASN", "end_label_asym_id": "A",
            "end_label_seq_id": "13", "pdbx_end_PDB_ins_code": "?",
            "beg_auth_comp_id": "VAL", "beg_auth_asym_id": "A",
            "beg_auth_seq_id": "11", "end_auth_comp_id": "ASN",
            "end_auth_asym_id": "A", "end_auth_seq_id": "13",
            "pdbx_PDB_helix_class": "5", "details": "?",
            "pdbx_PDB_helix_length": "3",
        })
        self.assertEqual(d["struct_conf"][-1], {
            "conf_type_id": "HELX_P", "id": "HELX_P22", "pdbx_PDB_helix_id": "22",
            "beg_label_comp_id": "ASN", "beg_label_asym_id": "B",
            "beg_label_seq_id": "210", "pdbx_beg_PDB_ins_code": "?",
            "end_label_comp_id": "LEU", "end_label_asym_id": "B",
            "end_label_seq_id": "225", "pdbx_end_PDB_ins_code": "?",
            "beg_auth_comp_id": "ASN", "beg_auth_asym_id": "B",
            "beg_auth_seq_id": "1210", "end_auth_comp_id": "LEU",
            "end_auth_asym_id": "B", "end_auth_seq_id": "1225",
            "pdbx_PDB_helix_class": "1", "details": "?",
            "pdbx_PDB_helix_length": "16",
        })
        self.assertEqual(d["struct_sheet"], [
            {"id": "A", "type": "?", "number_strands": "9", "details": "?"},
            {"id": "B", "type": "?", "number_strands": "9", "details": "?"},
        ])
        self.assertEqual(len(d["struct_sheet_range"]), 18)
        self.assertEqual(d["struct_sheet_range"][0], {
            "sheet_id": "A", "id": "1", "beg_label_comp_id": "LEU",
            "beg_label_asym_id": "A", "beg_label_seq_id": "15",
            "pdbx_beg_PDB_ins_code": "?", "end_label_comp_id": "MET",
            "end_label_asym_id": "A", "end_label_seq_id": "19",
            "pdbx_end_PDB_ins_code": "?", "beg_auth_comp_id": "LEU",
            "beg_auth_asym_id": "A", "beg_auth_seq_id": "15",
            "end_auth_comp_id": "MET", "end_auth_asym_id": "A",
            "end_auth_seq_id": "19",
        })
        self.assertEqual(d["struct_sheet_range"][-1], {
            "sheet_id": "B", "id": "9", "beg_label_comp_id": "LEU",
            "beg_label_asym_id": "B", "beg_label_seq_id": "15",
            "pdbx_beg_PDB_ins_code": "?", "end_label_comp_id": "ALA",
            "end_label_asym_id": "B", "end_label_seq_id": "18",
            "pdbx_end_PDB_ins_code": "?", "beg_auth_comp_id": "LEU",
            "beg_auth_asym_id": "B", "beg_auth_seq_id": "1015",
            "end_auth_comp_id": "ALA", "end_auth_asym_id": "B",
            "end_auth_seq_id": "1018",
        })
        self.assertEqual(len(d["struct_sheet_order"]), 16)
        self.assertEqual(d["struct_sheet_order"][0], {
            "sheet_id": "A", "range_id_1": "1", "range_id_2": "2", "offset": "?",
            "sense": "parallel",
        })
        self.assertEqual(d["struct_sheet_order"][-1], {
            "sheet_id": "B", "range_id_1": "8", "range_id_2": "9", "offset": "?",
            "sense": "parallel",
        })
        self.assertEqual(len(d["pdbx_struct_sheet_hbond"]), 16)
        self.assertEqual(d["pdbx_struct_sheet_hbond"][0], {
            "sheet_id": "A", "range_id_1": "1", "range_id_2": "2",
            "range_1_label_atom_id": "N", "range_1_label_comp_id": "LEU",
            "range_1_label_asym_id": "A", "range_1_label_seq_id": "17",
            "range_1_PDB_ins_code": "?", "range_1_auth_atom_id": "N",
            "range_1_auth_comp_id": "LEU", "range_1_auth_asym_id": "A",
            "range_1_auth_seq_id": "17", "range_2_label_atom_id": "O",
            "range_2_label_comp_id": "LYS", "range_2_label_asym_id": "A",
            "range_2_label_seq_id": "44", "range_2_PDB_ins_code": "?",
            "range_2_auth_atom_id": "O", "range_2_auth_comp_id": "LYS",
            "range_2_auth_asym_id": "A", "range_2_auth_seq_id": "44",
        })
        self.assertEqual(d["pdbx_struct_sheet_hbond"][-1], {
            "sheet_id": "B", "range_id_1": "8", "range_id_2": "9",
            "range_1_label_atom_id": "O", "range_1_label_comp_id": "VAL",
            "range_1_label_asym_id": "B", "range_1_label_seq_id": "201",
            "range_1_PDB_ins_code": "?", "range_1_auth_atom_id": "O",
            "range_1_auth_comp_id": "VAL", "range_1_auth_asym_id": "B",
            "range_1_auth_seq_id": "1201", "range_2_label_atom_id": "N",
            "range_2_label_comp_id": "ILE", "range_2_label_asym_id": "B",
            "range_2_label_seq_id": "18", "range_2_PDB_ins_code": "?",
            "range_2_auth_atom_id": "N", "range_2_auth_comp_id": "ILE",
            "range_2_auth_asym_id": "B", "range_2_auth_seq_id": "1018",
        })
        self.assertNotIn("struct_conn", d)
        self.assertEqual(d["struct_mon_prot_cis"], [{
            "pdbx_id": "1", "label_comp_id": "ASP", "label_seq_id": "188",
            "label_asym_id": "B", "label_alt_id": ".", "pdbx_PDB_ins_code": "?",
            "auth_comp_id": "ASP", "auth_seq_id": "1188", "auth_asym_id": "B",
            "pdbx_label_comp_id_2": "PRO", "pdbx_label_seq_id_2": "189",
            "pdbx_label_asym_id_2": "B", "pdbx_PDB_ins_code_2": "?",
            "pdbx_auth_comp_id_2": "PRO", "pdbx_auth_seq_id_2": "1189",
            "pdbx_auth_asym_id_2": "B", "pdbx_PDB_model_num": "1",
            "pdbx_omega_angle": "0.35",
        }])
        self.assertEqual(d["struct_site"], [{
            "id": "AC1", "pdbx_evidence_code": "Software", "pdbx_auth_asym_id": "?",
            "pdbx_auth_comp_id": "?", "pdbx_auth_seq_id": "?",
            "pdbx_auth_ins_code": "?", "pdbx_num_residues": "?",
            "details": "BINDING SITE FOR RESIDUE BU2 A 5001",
        }, {
            "id": "AC2", "pdbx_evidence_code": "Software", "pdbx_auth_asym_id": "?",
            "pdbx_auth_comp_id": "?", "pdbx_auth_seq_id": "?",
            "pdbx_auth_ins_code": "?", "pdbx_num_residues": "?",
            "details": "BINDING SITE FOR RESIDUE BU2 B 5002",
        }, {
            "id": "AC3", "pdbx_evidence_code": "Software", "pdbx_auth_asym_id": "?",
            "pdbx_auth_comp_id": "?", "pdbx_auth_seq_id": "?",
            "pdbx_auth_ins_code": "?", "pdbx_num_residues": "?",
            "details": "BINDING SITE FOR RESIDUE XMP A 2001",
        }, {
            "id": "AC4", "pdbx_evidence_code": "Software", "pdbx_auth_asym_id": "?",
            "pdbx_auth_comp_id": "?", "pdbx_auth_seq_id": "?",
            "pdbx_auth_ins_code": "?", "pdbx_num_residues": "?",
            "details": "BINDING SITE FOR RESIDUE XMP B 2002",
        }])
        self.assertEqual(len(d["struct_site_gen"]), 46)
        self.assertEqual(d["struct_site_gen"][0], {
            "id": "1", "site_id": "AC1", "pdbx_num_res": "6",
            "label_comp_id": "ASP", "label_asym_id": "A", "label_seq_id": "70",
            "pdbx_auth_ins_code": "?", "auth_comp_id": "ASP", "auth_asym_id": "A",
            "auth_seq_id": "70", "label_atom_id": ".", "label_alt_id": "?",
            "symmetry": "1_555", "details": "?",
        })
        self.assertEqual(d["struct_site_gen"][-1], {
            "id": "46", "site_id": "AC4", "pdbx_num_res": "16",
            "label_comp_id": "BU2", "label_asym_id": "E", "label_seq_id": ".",
            "pdbx_auth_ins_code": "?", "auth_comp_id": "BU2", "auth_asym_id": "B",
            "auth_seq_id": "5002", "label_atom_id": ".", "label_alt_id": "?",
            "symmetry": "1_555", "details": "?",
        })

        # Crystal categories
        self.assertEqual(d["cell"], [{
            "entry_id": "1LOL", "length_a": "57.570", "length_b": "55.482",
            "length_c": "66.129", "angle_alpha": "90.00", "angle_beta": "94.28",
            "angle_gamma": "90.00", "Z_pdb": "4", "pdbx_unique_axis": "?",
        }])
        self.assertEqual(d["symmetry"], [{
            "entry_id": "1LOL", "space_group_name_H-M": "P 1 21 1",
            "pdbx_full_space_group_name_H-M": "?", "cell_setting": "?",
            "Int_Tables_number": "?",
        }])
        self.assertEqual(d["atom_sites"], [{
            "entry_id": "1LOL", "fract_transf_matrix[1][1]": "0.017370",
            "fract_transf_matrix[1][2]": "0.000000",
            "fract_transf_matrix[1][3]": "0.001301",
            "fract_transf_matrix[2][1]": "0.000000",
            "fract_transf_matrix[2][2]": "0.018024",
            "fract_transf_matrix[2][3]": "0.000000",
            "fract_transf_matrix[3][1]": "0.000000",
            "fract_transf_matrix[3][2]": "0.000000",
            "fract_transf_matrix[3][3]": "0.015164",
            "fract_transf_vector[1]": "0.00000",
            "fract_transf_vector[2]": "0.00000",
            "fract_transf_vector[3]": "0.00000",
        }])
        self.assertEqual(d["database_PDB_matrix"], [{
            "entry_id": "1LOL", "origx[1][1]": "1.000000",
            "origx[1][2]": "0.000000", "origx[1][3]": "0.000000",
            "origx[2][1]": "0.000000", "origx[2][2]": "1.000000",
            "origx[2][3]": "0.000000", "origx[3][1]": "0.000000",
            "origx[3][2]": "0.000000", "origx[3][3]": "1.000000",
            "origx_vector[1]": "0.00000", "origx_vector[2]": "0.00000",
            "origx_vector[3]": "0.00000",
        }])
        self.assertNotIn("struct_ncs_oper", d)

        # Citation categories
        self.assertEqual(d["audit_author"], [
            {"name": "Wu, N", "pdbx_ordinal": "1"},
            {"name": "Pai, E.F", "pdbx_ordinal": "2"},
        ])
        self.assertEqual(d["citation_author"], [
            {"citation_id": "primary", "name": "Wu, N", "pdbx_ordinal": "1"},
            {"citation_id": "primary", "name": "Pai, E.F", "pdbx_ordinal": "2"},
        ])
        self.assertNotIn("citation_editor", d)
        self.assertEqual(d["citation"], [{
            "id": "primary",
            "title": "CRYSTAL STRUCTURES OF INHIBITOR COMPLEXES REVEAL AN ALTERNATE BINDING MODE IN OROTIDINE-5'-MONOPHOSPHATE DECARBOXYLASE.",
            "journal_abbrev": "J.Biol.Chem.", "journal_volume": "277",
            "page_first": "28080", "page_last": "?", "year": "2002",
            "journal_id_ASTM": "?", "country": "?", "journal_id_ISSN": "0021-9258",
            "journal_id_CSD": "?", "book_publisher": "?",
            "pdbx_database_id_PubMed": "12011084",
            "pdbx_database_id_DOI": "10.1074/jbc.m202362200",
        }])

        # Experimental categories
        self.assertEqual(d["reflns"], [{
            "entry_id": "1LOL", "observed_criterion_sigma_I": "?",
            "observed_criterion_sigma_F": "?", "d_resolution_low": "?",
            "d_resolution_high": "1.90", "number_obs": "?", "number_all": "?",
            "percent_possible_obs": "?", "pdbx_Rmerge_I_obs": "?",
            "pdbx_Rsym_value": "?", "pdbx_netI_over_sigmaI": "?",
            "B_iso_Wilson_estimate": "11.20", "pdbx_redundancy": "?",
            "R_free_details": "?", "limit_h_max": "?", "limit_h_min": "?",
            "limit_k_max": "?", "limit_k_min": "?", "limit_l_max": "?",
            "limit_l_min": "?", "observed_criterion_F_max": "?",
            "observed_criterion_F_min": "?", "pdbx_chi_squared": "?",
            "pdbx_scaling_rejects": "?", "pdbx_diffrn_id": "?", "pdbx_ordinal": "?",
        }])
        self.assertEqual(d["refine"], [{
            "entry_id": "1LOL", "ls_number_reflns_obs": "32092",
            "ls_number_reflns_all": "?", "pdbx_ls_sigma_I": "?",
            "pdbx_ls_sigma_F": "?", "pdbx_data_cutoff_high_absF": "679650.230",
            "pdbx_data_cutoff_low_absF": "0.0000", "ls_d_res_low": "27.07",
            "ls_d_res_high": "1.90", "ls_percent_reflns_obs": "97.2",
            "ls_R_factor_obs": "0.193", "ls_R_factor_all": "0.193",
            "ls_R_factor_R_work": "0.193", "ls_R_factor_R_free": "0.229",
            "ls_R_factor_R_free_error": "0.006",
            "ls_R_factor_R_free_error_details": "?",
            "ls_percent_reflns_R_free": "4.900", "ls_number_reflns_R_free": "1583",
            "ls_number_parameters": "?", "ls_number_restraints": "?",
            "occupancy_min": "?", "occupancy_max": "?",
            "correlation_coeff_Fo_to_Fc": "?",
            "correlation_coeff_Fo_to_Fc_free": "?", "B_iso_mean": "24.20",
            "aniso_B[1][1]": "7.19000", "aniso_B[2][2]": "-3.85000",
            "aniso_B[3][3]": "-3.34000", "aniso_B[1][2]": "0.00000",
            "aniso_B[1][3]": "3.48000", "aniso_B[2][3]": "0.00000",
            "solvent_model_details": "FLAT MODEL",
            "solvent_model_param_ksol": "0.39", "solvent_model_param_bsol": "54.02",
            "pdbx_solvent_vdw_probe_radii": "?",
            "pdbx_solvent_ion_probe_radii": "?",
            "pdbx_solvent_shrinkage_radii": "?",
            "pdbx_ls_cross_valid_method": "THROUGHOUT", "details": "?",
            "pdbx_starting_model": "?", "pdbx_method_to_determine_struct": "?",
            "pdbx_isotropic_thermal_model": "RESTRAINED",
            "pdbx_stereochemistry_target_values": "ENGH & HUBER",
            "pdbx_stereochem_target_val_spec_case": "?",
            "pdbx_R_Free_selection_details": "RANDOM",
            "pdbx_overall_ESU_R_Free": "?", "overall_SU_B": "?",
            "ls_redundancy_reflns_obs": "?", "B_iso_min": "?", "B_iso_max": "?",
            "overall_SU_R_Cruickshank_DPI": "?", "overall_SU_R_free": "?",
            "overall_SU_ML": "?", "pdbx_overall_ESU_R": "?",
            "pdbx_data_cutoff_high_rms_absF": "679650.230", "pdbx_refine_id": "?",
            "pdbx_overall_phase_error": "?", "ls_wR_factor_R_free": "?",
            "ls_wR_factor_R_work": "?", "overall_FOM_free_R_set": "?",
            "overall_FOM_work_R_set": "?", "pdbx_diffrn_id": "?",
            "pdbx_TLS_residual_ADP_flag": "?",
            "pdbx_overall_SU_R_free_Cruickshank_DPI": "?",
            "pdbx_overall_SU_R_Blow_DPI": "?",
            "pdbx_overall_SU_R_free_Blow_DPI": "?",
        }])

        # Missing items categories
        self.assertEqual(len(d["pdbx_unobs_or_zero_occ_residues"]), 40)
        self.assertEqual(d["pdbx_unobs_or_zero_occ_residues"][0], {
            "id": "1", "PDB_model_num": "1", "polymer_flag": "Y",
            "occupancy_flag": "1", "auth_asym_id": "A", "auth_comp_id": "LEU",
            "auth_seq_id": "1", "PDB_ins_code": "?", "label_asym_id": "A",
            "label_comp_id": "LEU", "label_seq_id": "1",
        })
        self.assertEqual(d["pdbx_unobs_or_zero_occ_residues"][-1], {
            "id": "40", "PDB_model_num": "1", "polymer_flag": "Y",
            "occupancy_flag": "1", "auth_asym_id": "B", "auth_comp_id": "GLY",
            "auth_seq_id": "1186", "PDB_ins_code": "?", "label_asym_id": "B",
            "label_comp_id": "GLY", "label_seq_id": "1186",
        })

        # Assembly categories
        self.assertEqual(d["pdbx_struct_assembly"], [
            {"id": "1", "details": "?", "method_details": "PISA", "oligomeric_details": "?", "oligomeric_count": "?"},
        ])
        self.assertEqual(d["pdbx_struct_assembly_gen"], [
            {"assembly_id": "1", "oper_expression": "1", "asym_id_list": "A,B,C,D,E,F,G,H"},
        ])
        self.assertEqual(d["pdbx_struct_assembly_prop"], [
            {"biol_id": "1", "type": "ABSA (A^2)", "value": "5230", "details": "?"},
            {"biol_id": "1", "type": "MORE", "value": "-31.0", "details": "?"},
            {"biol_id": "1", "type": "SSA (A^2)", "value": "16550", "details": "?"},
        ])
        self.assertEqual(d["pdbx_struct_oper_list"], [{
            "id": "1", "type": "?", "name": "?", "symmetry_operation": "x,y,z",
            "matrix[1][1]": "1.0", "matrix[1][2]": "0.0", "matrix[1][3]": "0.0",
            "vector[1]": "0.0", "matrix[2][1]": "0.0", "matrix[2][2]": "1.0",
            "matrix[2][3]": "0.0", "vector[2]": "0.0", "matrix[3][1]": "0.0",
            "matrix[3][2]": "0.0", "matrix[3][3]": "1.0", "vector[3]": "0.0",
        }])


    def test_5xme(self):
        d = atomium.open("tests/integration/files/5xme.pdb", dictionary=True)
        self.assertEqual(len(d.keys()), 33)

        # Metadata categories
        self.assertEqual(d["entry"], [{"id": "5XME"}])
        self.assertEqual(d["pdbx_database_status"], [
            {"status_code": "REL", "entry_id": "5XME", "recvd_initial_deposition_date": "2017-05-15"},
        ])
        self.assertEqual(d["struct_keywords"], [
            {"entry_id": "5XME", "pdbx_keywords": "APOPTOSIS", "text": "TRADD, DEATH DOMAIN, APOPTOSIS"},
        ])
        self.assertNotIn("pdbx_database_PDB_obs_spr", d)
        self.assertEqual(d["struct"], [{
            "entry_id": "5XME",
            "title": "SOLUTION STRUCTURE OF C-TERMINAL DOMAIN OF TRADD",
            "pdbx_descriptor": "?", "pdbx_model_details": "?",
            "pdbx_CASP_flag": "?", "pdbx_model_type_details": "?",
        }])
        self.assertNotIn("pdbx_database_related", d)
        self.assertNotIn("database_PDB_caveat", d)
        self.assertEqual(d["exptl"], [{"entry_id": "5XME", "method": "SOLUTION NMR"}])
        self.assertNotIn("pdbx_coordinate_model", d)
        self.assertEqual(d["pdbx_audit_revision_history"], [{
            "ordinal": "1", "data_content_type": "Structure model",
            "major_revision": "1", "minor_revision": "0",
            "revision_date": "2017-09-06",
        }])

        # Entity categories
        self.assertEqual(d["entity"], [{
            "id": "1", "type": "polymer", "src_method": "man",
            "pdbx_description": "TUMOR NECROSIS FACTOR RECEPTOR TYPE 1-ASSOCIATED DEATH DOMAIN PROTEIN",
            "formula_weight": "?", "pdbx_number_of_molecules": "1", "pdbx_ec": "?",
            "pdbx_mutation": "?", "pdbx_fragment": "UNP RESIDUES 199-312",
            "details": "?",
        }])
        self.assertEqual(d["entity_name_com"], [
            {"entity_id": "1", "name": "TNFR1-ASSOCIATED DEATH DOMAIN PROTEIN,TNFRSF1A-ASSOCIATED VIA DEATH DOMAIN"},
        ])
        self.assertEqual(d["entity_poly"], [{
            "entity_id": "1", "type": "polypeptide(L)", "nstd_linkage": "no",
            "nstd_monomer": "no",
            "pdbx_seq_one_letter_code": "MHHHHHHSSGRGSAQTFLFQGQPVVNRPLSLKDQQTFARSVGLKWRKVGRSLQRGCRALRDPALDSLAYEYEREGLYEQAFQLLRRFVQAEGRRATLQRLVEALEENELTSLAEDLLGLTDPNGGLA",
            "pdbx_seq_one_letter_code_can": "MHHHHHHSSGRGSAQTFLFQGQPVVNRPLSLKDQQTFARSVGLKWRKVGRSLQRGCRALRDPALDSLAYEYEREGLYEQAFQLLRRFVQAEGRRATLQRLVEALEENELTSLAEDLLGLTDPNGGLA",
            "pdbx_strand_id": "A", "pdbx_target_identifier": "?",
        }])
        self.assertEqual(len(d["entity_poly_seq"]), 127)
        self.assertEqual(d["entity_poly_seq"][0], {
            "entity_id": "1", "num": "1", "mon_id": "MET", "hetero": "n",
        })
        self.assertEqual(d["entity_poly_seq"][-1], {
            "entity_id": "1", "num": "127", "mon_id": "ALA", "hetero": "n",
        })
        self.assertEqual(d["struct_ref"], [{
            "id": "1", "db_name": "UNP", "db_code": "TRADD_HUMAN",
            "pdbx_db_accession": "Q15628", "pdbx_db_isoform": "?", "entity_id": "1",
            "pdbx_seq_one_letter_code": "?", "pdbx_align_begin": "?",
        }])
        self.assertEqual(d["struct_ref_seq"], [{
            "align_id": "1", "ref_id": "1", "pdbx_PDB_id_code": "5XME",
            "pdbx_strand_id": "A", "seq_align_beg": "?",
            "pdbx_seq_align_beg_ins_code": "?", "seq_align_end": "?",
            "pdbx_seq_align_end_ins_code": "?", "pdbx_db_accession": "Q15628",
            "db_align_beg": "199", "pdbx_db_align_beg_ins_code": "?",
            "db_align_end": "312", "pdbx_db_align_end_ins_code": "?",
            "pdbx_auth_seq_align_beg": "199", "pdbx_auth_seq_align_end": "312",
        }])
        self.assertEqual(len(d["struct_ref_seq_dif"]), 13)
        self.assertEqual(d["struct_ref_seq_dif"][0], {
            "align_id": "1", "pdbx_pdb_id_code": "5XME", "mon_id": "MET",
            "pdbx_pdb_strand_id": "A", "seq_num": "?", "pdbx_pdb_ins_code": "?",
            "pdbx_seq_db_name": "UNP", "pdbx_seq_db_accession_code": "Q15628",
            "db_mon_id": "?", "pdbx_seq_db_seq_num": "?",
            "details": "EXPRESSION TAG", "pdbx_auth_seq_num": "186",
            "pdbx_ordinal": "1",
        })
        self.assertEqual(d["struct_ref_seq_dif"][-1], {
            "align_id": "1", "pdbx_pdb_id_code": "5XME", "mon_id": "SER",
            "pdbx_pdb_strand_id": "A", "seq_num": "?", "pdbx_pdb_ins_code": "?",
            "pdbx_seq_db_name": "UNP", "pdbx_seq_db_accession_code": "Q15628",
            "db_mon_id": "?", "pdbx_seq_db_seq_num": "?",
            "details": "EXPRESSION TAG", "pdbx_auth_seq_num": "198",
            "pdbx_ordinal": "13",
        })
        self.assertNotIn("pdbx_struct_mod_residue", d)
        self.assertEqual(len(d["chem_comp"]), 19)
        self.assertEqual(d["chem_comp"][0], {
            "id": "ALA", "type": "L-peptide linking", "mon_nstd_flag": "y",
            "name": "ALANINE", "pdbx_synonyms": "?", "formula": "C3 H7 N O2",
            "formula_weight": "89.093",
        })
        self.assertEqual(d["chem_comp"][-1], {
            "id": "VAL", "type": "L-peptide linking", "mon_nstd_flag": "y",
            "name": "VALINE", "pdbx_synonyms": "?", "formula": "C5 H11 N O2",
            "formula_weight": "117.146",
        })
        self.assertNotIn("pdbx_entity_nonpoly", d)
        self.assertEqual(d["struct_asym"], [
            {"id": "A", "pdbx_blank_PDB_chainid_flag": "N", "pdbx_modified": "N", "entity_id": "1", "details": "?"},
        ])

        # Atom categories
        self.assertEqual(d["atom_type"], [
            {"symbol": "C"},
            {"symbol": "H"},
            {"symbol": "N"},
            {"symbol": "O"},
            {"symbol": "S"},
        ])
        self.assertEqual(len(d["atom_site"]), 18270)
        self.assertEqual(d["atom_site"][0], {
            "group_PDB": "ATOM", "id": "1", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "ALA",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "14",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "33.969", "Cartn_y": "-8.430",
            "Cartn_z": "-0.271", "occupancy": "1.00", "B_iso_or_equiv": "0.00",
            "pdbx_formal_charge": "?", "auth_seq_id": "199", "auth_comp_id": "ALA",
            "auth_asym_id": "A", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1827], {
            "group_PDB": "ATOM", "id": "1828", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "ALA",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "14",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "34.064", "Cartn_y": "-8.092",
            "Cartn_z": "-0.062", "occupancy": "1.00", "B_iso_or_equiv": "0.00",
            "pdbx_formal_charge": "?", "auth_seq_id": "199", "auth_comp_id": "ALA",
            "auth_asym_id": "A", "auth_atom_id": "N", "pdbx_PDB_model_num": "2",
        })
        self.assertEqual(d["atom_site"][3654], {
            "group_PDB": "ATOM", "id": "3655", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "ALA",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "14",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "37.369", "Cartn_y": "-7.973",
            "Cartn_z": "-0.242", "occupancy": "1.00", "B_iso_or_equiv": "0.00",
            "pdbx_formal_charge": "?", "auth_seq_id": "199", "auth_comp_id": "ALA",
            "auth_asym_id": "A", "auth_atom_id": "N", "pdbx_PDB_model_num": "3",
        })
        self.assertEqual(d["atom_site"][5481], {
            "group_PDB": "ATOM", "id": "5482", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "ALA",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "14",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "36.023", "Cartn_y": "-7.429",
            "Cartn_z": "-0.637", "occupancy": "1.00", "B_iso_or_equiv": "0.00",
            "pdbx_formal_charge": "?", "auth_seq_id": "199", "auth_comp_id": "ALA",
            "auth_asym_id": "A", "auth_atom_id": "N", "pdbx_PDB_model_num": "4",
        })
        self.assertEqual(d["atom_site"][7308], {
            "group_PDB": "ATOM", "id": "7309", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "ALA",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "14",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "35.245", "Cartn_y": "-9.094",
            "Cartn_z": "0.245", "occupancy": "1.00", "B_iso_or_equiv": "0.00",
            "pdbx_formal_charge": "?", "auth_seq_id": "199", "auth_comp_id": "ALA",
            "auth_asym_id": "A", "auth_atom_id": "N", "pdbx_PDB_model_num": "5",
        })
        self.assertEqual(d["atom_site"][9135], {
            "group_PDB": "ATOM", "id": "9136", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "ALA",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "14",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "35.835", "Cartn_y": "-7.648",
            "Cartn_z": "-0.888", "occupancy": "1.00", "B_iso_or_equiv": "0.00",
            "pdbx_formal_charge": "?", "auth_seq_id": "199", "auth_comp_id": "ALA",
            "auth_asym_id": "A", "auth_atom_id": "N", "pdbx_PDB_model_num": "6",
        })
        self.assertEqual(d["atom_site"][10962], {
            "group_PDB": "ATOM", "id": "10963", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "ALA",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "14",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "37.525", "Cartn_y": "-7.759",
            "Cartn_z": "-0.297", "occupancy": "1.00", "B_iso_or_equiv": "0.00",
            "pdbx_formal_charge": "?", "auth_seq_id": "199", "auth_comp_id": "ALA",
            "auth_asym_id": "A", "auth_atom_id": "N", "pdbx_PDB_model_num": "7",
        })
        self.assertEqual(d["atom_site"][12789], {
            "group_PDB": "ATOM", "id": "12790", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "ALA",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "14",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "35.062", "Cartn_y": "-9.220",
            "Cartn_z": "0.031", "occupancy": "1.00", "B_iso_or_equiv": "0.00",
            "pdbx_formal_charge": "?", "auth_seq_id": "199", "auth_comp_id": "ALA",
            "auth_asym_id": "A", "auth_atom_id": "N", "pdbx_PDB_model_num": "8",
        })
        self.assertEqual(d["atom_site"][14616], {
            "group_PDB": "ATOM", "id": "14617", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "ALA",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "14",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "36.244", "Cartn_y": "-7.381",
            "Cartn_z": "-0.745", "occupancy": "1.00", "B_iso_or_equiv": "0.00",
            "pdbx_formal_charge": "?", "auth_seq_id": "199", "auth_comp_id": "ALA",
            "auth_asym_id": "A", "auth_atom_id": "N", "pdbx_PDB_model_num": "9",
        })
        self.assertEqual(d["atom_site"][16443], {
            "group_PDB": "ATOM", "id": "16444", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "ALA",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "14",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "37.677", "Cartn_y": "-7.651",
            "Cartn_z": "-0.218", "occupancy": "1.00", "B_iso_or_equiv": "0.00",
            "pdbx_formal_charge": "?", "auth_seq_id": "199", "auth_comp_id": "ALA",
            "auth_asym_id": "A", "auth_atom_id": "N", "pdbx_PDB_model_num": "10",
        })
        self.assertEqual(d["atom_site"][-1], {
            "group_PDB": "ATOM", "id": "18270", "type_symbol": "H",
            "label_atom_id": "HB3", "label_alt_id": ".", "label_comp_id": "ALA",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "127",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "22.641", "Cartn_y": "-15.379",
            "Cartn_z": "-10.884", "occupancy": "1.00", "B_iso_or_equiv": "0.00",
            "pdbx_formal_charge": "?", "auth_seq_id": "312", "auth_comp_id": "ALA",
            "auth_asym_id": "A", "auth_atom_id": "HB3", "pdbx_PDB_model_num": "10",
        })
        self.assertNotIn("atom_site_anisotrop", d)

        # Annotation categories
        self.assertEqual(len(d["struct_conf"]), 7)
        self.assertEqual(d["struct_conf"][0], {
            "conf_type_id": "HELX_P", "id": "HELX_P1", "pdbx_PDB_helix_id": "1",
            "beg_label_comp_id": "SER", "beg_label_asym_id": "A",
            "beg_label_seq_id": "30", "pdbx_beg_PDB_ins_code": "?",
            "end_label_comp_id": "GLY", "end_label_asym_id": "A",
            "end_label_seq_id": "42", "pdbx_end_PDB_ins_code": "?",
            "beg_auth_comp_id": "SER", "beg_auth_asym_id": "A",
            "beg_auth_seq_id": "215", "end_auth_comp_id": "GLY",
            "end_auth_asym_id": "A", "end_auth_seq_id": "227",
            "pdbx_PDB_helix_class": "1", "details": "?",
            "pdbx_PDB_helix_length": "13",
        })
        self.assertEqual(d["struct_conf"][-1], {
            "conf_type_id": "HELX_P", "id": "HELX_P7", "pdbx_PDB_helix_id": "7",
            "beg_label_comp_id": "LEU", "beg_label_asym_id": "A",
            "beg_label_seq_id": "109", "pdbx_beg_PDB_ins_code": "?",
            "end_label_comp_id": "LEU", "end_label_asym_id": "A",
            "end_label_seq_id": "117", "pdbx_end_PDB_ins_code": "?",
            "beg_auth_comp_id": "LEU", "beg_auth_asym_id": "A",
            "beg_auth_seq_id": "294", "end_auth_comp_id": "LEU",
            "end_auth_asym_id": "A", "end_auth_seq_id": "302",
            "pdbx_PDB_helix_class": "1", "details": "?",
            "pdbx_PDB_helix_length": "9",
        })
        self.assertEqual(d["struct_sheet"], [{"id": "AA1", "type": "?", "number_strands": "2", "details": "?"}])
        self.assertEqual(d["struct_sheet_range"], [{
            "sheet_id": "AA1", "id": "1", "beg_label_comp_id": "THR",
            "beg_label_asym_id": "A", "beg_label_seq_id": "16",
            "pdbx_beg_PDB_ins_code": "?", "end_label_comp_id": "PHE",
            "end_label_asym_id": "A", "end_label_seq_id": "19",
            "pdbx_end_PDB_ins_code": "?", "beg_auth_comp_id": "THR",
            "beg_auth_asym_id": "A", "beg_auth_seq_id": "201",
            "end_auth_comp_id": "PHE", "end_auth_asym_id": "A",
            "end_auth_seq_id": "204",
        }, {
            "sheet_id": "AA1", "id": "2", "beg_label_comp_id": "GLN",
            "beg_label_asym_id": "A", "beg_label_seq_id": "22",
            "pdbx_beg_PDB_ins_code": "?", "end_label_comp_id": "VAL",
            "end_label_asym_id": "A", "end_label_seq_id": "25",
            "pdbx_end_PDB_ins_code": "?", "beg_auth_comp_id": "GLN",
            "beg_auth_asym_id": "A", "beg_auth_seq_id": "207",
            "end_auth_comp_id": "VAL", "end_auth_asym_id": "A",
            "end_auth_seq_id": "210",
        }])
        self.assertEqual(d["struct_sheet_order"], [
            {"sheet_id": "AA1", "range_id_1": "1", "range_id_2": "2", "offset": "?", "sense": "anti-parallel"},
        ])
        self.assertEqual(d["pdbx_struct_sheet_hbond"], [{
            "sheet_id": "AA1", "range_id_1": "1", "range_id_2": "2",
            "range_1_label_atom_id": "N", "range_1_label_comp_id": "PHE",
            "range_1_label_asym_id": "A", "range_1_label_seq_id": "19",
            "range_1_PDB_ins_code": "?", "range_1_auth_atom_id": "N",
            "range_1_auth_comp_id": "PHE", "range_1_auth_asym_id": "A",
            "range_1_auth_seq_id": "204", "range_2_label_atom_id": "O",
            "range_2_label_comp_id": "GLN", "range_2_label_asym_id": "A",
            "range_2_label_seq_id": "25", "range_2_PDB_ins_code": "?",
            "range_2_auth_atom_id": "O", "range_2_auth_comp_id": "GLN",
            "range_2_auth_asym_id": "A", "range_2_auth_seq_id": "210",
        }])
        self.assertNotIn("struct_conn", d)
        self.assertNotIn("struct_mon_prot_cis", d)
        self.assertNotIn("struct_site", d)
        self.assertNotIn("struct_site_gen", d)

        # Crystal categories
        self.assertEqual(d["cell"], [{
            "entry_id": "5XME", "length_a": "1.000", "length_b": "1.000",
            "length_c": "1.000", "angle_alpha": "90.00", "angle_beta": "90.00",
            "angle_gamma": "90.00", "Z_pdb": "1", "pdbx_unique_axis": "?",
        }])
        self.assertEqual(d["symmetry"], [{
            "entry_id": "5XME", "space_group_name_H-M": "P 1",
            "pdbx_full_space_group_name_H-M": "?", "cell_setting": "?",
            "Int_Tables_number": "?",
        }])
        self.assertEqual(d["atom_sites"], [{
            "entry_id": "5XME", "fract_transf_matrix[1][1]": "1.000000",
            "fract_transf_matrix[1][2]": "0.000000",
            "fract_transf_matrix[1][3]": "0.000000",
            "fract_transf_matrix[2][1]": "0.000000",
            "fract_transf_matrix[2][2]": "1.000000",
            "fract_transf_matrix[2][3]": "0.000000",
            "fract_transf_matrix[3][1]": "0.000000",
            "fract_transf_matrix[3][2]": "0.000000",
            "fract_transf_matrix[3][3]": "1.000000",
            "fract_transf_vector[1]": "0.00000",
            "fract_transf_vector[2]": "0.00000",
            "fract_transf_vector[3]": "0.00000",
        }])
        self.assertEqual(d["database_PDB_matrix"], [{
            "entry_id": "5XME", "origx[1][1]": "1.000000",
            "origx[1][2]": "0.000000", "origx[1][3]": "0.000000",
            "origx[2][1]": "0.000000", "origx[2][2]": "1.000000",
            "origx[2][3]": "0.000000", "origx[3][1]": "0.000000",
            "origx[3][2]": "0.000000", "origx[3][3]": "1.000000",
            "origx_vector[1]": "0.00000", "origx_vector[2]": "0.00000",
            "origx_vector[3]": "0.00000",
        }])
        self.assertNotIn("struct_ncs_oper", d)

        # Citation categories
        self.assertEqual(d["audit_author"], [
            {"name": "Lin, Z", "pdbx_ordinal": "1"},
            {"name": "Zhang, N", "pdbx_ordinal": "2"},
        ])
        self.assertEqual(d["citation_author"], [
            {"citation_id": "primary", "name": "Zhang, N", "pdbx_ordinal": "1"},
            {"citation_id": "primary", "name": "Yuan, W", "pdbx_ordinal": "2"},
            {"citation_id": "primary", "name": "Fan, J.S", "pdbx_ordinal": "3"},
            {"citation_id": "primary", "name": "Lin, Z", "pdbx_ordinal": "4"},
        ])
        self.assertNotIn("citation_editor", d)
        self.assertEqual(d["citation"], [{
            "id": "primary",
            "title": "STRUCTURE OF THE C-TERMINAL DOMAIN OF TRADD REVEALS A NOVEL FOLD IN THE DEATH DOMAIN SUPERFAMILY.",
            "journal_abbrev": "Sci Rep", "journal_volume": "7",
            "page_first": "7073", "page_last": "?", "year": "2017",
            "journal_id_ASTM": "?", "country": "?", "journal_id_ISSN": "2045-2322",
            "journal_id_CSD": "?", "book_publisher": "?",
            "pdbx_database_id_PubMed": "28765645",
            "pdbx_database_id_DOI": "10.1038/s41598-017-07348-9",
        }])

        # Experimental categories
        self.assertNotIn("reflns", d)
        self.assertNotIn("refine", d)

        # Missing items categories
        self.assertEqual(len(d["pdbx_unobs_or_zero_occ_residues"]), 13)
        self.assertEqual(d["pdbx_unobs_or_zero_occ_residues"][0], {
            "id": "1", "PDB_model_num": "1", "polymer_flag": "Y",
            "occupancy_flag": "1", "auth_asym_id": "A", "auth_comp_id": "MET",
            "auth_seq_id": "186", "PDB_ins_code": "?", "label_asym_id": "A",
            "label_comp_id": "MET", "label_seq_id": "186",
        })
        self.assertEqual(d["pdbx_unobs_or_zero_occ_residues"][-1], {
            "id": "13", "PDB_model_num": "1", "polymer_flag": "Y",
            "occupancy_flag": "1", "auth_asym_id": "A", "auth_comp_id": "SER",
            "auth_seq_id": "198", "PDB_ins_code": "?", "label_asym_id": "A",
            "label_comp_id": "SER", "label_seq_id": "198",
        })

        # Assembly categories
        self.assertEqual(d["pdbx_struct_assembly"], [
            {"id": "1", "details": "?", "method_details": "?", "oligomeric_details": "?", "oligomeric_count": "?"},
        ])
        self.assertEqual(d["pdbx_struct_assembly_gen"], [{"assembly_id": "1", "oper_expression": "1", "asym_id_list": "A"}])
        self.assertNotIn("pdbx_struct_assembly_prop", d)
        self.assertEqual(d["pdbx_struct_oper_list"], [{
            "id": "1", "type": "?", "name": "?", "symmetry_operation": "x,y,z",
            "matrix[1][1]": "1.0", "matrix[1][2]": "0.0", "matrix[1][3]": "0.0",
            "vector[1]": "0.0", "matrix[2][1]": "0.0", "matrix[2][2]": "1.0",
            "matrix[2][3]": "0.0", "vector[2]": "0.0", "matrix[3][1]": "0.0",
            "matrix[3][2]": "0.0", "matrix[3][3]": "1.0", "vector[3]": "0.0",
        }])


    def test_1gsg(self):
        d = atomium.open("tests/integration/files/1gsg.pdb", dictionary=True)
        self.assertEqual(len(d.keys()), 30)

        # Metadata categories
        self.assertEqual(d["entry"], [{"id": "1GSG"}])
        self.assertEqual(d["pdbx_database_status"], [
            {"status_code": "REL", "entry_id": "1GSG", "recvd_initial_deposition_date": "1990-04-03"},
        ])
        self.assertEqual(d["struct_keywords"], [{
            "entry_id": "1GSG", "pdbx_keywords": "LIGASE/RNA",
            "text": "PROTEIN-T-RNA COMPLEX, SINGLE STRAND, PROTEIN/RNA, LIGASE-RNA COMPLEX",
        }])
        self.assertNotIn("pdbx_database_PDB_obs_spr", d)
        self.assertEqual(d["struct"], [{
            "entry_id": "1GSG",
            "title": "STRUCTURE OF E.COLI GLUTAMINYL-TRNA SYNTHETASE COMPLEXED WITH TRNAGLN AND ATP AT 2.8 ANGSTROMS RESOLUTION",
            "pdbx_descriptor": "?", "pdbx_model_details": "?",
            "pdbx_CASP_flag": "?", "pdbx_model_type_details": "?",
        }])
        self.assertNotIn("pdbx_database_related", d)
        self.assertNotIn("database_PDB_caveat", d)
        self.assertEqual(d["exptl"], [{"entry_id": "1GSG", "method": "X-RAY DIFFRACTION"}])
        self.assertEqual(d["pdbx_coordinate_model"], [
            {"asym_id": "P", "type": "CA ATOMS ONLY"},
            {"asym_id": "T", "type": "P ATOMS ONLY"},
        ])
        self.assertEqual(d["pdbx_audit_revision_history"], [{
            "ordinal": "1", "data_content_type": "Structure model",
            "major_revision": "1", "minor_revision": "0",
            "revision_date": "1992-02-24",
        }, {
            "ordinal": "2", "data_content_type": "Structure model",
            "major_revision": "2", "minor_revision": "0",
            "revision_date": "2001-09-21",
        }, {
            "ordinal": "3", "data_content_type": "Structure model",
            "major_revision": "3", "minor_revision": "0",
            "revision_date": "2003-04-01",
        }, {
            "ordinal": "4", "data_content_type": "Structure model",
            "major_revision": "4", "minor_revision": "0",
            "revision_date": "2009-02-24",
        }, {
            "ordinal": "5", "data_content_type": "Structure model",
            "major_revision": "5", "minor_revision": "0",
            "revision_date": "2011-12-07",
        }, {
            "ordinal": "6", "data_content_type": "Structure model",
            "major_revision": "6", "minor_revision": "0",
            "revision_date": "2012-05-09",
        }])

        # Entity categories
        self.assertEqual(d["entity"], [{
            "id": "1", "type": "polymer", "src_method": "syn",
            "pdbx_description": "TRNAGLN", "formula_weight": "?",
            "pdbx_number_of_molecules": "1", "pdbx_ec": "?", "pdbx_mutation": "?",
            "pdbx_fragment": "?", "details": "?",
        }, {
            "id": "2", "type": "polymer", "src_method": "nat",
            "pdbx_description": "GLUTAMINYL-TRNA SYNTHETASE", "formula_weight": "?",
            "pdbx_number_of_molecules": "1", "pdbx_ec": "6.1.1.18",
            "pdbx_mutation": "?", "pdbx_fragment": "?", "details": "?",
        }])
        self.assertNotIn("entity_name_com", d)
        self.assertEqual(d["entity_poly"], [{
            "entity_id": "1", "type": "polyribonucleotide", "nstd_linkage": "no",
            "nstd_monomer": "yes",
            "pdbx_seq_one_letter_code": "UGGGGUA(4SU)CGCCAAGC(OMG)G(H2U)AAGGCACCGGA(OMU)UCUG(2MA)(PSU)(PSU)CCGGCAUUCCGAGG(5MU)(PSU)CGAAUCCUCGUACCCCAGCCA",
            "pdbx_seq_one_letter_code_can": "UGGGGUAUCGCCAAGCGGUAAGGCACCGGAUUCUGAUUCCGGCAUUCCGAGGUUCGAAUCCUCGUACCCCAGCCA",
            "pdbx_strand_id": "T", "pdbx_target_identifier": "?",
        }, {
            "entity_id": "2", "type": "polypeptide(L)", "nstd_linkage": "no",
            "nstd_monomer": "no",
            "pdbx_seq_one_letter_code": "SEAEARPTNFIRQIIDEDLASGKHTTVHTRFPPEPNGYLHIGHAKSICLNFGIAQDYKGQCNLRFDDTNPVKEDIEYVESIKNDVEWLGFHWSGNVRYSSDYFDQLHAYAIELINKGLAYVDELTPEQIREYRGTLTQPGKNSPYRDRSVEENLALFEKMRAGGFEEGKACLRAKIDMASPFIVMRDPVLYRIKFAEHHQTGNKWCIYPMYDFTHCISDALEGITHSLCTLEFQDNRRLYDWVLDNITIPVHPRQYEFSRLNLEYTVMSKRKLNLLVTDKHVEGWDDPRMPTISGLRRRGYTAASIREFCKRIGVTKQDNTIEMASLESCIREDLNENAPRAMAVIDPVKLVIENYQGEGEMVTMPNHPNKPEMGSRQVPFSGEIWIDRADFREEANKQYKRLVLGKEVRLRNAYVIKAERVEKDAEGNITTIFCTYDADTLSKDPADGRKVKGVIHWVSAAHALPVEIRLYDRLFSVPNPGAADDFLSVINPESLVIKQGFAEPSLKDAVAGKAFQFEREGYFCLDSRHSTAEKPVFNRTVGLRDTWAKVGE",
            "pdbx_seq_one_letter_code_can": "SEAEARPTNFIRQIIDEDLASGKHTTVHTRFPPEPNGYLHIGHAKSICLNFGIAQDYKGQCNLRFDDTNPVKEDIEYVESIKNDVEWLGFHWSGNVRYSSDYFDQLHAYAIELINKGLAYVDELTPEQIREYRGTLTQPGKNSPYRDRSVEENLALFEKMRAGGFEEGKACLRAKIDMASPFIVMRDPVLYRIKFAEHHQTGNKWCIYPMYDFTHCISDALEGITHSLCTLEFQDNRRLYDWVLDNITIPVHPRQYEFSRLNLEYTVMSKRKLNLLVTDKHVEGWDDPRMPTISGLRRRGYTAASIREFCKRIGVTKQDNTIEMASLESCIREDLNENAPRAMAVIDPVKLVIENYQGEGEMVTMPNHPNKPEMGSRQVPFSGEIWIDRADFREEANKQYKRLVLGKEVRLRNAYVIKAERVEKDAEGNITTIFCTYDADTLSKDPADGRKVKGVIHWVSAAHALPVEIRLYDRLFSVPNPGAADDFLSVINPESLVIKQGFAEPSLKDAVAGKAFQFEREGYFCLDSRHSTAEKPVFNRTVGLRDTWAKVGE",
            "pdbx_strand_id": "P", "pdbx_target_identifier": "?",
        }])
        self.assertEqual(len(d["entity_poly_seq"]), 628)
        self.assertEqual(d["entity_poly_seq"][0], {
            "entity_id": "1", "num": "1", "mon_id": "U", "hetero": "n",
        })
        self.assertEqual(d["entity_poly_seq"][-1], {
            "entity_id": "2", "num": "553", "mon_id": "GLU", "hetero": "n",
        })
        self.assertEqual(d["struct_ref"], [{
            "id": "1", "db_name": "PDB", "db_code": "1GSG",
            "pdbx_db_accession": "1GSG", "pdbx_db_isoform": "?", "entity_id": "1",
            "pdbx_seq_one_letter_code": "?", "pdbx_align_begin": "?",
        }, {
            "id": "2", "db_name": "UNP", "db_code": "SYQ_ECOLI",
            "pdbx_db_accession": "P00962", "pdbx_db_isoform": "?", "entity_id": "2",
            "pdbx_seq_one_letter_code": "?", "pdbx_align_begin": "?",
        }])
        self.assertEqual(d["struct_ref_seq"], [{
            "align_id": "1", "ref_id": "2", "pdbx_PDB_id_code": "1GSG",
            "pdbx_strand_id": "P", "seq_align_beg": "?",
            "pdbx_seq_align_beg_ins_code": "?", "seq_align_end": "?",
            "pdbx_seq_align_end_ins_code": "?", "pdbx_db_accession": "P00962",
            "db_align_beg": "1", "pdbx_db_align_beg_ins_code": "?",
            "db_align_end": "553", "pdbx_db_align_end_ins_code": "?",
            "pdbx_auth_seq_align_beg": "1", "pdbx_auth_seq_align_end": "553",
        }, {
            "align_id": "2", "ref_id": "1", "pdbx_PDB_id_code": "1GSG",
            "pdbx_strand_id": "T", "seq_align_beg": "?",
            "pdbx_seq_align_beg_ins_code": "?", "seq_align_end": "?",
            "pdbx_seq_align_end_ins_code": "?", "pdbx_db_accession": "1GSG",
            "db_align_beg": "1", "pdbx_db_align_beg_ins_code": "?",
            "db_align_end": "76", "pdbx_db_align_end_ins_code": "?",
            "pdbx_auth_seq_align_beg": "1", "pdbx_auth_seq_align_end": "76",
        }])
        self.assertNotIn("struct_ref_seq_dif", d)
        self.assertEqual(len(d["pdbx_struct_mod_residue"]), 9)
        self.assertEqual(d["pdbx_struct_mod_residue"][0], {
            "id": "1", "label_asym_id": "A", "label_comp_id": "4SU",
            "label_seq_id": "8", "auth_asym_id": "T", "auth_comp_id": "4SU",
            "auth_seq_id": "8", "PDB_ins_code": "?", "parent_comp_id": "U",
            "details": "4-THIOURIDINE-5'-MONOPHOSPHATE",
        })
        self.assertEqual(d["pdbx_struct_mod_residue"][-1], {
            "id": "9", "label_asym_id": "A", "label_comp_id": "PSU",
            "label_seq_id": "54", "auth_asym_id": "T", "auth_comp_id": "PSU",
            "auth_seq_id": "55", "PDB_ins_code": "?", "parent_comp_id": "U",
            "details": "PSEUDOURIDINE-5'-MONOPHOSPHATE",
        })
        self.assertEqual(len(d["chem_comp"]), 31)
        self.assertEqual(d["chem_comp"][0], {
            "id": "2MA", "type": "L-peptide linking", "mon_nstd_flag": "n",
            "name": "2-METHYLADENOSINE-5'-MONOPHOSPHATE", "pdbx_synonyms": "?",
            "formula": "C11 H16 N5 O7 P", "formula_weight": "361.247",
        })
        self.assertEqual(d["chem_comp"][1], {
            "id": "4SU", "type": "L-peptide linking", "mon_nstd_flag": "n",
            "name": "4-THIOURIDINE-5'-MONOPHOSPHATE", "pdbx_synonyms": "?",
            "formula": "C9 H13 N2 O8 P S", "formula_weight": "340.246",
        })
        self.assertEqual(d["chem_comp"][2], {
            "id": "5MU", "type": "L-peptide linking", "mon_nstd_flag": "n",
            "name": "5-METHYLURIDINE 5'-MONOPHOSPHATE", "pdbx_synonyms": "?",
            "formula": "C10 H15 N2 O9 P", "formula_weight": "338.207",
        })
        self.assertEqual(d["chem_comp"][14], {
            "id": "H2U", "type": "L-peptide linking", "mon_nstd_flag": "n",
            "name": "5,6-DIHYDROURIDINE-5'-MONOPHOSPHATE", "pdbx_synonyms": "?",
            "formula": "C9 H15 N2 O9 P", "formula_weight": "326.197",
        })
        self.assertEqual(d["chem_comp"][20], {
            "id": "OMG", "type": "L-peptide linking", "mon_nstd_flag": "n",
            "name": "O2'-METHYLGUANOSINE-5'-MONOPHOSPHATE", "pdbx_synonyms": "?",
            "formula": "C11 H16 N5 O8 P", "formula_weight": "377.247",
        })
        self.assertEqual(d["chem_comp"][21], {
            "id": "OMU", "type": "L-peptide linking", "mon_nstd_flag": "n",
            "name": "O2'-METHYLURIDINE 5'-MONOPHOSPHATE", "pdbx_synonyms": "?",
            "formula": "C10 H15 N2 O9 P", "formula_weight": "338.207",
        })
        self.assertEqual(d["chem_comp"][24], {
            "id": "PSU", "type": "L-peptide linking", "mon_nstd_flag": "n",
            "name": "PSEUDOURIDINE-5'-MONOPHOSPHATE", "pdbx_synonyms": "?",
            "formula": "C9 H13 N2 O9 P", "formula_weight": "324.181",
        })
        self.assertEqual(d["chem_comp"][-1], {
            "id": "VAL", "type": "L-peptide linking", "mon_nstd_flag": "y",
            "name": "VALINE", "pdbx_synonyms": "?", "formula": "C5 H11 N O2",
            "formula_weight": "117.146",
        })
        self.assertNotIn("pdbx_entity_nonpoly", d)
        self.assertEqual(d["struct_asym"], [
            {"id": "A", "pdbx_blank_PDB_chainid_flag": "N", "pdbx_modified": "N", "entity_id": "1", "details": "?"},
            {"id": "B", "pdbx_blank_PDB_chainid_flag": "N", "pdbx_modified": "N", "entity_id": "2", "details": "?"},
        ])

        # Atom categories
        self.assertEqual(d["atom_type"], [{"symbol": "C"}, {"symbol": "P"}])
        self.assertEqual(len(d["atom_site"]), 601)
        self.assertEqual(d["atom_site"][0], {
            "group_PDB": "ATOM", "id": "1", "type_symbol": "P",
            "label_atom_id": "P", "label_alt_id": ".", "label_comp_id": "U",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "41.703", "Cartn_y": "40.142",
            "Cartn_z": "39.963", "occupancy": "1.00", "B_iso_or_equiv": "71.49",
            "pdbx_formal_charge": "?", "auth_seq_id": "1", "auth_comp_id": "U",
            "auth_asym_id": "T", "auth_atom_id": "P", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][75], {
            "group_PDB": "ATOM", "id": "76", "type_symbol": "C",
            "label_atom_id": "CA", "label_alt_id": ".", "label_comp_id": "THR",
            "label_asym_id": "B", "label_entity_id": "2", "label_seq_id": "8",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "22.968", "Cartn_y": "28.414",
            "Cartn_z": "13.890", "occupancy": "1.00", "B_iso_or_equiv": "34.56",
            "pdbx_formal_charge": "?", "auth_seq_id": "8", "auth_comp_id": "THR",
            "auth_asym_id": "P", "auth_atom_id": "CA", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][-1], {
            "group_PDB": "ATOM", "id": "601", "type_symbol": "C",
            "label_atom_id": "CA", "label_alt_id": ".", "label_comp_id": "THR",
            "label_asym_id": "B", "label_entity_id": "2", "label_seq_id": "547",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "26.659", "Cartn_y": "-2.959",
            "Cartn_z": "28.778", "occupancy": "1.00", "B_iso_or_equiv": "49.48",
            "pdbx_formal_charge": "?", "auth_seq_id": "547", "auth_comp_id": "THR",
            "auth_asym_id": "P", "auth_atom_id": "CA", "pdbx_PDB_model_num": "1",
        })
        self.assertNotIn("atom_site_anisotrop", d)

        # Annotation categories
        self.assertNotIn("struct_conf", d)
        self.assertNotIn("struct_sheet", d)
        self.assertNotIn("struct_sheet_range", d)
        self.assertNotIn("struct_sheet_order", d)
        self.assertNotIn("pdbx_struct_sheet_hbond", d)
        self.assertNotIn("struct_conn", d)
        self.assertNotIn("struct_mon_prot_cis", d)
        self.assertNotIn("struct_site", d)
        self.assertNotIn("struct_site_gen", d)

        # Crystal categories
        self.assertEqual(d["cell"], [{
            "entry_id": "1GSG", "length_a": "242.800", "length_b": "93.600",
            "length_c": "115.700", "angle_alpha": "90.00", "angle_beta": "90.00",
            "angle_gamma": "90.00", "Z_pdb": "8", "pdbx_unique_axis": "?",
        }])
        self.assertEqual(d["symmetry"], [{
            "entry_id": "1GSG", "space_group_name_H-M": "C 2 2 21",
            "pdbx_full_space_group_name_H-M": "?", "cell_setting": "?",
            "Int_Tables_number": "?",
        }])
        self.assertEqual(d["atom_sites"], [{
            "entry_id": "1GSG", "fract_transf_matrix[1][1]": "0.004119",
            "fract_transf_matrix[1][2]": "0.000000",
            "fract_transf_matrix[1][3]": "0.000000",
            "fract_transf_matrix[2][1]": "0.000000",
            "fract_transf_matrix[2][2]": "0.010684",
            "fract_transf_matrix[2][3]": "0.000000",
            "fract_transf_matrix[3][1]": "0.000000",
            "fract_transf_matrix[3][2]": "0.000000",
            "fract_transf_matrix[3][3]": "0.008643",
            "fract_transf_vector[1]": "0.00000",
            "fract_transf_vector[2]": "0.00000",
            "fract_transf_vector[3]": "0.00000",
        }])
        self.assertEqual(d["database_PDB_matrix"], [{
            "entry_id": "1GSG", "origx[1][1]": "1.000000",
            "origx[1][2]": "0.000000", "origx[1][3]": "0.000000",
            "origx[2][1]": "0.000000", "origx[2][2]": "1.000000",
            "origx[2][3]": "0.000000", "origx[3][1]": "0.000000",
            "origx[3][2]": "0.000000", "origx[3][3]": "1.000000",
            "origx_vector[1]": "0.00000", "origx_vector[2]": "0.00000",
            "origx_vector[3]": "0.00000",
        }])
        self.assertNotIn("struct_ncs_oper", d)

        # Citation categories
        self.assertEqual(d["audit_author"], [
            {"name": "Rould, M.A", "pdbx_ordinal": "1"},
            {"name": "Perona, J.J", "pdbx_ordinal": "2"},
            {"name": "Soell, D", "pdbx_ordinal": "3"},
            {"name": "Steitz, T.A", "pdbx_ordinal": "4"},
        ])
        self.assertEqual(d["citation_author"], [
            {"citation_id": "primary", "name": "Rould, M.A", "pdbx_ordinal": "1"},
            {"citation_id": "primary", "name": "Perona, J.J", "pdbx_ordinal": "2"},
            {"citation_id": "primary", "name": "Soll, D", "pdbx_ordinal": "3"},
            {"citation_id": "primary", "name": "Steitz, T.A", "pdbx_ordinal": "4"},
        ])
        self.assertNotIn("citation_editor", d)
        self.assertEqual(d["citation"], [{
            "id": "primary",
            "title": "STRUCTURE OF E. COLI GLUTAMINYL-TRNA SYNTHETASE COMPLEXED WITH TRNA(GLN) AND ATP AT 2.8 A RESOLUTION.",
            "journal_abbrev": "Science", "journal_volume": "246",
            "page_first": "1135", "page_last": "?", "year": "1989",
            "journal_id_ASTM": "?", "country": "?", "journal_id_ISSN": "0036-8075",
            "journal_id_CSD": "?", "book_publisher": "?",
            "pdbx_database_id_PubMed": "2479982", "pdbx_database_id_DOI": "?",
        }])

        # Experimental categories
        self.assertEqual(d["reflns"], [{
            "entry_id": "1GSG", "observed_criterion_sigma_I": "?",
            "observed_criterion_sigma_F": "?", "d_resolution_low": "?",
            "d_resolution_high": "2.80", "number_obs": "?", "number_all": "?",
            "percent_possible_obs": "?", "pdbx_Rmerge_I_obs": "?",
            "pdbx_Rsym_value": "?", "pdbx_netI_over_sigmaI": "?",
            "B_iso_Wilson_estimate": "?", "pdbx_redundancy": "?",
            "R_free_details": "?", "limit_h_max": "?", "limit_h_min": "?",
            "limit_k_max": "?", "limit_k_min": "?", "limit_l_max": "?",
            "limit_l_min": "?", "observed_criterion_F_max": "?",
            "observed_criterion_F_min": "?", "pdbx_chi_squared": "?",
            "pdbx_scaling_rejects": "?", "pdbx_diffrn_id": "?", "pdbx_ordinal": "?",
        }])
        self.assertEqual(d["refine"], [{
            "entry_id": "1GSG", "ls_number_reflns_obs": "?",
            "ls_number_reflns_all": "?", "pdbx_ls_sigma_I": "?",
            "pdbx_ls_sigma_F": "?", "pdbx_data_cutoff_high_absF": "?",
            "pdbx_data_cutoff_low_absF": "?", "ls_d_res_low": "7.00",
            "ls_d_res_high": "2.80", "ls_percent_reflns_obs": "?",
            "ls_R_factor_obs": "0.279", "ls_R_factor_all": "0.279",
            "ls_R_factor_R_work": "0.279", "ls_R_factor_R_free": "?",
            "ls_R_factor_R_free_error": "?",
            "ls_R_factor_R_free_error_details": "?",
            "ls_percent_reflns_R_free": "?", "ls_number_reflns_R_free": "?",
            "ls_number_parameters": "?", "ls_number_restraints": "?",
            "occupancy_min": "?", "occupancy_max": "?",
            "correlation_coeff_Fo_to_Fc": "?",
            "correlation_coeff_Fo_to_Fc_free": "?", "B_iso_mean": "?",
            "aniso_B[1][1]": "?", "aniso_B[2][2]": "?", "aniso_B[3][3]": "?",
            "aniso_B[1][2]": "?", "aniso_B[1][3]": "?", "aniso_B[2][3]": "?",
            "solvent_model_details": "?", "solvent_model_param_ksol": "?",
            "solvent_model_param_bsol": "?", "pdbx_solvent_vdw_probe_radii": "?",
            "pdbx_solvent_ion_probe_radii": "?",
            "pdbx_solvent_shrinkage_radii": "?", "pdbx_ls_cross_valid_method": "?",
            "details": "?", "pdbx_starting_model": "?",
            "pdbx_method_to_determine_struct": "?",
            "pdbx_isotropic_thermal_model": "?",
            "pdbx_stereochemistry_target_values": "?",
            "pdbx_stereochem_target_val_spec_case": "?",
            "pdbx_R_Free_selection_details": "?", "pdbx_overall_ESU_R_Free": "?",
            "overall_SU_B": "?", "ls_redundancy_reflns_obs": "?", "B_iso_min": "?",
            "B_iso_max": "?", "overall_SU_R_Cruickshank_DPI": "?",
            "overall_SU_R_free": "?", "overall_SU_ML": "?",
            "pdbx_overall_ESU_R": "?", "pdbx_data_cutoff_high_rms_absF": "?",
            "pdbx_refine_id": "?", "pdbx_overall_phase_error": "?",
            "ls_wR_factor_R_free": "?", "ls_wR_factor_R_work": "?",
            "overall_FOM_free_R_set": "?", "overall_FOM_work_R_set": "?",
            "pdbx_diffrn_id": "?", "pdbx_TLS_residual_ADP_flag": "?",
            "pdbx_overall_SU_R_free_Cruickshank_DPI": "?",
            "pdbx_overall_SU_R_Blow_DPI": "?",
            "pdbx_overall_SU_R_free_Blow_DPI": "?",
        }])

        # Missing items categories
        self.assertEqual(len(d["pdbx_unobs_or_zero_occ_residues"]), 27)
        self.assertEqual(d["pdbx_unobs_or_zero_occ_residues"][0], {
            "id": "1", "PDB_model_num": "1", "polymer_flag": "Y",
            "occupancy_flag": "1", "auth_asym_id": "P", "auth_comp_id": "SER",
            "auth_seq_id": "1", "PDB_ins_code": "?", "label_asym_id": "P",
            "label_comp_id": "SER", "label_seq_id": "1",
        })
        self.assertEqual(d["pdbx_unobs_or_zero_occ_residues"][-1], {
            "id": "27", "PDB_model_num": "1", "polymer_flag": "Y",
            "occupancy_flag": "1", "auth_asym_id": "P", "auth_comp_id": "GLU",
            "auth_seq_id": "553", "PDB_ins_code": "?", "label_asym_id": "P",
            "label_comp_id": "GLU", "label_seq_id": "553",
        })

        # Assembly categories
        self.assertEqual(d["pdbx_struct_assembly"], [
            {"id": "1", "details": "?", "method_details": "?", "oligomeric_details": "?", "oligomeric_count": "?"},
        ])
        self.assertEqual(d["pdbx_struct_assembly_gen"], [{"assembly_id": "1", "oper_expression": "1", "asym_id_list": "A,B"}])
        self.assertNotIn("pdbx_struct_assembly_prop", d)
        self.assertEqual(d["pdbx_struct_oper_list"], [{
            "id": "1", "type": "?", "name": "?", "symmetry_operation": "x,y,z",
            "matrix[1][1]": "1.0", "matrix[1][2]": "0.0", "matrix[1][3]": "0.0",
            "vector[1]": "0.0", "matrix[2][1]": "0.0", "matrix[2][2]": "1.0",
            "matrix[2][3]": "0.0", "vector[2]": "0.0", "matrix[3][1]": "0.0",
            "matrix[3][2]": "0.0", "matrix[3][3]": "1.0", "vector[3]": "0.0",
        }])


    def test_1xda(self):
        d = atomium.open("tests/integration/files/1xda.pdb", dictionary=True)
        self.assertEqual(len(d.keys()), 38)

        # Metadata categories
        self.assertEqual(d["entry"], [{"id": "1XDA"}])
        self.assertEqual(d["pdbx_database_status"], [
            {"status_code": "REL", "entry_id": "1XDA", "recvd_initial_deposition_date": "1996-12-18"},
        ])
        self.assertEqual(d["struct_keywords"], [{
            "entry_id": "1XDA", "pdbx_keywords": "HORMONE",
            "text": "HORMONE, METABOLIC ROLE, CHEMICAL ACTIVITY, INSULIN ALBUMIN, FATTY ACID, GLUCOSE METABOLISM, DIABETES",
        }])
        self.assertNotIn("pdbx_database_PDB_obs_spr", d)
        self.assertEqual(d["struct"], [{
            "entry_id": "1XDA", "title": "STRUCTURE OF INSULIN",
            "pdbx_descriptor": "?", "pdbx_model_details": "?",
            "pdbx_CASP_flag": "?", "pdbx_model_type_details": "?",
        }])
        self.assertNotIn("pdbx_database_related", d)
        self.assertNotIn("database_PDB_caveat", d)
        self.assertEqual(d["exptl"], [{"entry_id": "1XDA", "method": "X-RAY DIFFRACTION"}])
        self.assertNotIn("pdbx_coordinate_model", d)
        self.assertEqual(d["pdbx_audit_revision_history"], [{
            "ordinal": "1", "data_content_type": "Structure model",
            "major_revision": "1", "minor_revision": "0",
            "revision_date": "1997-07-07",
        }, {
            "ordinal": "2", "data_content_type": "Structure model",
            "major_revision": "2", "minor_revision": "0",
            "revision_date": "2009-02-24",
        }])

        # Entity categories
        self.assertEqual(len(d["entity"]), 7)
        self.assertEqual(d["entity"][0], {
            "id": "1", "type": "polymer", "src_method": "man",
            "pdbx_description": "FATTY ACID ACYLATED INSULIN",
            "formula_weight": "?", "pdbx_number_of_molecules": "4", "pdbx_ec": "?",
            "pdbx_mutation": "?", "pdbx_fragment": "?", "details": "?",
        })
        self.assertEqual(d["entity"][-1], {
            "id": "7", "type": "water", "src_method": "nat",
            "pdbx_description": "water", "formula_weight": "?",
            "pdbx_number_of_molecules": "154", "pdbx_ec": "?", "pdbx_mutation": "?",
            "pdbx_fragment": "?", "details": "?",
        })
        self.assertEqual(d["entity_name_com"], [
            {"entity_id": "1", "name": "NN304 INSULIN"},
            {"entity_id": "2", "name": "NN304 INSULIN"},
        ])
        self.assertEqual(d["entity_poly"], [{
            "entity_id": "1", "type": "polypeptide(L)", "nstd_linkage": "no",
            "nstd_monomer": "no",
            "pdbx_seq_one_letter_code": "GIVEQCCTSICSLYQLENYCN",
            "pdbx_seq_one_letter_code_can": "GIVEQCCTSICSLYQLENYCN",
            "pdbx_strand_id": "A,C,E,G", "pdbx_target_identifier": "?",
        }, {
            "entity_id": "2", "type": "polypeptide(L)", "nstd_linkage": "no",
            "nstd_monomer": "no",
            "pdbx_seq_one_letter_code": "FVNQHLCGSHLVEALYLVCGERGFFYTPK",
            "pdbx_seq_one_letter_code_can": "FVNQHLCGSHLVEALYLVCGERGFFYTPK",
            "pdbx_strand_id": "B,D,F,H", "pdbx_target_identifier": "?",
        }])
        self.assertEqual(len(d["entity_poly_seq"]), 50)
        self.assertEqual(d["entity_poly_seq"][0], {
            "entity_id": "1", "num": "1", "mon_id": "GLY", "hetero": "n",
        })
        self.assertEqual(d["entity_poly_seq"][-1], {
            "entity_id": "2", "num": "29", "mon_id": "LYS", "hetero": "n",
        })
        self.assertEqual(d["struct_ref"], [{
            "id": "1", "db_name": "UNP", "db_code": "INS_HUMAN",
            "pdbx_db_accession": "P01308", "pdbx_db_isoform": "?", "entity_id": "1",
            "pdbx_seq_one_letter_code": "?", "pdbx_align_begin": "?",
        }, {
            "id": "2", "db_name": "UNP", "db_code": "INS_HUMAN",
            "pdbx_db_accession": "P01308", "pdbx_db_isoform": "?", "entity_id": "2",
            "pdbx_seq_one_letter_code": "?", "pdbx_align_begin": "?",
        }])
        self.assertEqual(len(d["struct_ref_seq"]), 8)
        self.assertEqual(d["struct_ref_seq"][0], {
            "align_id": "1", "ref_id": "1", "pdbx_PDB_id_code": "1XDA",
            "pdbx_strand_id": "A", "seq_align_beg": "?",
            "pdbx_seq_align_beg_ins_code": "?", "seq_align_end": "?",
            "pdbx_seq_align_end_ins_code": "?", "pdbx_db_accession": "P01308",
            "db_align_beg": "90", "pdbx_db_align_beg_ins_code": "?",
            "db_align_end": "110", "pdbx_db_align_end_ins_code": "?",
            "pdbx_auth_seq_align_beg": "1", "pdbx_auth_seq_align_end": "21",
        })
        self.assertEqual(d["struct_ref_seq"][-1], {
            "align_id": "8", "ref_id": "1", "pdbx_PDB_id_code": "1XDA",
            "pdbx_strand_id": "H", "seq_align_beg": "?",
            "pdbx_seq_align_beg_ins_code": "?", "seq_align_end": "?",
            "pdbx_seq_align_end_ins_code": "?", "pdbx_db_accession": "P01308",
            "db_align_beg": "25", "pdbx_db_align_beg_ins_code": "?",
            "db_align_end": "53", "pdbx_db_align_end_ins_code": "?",
            "pdbx_auth_seq_align_beg": "1", "pdbx_auth_seq_align_end": "29",
        })
        self.assertNotIn("struct_ref_seq_dif", d)
        self.assertNotIn("pdbx_struct_mod_residue", d)
        self.assertEqual(len(d["chem_comp"]), 22)
        self.assertEqual(d["chem_comp"][0], {
            "id": "ALA", "type": "L-peptide linking", "mon_nstd_flag": "y",
            "name": "ALANINE", "pdbx_synonyms": "?", "formula": "C3 H7 N O2",
            "formula_weight": "89.093",
        })
        self.assertEqual(d["chem_comp"][3], {
            "id": "CL", "type": "non-polymer", "mon_nstd_flag": ".",
            "name": "CHLORIDE ION", "pdbx_synonyms": "?", "formula": "CL 1-",
            "formula_weight": "35.453",
        })
        self.assertEqual(d["chem_comp"][9], {
            "id": "HOH", "type": "non-polymer", "mon_nstd_flag": ".",
            "name": "WATER", "pdbx_synonyms": "?", "formula": "H2 O",
            "formula_weight": "18.015",
        })
        self.assertEqual(d["chem_comp"][11], {
            "id": "IPH", "type": "non-polymer", "mon_nstd_flag": ".",
            "name": "PHENOL", "pdbx_synonyms": "?", "formula": "C6 H6 O",
            "formula_weight": "94.111",
        })
        self.assertEqual(d["chem_comp"][14], {
            "id": "MYR", "type": "non-polymer", "mon_nstd_flag": ".",
            "name": "MYRISTIC ACID", "pdbx_synonyms": "?", "formula": "C14 H28 O2",
            "formula_weight": "228.370",
        })
        self.assertEqual(d["chem_comp"][-1], {
            "id": "ZN", "type": "non-polymer", "mon_nstd_flag": ".",
            "name": "ZINC ION", "pdbx_synonyms": "?", "formula": "ZN 2+",
            "formula_weight": "65.390",
        })
        self.assertEqual(d["pdbx_entity_nonpoly"], [
            {"entity_id": "3", "name": "ZINC ION", "comp_id": "ZN"},
            {"entity_id": "4", "name": "CHLORIDE ION", "comp_id": "CL"},
            {"entity_id": "5", "name": "PHENOL", "comp_id": "IPH"},
            {"entity_id": "6", "name": "MYRISTIC ACID", "comp_id": "MYR"},
            {"entity_id": "7", "name": "water", "comp_id": "HOH"},
        ])
        self.assertEqual(len(d["struct_asym"]), 32)
        self.assertEqual(d["struct_asym"][0], {
            "id": "A", "pdbx_blank_PDB_chainid_flag": "N", "pdbx_modified": "N",
            "entity_id": "1", "details": "?",
        })
        self.assertEqual(d["struct_asym"][-1], {
            "id": "FA", "pdbx_blank_PDB_chainid_flag": "N", "pdbx_modified": "N",
            "entity_id": "7", "details": "?",
        })

        # Atom categories
        self.assertEqual(d["atom_type"], [
            {"symbol": "C"},
            {"symbol": "CL"},
            {"symbol": "N"},
            {"symbol": "O"},
            {"symbol": "S"},
            {"symbol": "ZN"},
        ])
        self.assertEqual(len(d["atom_site"]), 1860)
        self.assertEqual(d["atom_site"][0], {
            "group_PDB": "ATOM", "id": "1", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "GLY",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "16.136", "Cartn_y": "11.725",
            "Cartn_z": "12.834", "occupancy": "1.00", "B_iso_or_equiv": "17.65",
            "pdbx_formal_charge": "?", "auth_seq_id": "1", "auth_comp_id": "GLY",
            "auth_asym_id": "A", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][163], {
            "group_PDB": "ATOM", "id": "164", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "PHE",
            "label_asym_id": "B", "label_entity_id": "2", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-0.033", "Cartn_y": "11.515",
            "Cartn_z": "9.458", "occupancy": "1.00", "B_iso_or_equiv": "19.34",
            "pdbx_formal_charge": "?", "auth_seq_id": "1", "auth_comp_id": "PHE",
            "auth_asym_id": "B", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][405], {
            "group_PDB": "ATOM", "id": "406", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "GLY",
            "label_asym_id": "C", "label_entity_id": "1", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "9.561", "Cartn_y": "16.902",
            "Cartn_z": "40.400", "occupancy": "1.00", "B_iso_or_equiv": "25.15",
            "pdbx_formal_charge": "?", "auth_seq_id": "1", "auth_comp_id": "GLY",
            "auth_asym_id": "C", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][568], {
            "group_PDB": "ATOM", "id": "569", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "PHE",
            "label_asym_id": "D", "label_entity_id": "2", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "11.900", "Cartn_y": "0.573",
            "Cartn_z": "42.957", "occupancy": "1.00", "B_iso_or_equiv": "26.70",
            "pdbx_formal_charge": "?", "auth_seq_id": "1", "auth_comp_id": "PHE",
            "auth_asym_id": "D", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][806], {
            "group_PDB": "ATOM", "id": "807", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "GLY",
            "label_asym_id": "E", "label_entity_id": "1", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "17.335", "Cartn_y": "9.341",
            "Cartn_z": "52.401", "occupancy": "1.00", "B_iso_or_equiv": "19.16",
            "pdbx_formal_charge": "?", "auth_seq_id": "1", "auth_comp_id": "GLY",
            "auth_asym_id": "E", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][969], {
            "group_PDB": "ATOM", "id": "970", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "PHE",
            "label_asym_id": "F", "label_entity_id": "2", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "1.175", "Cartn_y": "11.882",
            "Cartn_z": "49.858", "occupancy": "1.00", "B_iso_or_equiv": "28.70",
            "pdbx_formal_charge": "?", "auth_seq_id": "1", "auth_comp_id": "PHE",
            "auth_asym_id": "F", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1207], {
            "group_PDB": "ATOM", "id": "1208", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "GLY",
            "label_asym_id": "G", "label_entity_id": "1", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "12.506", "Cartn_y": "15.301",
            "Cartn_z": "80.111", "occupancy": "1.00", "B_iso_or_equiv": "22.75",
            "pdbx_formal_charge": "?", "auth_seq_id": "1", "auth_comp_id": "GLY",
            "auth_asym_id": "G", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1373], {
            "group_PDB": "ATOM", "id": "1374", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "PHE",
            "label_asym_id": "H", "label_entity_id": "2", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "11.577", "Cartn_y": "-0.737",
            "Cartn_z": "83.108", "occupancy": "1.00", "B_iso_or_equiv": "17.85",
            "pdbx_formal_charge": "?", "auth_seq_id": "1", "auth_comp_id": "PHE",
            "auth_asym_id": "H", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1610], {
            "group_PDB": "HETATM", "id": "1611", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "IPH",
            "label_asym_id": "I", "label_entity_id": "5", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-5.661", "Cartn_y": "9.282",
            "Cartn_z": "16.887", "occupancy": "1.00", "B_iso_or_equiv": "13.71",
            "pdbx_formal_charge": "?", "auth_seq_id": "22", "auth_comp_id": "IPH",
            "auth_asym_id": "A", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1617], {
            "group_PDB": "HETATM", "id": "1618", "type_symbol": "ZN",
            "label_atom_id": "ZN", "label_alt_id": ".", "label_comp_id": "ZN",
            "label_asym_id": "J", "label_entity_id": "3", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "0.000", "Cartn_y": "0.000",
            "Cartn_z": "18.779", "occupancy": "0.33", "B_iso_or_equiv": "10.69",
            "pdbx_formal_charge": "?", "auth_seq_id": "30", "auth_comp_id": "ZN",
            "auth_asym_id": "B", "auth_atom_id": "ZN", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1618], {
            "group_PDB": "HETATM", "id": "1619", "type_symbol": "CL",
            "label_atom_id": "CL", "label_alt_id": ".", "label_comp_id": "CL",
            "label_asym_id": "K", "label_entity_id": "4", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "0.000", "Cartn_y": "0.000",
            "Cartn_z": "16.611", "occupancy": "0.33", "B_iso_or_equiv": "12.22",
            "pdbx_formal_charge": "?", "auth_seq_id": "31", "auth_comp_id": "CL",
            "auth_asym_id": "B", "auth_atom_id": "CL", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1619], {
            "group_PDB": "HETATM", "id": "1620", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "MYR",
            "label_asym_id": "L", "label_entity_id": "6", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "7.835", "Cartn_y": "13.475",
            "Cartn_z": "11.645", "occupancy": "0.00", "B_iso_or_equiv": "20.00",
            "pdbx_formal_charge": "?", "auth_seq_id": "39", "auth_comp_id": "MYR",
            "auth_asym_id": "B", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1634], {
            "group_PDB": "HETATM", "id": "1635", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "IPH",
            "label_asym_id": "M", "label_entity_id": "5", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "9.806", "Cartn_y": "-4.594",
            "Cartn_z": "35.946", "occupancy": "1.00", "B_iso_or_equiv": "18.13",
            "pdbx_formal_charge": "?", "auth_seq_id": "22", "auth_comp_id": "IPH",
            "auth_asym_id": "C", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1641], {
            "group_PDB": "HETATM", "id": "1642", "type_symbol": "ZN",
            "label_atom_id": "ZN", "label_alt_id": ".", "label_comp_id": "ZN",
            "label_asym_id": "N", "label_entity_id": "3", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "0.000", "Cartn_y": "0.000",
            "Cartn_z": "34.196", "occupancy": "0.33", "B_iso_or_equiv": "12.34",
            "pdbx_formal_charge": "?", "auth_seq_id": "30", "auth_comp_id": "ZN",
            "auth_asym_id": "D", "auth_atom_id": "ZN", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1642], {
            "group_PDB": "HETATM", "id": "1643", "type_symbol": "CL",
            "label_atom_id": "CL", "label_alt_id": ".", "label_comp_id": "CL",
            "label_asym_id": "O", "label_entity_id": "4", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "0.000", "Cartn_y": "0.000",
            "Cartn_z": "36.366", "occupancy": "0.33", "B_iso_or_equiv": "11.37",
            "pdbx_formal_charge": "?", "auth_seq_id": "31", "auth_comp_id": "CL",
            "auth_asym_id": "D", "auth_atom_id": "CL", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1643], {
            "group_PDB": "HETATM", "id": "1644", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "MYR",
            "label_asym_id": "P", "label_entity_id": "6", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "13.145", "Cartn_y": "5.879",
            "Cartn_z": "44.308", "occupancy": "1.00", "B_iso_or_equiv": "43.86",
            "pdbx_formal_charge": "?", "auth_seq_id": "39", "auth_comp_id": "MYR",
            "auth_asym_id": "D", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1658], {
            "group_PDB": "HETATM", "id": "1659", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "IPH",
            "label_asym_id": "Q", "label_entity_id": "5", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-4.433", "Cartn_y": "9.971",
            "Cartn_z": "56.719", "occupancy": "1.00", "B_iso_or_equiv": "15.57",
            "pdbx_formal_charge": "?", "auth_seq_id": "22", "auth_comp_id": "IPH",
            "auth_asym_id": "E", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1665], {
            "group_PDB": "HETATM", "id": "1666", "type_symbol": "ZN",
            "label_atom_id": "ZN", "label_alt_id": ".", "label_comp_id": "ZN",
            "label_asym_id": "R", "label_entity_id": "3", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "0.000", "Cartn_y": "0.000",
            "Cartn_z": "58.639", "occupancy": "0.33", "B_iso_or_equiv": "11.99",
            "pdbx_formal_charge": "?", "auth_seq_id": "30", "auth_comp_id": "ZN",
            "auth_asym_id": "F", "auth_atom_id": "ZN", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1666], {
            "group_PDB": "HETATM", "id": "1667", "type_symbol": "CL",
            "label_atom_id": "CL", "label_alt_id": ".", "label_comp_id": "CL",
            "label_asym_id": "S", "label_entity_id": "4", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "0.000", "Cartn_y": "0.000",
            "Cartn_z": "56.498", "occupancy": "0.33", "B_iso_or_equiv": "13.93",
            "pdbx_formal_charge": "?", "auth_seq_id": "31", "auth_comp_id": "CL",
            "auth_asym_id": "F", "auth_atom_id": "CL", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1667], {
            "group_PDB": "HETATM", "id": "1668", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "MYR",
            "label_asym_id": "T", "label_entity_id": "6", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "12.146", "Cartn_y": "13.183",
            "Cartn_z": "45.515", "occupancy": "0.00", "B_iso_or_equiv": "20.00",
            "pdbx_formal_charge": "?", "auth_seq_id": "39", "auth_comp_id": "MYR",
            "auth_asym_id": "F", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1682], {
            "group_PDB": "HETATM", "id": "1683", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "IPH",
            "label_asym_id": "U", "label_entity_id": "5", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "8.823", "Cartn_y": "-6.017",
            "Cartn_z": "75.855", "occupancy": "1.00", "B_iso_or_equiv": "15.74",
            "pdbx_formal_charge": "?", "auth_seq_id": "22", "auth_comp_id": "IPH",
            "auth_asym_id": "G", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1689], {
            "group_PDB": "HETATM", "id": "1690", "type_symbol": "ZN",
            "label_atom_id": "ZN", "label_alt_id": ".", "label_comp_id": "ZN",
            "label_asym_id": "V", "label_entity_id": "3", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "0.000", "Cartn_y": "0.000",
            "Cartn_z": "74.024", "occupancy": "0.33", "B_iso_or_equiv": "14.68",
            "pdbx_formal_charge": "?", "auth_seq_id": "30", "auth_comp_id": "ZN",
            "auth_asym_id": "H", "auth_atom_id": "ZN", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1690], {
            "group_PDB": "HETATM", "id": "1691", "type_symbol": "CL",
            "label_atom_id": "CL", "label_alt_id": ".", "label_comp_id": "CL",
            "label_asym_id": "W", "label_entity_id": "4", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "0.000", "Cartn_y": "0.000",
            "Cartn_z": "76.191", "occupancy": "0.33", "B_iso_or_equiv": "15.86",
            "pdbx_formal_charge": "?", "auth_seq_id": "31", "auth_comp_id": "CL",
            "auth_asym_id": "H", "auth_atom_id": "CL", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1691], {
            "group_PDB": "HETATM", "id": "1692", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "MYR",
            "label_asym_id": "X", "label_entity_id": "6", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "14.782", "Cartn_y": "6.834",
            "Cartn_z": "84.170", "occupancy": "1.00", "B_iso_or_equiv": "43.17",
            "pdbx_formal_charge": "?", "auth_seq_id": "39", "auth_comp_id": "MYR",
            "auth_asym_id": "H", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1706], {
            "group_PDB": "HETATM", "id": "1707", "type_symbol": "O",
            "label_atom_id": "O", "label_alt_id": ".", "label_comp_id": "HOH",
            "label_asym_id": "Y", "label_entity_id": "7", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "17.639", "Cartn_y": "-5.979",
            "Cartn_z": "20.330", "occupancy": "1.00", "B_iso_or_equiv": "27.85",
            "pdbx_formal_charge": "?", "auth_seq_id": "23", "auth_comp_id": "HOH",
            "auth_asym_id": "A", "auth_atom_id": "O", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1721], {
            "group_PDB": "HETATM", "id": "1722", "type_symbol": "O",
            "label_atom_id": "O", "label_alt_id": ".", "label_comp_id": "HOH",
            "label_asym_id": "Z", "label_entity_id": "7", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "13.817", "Cartn_y": "11.255",
            "Cartn_z": "20.074", "occupancy": "1.00", "B_iso_or_equiv": "18.80",
            "pdbx_formal_charge": "?", "auth_seq_id": "40", "auth_comp_id": "HOH",
            "auth_asym_id": "B", "auth_atom_id": "O", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1746], {
            "group_PDB": "HETATM", "id": "1747", "type_symbol": "O",
            "label_atom_id": "O", "label_alt_id": ".", "label_comp_id": "HOH",
            "label_asym_id": "AA", "label_entity_id": "7", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "10.432", "Cartn_y": "18.286",
            "Cartn_z": "28.169", "occupancy": "1.00", "B_iso_or_equiv": "17.65",
            "pdbx_formal_charge": "?", "auth_seq_id": "42", "auth_comp_id": "HOH",
            "auth_asym_id": "C", "auth_atom_id": "O", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1762], {
            "group_PDB": "HETATM", "id": "1763", "type_symbol": "O",
            "label_atom_id": "O", "label_alt_id": ".", "label_comp_id": "HOH",
            "label_asym_id": "BA", "label_entity_id": "7", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-1.686", "Cartn_y": "4.748",
            "Cartn_z": "29.405", "occupancy": "1.00", "B_iso_or_equiv": "28.76",
            "pdbx_formal_charge": "?", "auth_seq_id": "40", "auth_comp_id": "HOH",
            "auth_asym_id": "D", "auth_atom_id": "O", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1779], {
            "group_PDB": "HETATM", "id": "1780", "type_symbol": "O",
            "label_atom_id": "O", "label_alt_id": ".", "label_comp_id": "HOH",
            "label_asym_id": "CA", "label_entity_id": "7", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "19.132", "Cartn_y": "10.418",
            "Cartn_z": "56.864", "occupancy": "1.00", "B_iso_or_equiv": "17.46",
            "pdbx_formal_charge": "?", "auth_seq_id": "23", "auth_comp_id": "HOH",
            "auth_asym_id": "E", "auth_atom_id": "O", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1794], {
            "group_PDB": "HETATM", "id": "1795", "type_symbol": "O",
            "label_atom_id": "O", "label_alt_id": ".", "label_comp_id": "HOH",
            "label_asym_id": "DA", "label_entity_id": "7", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "14.151", "Cartn_y": "4.653",
            "Cartn_z": "70.826", "occupancy": "1.00", "B_iso_or_equiv": "18.08",
            "pdbx_formal_charge": "?", "auth_seq_id": "40", "auth_comp_id": "HOH",
            "auth_asym_id": "F", "auth_atom_id": "O", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1810], {
            "group_PDB": "HETATM", "id": "1811", "type_symbol": "O",
            "label_atom_id": "O", "label_alt_id": ".", "label_comp_id": "HOH",
            "label_asym_id": "EA", "label_entity_id": "7", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "5.232", "Cartn_y": "20.030",
            "Cartn_z": "75.975", "occupancy": "1.00", "B_iso_or_equiv": "25.28",
            "pdbx_formal_charge": "?", "auth_seq_id": "23", "auth_comp_id": "HOH",
            "auth_asym_id": "G", "auth_atom_id": "O", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1840], {
            "group_PDB": "HETATM", "id": "1841", "type_symbol": "O",
            "label_atom_id": "O", "label_alt_id": ".", "label_comp_id": "HOH",
            "label_asym_id": "FA", "label_entity_id": "7", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "12.513", "Cartn_y": "14.380",
            "Cartn_z": "75.003", "occupancy": "1.00", "B_iso_or_equiv": "15.41",
            "pdbx_formal_charge": "?", "auth_seq_id": "40", "auth_comp_id": "HOH",
            "auth_asym_id": "H", "auth_atom_id": "O", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][-1], {
            "group_PDB": "HETATM", "id": "1860", "type_symbol": "O",
            "label_atom_id": "O", "label_alt_id": ".", "label_comp_id": "HOH",
            "label_asym_id": "FA", "label_entity_id": "7", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "11.353", "Cartn_y": "5.616",
            "Cartn_z": "75.248", "occupancy": "1.00", "B_iso_or_equiv": "33.51",
            "pdbx_formal_charge": "?", "auth_seq_id": "59", "auth_comp_id": "HOH",
            "auth_asym_id": "H", "auth_atom_id": "O", "pdbx_PDB_model_num": "1",
        })
        self.assertNotIn("atom_site_anisotrop", d)

        # Annotation categories
        self.assertEqual(len(d["struct_conf"]), 12)
        self.assertEqual(d["struct_conf"][0], {
            "conf_type_id": "HELX_P", "id": "HELX_P1", "pdbx_PDB_helix_id": "1",
            "beg_label_comp_id": "ILE", "beg_label_asym_id": "A",
            "beg_label_seq_id": "2", "pdbx_beg_PDB_ins_code": "?",
            "end_label_comp_id": "CYS", "end_label_asym_id": "A",
            "end_label_seq_id": "6", "pdbx_end_PDB_ins_code": "?",
            "beg_auth_comp_id": "ILE", "beg_auth_asym_id": "A",
            "beg_auth_seq_id": "2", "end_auth_comp_id": "CYS",
            "end_auth_asym_id": "A", "end_auth_seq_id": "6",
            "pdbx_PDB_helix_class": "1", "details": "?",
            "pdbx_PDB_helix_length": "5",
        })
        self.assertEqual(d["struct_conf"][-1], {
            "conf_type_id": "HELX_P", "id": "HELX_P12", "pdbx_PDB_helix_id": "12",
            "beg_label_comp_id": "VAL", "beg_label_asym_id": "H",
            "beg_label_seq_id": "2", "pdbx_beg_PDB_ins_code": "?",
            "end_label_comp_id": "ARG", "end_label_asym_id": "H",
            "end_label_seq_id": "22", "pdbx_end_PDB_ins_code": "?",
            "beg_auth_comp_id": "VAL", "beg_auth_asym_id": "H",
            "beg_auth_seq_id": "2", "end_auth_comp_id": "ARG",
            "end_auth_asym_id": "H", "end_auth_seq_id": "22",
            "pdbx_PDB_helix_class": "1", "details": "?",
            "pdbx_PDB_helix_length": "21",
        })
        self.assertEqual(d["struct_sheet"], [
            {"id": "A", "type": "?", "number_strands": "2", "details": "?"},
            {"id": "B", "type": "?", "number_strands": "2", "details": "?"},
        ])
        self.assertEqual(len(d["struct_sheet_range"]), 4)
        self.assertEqual(d["struct_sheet_range"][0], {
            "sheet_id": "A", "id": "1", "beg_label_comp_id": "PHE",
            "beg_label_asym_id": "B", "beg_label_seq_id": "24",
            "pdbx_beg_PDB_ins_code": "?", "end_label_comp_id": "TYR",
            "end_label_asym_id": "B", "end_label_seq_id": "26",
            "pdbx_end_PDB_ins_code": "?", "beg_auth_comp_id": "PHE",
            "beg_auth_asym_id": "B", "beg_auth_seq_id": "24",
            "end_auth_comp_id": "TYR", "end_auth_asym_id": "B",
            "end_auth_seq_id": "26",
        })
        self.assertEqual(d["struct_sheet_range"][-1], {
            "sheet_id": "B", "id": "2", "beg_label_comp_id": "PHE",
            "beg_label_asym_id": "H", "beg_label_seq_id": "24",
            "pdbx_beg_PDB_ins_code": "?", "end_label_comp_id": "TYR",
            "end_label_asym_id": "H", "end_label_seq_id": "26",
            "pdbx_end_PDB_ins_code": "?", "beg_auth_comp_id": "PHE",
            "beg_auth_asym_id": "H", "beg_auth_seq_id": "24",
            "end_auth_comp_id": "TYR", "end_auth_asym_id": "H",
            "end_auth_seq_id": "26",
        })
        self.assertEqual(d["struct_sheet_order"], [
            {"sheet_id": "A", "range_id_1": "1", "range_id_2": "2", "offset": "?", "sense": "anti-parallel"},
            {"sheet_id": "B", "range_id_1": "1", "range_id_2": "2", "offset": "?", "sense": "anti-parallel"},
        ])
        self.assertEqual(d["pdbx_struct_sheet_hbond"], [{
            "sheet_id": "A", "range_id_1": "1", "range_id_2": "2",
            "range_1_label_atom_id": "O", "range_1_label_comp_id": "PHE",
            "range_1_label_asym_id": "B", "range_1_label_seq_id": "24",
            "range_1_PDB_ins_code": "?", "range_1_auth_atom_id": "O",
            "range_1_auth_comp_id": "PHE", "range_1_auth_asym_id": "B",
            "range_1_auth_seq_id": "24", "range_2_label_atom_id": "N",
            "range_2_label_comp_id": "TYR", "range_2_label_asym_id": "D",
            "range_2_label_seq_id": "26", "range_2_PDB_ins_code": "?",
            "range_2_auth_atom_id": "N", "range_2_auth_comp_id": "TYR",
            "range_2_auth_asym_id": "D", "range_2_auth_seq_id": "26",
        }, {
            "sheet_id": "B", "range_id_1": "1", "range_id_2": "2",
            "range_1_label_atom_id": "O", "range_1_label_comp_id": "PHE",
            "range_1_label_asym_id": "F", "range_1_label_seq_id": "24",
            "range_1_PDB_ins_code": "?", "range_1_auth_atom_id": "O",
            "range_1_auth_comp_id": "PHE", "range_1_auth_asym_id": "F",
            "range_1_auth_seq_id": "24", "range_2_label_atom_id": "N",
            "range_2_label_comp_id": "TYR", "range_2_label_asym_id": "H",
            "range_2_label_seq_id": "26", "range_2_PDB_ins_code": "?",
            "range_2_auth_atom_id": "N", "range_2_auth_comp_id": "TYR",
            "range_2_auth_asym_id": "H", "range_2_auth_seq_id": "26",
        }])
        self.assertEqual(len(d["struct_conn"]), 41)
        self.assertEqual(d["struct_conn"][0], {
            "id": "disulf1", "conn_type_id": "disulf",
            "pdbx_leaving_atom_flag": "?", "pdbx_PDB_id": "?",
            "ptnr1_label_asym_id": "A", "ptnr1_label_comp_id": "CYS",
            "ptnr1_label_seq_id": "6", "ptnr1_label_atom_id": "?",
            "pdbx_ptnr1_label_alt_id": "?", "pdbx_ptnr1_PDB_ins_code": "?",
            "pdbx_ptnr1_standard_comp_id": "?", "ptnr1_symmetry": "1555",
            "ptnr2_label_asym_id": "A", "ptnr2_label_comp_id": "CYS",
            "ptnr2_label_seq_id": "11", "ptnr2_label_atom_id": "?",
            "pdbx_ptnr2_label_alt_id": "?", "pdbx_ptnr2_PDB_ins_code": "?",
            "ptnr1_auth_asym_id": "A", "ptnr1_auth_comp_id": "CYS",
            "ptnr1_auth_seq_id": "6", "ptnr2_auth_asym_id": "A",
            "ptnr2_auth_comp_id": "CYS", "ptnr2_auth_seq_id": "11",
            "ptnr2_symmetry": "1555", "pdbx_ptnr3_label_atom_id": "?",
            "pdbx_ptnr3_label_seq_id": "?", "pdbx_ptnr3_label_comp_id": "?",
            "pdbx_ptnr3_label_asym_id": "?", "pdbx_ptnr3_label_alt_id": "?",
            "pdbx_ptnr3_PDB_ins_code": "?", "details": "?",
            "pdbx_dist_value": "2.02", "pdbx_value_order": "?",
        })
        self.assertEqual(d["struct_conn"][-1], {
            "id": "covale29", "conn_type_id": "covale",
            "pdbx_leaving_atom_flag": "?", "pdbx_PDB_id": "?",
            "ptnr1_label_asym_id": "V", "ptnr1_label_comp_id": "ZN",
            "ptnr1_label_seq_id": ".", "ptnr1_label_atom_id": "ZN",
            "pdbx_ptnr1_label_alt_id": "?", "pdbx_ptnr1_PDB_ins_code": "?",
            "pdbx_ptnr1_standard_comp_id": "?", "ptnr1_symmetry": "1555",
            "ptnr2_label_asym_id": "H", "ptnr2_label_comp_id": "HIS",
            "ptnr2_label_seq_id": "10", "ptnr2_label_atom_id": "NE2",
            "pdbx_ptnr2_label_alt_id": "?", "pdbx_ptnr2_PDB_ins_code": "?",
            "ptnr1_auth_asym_id": "H", "ptnr1_auth_comp_id": "ZN",
            "ptnr1_auth_seq_id": "30", "ptnr2_auth_asym_id": "H",
            "ptnr2_auth_comp_id": "HIS", "ptnr2_auth_seq_id": "10",
            "ptnr2_symmetry": "2555", "pdbx_ptnr3_label_atom_id": "?",
            "pdbx_ptnr3_label_seq_id": "?", "pdbx_ptnr3_label_comp_id": "?",
            "pdbx_ptnr3_label_asym_id": "?", "pdbx_ptnr3_label_alt_id": "?",
            "pdbx_ptnr3_PDB_ins_code": "?", "details": "?",
            "pdbx_dist_value": "1.99", "pdbx_value_order": "?",
        })
        self.assertNotIn("struct_mon_prot_cis", d)
        self.assertEqual(len(d["struct_site"]), 16)
        self.assertEqual(d["struct_site"][0], {
            "id": "AC1", "pdbx_evidence_code": "Software", "pdbx_auth_asym_id": "?",
            "pdbx_auth_comp_id": "?", "pdbx_auth_seq_id": "?",
            "pdbx_auth_ins_code": "?", "pdbx_num_residues": "?",
            "details": "BINDING SITE FOR RESIDUE ZN B 30",
        })
        self.assertEqual(d["struct_site"][-1], {
            "id": "BC7", "pdbx_evidence_code": "Software", "pdbx_auth_asym_id": "?",
            "pdbx_auth_comp_id": "?", "pdbx_auth_seq_id": "?",
            "pdbx_auth_ins_code": "?", "pdbx_num_residues": "?",
            "details": "BINDING SITE FOR RESIDUE MYR H 39",
        })
        self.assertEqual(len(d["struct_site_gen"]), 55)
        self.assertEqual(d["struct_site_gen"][0], {
            "id": "1", "site_id": "AC1", "pdbx_num_res": "2",
            "label_comp_id": "HIS", "label_asym_id": "B", "label_seq_id": "10",
            "pdbx_auth_ins_code": "?", "auth_comp_id": "HIS", "auth_asym_id": "B",
            "auth_seq_id": "10", "label_atom_id": ".", "label_alt_id": "?",
            "symmetry": "1_555", "details": "?",
        })
        self.assertEqual(d["struct_site_gen"][-1], {
            "id": "55", "site_id": "BC7", "pdbx_num_res": "5",
            "label_comp_id": "GLN", "label_asym_id": "H", "label_seq_id": "4",
            "pdbx_auth_ins_code": "?", "auth_comp_id": "GLN", "auth_asym_id": "H",
            "auth_seq_id": "4", "label_atom_id": ".", "label_alt_id": "?",
            "symmetry": "1_555", "details": "?",
        })

        # Crystal categories
        self.assertEqual(d["cell"], [{
            "entry_id": "1XDA", "length_a": "78.752", "length_b": "78.752",
            "length_c": "79.199", "angle_alpha": "90.00", "angle_beta": "90.00",
            "angle_gamma": "120.00", "Z_pdb": "36", "pdbx_unique_axis": "?",
        }])
        self.assertEqual(d["symmetry"], [{
            "entry_id": "1XDA", "space_group_name_H-M": "H 3",
            "pdbx_full_space_group_name_H-M": "?", "cell_setting": "?",
            "Int_Tables_number": "?",
        }])
        self.assertEqual(d["atom_sites"], [{
            "entry_id": "1XDA", "fract_transf_matrix[1][1]": "0.012698",
            "fract_transf_matrix[1][2]": "0.007331",
            "fract_transf_matrix[1][3]": "0.000000",
            "fract_transf_matrix[2][1]": "0.000000",
            "fract_transf_matrix[2][2]": "0.014662",
            "fract_transf_matrix[2][3]": "0.000000",
            "fract_transf_matrix[3][1]": "0.000000",
            "fract_transf_matrix[3][2]": "0.000000",
            "fract_transf_matrix[3][3]": "0.012626",
            "fract_transf_vector[1]": "0.00000",
            "fract_transf_vector[2]": "0.00000",
            "fract_transf_vector[3]": "0.00000",
        }])
        self.assertEqual(d["database_PDB_matrix"], [{
            "entry_id": "1XDA", "origx[1][1]": "1.000000",
            "origx[1][2]": "0.000000", "origx[1][3]": "0.000000",
            "origx[2][1]": "0.000000", "origx[2][2]": "1.000000",
            "origx[2][3]": "0.000000", "origx[3][1]": "0.000000",
            "origx[3][2]": "0.000000", "origx[3][3]": "1.000000",
            "origx_vector[1]": "0.00000", "origx_vector[2]": "0.00000",
            "origx_vector[3]": "0.00000",
        }])
        self.assertNotIn("struct_ncs_oper", d)

        # Citation categories
        self.assertEqual(d["audit_author"], [
            {"name": "Whittingham, J.L", "pdbx_ordinal": "1"},
            {"name": "Havelund, S", "pdbx_ordinal": "2"},
            {"name": "Jonassen, I", "pdbx_ordinal": "3"},
        ])
        self.assertEqual(d["citation_author"], [
            {"citation_id": "primary", "name": "Whittingham, J.L", "pdbx_ordinal": "1"},
            {"citation_id": "primary", "name": "Havelund, S", "pdbx_ordinal": "2"},
            {"citation_id": "primary", "name": "Jonassen, I", "pdbx_ordinal": "3"},
        ])
        self.assertNotIn("citation_editor", d)
        self.assertEqual(d["citation"], [{
            "id": "primary",
            "title": "CRYSTAL STRUCTURE OF A PROLONGED-ACTING INSULIN WITH ALBUMIN-BINDING PROPERTIES.",
            "journal_abbrev": "Biochemistry", "journal_volume": "36",
            "page_first": "2826", "page_last": "?", "year": "1997",
            "journal_id_ASTM": "?", "country": "?", "journal_id_ISSN": "0006-2960",
            "journal_id_CSD": "?", "book_publisher": "?",
            "pdbx_database_id_PubMed": "9062110",
            "pdbx_database_id_DOI": "10.1021/bi9625105",
        }])

        # Experimental categories
        self.assertEqual(d["reflns"], [{
            "entry_id": "1XDA", "observed_criterion_sigma_I": "?",
            "observed_criterion_sigma_F": "?", "d_resolution_low": "?",
            "d_resolution_high": "1.80", "number_obs": "?", "number_all": "?",
            "percent_possible_obs": "?", "pdbx_Rmerge_I_obs": "?",
            "pdbx_Rsym_value": "?", "pdbx_netI_over_sigmaI": "?",
            "B_iso_Wilson_estimate": "17.30", "pdbx_redundancy": "?",
            "R_free_details": "?", "limit_h_max": "?", "limit_h_min": "?",
            "limit_k_max": "?", "limit_k_min": "?", "limit_l_max": "?",
            "limit_l_min": "?", "observed_criterion_F_max": "?",
            "observed_criterion_F_min": "?", "pdbx_chi_squared": "?",
            "pdbx_scaling_rejects": "?", "pdbx_diffrn_id": "?", "pdbx_ordinal": "?",
        }])
        self.assertEqual(d["refine"], [{
            "entry_id": "1XDA", "ls_number_reflns_obs": "16624",
            "ls_number_reflns_all": "?", "pdbx_ls_sigma_I": "?",
            "pdbx_ls_sigma_F": "?", "pdbx_data_cutoff_high_absF": "?",
            "pdbx_data_cutoff_low_absF": "?", "ls_d_res_low": "15.00",
            "ls_d_res_high": "1.80", "ls_percent_reflns_obs": "?",
            "ls_R_factor_obs": "0.174", "ls_R_factor_all": "0.174",
            "ls_R_factor_R_work": "0.174", "ls_R_factor_R_free": "?",
            "ls_R_factor_R_free_error": "?",
            "ls_R_factor_R_free_error_details": "?",
            "ls_percent_reflns_R_free": "?", "ls_number_reflns_R_free": "?",
            "ls_number_parameters": "?", "ls_number_restraints": "?",
            "occupancy_min": "?", "occupancy_max": "?",
            "correlation_coeff_Fo_to_Fc": "?",
            "correlation_coeff_Fo_to_Fc_free": "?", "B_iso_mean": "19.40",
            "aniso_B[1][1]": "?", "aniso_B[2][2]": "?", "aniso_B[3][3]": "?",
            "aniso_B[1][2]": "?", "aniso_B[1][3]": "?", "aniso_B[2][3]": "?",
            "solvent_model_details": "?", "solvent_model_param_ksol": "?",
            "solvent_model_param_bsol": "?", "pdbx_solvent_vdw_probe_radii": "?",
            "pdbx_solvent_ion_probe_radii": "?",
            "pdbx_solvent_shrinkage_radii": "?", "pdbx_ls_cross_valid_method": "?",
            "details": "?", "pdbx_starting_model": "?",
            "pdbx_method_to_determine_struct": "?",
            "pdbx_isotropic_thermal_model": "?",
            "pdbx_stereochemistry_target_values": "?",
            "pdbx_stereochem_target_val_spec_case": "?",
            "pdbx_R_Free_selection_details": "?", "pdbx_overall_ESU_R_Free": "?",
            "overall_SU_B": "?", "ls_redundancy_reflns_obs": "?", "B_iso_min": "?",
            "B_iso_max": "?", "overall_SU_R_Cruickshank_DPI": "?",
            "overall_SU_R_free": "?", "overall_SU_ML": "?",
            "pdbx_overall_ESU_R": "?", "pdbx_data_cutoff_high_rms_absF": "?",
            "pdbx_refine_id": "?", "pdbx_overall_phase_error": "?",
            "ls_wR_factor_R_free": "?", "ls_wR_factor_R_work": "?",
            "overall_FOM_free_R_set": "?", "overall_FOM_work_R_set": "?",
            "pdbx_diffrn_id": "?", "pdbx_TLS_residual_ADP_flag": "?",
            "pdbx_overall_SU_R_free_Cruickshank_DPI": "?",
            "pdbx_overall_SU_R_Blow_DPI": "?",
            "pdbx_overall_SU_R_free_Blow_DPI": "?",
        }])

        # Missing items categories
        self.assertNotIn("pdbx_unobs_or_zero_occ_residues", d)

        # Assembly categories
        self.assertEqual(len(d["pdbx_struct_assembly"]), 12)
        self.assertEqual(d["pdbx_struct_assembly"][0], {
            "id": "1", "details": "?", "method_details": "PISA",
            "oligomeric_details": "?", "oligomeric_count": "?",
        })
        self.assertEqual(d["pdbx_struct_assembly"][-1], {
            "id": "12", "details": "?", "method_details": "PISA",
            "oligomeric_details": "?", "oligomeric_count": "?",
        })
        self.assertEqual(len(d["pdbx_struct_assembly_gen"]), 12)
        self.assertEqual(d["pdbx_struct_assembly_gen"][0], {
            "assembly_id": "1", "oper_expression": "1",
            "asym_id_list": "A,B,I,J,K,L,Y,Z",
        })
        self.assertEqual(d["pdbx_struct_assembly_gen"][-1], {
            "assembly_id": "12", "oper_expression": "1",
            "asym_id_list": "A,B,C,D,I,J,K,L,M,N,O,P,Y,Z,AA,BA",
        })
        self.assertEqual(len(d["pdbx_struct_assembly_prop"]), 36)
        self.assertEqual(d["pdbx_struct_assembly_prop"][0], {
            "biol_id": "1", "type": "ABSA (A^2)", "value": "1720", "details": "?",
        })
        self.assertEqual(d["pdbx_struct_assembly_prop"][-1], {
            "biol_id": "12", "type": "SSA (A^2)", "value": "6770", "details": "?",
        })
        self.assertEqual(d["pdbx_struct_oper_list"], [{
            "id": "1", "type": "?", "name": "?", "symmetry_operation": "x,y,z",
            "matrix[1][1]": "1.0", "matrix[1][2]": "0.0", "matrix[1][3]": "0.0",
            "vector[1]": "0.0", "matrix[2][1]": "0.0", "matrix[2][2]": "1.0",
            "matrix[2][3]": "0.0", "vector[2]": "0.0", "matrix[3][1]": "0.0",
            "matrix[3][2]": "0.0", "matrix[3][3]": "1.0", "vector[3]": "0.0",
        }, {
            "id": "2", "type": "?", "name": "?", "symmetry_operation": "x,y,z",
            "matrix[1][1]": "-0.5", "matrix[1][2]": "-0.866025",
            "matrix[1][3]": "0.0", "vector[1]": "0.0", "matrix[2][1]": "0.866025",
            "matrix[2][2]": "-0.5", "matrix[2][3]": "0.0", "vector[2]": "0.0",
            "matrix[3][1]": "0.0", "matrix[3][2]": "0.0", "matrix[3][3]": "1.0",
            "vector[3]": "0.0",
        }, {
            "id": "3", "type": "?", "name": "?", "symmetry_operation": "x,y,z",
            "matrix[1][1]": "-0.5", "matrix[1][2]": "0.866025",
            "matrix[1][3]": "0.0", "vector[1]": "0.0", "matrix[2][1]": "-0.866025",
            "matrix[2][2]": "-0.5", "matrix[2][3]": "0.0", "vector[2]": "0.0",
            "matrix[3][1]": "0.0", "matrix[3][2]": "0.0", "matrix[3][3]": "1.0",
            "vector[3]": "0.0",
        }])


    def test_1m4x(self):
        d = atomium.open("tests/integration/files/1m4x.pdb", dictionary=True)
        self.assertEqual(len(d.keys()), 32)

        # Metadata categories
        self.assertEqual(d["entry"], [{"id": "1M4X"}])
        self.assertEqual(d["pdbx_database_status"], [
            {"status_code": "REL", "entry_id": "1M4X", "recvd_initial_deposition_date": "2002-07-05"},
        ])
        self.assertEqual(d["struct_keywords"], [{
            "entry_id": "1M4X", "pdbx_keywords": "VIRUS",
            "text": "ICOSAHEDRAL VIRUS CAPSID, BETA BARREL, ICOSAHEDRAL VIRUS, VIRUS",
        }])
        self.assertNotIn("pdbx_database_PDB_obs_spr", d)
        self.assertEqual(d["struct"], [{
            "entry_id": "1M4X", "title": "PBCV-1 VIRUS CAPSID, QUASI-ATOMIC MODEL",
            "pdbx_descriptor": "?", "pdbx_model_details": "?",
            "pdbx_CASP_flag": "?", "pdbx_model_type_details": "?",
        }])
        self.assertNotIn("pdbx_database_related", d)
        self.assertNotIn("database_PDB_caveat", d)
        self.assertEqual(d["exptl"], [{"entry_id": "1M4X", "method": "ELECTRON MICROSCOPY"}])
        self.assertNotIn("pdbx_coordinate_model", d)
        self.assertEqual(d["pdbx_audit_revision_history"], [{
            "ordinal": "1", "data_content_type": "Structure model",
            "major_revision": "1", "minor_revision": "0",
            "revision_date": "2002-12-04",
        }, {
            "ordinal": "2", "data_content_type": "Structure model",
            "major_revision": "2", "minor_revision": "0",
            "revision_date": "2003-04-22",
        }, {
            "ordinal": "3", "data_content_type": "Structure model",
            "major_revision": "3", "minor_revision": "0",
            "revision_date": "2003-10-21",
        }, {
            "ordinal": "4", "data_content_type": "Structure model",
            "major_revision": "4", "minor_revision": "0",
            "revision_date": "2009-02-24",
        }, {
            "ordinal": "5", "data_content_type": "Structure model",
            "major_revision": "5", "minor_revision": "0",
            "revision_date": "2018-07-18",
        }])

        # Entity categories
        self.assertEqual(d["entity"], [{
            "id": "1", "type": "polymer", "src_method": "nat",
            "pdbx_description": "PBCV-1 VIRUS CAPSID", "formula_weight": "?",
            "pdbx_number_of_molecules": "3", "pdbx_ec": "?", "pdbx_mutation": "?",
            "pdbx_fragment": "VIRUS CAPSID",
            "details": "VIRUS INFECTS CHLORELLA ALGAE, WHICH ARE SYMBIONTS WITH PARAMECIUM BURSARIA",
        }])
        self.assertNotIn("entity_name_com", d)
        self.assertEqual(d["entity_poly"], [{
            "entity_id": "1", "type": "polypeptide(L)", "nstd_linkage": "no",
            "nstd_monomer": "no",
            "pdbx_seq_one_letter_code": "TFFKTVYRRYTNFAIESIQQTINGSVGFGNKVSTQISRNGDLITDIVVEFVLTKGGNGGTTYYPAEELLQDVELEIGGQRIDKHYNDWFRTYDALFRMNDDRYNYRRMTDWVNNELVGAQKRFYVPLIFFFNQTPGLALPLIALQYHEVKLYFTLASQVQGVNYNGSSAIAGAAQPTMSVWVDYIFLDTQERTRFAQLPHEYLIEQLQFTGSETATPSATTQAAQNIRLNFNHPTKYLAWNFNNPTNYGQYTALANIPGACSGAGTAAATVTTPDYGNTGTYNEQLAVLDSAKIQLNGQDRFATRKGSYFNKVQPYQSIGGVTPAGVYLYSFALKPAGRQPSGTCNFSRIDNATLSLTYKTCSIDATSPAAVLGNTETVTANTATLLTALNIYAKNYNVLRIMSGMGGLAYAN",
            "pdbx_seq_one_letter_code_can": "TFFKTVYRRYTNFAIESIQQTINGSVGFGNKVSTQISRNGDLITDIVVEFVLTKGGNGGTTYYPAEELLQDVELEIGGQRIDKHYNDWFRTYDALFRMNDDRYNYRRMTDWVNNELVGAQKRFYVPLIFFFNQTPGLALPLIALQYHEVKLYFTLASQVQGVNYNGSSAIAGAAQPTMSVWVDYIFLDTQERTRFAQLPHEYLIEQLQFTGSETATPSATTQAAQNIRLNFNHPTKYLAWNFNNPTNYGQYTALANIPGACSGAGTAAATVTTPDYGNTGTYNEQLAVLDSAKIQLNGQDRFATRKGSYFNKVQPYQSIGGVTPAGVYLYSFALKPAGRQPSGTCNFSRIDNATLSLTYKTCSIDATSPAAVLGNTETVTANTATLLTALNIYAKNYNVLRIMSGMGGLAYAN",
            "pdbx_strand_id": "A,B,C", "pdbx_target_identifier": "?",
        }])
        self.assertEqual(len(d["entity_poly_seq"]), 413)
        self.assertEqual(d["entity_poly_seq"][0], {
            "entity_id": "1", "num": "1", "mon_id": "THR", "hetero": "n",
        })
        self.assertEqual(d["entity_poly_seq"][-1], {
            "entity_id": "1", "num": "413", "mon_id": "ASN", "hetero": "n",
        })
        self.assertEqual(d["struct_ref"], [{
            "id": "1", "db_name": "GB", "db_code": "AAA88828",
            "pdbx_db_accession": "323324", "pdbx_db_isoform": "?", "entity_id": "1",
            "pdbx_seq_one_letter_code": "?", "pdbx_align_begin": "?",
        }])
        self.assertEqual(d["struct_ref_seq"], [{
            "align_id": "1", "ref_id": "1", "pdbx_PDB_id_code": "1M4X",
            "pdbx_strand_id": "A", "seq_align_beg": "?",
            "pdbx_seq_align_beg_ins_code": "?", "seq_align_end": "?",
            "pdbx_seq_align_end_ins_code": "?", "pdbx_db_accession": "323324",
            "db_align_beg": "25", "pdbx_db_align_beg_ins_code": "?",
            "db_align_end": "437", "pdbx_db_align_end_ins_code": "?",
            "pdbx_auth_seq_align_beg": "25", "pdbx_auth_seq_align_end": "437",
        }, {
            "align_id": "2", "ref_id": "1", "pdbx_PDB_id_code": "1M4X",
            "pdbx_strand_id": "B", "seq_align_beg": "?",
            "pdbx_seq_align_beg_ins_code": "?", "seq_align_end": "?",
            "pdbx_seq_align_end_ins_code": "?", "pdbx_db_accession": "323324",
            "db_align_beg": "25", "pdbx_db_align_beg_ins_code": "?",
            "db_align_end": "437", "pdbx_db_align_end_ins_code": "?",
            "pdbx_auth_seq_align_beg": "25", "pdbx_auth_seq_align_end": "437",
        }, {
            "align_id": "3", "ref_id": "1", "pdbx_PDB_id_code": "1M4X",
            "pdbx_strand_id": "C", "seq_align_beg": "?",
            "pdbx_seq_align_beg_ins_code": "?", "seq_align_end": "?",
            "pdbx_seq_align_end_ins_code": "?", "pdbx_db_accession": "323324",
            "db_align_beg": "25", "pdbx_db_align_beg_ins_code": "?",
            "db_align_end": "437", "pdbx_db_align_end_ins_code": "?",
            "pdbx_auth_seq_align_beg": "25", "pdbx_auth_seq_align_end": "437",
        }])
        self.assertEqual(d["struct_ref_seq_dif"], [{
            "align_id": "1", "pdbx_pdb_id_code": "1M4X", "mon_id": "ALA",
            "pdbx_pdb_strand_id": "A", "seq_num": "?", "pdbx_pdb_ins_code": "?",
            "pdbx_seq_db_name": "GB", "pdbx_seq_db_accession_code": "323324",
            "db_mon_id": "SER", "pdbx_seq_db_seq_num": "248", "details": "CONFLICT",
            "pdbx_auth_seq_num": "248", "pdbx_ordinal": "1",
        }, {
            "align_id": "2", "pdbx_pdb_id_code": "1M4X", "mon_id": "ALA",
            "pdbx_pdb_strand_id": "B", "seq_num": "?", "pdbx_pdb_ins_code": "?",
            "pdbx_seq_db_name": "GB", "pdbx_seq_db_accession_code": "323324",
            "db_mon_id": "SER", "pdbx_seq_db_seq_num": "248", "details": "CONFLICT",
            "pdbx_auth_seq_num": "248", "pdbx_ordinal": "2",
        }, {
            "align_id": "3", "pdbx_pdb_id_code": "1M4X", "mon_id": "ALA",
            "pdbx_pdb_strand_id": "C", "seq_num": "?", "pdbx_pdb_ins_code": "?",
            "pdbx_seq_db_name": "GB", "pdbx_seq_db_accession_code": "323324",
            "db_mon_id": "SER", "pdbx_seq_db_seq_num": "248", "details": "CONFLICT",
            "pdbx_auth_seq_num": "248", "pdbx_ordinal": "3",
        }])
        self.assertNotIn("pdbx_struct_mod_residue", d)
        self.assertEqual(len(d["chem_comp"]), 20)
        self.assertEqual(d["chem_comp"][0], {
            "id": "ALA", "type": "L-peptide linking", "mon_nstd_flag": "y",
            "name": "ALANINE", "pdbx_synonyms": "?", "formula": "C3 H7 N O2",
            "formula_weight": "89.093",
        })
        self.assertEqual(d["chem_comp"][-1], {
            "id": "VAL", "type": "L-peptide linking", "mon_nstd_flag": "y",
            "name": "VALINE", "pdbx_synonyms": "?", "formula": "C5 H11 N O2",
            "formula_weight": "117.146",
        })
        self.assertNotIn("pdbx_entity_nonpoly", d)
        self.assertEqual(d["struct_asym"], [
            {"id": "A", "pdbx_blank_PDB_chainid_flag": "N", "pdbx_modified": "N", "entity_id": "1", "details": "?"},
            {"id": "B", "pdbx_blank_PDB_chainid_flag": "N", "pdbx_modified": "N", "entity_id": "1", "details": "?"},
            {"id": "C", "pdbx_blank_PDB_chainid_flag": "N", "pdbx_modified": "N", "entity_id": "1", "details": "?"},
        ])

        # Atom categories
        self.assertEqual(d["atom_type"], [{"symbol": "C"}, {"symbol": "N"}, {"symbol": "O"}, {"symbol": "S"}])
        self.assertEqual(len(d["atom_site"]), 9693)
        self.assertEqual(d["atom_site"][0], {
            "group_PDB": "ATOM", "id": "1", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "THR",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "529.720", "Cartn_y": "577.569",
            "Cartn_z": "66.498", "occupancy": "1.00", "B_iso_or_equiv": "32.89",
            "pdbx_formal_charge": "?", "auth_seq_id": "25", "auth_comp_id": "THR",
            "auth_asym_id": "A", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][3231], {
            "group_PDB": "ATOM", "id": "3232", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "THR",
            "label_asym_id": "B", "label_entity_id": "1", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "528.258", "Cartn_y": "582.202",
            "Cartn_z": "44.416", "occupancy": "1.00", "B_iso_or_equiv": "32.89",
            "pdbx_formal_charge": "?", "auth_seq_id": "25", "auth_comp_id": "THR",
            "auth_asym_id": "B", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][6462], {
            "group_PDB": "ATOM", "id": "6463", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "THR",
            "label_asym_id": "C", "label_entity_id": "1", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "546.577", "Cartn_y": "571.770",
            "Cartn_z": "52.590", "occupancy": "1.00", "B_iso_or_equiv": "32.89",
            "pdbx_formal_charge": "?", "auth_seq_id": "25", "auth_comp_id": "THR",
            "auth_asym_id": "C", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][-1], {
            "group_PDB": "ATOM", "id": "9693", "type_symbol": "O",
            "label_atom_id": "OXT", "label_alt_id": ".", "label_comp_id": "ASN",
            "label_asym_id": "C", "label_entity_id": "1", "label_seq_id": "413",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "511.742", "Cartn_y": "599.344",
            "Cartn_z": "66.077", "occupancy": "1.00", "B_iso_or_equiv": "52.28",
            "pdbx_formal_charge": "?", "auth_seq_id": "437", "auth_comp_id": "ASN",
            "auth_asym_id": "C", "auth_atom_id": "OXT", "pdbx_PDB_model_num": "1",
        })
        self.assertNotIn("atom_site_anisotrop", d)

        # Annotation categories
        self.assertEqual(len(d["struct_conf"]), 42)
        self.assertEqual(d["struct_conf"][0], {
            "conf_type_id": "HELX_P", "id": "HELX_P1", "pdbx_PDB_helix_id": "1",
            "beg_label_comp_id": "TYR", "beg_label_asym_id": "A",
            "beg_label_seq_id": "63", "pdbx_beg_PDB_ins_code": "?",
            "end_label_comp_id": "LEU", "end_label_asym_id": "A",
            "end_label_seq_id": "68", "pdbx_end_PDB_ins_code": "?",
            "beg_auth_comp_id": "TYR", "beg_auth_asym_id": "A",
            "beg_auth_seq_id": "87", "end_auth_comp_id": "LEU",
            "end_auth_asym_id": "A", "end_auth_seq_id": "92",
            "pdbx_PDB_helix_class": "1", "details": "?",
            "pdbx_PDB_helix_length": "6",
        })
        self.assertEqual(d["struct_conf"][-1], {
            "conf_type_id": "HELX_P", "id": "HELX_P42", "pdbx_PDB_helix_id": "42",
            "beg_label_comp_id": "THR", "beg_label_asym_id": "C",
            "beg_label_seq_id": "383", "pdbx_beg_PDB_ins_code": "?",
            "end_label_comp_id": "LEU", "end_label_asym_id": "C",
            "end_label_seq_id": "387", "pdbx_end_PDB_ins_code": "?",
            "beg_auth_comp_id": "THR", "beg_auth_asym_id": "C",
            "beg_auth_seq_id": "407", "end_auth_comp_id": "LEU",
            "end_auth_asym_id": "C", "end_auth_seq_id": "411",
            "pdbx_PDB_helix_class": "5", "details": "?",
            "pdbx_PDB_helix_length": "5",
        })
        self.assertEqual(len(d["struct_sheet"]), 18)
        self.assertEqual(d["struct_sheet"][0], {
            "id": "A", "type": "?", "number_strands": "6", "details": "?",
        })
        self.assertEqual(d["struct_sheet"][-1], {
            "id": "R", "type": "?", "number_strands": "2", "details": "?",
        })
        self.assertEqual(len(d["struct_sheet_range"]), 84)
        self.assertEqual(d["struct_sheet_range"][0], {
            "sheet_id": "A", "id": "1", "beg_label_comp_id": "PHE",
            "beg_label_asym_id": "A", "beg_label_seq_id": "13",
            "pdbx_beg_PDB_ins_code": "?", "end_label_comp_id": "GLN",
            "end_label_asym_id": "A", "end_label_seq_id": "19",
            "pdbx_end_PDB_ins_code": "?", "beg_auth_comp_id": "PHE",
            "beg_auth_asym_id": "A", "beg_auth_seq_id": "37",
            "end_auth_comp_id": "GLN", "end_auth_asym_id": "A",
            "end_auth_seq_id": "43",
        })
        self.assertEqual(d["struct_sheet_range"][-1], {
            "sheet_id": "R", "id": "2", "beg_label_comp_id": "VAL",
            "beg_label_asym_id": "C", "beg_label_seq_id": "379",
            "pdbx_beg_PDB_ins_code": "?", "end_label_comp_id": "ALA",
            "end_label_asym_id": "C", "end_label_seq_id": "381",
            "pdbx_end_PDB_ins_code": "?", "beg_auth_comp_id": "VAL",
            "beg_auth_asym_id": "C", "beg_auth_seq_id": "403",
            "end_auth_comp_id": "ALA", "end_auth_asym_id": "C",
            "end_auth_seq_id": "405",
        })
        self.assertEqual(len(d["struct_sheet_order"]), 66)
        self.assertEqual(d["struct_sheet_order"][0], {
            "sheet_id": "A", "range_id_1": "1", "range_id_2": "2", "offset": "?",
            "sense": "anti-parallel",
        })
        self.assertEqual(d["struct_sheet_order"][-1], {
            "sheet_id": "R", "range_id_1": "1", "range_id_2": "2", "offset": "?",
            "sense": "parallel",
        })
        self.assertEqual(len(d["pdbx_struct_sheet_hbond"]), 66)
        self.assertEqual(d["pdbx_struct_sheet_hbond"][0], {
            "sheet_id": "A", "range_id_1": "1", "range_id_2": "2",
            "range_1_label_atom_id": "N", "range_1_label_comp_id": "ILE",
            "range_1_label_asym_id": "A", "range_1_label_seq_id": "18",
            "range_1_PDB_ins_code": "?", "range_1_auth_atom_id": "N",
            "range_1_auth_comp_id": "ILE", "range_1_auth_asym_id": "A",
            "range_1_auth_seq_id": "42", "range_2_label_atom_id": "O",
            "range_2_label_comp_id": "VAL", "range_2_label_asym_id": "A",
            "range_2_label_seq_id": "187", "range_2_PDB_ins_code": "?",
            "range_2_auth_atom_id": "O", "range_2_auth_comp_id": "VAL",
            "range_2_auth_asym_id": "A", "range_2_auth_seq_id": "211",
        })
        self.assertEqual(d["pdbx_struct_sheet_hbond"][-1], {
            "sheet_id": "R", "range_id_1": "1", "range_id_2": "2",
            "range_1_label_atom_id": "N", "range_1_label_comp_id": "ILE",
            "range_1_label_asym_id": "C", "range_1_label_seq_id": "364",
            "range_1_PDB_ins_code": "?", "range_1_auth_atom_id": "N",
            "range_1_auth_comp_id": "ILE", "range_1_auth_asym_id": "C",
            "range_1_auth_seq_id": "388", "range_2_label_atom_id": "O",
            "range_2_label_comp_id": "THR", "range_2_label_asym_id": "C",
            "range_2_label_seq_id": "381", "range_2_PDB_ins_code": "?",
            "range_2_auth_atom_id": "O", "range_2_auth_comp_id": "THR",
            "range_2_auth_asym_id": "C", "range_2_auth_seq_id": "405",
        })
        self.assertNotIn("struct_conn", d)
        self.assertNotIn("struct_mon_prot_cis", d)
        self.assertNotIn("struct_site", d)
        self.assertNotIn("struct_site_gen", d)

        # Crystal categories
        self.assertEqual(d["cell"], [{
            "entry_id": "1M4X", "length_a": "1927.000", "length_b": "1927.000",
            "length_c": "1927.000", "angle_alpha": "90.00", "angle_beta": "90.00",
            "angle_gamma": "90.00", "Z_pdb": "36", "pdbx_unique_axis": "?",
        }])
        self.assertEqual(d["symmetry"], [{
            "entry_id": "1M4X", "space_group_name_H-M": "P 2 3",
            "pdbx_full_space_group_name_H-M": "?", "cell_setting": "?",
            "Int_Tables_number": "?",
        }])
        self.assertEqual(d["atom_sites"], [{
            "entry_id": "1M4X", "fract_transf_matrix[1][1]": "0.000519",
            "fract_transf_matrix[1][2]": "0.000000",
            "fract_transf_matrix[1][3]": "0.000000",
            "fract_transf_matrix[2][1]": "0.000000",
            "fract_transf_matrix[2][2]": "0.000519",
            "fract_transf_matrix[2][3]": "0.000000",
            "fract_transf_matrix[3][1]": "0.000000",
            "fract_transf_matrix[3][2]": "0.000000",
            "fract_transf_matrix[3][3]": "0.000519",
            "fract_transf_vector[1]": "0.00000",
            "fract_transf_vector[2]": "0.00000",
            "fract_transf_vector[3]": "0.00000",
        }])
        self.assertEqual(d["database_PDB_matrix"], [{
            "entry_id": "1M4X", "origx[1][1]": "1.000000",
            "origx[1][2]": "0.000000", "origx[1][3]": "0.000000",
            "origx[2][1]": "0.000000", "origx[2][2]": "1.000000",
            "origx[2][3]": "0.000000", "origx[3][1]": "0.000000",
            "origx[3][2]": "0.000000", "origx[3][3]": "1.000000",
            "origx_vector[1]": "0.00000", "origx_vector[2]": "0.00000",
            "origx_vector[3]": "0.00000",
        }])
        self.assertNotIn("struct_ncs_oper", d)

        # Citation categories
        self.assertEqual(d["audit_author"], [
            {"name": "Nandhagopal, N", "pdbx_ordinal": "1"},
            {"name": "Simpson, A.A", "pdbx_ordinal": "2"},
            {"name": "Gurnon, J.R", "pdbx_ordinal": "3"},
            {"name": "Yan, X", "pdbx_ordinal": "4"},
            {"name": "Baker, T.S", "pdbx_ordinal": "5"},
            {"name": "Graves, M.V", "pdbx_ordinal": "6"},
            {"name": "Van Etten, J.L", "pdbx_ordinal": "7"},
            {"name": "Rossmann, M.G", "pdbx_ordinal": "8"},
        ])
        self.assertEqual(d["citation_author"], [
            {"citation_id": "primary", "name": "Nandhagopal, N", "pdbx_ordinal": "1"},
            {"citation_id": "primary", "name": "Simpson, A", "pdbx_ordinal": "2"},
            {"citation_id": "primary", "name": "Gurnon, J.R", "pdbx_ordinal": "3"},
            {"citation_id": "primary", "name": "Yan, X", "pdbx_ordinal": "4"},
            {"citation_id": "primary", "name": "Baker, T.S", "pdbx_ordinal": "5"},
            {"citation_id": "primary", "name": "Graves, M.V", "pdbx_ordinal": "6"},
            {"citation_id": "primary", "name": "Van Etten, J.L", "pdbx_ordinal": "7"},
            {"citation_id": "primary", "name": "Rossmann, M.G", "pdbx_ordinal": "8"},
        ])
        self.assertNotIn("citation_editor", d)
        self.assertEqual(d["citation"], [{
            "id": "primary",
            "title": "THE STRUCTURE AND EVOLUTION OF THE MAJOR CAPSID PROTEIN OF A LARGE, LIPID CONTAINING, DNA VIRUS.",
            "journal_abbrev": "Proc.Natl.Acad.Sci.Usa", "journal_volume": "99",
            "page_first": "14758", "page_last": "?", "year": "2002",
            "journal_id_ASTM": "?", "country": "?", "journal_id_ISSN": "0027-8424",
            "journal_id_CSD": "?", "book_publisher": "?",
            "pdbx_database_id_PubMed": "12411581",
            "pdbx_database_id_DOI": "10.1073/pnas.232580699",
        }])

        # Experimental categories
        self.assertNotIn("reflns", d)
        self.assertEqual(d["refine"], [{
            "entry_id": "1M4X", "ls_number_reflns_obs": "?",
            "ls_number_reflns_all": "?", "pdbx_ls_sigma_I": "?",
            "pdbx_ls_sigma_F": "?", "pdbx_data_cutoff_high_absF": "?",
            "pdbx_data_cutoff_low_absF": "?", "ls_d_res_low": "?",
            "ls_d_res_high": "?", "ls_percent_reflns_obs": "?",
            "ls_R_factor_obs": "?", "ls_R_factor_all": "?",
            "ls_R_factor_R_work": "?", "ls_R_factor_R_free": "?",
            "ls_R_factor_R_free_error": "?",
            "ls_R_factor_R_free_error_details": "?",
            "ls_percent_reflns_R_free": "?", "ls_number_reflns_R_free": "?",
            "ls_number_parameters": "?", "ls_number_restraints": "?",
            "occupancy_min": "?", "occupancy_max": "?",
            "correlation_coeff_Fo_to_Fc": "?",
            "correlation_coeff_Fo_to_Fc_free": "?", "B_iso_mean": "?",
            "aniso_B[1][1]": "?", "aniso_B[2][2]": "?", "aniso_B[3][3]": "?",
            "aniso_B[1][2]": "?", "aniso_B[1][3]": "?", "aniso_B[2][3]": "?",
            "solvent_model_details": "?", "solvent_model_param_ksol": "?",
            "solvent_model_param_bsol": "?", "pdbx_solvent_vdw_probe_radii": "?",
            "pdbx_solvent_ion_probe_radii": "?",
            "pdbx_solvent_shrinkage_radii": "?", "pdbx_ls_cross_valid_method": "?",
            "details": "?", "pdbx_starting_model": "?",
            "pdbx_method_to_determine_struct": "?",
            "pdbx_isotropic_thermal_model": "?",
            "pdbx_stereochemistry_target_values": "RIGID BODY REFINEMENT IN REAL",
            "pdbx_stereochem_target_val_spec_case": "?",
            "pdbx_R_Free_selection_details": "?", "pdbx_overall_ESU_R_Free": "?",
            "overall_SU_B": "?", "ls_redundancy_reflns_obs": "?", "B_iso_min": "?",
            "B_iso_max": "?", "overall_SU_R_Cruickshank_DPI": "?",
            "overall_SU_R_free": "?", "overall_SU_ML": "?",
            "pdbx_overall_ESU_R": "?", "pdbx_data_cutoff_high_rms_absF": "?",
            "pdbx_refine_id": "?", "pdbx_overall_phase_error": "?",
            "ls_wR_factor_R_free": "?", "ls_wR_factor_R_work": "?",
            "overall_FOM_free_R_set": "?", "overall_FOM_work_R_set": "?",
            "pdbx_diffrn_id": "?", "pdbx_TLS_residual_ADP_flag": "?",
            "pdbx_overall_SU_R_free_Cruickshank_DPI": "?",
            "pdbx_overall_SU_R_Blow_DPI": "?",
            "pdbx_overall_SU_R_free_Blow_DPI": "?",
        }])

        # Missing items categories
        self.assertNotIn("pdbx_unobs_or_zero_occ_residues", d)

        # Assembly categories
        self.assertEqual(d["pdbx_struct_assembly"], [
            {"id": "1", "details": "?", "method_details": "?", "oligomeric_details": "?", "oligomeric_count": "?"},
        ])
        self.assertEqual(d["pdbx_struct_assembly_gen"], [{
            "assembly_id": "1",
            "oper_expression": "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,327,328,329,330,331,332,333,334,335,336,337,338,339,340,341,342,343,344,345,346,347,348,349,350,351,352,353,354,355,356,357,358,359,360,361,362,363,364,365,366,367,368,369,370,371,372,373,374,375,376,377,378,379,380,381,382,383,384,385,386,387,388,389,390,391,392,393,394,395,396,397,398,399,400,401,402,403,404,405,406,407,408,409,410,411,412,413,414,415,416,417,418,419,420,421,422,423,424,425,426,427,428,429,430,431,432,433,434,435,436,437,438,439,440,441,442,443,444,445,446,447,448,449,450,451,452,453,454,455,456,457,458,459,460,461,462,463,464,465,466,467,468,469,470,471,472,473,474,475,476,477,478,479,480,481,482,483,484,485,486,487,488,489,490,491,492,493,494,495,496,497,498,499,500,501,502,503,504,505,506,507,508,509,510,511,512,513,514,515,516,517,518,519,520,521,522,523,524,525,526,527,528,529,530,531,532,533,534,535,536,537,538,539,540,541,542,543,544,545,546,547,548,549,550,551,552,553,554,555,556,557,558,559,560,561,562,563,564,565,566,567,568,569,570,571,572,573,574,575,576,577,578,579,580,581,582,583,584,585,586,587,588,589,590,591,592,593,594,595,596,597,598,599,600,601,602,603,604,605,606,607,608,609,610,611,612,613,614,615,616,617,618,619,620,621,622,623,624,625,626,627,628,629,630,631,632,633,634,635,636,637,638,639,640,641,642,643,644,645,646,647,648,649,650,651,652,653,654,655,656,657,658,659,660,661,662,663,664,665,666,667,668,669,670,671,672,673,674,675,676,677,678,679,680,681,682,683,684,685,686,687,688,689,690,691,692,693,694,695,696,697,698,699,700,701,702,703,704,705,706,707,708,709,710,711,712,713,714,715,716,717,718,719,720,721,722,723,724,725,726,727,728,729,730,731,732,733,734,735,736,737,738,739,740,741,742,743,744,745,746,747,748,749,750,751,752,753,754,755,756,757,758,759,760,761,762,763,764,765,766,767,768,769,770,771,772,773,774,775,776,777,778,779,780,781,782,783,784,785,786,787,788,789,790,791,792,793,794,795,796,797,798,799,800,801,802,803,804,805,806,807,808,809,810,811,812,813,814,815,816,817,818,819,820,821,822,823,824,825,826,827,828,829,830,831,832,833,834,835,836,837,838,839,840,841,842,843,844,845,846,847,848,849,850,851,852,853,854,855,856,857,858,859,860,861,862,863,864,865,866,867,868,869,870,871,872,873,874,875,876,877,878,879,880,881,882,883,884,885,886,887,888,889,890,891,892,893,894,895,896,897,898,899,900,901,902,903,904,905,906,907,908,909,910,911,912,913,914,915,916,917,918,919,920,921,922,923,924,925,926,927,928,929,930,931,932,933,934,935,936,937,938,939,940,941,942,943,944,945,946,947,948,949,950,951,952,953,954,955,956,957,958,959,960,961,962,963,964,965,966,967,968,969,970,971,972,973,974,975,976,977,978,979,980,981,982,983,984,985,986,987,988,989,990,991,992,993,994,995,996,997,998,999,1000,1001,1002,1003,1004,1005,1006,1007,1008,1009,1010,1011,1012,1013,1014,1015,1016,1017,1018,1019,1020,1021,1022,1023,1024,1025,1026,1027,1028,1029,1030,1031,1032,1033,1034,1035,1036,1037,1038,1039,1040,1041,1042,1043,1044,1045,1046,1047,1048,1049,1050,1051,1052,1053,1054,1055,1056,1057,1058,1059,1060,1061,1062,1063,1064,1065,1066,1067,1068,1069,1070,1071,1072,1073,1074,1075,1076,1077,1078,1079,1080,1081,1082,1083,1084,1085,1086,1087,1088,1089,1090,1091,1092,1093,1094,1095,1096,1097,1098,1099,1100,1101,1102,1103,1104,1105,1106,1107,1108,1109,1110,1111,1112,1113,1114,1115,1116,1117,1118,1119,1120,1121,1122,1123,1124,1125,1126,1127,1128,1129,1130,1131,1132,1133,1134,1135,1136,1137,1138,1139,1140,1141,1142,1143,1144,1145,1146,1147,1148,1149,1150,1151,1152,1153,1154,1155,1156,1157,1158,1159,1160,1161,1162,1163,1164,1165,1166,1167,1168,1169,1170,1171,1172,1173,1174,1175,1176,1177,1178,1179,1180,1181,1182,1183,1184,1185,1186,1187,1188,1189,1190,1191,1192,1193,1194,1195,1196,1197,1198,1199,1200,1201,1202,1203,1204,1205,1206,1207,1208,1209,1210,1211,1212,1213,1214,1215,1216,1217,1218,1219,1220,1221,1222,1223,1224,1225,1226,1227,1228,1229,1230,1231,1232,1233,1234,1235,1236,1237,1238,1239,1240,1241,1242,1243,1244,1245,1246,1247,1248,1249,1250,1251,1252,1253,1254,1255,1256,1257,1258,1259,1260,1261,1262,1263,1264,1265,1266,1267,1268,1269,1270,1271,1272,1273,1274,1275,1276,1277,1278,1279,1280,1281,1282,1283,1284,1285,1286,1287,1288,1289,1290,1291,1292,1293,1294,1295,1296,1297,1298,1299,1300,1301,1302,1303,1304,1305,1306,1307,1308,1309,1310,1311,1312,1313,1314,1315,1316,1317,1318,1319,1320,1321,1322,1323,1324,1325,1326,1327,1328,1329,1330,1331,1332,1333,1334,1335,1336,1337,1338,1339,1340,1341,1342,1343,1344,1345,1346,1347,1348,1349,1350,1351,1352,1353,1354,1355,1356,1357,1358,1359,1360,1361,1362,1363,1364,1365,1366,1367,1368,1369,1370,1371,1372,1373,1374,1375,1376,1377,1378,1379,1380,1381,1382,1383,1384,1385,1386,1387,1388,1389,1390,1391,1392,1393,1394,1395,1396,1397,1398,1399,1400,1401,1402,1403,1404,1405,1406,1407,1408,1409,1410,1411,1412,1413,1414,1415,1416,1417,1418,1419,1420,1421,1422,1423,1424,1425,1426,1427,1428,1429,1430,1431,1432,1433,1434,1435,1436,1437,1438,1439,1440,1441,1442,1443,1444,1445,1446,1447,1448,1449,1450,1451,1452,1453,1454,1455,1456,1457,1458,1459,1460,1461,1462,1463,1464,1465,1466,1467,1468,1469,1470,1471,1472,1473,1474,1475,1476,1477,1478,1479,1480,1481,1482,1483,1484,1485,1486,1487,1488,1489,1490,1491,1492,1493,1494,1495,1496,1497,1498,1499,1500,1501,1502,1503,1504,1505,1506,1507,1508,1509,1510,1511,1512,1513,1514,1515,1516,1517,1518,1519,1520,1521,1522,1523,1524,1525,1526,1527,1528,1529,1530,1531,1532,1533,1534,1535,1536,1537,1538,1539,1540,1541,1542,1543,1544,1545,1546,1547,1548,1549,1550,1551,1552,1553,1554,1555,1556,1557,1558,1559,1560,1561,1562,1563,1564,1565,1566,1567,1568,1569,1570,1571,1572,1573,1574,1575,1576,1577,1578,1579,1580,1581,1582,1583,1584,1585,1586,1587,1588,1589,1590,1591,1592,1593,1594,1595,1596,1597,1598,1599,1600,1601,1602,1603,1604,1605,1606,1607,1608,1609,1610,1611,1612,1613,1614,1615,1616,1617,1618,1619,1620,1621,1622,1623,1624,1625,1626,1627,1628,1629,1630,1631,1632,1633,1634,1635,1636,1637,1638,1639,1640,1641,1642,1643,1644,1645,1646,1647,1648,1649,1650,1651,1652,1653,1654,1655,1656,1657,1658,1659,1660,1661,1662,1663,1664,1665,1666,1667,1668,1669,1670,1671,1672,1673,1674,1675,1676,1677,1678,1679,1680",
            "asym_id_list": "A,B,C",
        }])
        self.assertNotIn("pdbx_struct_assembly_prop", d)
        self.assertEqual(len(d["pdbx_struct_oper_list"]), 1680)
        self.assertEqual(d["pdbx_struct_oper_list"][0], {
            "id": "1", "type": "?", "name": "?", "symmetry_operation": "x,y,z",
            "matrix[1][1]": "1.0", "matrix[1][2]": "0.0", "matrix[1][3]": "0.0",
            "vector[1]": "0.0", "matrix[2][1]": "0.0", "matrix[2][2]": "1.0",
            "matrix[2][3]": "0.0", "vector[2]": "0.0", "matrix[3][1]": "0.0",
            "matrix[3][2]": "0.0", "matrix[3][3]": "1.0", "vector[3]": "0.0",
        })
        self.assertEqual(d["pdbx_struct_oper_list"][-1], {
            "id": "1680", "type": "?", "name": "?", "symmetry_operation": "x,y,z",
            "matrix[1][1]": "-0.384039", "matrix[1][2]": "-0.917901",
            "matrix[1][3]": "-0.099862", "vector[1]": "-11.9906",
            "matrix[2][1]": "0.755116", "matrix[2][2]": "-0.249996",
            "matrix[2][3]": "-0.606053", "vector[2]": "-134.0685",
            "matrix[3][1]": "0.531332", "matrix[3][2]": "-0.308156",
            "matrix[3][3]": "0.78913", "vector[3]": "-89.3239",
        })


    def test_4opj(self):
        d = atomium.open("tests/integration/files/4opj.pdb", dictionary=True)
        self.assertEqual(len(d.keys()), 41)

        # Metadata categories
        self.assertEqual(d["entry"], [{"id": "4OPJ"}])
        self.assertEqual(d["pdbx_database_status"], [
            {"status_code": "REL", "entry_id": "4OPJ", "recvd_initial_deposition_date": "2014-02-05"},
        ])
        self.assertEqual(d["struct_keywords"], [{
            "entry_id": "4OPJ", "pdbx_keywords": "HYDROLASE/DNA",
            "text": "BH RNASE-H:DNA COMPLEX, PROTEIN-DNA COMPLEX, RNASE H, RIBONUCLEASE H, TRICYCLO DNA, HYDROLASE-DNA COMPLEX",
        }])
        self.assertNotIn("pdbx_database_PDB_obs_spr", d)
        self.assertEqual(d["struct"], [{
            "entry_id": "4OPJ", "title": "BH-RNASEH:TCDA-DNA COMPLEX",
            "pdbx_descriptor": "?", "pdbx_model_details": "?",
            "pdbx_CASP_flag": "?", "pdbx_model_type_details": "?",
        }])
        self.assertNotIn("pdbx_database_related", d)
        self.assertNotIn("database_PDB_caveat", d)
        self.assertEqual(d["exptl"], [{"entry_id": "4OPJ", "method": "X-RAY DIFFRACTION"}])
        self.assertNotIn("pdbx_coordinate_model", d)
        self.assertEqual(d["pdbx_audit_revision_history"], [{
            "ordinal": "1", "data_content_type": "Structure model",
            "major_revision": "1", "minor_revision": "0",
            "revision_date": "2015-02-11",
        }, {
            "ordinal": "2", "data_content_type": "Structure model",
            "major_revision": "2", "minor_revision": "0",
            "revision_date": "2015-08-12",
        }])

        # Entity categories
        self.assertEqual(d["entity"], [{
            "id": "1", "type": "polymer", "src_method": "man",
            "pdbx_description": "RIBONUCLEASE H", "formula_weight": "?",
            "pdbx_number_of_molecules": "2", "pdbx_ec": "3.1.26.4",
            "pdbx_mutation": "?", "pdbx_fragment": "UNP RESIDUES 59-196",
            "details": "?",
        }, {
            "id": "2", "type": "polymer", "src_method": "syn",
            "pdbx_description": "5'-D(*CP*GP*CP*GP*AP*(TCY)P*TP*TP*CP*GP*CP*G)-3'",
            "formula_weight": "?", "pdbx_number_of_molecules": "2", "pdbx_ec": "?",
            "pdbx_mutation": "?", "pdbx_fragment": "?",
            "details": "TCD-ADENINE DNA",
        }, {
            "id": "3", "type": "non-polymer", "src_method": "syn",
            "pdbx_description": "GLYCEROL", "formula_weight": "?",
            "pdbx_number_of_molecules": "2", "pdbx_ec": "?", "pdbx_mutation": "?",
            "pdbx_fragment": "?", "details": "?",
        }, {
            "id": "4", "type": "water", "src_method": "nat",
            "pdbx_description": "water", "formula_weight": "?",
            "pdbx_number_of_molecules": "170", "pdbx_ec": "?", "pdbx_mutation": "?",
            "pdbx_fragment": "?", "details": "?",
        }])
        self.assertEqual(d["entity_name_com"], [{"entity_id": "1", "name": "RNASE H"}])
        self.assertEqual(d["entity_poly"], [{
            "entity_id": "1", "type": "polypeptide(L)", "nstd_linkage": "no",
            "nstd_monomer": "no",
            "pdbx_seq_one_letter_code": "GSHMAKEEIIWESLSVDVGSQGNPGIVEYKGVDTKTGEVLFEREPIPIGTNNMGEFLAIVHGLRYLKERNSRKPIYSNSQTAIKWVKDKKAKSTLVRNEETALIWKLVDEAEEWLNTHTYETPILKWQTDKWGEIKADYGRK",
            "pdbx_seq_one_letter_code_can": "GSHMAKEEIIWESLSVDVGSQGNPGIVEYKGVDTKTGEVLFEREPIPIGTNNMGEFLAIVHGLRYLKERNSRKPIYSNSQTAIKWVKDKKAKSTLVRNEETALIWKLVDEAEEWLNTHTYETPILKWQTDKWGEIKADYGRK",
            "pdbx_strand_id": "A,C", "pdbx_target_identifier": "?",
        }, {
            "entity_id": "2", "type": "polypeptide(L)", "nstd_linkage": "no",
            "nstd_monomer": "yes", "pdbx_seq_one_letter_code": "CGCGA(TCY)TTCGCG",
            "pdbx_seq_one_letter_code_can": "CGCGAATTCGCG", "pdbx_strand_id": "B,D",
            "pdbx_target_identifier": "?",
        }])
        self.assertEqual(len(d["entity_poly_seq"]), 154)
        self.assertEqual(d["entity_poly_seq"][0], {
            "entity_id": "1", "num": "1", "mon_id": "GLY", "hetero": "n",
        })
        self.assertEqual(d["entity_poly_seq"][-1], {
            "entity_id": "2", "num": "12", "mon_id": "DG", "hetero": "n",
        })
        self.assertEqual(d["struct_ref"], [{
            "id": "1", "db_name": "UNP", "db_code": "RNH1_BACHD",
            "pdbx_db_accession": "Q9KEI9", "pdbx_db_isoform": "?", "entity_id": "1",
            "pdbx_seq_one_letter_code": "?", "pdbx_align_begin": "?",
        }, {
            "id": "2", "db_name": "PDB", "db_code": "4OPJ",
            "pdbx_db_accession": "4OPJ", "pdbx_db_isoform": "?", "entity_id": "2",
            "pdbx_seq_one_letter_code": "?", "pdbx_align_begin": "?",
        }])
        self.assertEqual(len(d["struct_ref_seq"]), 4)
        self.assertEqual(d["struct_ref_seq"][0], {
            "align_id": "1", "ref_id": "1", "pdbx_PDB_id_code": "4OPJ",
            "pdbx_strand_id": "A", "seq_align_beg": "?",
            "pdbx_seq_align_beg_ins_code": "?", "seq_align_end": "?",
            "pdbx_seq_align_end_ins_code": "?", "pdbx_db_accession": "Q9KEI9",
            "db_align_beg": "59", "pdbx_db_align_beg_ins_code": "?",
            "db_align_end": "196", "pdbx_db_align_end_ins_code": "?",
            "pdbx_auth_seq_align_beg": "59", "pdbx_auth_seq_align_end": "196",
        })
        self.assertEqual(d["struct_ref_seq"][-1], {
            "align_id": "4", "ref_id": "2", "pdbx_PDB_id_code": "4OPJ",
            "pdbx_strand_id": "D", "seq_align_beg": "?",
            "pdbx_seq_align_beg_ins_code": "?", "seq_align_end": "?",
            "pdbx_seq_align_end_ins_code": "?", "pdbx_db_accession": "4OPJ",
            "db_align_beg": "1", "pdbx_db_align_beg_ins_code": "?",
            "db_align_end": "12", "pdbx_db_align_end_ins_code": "?",
            "pdbx_auth_seq_align_beg": "1", "pdbx_auth_seq_align_end": "12",
        })
        self.assertEqual(len(d["struct_ref_seq_dif"]), 10)
        self.assertEqual(d["struct_ref_seq_dif"][0], {
            "align_id": "1", "pdbx_pdb_id_code": "4OPJ", "mon_id": "GLY",
            "pdbx_pdb_strand_id": "A", "seq_num": "?", "pdbx_pdb_ins_code": "?",
            "pdbx_seq_db_name": "UNP", "pdbx_seq_db_accession_code": "Q9KEI9",
            "db_mon_id": "?", "pdbx_seq_db_seq_num": "?",
            "details": "EXPRESSION TAG", "pdbx_auth_seq_num": "55",
            "pdbx_ordinal": "1",
        })
        self.assertEqual(d["struct_ref_seq_dif"][-1], {
            "align_id": "2", "pdbx_pdb_id_code": "4OPJ", "mon_id": "ASN",
            "pdbx_pdb_strand_id": "C", "seq_num": "?", "pdbx_pdb_ins_code": "?",
            "pdbx_seq_db_name": "UNP", "pdbx_seq_db_accession_code": "Q9KEI9",
            "db_mon_id": "ASP", "pdbx_seq_db_seq_num": "132",
            "details": "ENGINEERED MUTATION", "pdbx_auth_seq_num": "132",
            "pdbx_ordinal": "10",
        })
        self.assertEqual(d["pdbx_struct_mod_residue"], [{
            "id": "1", "label_asym_id": "B", "label_comp_id": "TCY",
            "label_seq_id": "6", "auth_asym_id": "B", "auth_comp_id": "TCY",
            "auth_seq_id": "6", "PDB_ins_code": "?", "parent_comp_id": "DA",
            "details": "?",
        }, {
            "id": "2", "label_asym_id": "D", "label_comp_id": "TCY",
            "label_seq_id": "6", "auth_asym_id": "D", "auth_comp_id": "TCY",
            "auth_seq_id": "6", "PDB_ins_code": "?", "parent_comp_id": "DA",
            "details": "?",
        }])
        self.assertEqual(len(d["chem_comp"]), 26)
        self.assertEqual(d["chem_comp"][0], {
            "id": "ALA", "type": "L-peptide linking", "mon_nstd_flag": "y",
            "name": "ALANINE", "pdbx_synonyms": "?", "formula": "C3 H7 N O2",
            "formula_weight": "89.093",
        })
        self.assertEqual(d["chem_comp"][11], {
            "id": "GOL", "type": "non-polymer", "mon_nstd_flag": ".",
            "name": "GLYCEROL", "pdbx_synonyms": "GLYCERIN, PROPANE-1,2,3-TRIOL",
            "formula": "C3 H8 O3", "formula_weight": "92.094",
        })
        self.assertEqual(d["chem_comp"][13], {
            "id": "HOH", "type": "non-polymer", "mon_nstd_flag": ".",
            "name": "WATER", "pdbx_synonyms": "?", "formula": "H2 O",
            "formula_weight": "18.015",
        })
        self.assertEqual(d["chem_comp"][21], {
            "id": "TCY", "type": "L-peptide linking", "mon_nstd_flag": "n",
            "name": "(2R,3AS,4AR,5AR,5BS)-2-(6-AMINO-9H-PURIN-9-YL)-3A- HYDROXYHEXAHYDROCYCLOPROPA[4,5]CYCLOPENTA[1,2-B]FURAN- 5A(4H)-YL DIHYDROGEN PHOSPHATE",
            "pdbx_synonyms": "?", "formula": "C13 H16 N5 O6 P",
            "formula_weight": "369.269",
        })
        self.assertEqual(d["chem_comp"][-1], {
            "id": "VAL", "type": "L-peptide linking", "mon_nstd_flag": "y",
            "name": "VALINE", "pdbx_synonyms": "?", "formula": "C5 H11 N O2",
            "formula_weight": "117.146",
        })
        self.assertEqual(d["pdbx_entity_nonpoly"], [
            {"entity_id": "3", "name": "GLYCEROL", "comp_id": "GOL"},
            {"entity_id": "4", "name": "water", "comp_id": "HOH"},
        ])
        self.assertEqual(len(d["struct_asym"]), 10)
        self.assertEqual(d["struct_asym"][0], {
            "id": "A", "pdbx_blank_PDB_chainid_flag": "N", "pdbx_modified": "N",
            "entity_id": "1", "details": "?",
        })
        self.assertEqual(d["struct_asym"][-1], {
            "id": "J", "pdbx_blank_PDB_chainid_flag": "N", "pdbx_modified": "N",
            "entity_id": "4", "details": "?",
        })

        # Atom categories
        self.assertEqual(d["atom_type"], [
            {"symbol": "C"},
            {"symbol": "N"},
            {"symbol": "O"},
            {"symbol": "P"},
            {"symbol": "S"},
        ])
        self.assertEqual(len(d["atom_site"]), 2891)
        self.assertEqual(d["atom_site"][0], {
            "group_PDB": "ATOM", "id": "1", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "GLU",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "8",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-7.766", "Cartn_y": "-3.618",
            "Cartn_z": "-23.005", "occupancy": "1.00", "B_iso_or_equiv": "96.06",
            "pdbx_formal_charge": "?", "auth_seq_id": "62", "auth_comp_id": "GLU",
            "auth_asym_id": "A", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1079], {
            "group_PDB": "ATOM", "id": "1080", "type_symbol": "O",
            "label_atom_id": "O5'", "label_alt_id": ".", "label_comp_id": "DC",
            "label_asym_id": "B", "label_entity_id": "2", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "27.965", "Cartn_y": "12.580",
            "Cartn_z": "-15.037", "occupancy": "1.00", "B_iso_or_equiv": "42.59",
            "pdbx_formal_charge": "?", "auth_seq_id": "1", "auth_comp_id": "DC",
            "auth_asym_id": "B", "auth_atom_id": "O5'", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1325], {
            "group_PDB": "ATOM", "id": "1326", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "GLU",
            "label_asym_id": "C", "label_entity_id": "1", "label_seq_id": "8",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "8.483", "Cartn_y": "-14.291",
            "Cartn_z": "-3.295", "occupancy": "1.00", "B_iso_or_equiv": "63.33",
            "pdbx_formal_charge": "?", "auth_seq_id": "62", "auth_comp_id": "GLU",
            "auth_asym_id": "C", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][2401], {
            "group_PDB": "ATOM", "id": "2402", "type_symbol": "O",
            "label_atom_id": "O5'", "label_alt_id": ".", "label_comp_id": "DC",
            "label_asym_id": "D", "label_entity_id": "2", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-26.652", "Cartn_y": "-17.373",
            "Cartn_z": "3.861", "occupancy": "1.00", "B_iso_or_equiv": "109.81",
            "pdbx_formal_charge": "?", "auth_seq_id": "1", "auth_comp_id": "DC",
            "auth_asym_id": "D", "auth_atom_id": "O5'", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][2709], {
            "group_PDB": "HETATM", "id": "2710", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "GOL",
            "label_asym_id": "E", "label_entity_id": "3", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "9.641", "Cartn_y": "-0.356",
            "Cartn_z": "-9.592", "occupancy": "1.00", "B_iso_or_equiv": "43.89",
            "pdbx_formal_charge": "?", "auth_seq_id": "101", "auth_comp_id": "GOL",
            "auth_asym_id": "B", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][2715], {
            "group_PDB": "HETATM", "id": "2716", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "GOL",
            "label_asym_id": "F", "label_entity_id": "3", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-8.230", "Cartn_y": "-7.680",
            "Cartn_z": "20.558", "occupancy": "1.00", "B_iso_or_equiv": "26.15",
            "pdbx_formal_charge": "?", "auth_seq_id": "201", "auth_comp_id": "GOL",
            "auth_asym_id": "C", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][2721], {
            "group_PDB": "HETATM", "id": "2722", "type_symbol": "O",
            "label_atom_id": "O", "label_alt_id": ".", "label_comp_id": "HOH",
            "label_asym_id": "G", "label_entity_id": "4", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "18.033", "Cartn_y": "3.592",
            "Cartn_z": "-16.835", "occupancy": "1.00", "B_iso_or_equiv": "33.84",
            "pdbx_formal_charge": "?", "auth_seq_id": "201", "auth_comp_id": "HOH",
            "auth_asym_id": "B", "auth_atom_id": "O", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][2743], {
            "group_PDB": "HETATM", "id": "2744", "type_symbol": "O",
            "label_atom_id": "O", "label_alt_id": ".", "label_comp_id": "HOH",
            "label_asym_id": "H", "label_entity_id": "4", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-13.971", "Cartn_y": "-1.858",
            "Cartn_z": "3.636", "occupancy": "1.00", "B_iso_or_equiv": "47.28",
            "pdbx_formal_charge": "?", "auth_seq_id": "101", "auth_comp_id": "HOH",
            "auth_asym_id": "D", "auth_atom_id": "O", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][2756], {
            "group_PDB": "HETATM", "id": "2757", "type_symbol": "O",
            "label_atom_id": "O", "label_alt_id": ".", "label_comp_id": "HOH",
            "label_asym_id": "I", "label_entity_id": "4", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "11.111", "Cartn_y": "4.383",
            "Cartn_z": "-9.236", "occupancy": "1.00", "B_iso_or_equiv": "24.76",
            "pdbx_formal_charge": "?", "auth_seq_id": "201", "auth_comp_id": "HOH",
            "auth_asym_id": "A", "auth_atom_id": "O", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][2835], {
            "group_PDB": "HETATM", "id": "2836", "type_symbol": "O",
            "label_atom_id": "O", "label_alt_id": ".", "label_comp_id": "HOH",
            "label_asym_id": "J", "label_entity_id": "4", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-23.952", "Cartn_y": "-3.442",
            "Cartn_z": "13.303", "occupancy": "0.50", "B_iso_or_equiv": "31.76",
            "pdbx_formal_charge": "?", "auth_seq_id": "301", "auth_comp_id": "HOH",
            "auth_asym_id": "C", "auth_atom_id": "O", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][-1], {
            "group_PDB": "HETATM", "id": "2891", "type_symbol": "O",
            "label_atom_id": "O", "label_alt_id": ".", "label_comp_id": "HOH",
            "label_asym_id": "J", "label_entity_id": "4", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "1.573", "Cartn_y": "3.197",
            "Cartn_z": "10.904", "occupancy": "1.00", "B_iso_or_equiv": "51.11",
            "pdbx_formal_charge": "?", "auth_seq_id": "356", "auth_comp_id": "HOH",
            "auth_asym_id": "C", "auth_atom_id": "O", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(len(d["atom_site_anisotrop"]), 2891)
        self.assertEqual(d["atom_site_anisotrop"][0], {
            "id": "1", "type_symbol": "N", "pdbx_label_atom_id": "N",
            "pdbx_label_alt_id": ".", "pdbx_label_comp_id": "GLU",
            "pdbx_label_asym_id": "A", "pdbx_label_seq_id": "8",
            "pdbx_PDB_ins_code": "?", "U[1][1]": "1.777", "U[2][2]": "0.5569",
            "U[3][3]": "1.3158", "U[1][2]": "-0.3221", "U[1][3]": "0.5535",
            "U[2][3]": "-0.2064", "pdbx_auth_seq_id": "62",
            "pdbx_auth_comp_id": "GLU", "pdbx_auth_asym_id": "A",
            "pdbx_auth_atom_id ": "N",
        })
        self.assertEqual(d["atom_site_anisotrop"][-1], {
            "id": "2891", "type_symbol": "O", "pdbx_label_atom_id": "O",
            "pdbx_label_alt_id": ".", "pdbx_label_comp_id": "HOH",
            "pdbx_label_asym_id": "J", "pdbx_label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "U[1][1]": "0.8693", "U[2][2]": "0.3528",
            "U[3][3]": "0.7197", "U[1][2]": "0.1662", "U[1][3]": "0.1242",
            "U[2][3]": "0.4709", "pdbx_auth_seq_id": "356",
            "pdbx_auth_comp_id": "HOH", "pdbx_auth_asym_id": "C",
            "pdbx_auth_atom_id ": "O",
        })

        # Annotation categories
        self.assertEqual(len(d["struct_conf"]), 8)
        self.assertEqual(d["struct_conf"][0], {
            "conf_type_id": "HELX_P", "id": "HELX_P1", "pdbx_PDB_helix_id": "1",
            "beg_label_comp_id": "THR", "beg_label_asym_id": "A",
            "beg_label_seq_id": "50", "pdbx_beg_PDB_ins_code": "?",
            "end_label_comp_id": "ARG", "end_label_asym_id": "A",
            "end_label_seq_id": "69", "pdbx_end_PDB_ins_code": "?",
            "beg_auth_comp_id": "THR", "beg_auth_asym_id": "A",
            "beg_auth_seq_id": "104", "end_auth_comp_id": "ARG",
            "end_auth_asym_id": "A", "end_auth_seq_id": "123",
            "pdbx_PDB_helix_class": "1", "details": "?",
            "pdbx_PDB_helix_length": "20",
        })
        self.assertEqual(d["struct_conf"][-1], {
            "conf_type_id": "HELX_P", "id": "HELX_P8", "pdbx_PDB_helix_id": "8",
            "beg_label_comp_id": "GLN", "beg_label_asym_id": "C",
            "beg_label_seq_id": "128", "pdbx_beg_PDB_ins_code": "?",
            "end_label_comp_id": "GLY", "end_label_asym_id": "C",
            "end_label_seq_id": "133", "pdbx_end_PDB_ins_code": "?",
            "beg_auth_comp_id": "GLN", "beg_auth_asym_id": "C",
            "beg_auth_seq_id": "182", "end_auth_comp_id": "GLY",
            "end_auth_asym_id": "C", "end_auth_seq_id": "187",
            "pdbx_PDB_helix_class": "1", "details": "?",
            "pdbx_PDB_helix_length": "6",
        })
        self.assertEqual(d["struct_sheet"], [
            {"id": "A", "type": "?", "number_strands": "5", "details": "?"},
            {"id": "B", "type": "?", "number_strands": "5", "details": "?"},
        ])
        self.assertEqual(len(d["struct_sheet_range"]), 10)
        self.assertEqual(d["struct_sheet_range"][0], {
            "sheet_id": "A", "id": "1", "beg_label_comp_id": "VAL",
            "beg_label_asym_id": "A", "beg_label_seq_id": "39",
            "pdbx_beg_PDB_ins_code": "?", "end_label_comp_id": "GLY",
            "end_label_asym_id": "A", "end_label_seq_id": "49",
            "pdbx_end_PDB_ins_code": "?", "beg_auth_comp_id": "VAL",
            "beg_auth_asym_id": "A", "beg_auth_seq_id": "93",
            "end_auth_comp_id": "GLY", "end_auth_asym_id": "A",
            "end_auth_seq_id": "103",
        })
        self.assertEqual(d["struct_sheet_range"][-1], {
            "sheet_id": "B", "id": "5", "beg_label_comp_id": "ILE",
            "beg_label_asym_id": "C", "beg_label_seq_id": "124",
            "pdbx_beg_PDB_ins_code": "?", "end_label_comp_id": "LYS",
            "end_label_asym_id": "C", "end_label_seq_id": "126",
            "pdbx_end_PDB_ins_code": "?", "beg_auth_comp_id": "ILE",
            "beg_auth_asym_id": "C", "beg_auth_seq_id": "178",
            "end_auth_comp_id": "LYS", "end_auth_asym_id": "C",
            "end_auth_seq_id": "180",
        })
        self.assertEqual(d["struct_sheet_order"], [
            {"sheet_id": "A", "range_id_1": "1", "range_id_2": "2", "offset": "?", "sense": "anti-parallel"},
            {"sheet_id": "A", "range_id_1": "2", "range_id_2": "3", "offset": "?", "sense": "anti-parallel"},
            {"sheet_id": "A", "range_id_1": "3", "range_id_2": "4", "offset": "?", "sense": "parallel"},
            {"sheet_id": "A", "range_id_1": "4", "range_id_2": "5", "offset": "?", "sense": "parallel"},
            {"sheet_id": "B", "range_id_1": "1", "range_id_2": "2", "offset": "?", "sense": "anti-parallel"},
            {"sheet_id": "B", "range_id_1": "2", "range_id_2": "3", "offset": "?", "sense": "anti-parallel"},
            {"sheet_id": "B", "range_id_1": "3", "range_id_2": "4", "offset": "?", "sense": "parallel"},
            {"sheet_id": "B", "range_id_1": "4", "range_id_2": "5", "offset": "?", "sense": "parallel"},
        ])
        self.assertEqual(len(d["pdbx_struct_sheet_hbond"]), 8)
        self.assertEqual(d["pdbx_struct_sheet_hbond"][0], {
            "sheet_id": "A", "range_id_1": "1", "range_id_2": "2",
            "range_1_label_atom_id": "O", "range_1_label_comp_id": "ARG",
            "range_1_label_asym_id": "A", "range_1_label_seq_id": "43",
            "range_1_PDB_ins_code": "?", "range_1_auth_atom_id": "O",
            "range_1_auth_comp_id": "ARG", "range_1_auth_asym_id": "A",
            "range_1_auth_seq_id": "97", "range_2_label_atom_id": "N",
            "range_2_label_comp_id": "TYR", "range_2_label_asym_id": "A",
            "range_2_label_seq_id": "33", "range_2_PDB_ins_code": "?",
            "range_2_auth_atom_id": "N", "range_2_auth_comp_id": "TYR",
            "range_2_auth_asym_id": "A", "range_2_auth_seq_id": "87",
        })
        self.assertEqual(d["pdbx_struct_sheet_hbond"][-1], {
            "sheet_id": "B", "range_id_1": "4", "range_id_2": "5",
            "range_1_label_atom_id": "N", "range_1_label_comp_id": "ILE",
            "range_1_label_asym_id": "C", "range_1_label_seq_id": "75",
            "range_1_PDB_ins_code": "?", "range_1_auth_atom_id": "N",
            "range_1_auth_comp_id": "ILE", "range_1_auth_asym_id": "C",
            "range_1_auth_seq_id": "129", "range_2_label_atom_id": "O",
            "range_2_label_comp_id": "LEU", "range_2_label_asym_id": "C",
            "range_2_label_seq_id": "126", "range_2_PDB_ins_code": "?",
            "range_2_auth_atom_id": "O", "range_2_auth_comp_id": "LEU",
            "range_2_auth_asym_id": "C", "range_2_auth_seq_id": "180",
        })
        self.assertEqual(len(d["struct_conn"]), 5)
        self.assertEqual(d["struct_conn"][0], {
            "id": "covale1", "conn_type_id": "covale",
            "pdbx_leaving_atom_flag": "?", "pdbx_PDB_id": "?",
            "ptnr1_label_asym_id": "B", "ptnr1_label_comp_id": "DA",
            "ptnr1_label_seq_id": "5", "ptnr1_label_atom_id": "O3'",
            "pdbx_ptnr1_label_alt_id": "?", "pdbx_ptnr1_PDB_ins_code": "?",
            "pdbx_ptnr1_standard_comp_id": "?", "ptnr1_symmetry": "1555",
            "ptnr2_label_asym_id": "B", "ptnr2_label_comp_id": "TCY",
            "ptnr2_label_seq_id": "6", "ptnr2_label_atom_id": "P",
            "pdbx_ptnr2_label_alt_id": "?", "pdbx_ptnr2_PDB_ins_code": "?",
            "ptnr1_auth_asym_id": "B", "ptnr1_auth_comp_id": "DA",
            "ptnr1_auth_seq_id": "5", "ptnr2_auth_asym_id": "B",
            "ptnr2_auth_comp_id": "TCY", "ptnr2_auth_seq_id": "6",
            "ptnr2_symmetry": "1555", "pdbx_ptnr3_label_atom_id": "?",
            "pdbx_ptnr3_label_seq_id": "?", "pdbx_ptnr3_label_comp_id": "?",
            "pdbx_ptnr3_label_asym_id": "?", "pdbx_ptnr3_label_alt_id": "?",
            "pdbx_ptnr3_PDB_ins_code": "?", "details": "?",
            "pdbx_dist_value": "1.63", "pdbx_value_order": "?",
        })
        self.assertEqual(d["struct_conn"][-1], {
            "id": "covale5", "conn_type_id": "covale",
            "pdbx_leaving_atom_flag": "?", "pdbx_PDB_id": "?",
            "ptnr1_label_asym_id": "D", "ptnr1_label_comp_id": "TCY",
            "ptnr1_label_seq_id": "6", "ptnr1_label_atom_id": "O3'",
            "pdbx_ptnr1_label_alt_id": "B", "pdbx_ptnr1_PDB_ins_code": "?",
            "pdbx_ptnr1_standard_comp_id": "?", "ptnr1_symmetry": "1555",
            "ptnr2_label_asym_id": "D", "ptnr2_label_comp_id": "DT",
            "ptnr2_label_seq_id": "7", "ptnr2_label_atom_id": "P",
            "pdbx_ptnr2_label_alt_id": "B", "pdbx_ptnr2_PDB_ins_code": "?",
            "ptnr1_auth_asym_id": "D", "ptnr1_auth_comp_id": "TCY",
            "ptnr1_auth_seq_id": "6", "ptnr2_auth_asym_id": "D",
            "ptnr2_auth_comp_id": "DT", "ptnr2_auth_seq_id": "7",
            "ptnr2_symmetry": "1555", "pdbx_ptnr3_label_atom_id": "?",
            "pdbx_ptnr3_label_seq_id": "?", "pdbx_ptnr3_label_comp_id": "?",
            "pdbx_ptnr3_label_asym_id": "?", "pdbx_ptnr3_label_alt_id": "?",
            "pdbx_ptnr3_PDB_ins_code": "?", "details": "?",
            "pdbx_dist_value": "1.66", "pdbx_value_order": "?",
        })
        self.assertEqual(d["struct_mon_prot_cis"], [{
            "pdbx_id": "1", "label_comp_id": "ASN", "label_seq_id": "23",
            "label_asym_id": "A", "label_alt_id": ".", "pdbx_PDB_ins_code": "?",
            "auth_comp_id": "ASN", "auth_seq_id": "77", "auth_asym_id": "A",
            "pdbx_label_comp_id_2": "PRO", "pdbx_label_seq_id_2": "24",
            "pdbx_label_asym_id_2": "A", "pdbx_PDB_ins_code_2": "?",
            "pdbx_auth_comp_id_2": "PRO", "pdbx_auth_seq_id_2": "78",
            "pdbx_auth_asym_id_2": "A", "pdbx_PDB_model_num": "1",
            "pdbx_omega_angle": "1.58",
        }, {
            "pdbx_id": "2", "label_comp_id": "ASN", "label_seq_id": "23",
            "label_asym_id": "C", "label_alt_id": ".", "pdbx_PDB_ins_code": "?",
            "auth_comp_id": "ASN", "auth_seq_id": "77", "auth_asym_id": "C",
            "pdbx_label_comp_id_2": "PRO", "pdbx_label_seq_id_2": "24",
            "pdbx_label_asym_id_2": "C", "pdbx_PDB_ins_code_2": "?",
            "pdbx_auth_comp_id_2": "PRO", "pdbx_auth_seq_id_2": "78",
            "pdbx_auth_asym_id_2": "C", "pdbx_PDB_model_num": "1",
            "pdbx_omega_angle": "-1.32",
        }])
        self.assertNotIn("struct_site", d)
        self.assertNotIn("struct_site_gen", d)

        # Crystal categories
        self.assertEqual(d["cell"], [{
            "entry_id": "4OPJ", "length_a": "42.387", "length_b": "47.501",
            "length_c": "55.216", "angle_alpha": "100.87", "angle_beta": "101.77",
            "angle_gamma": "89.75", "Z_pdb": "2", "pdbx_unique_axis": "?",
        }])
        self.assertEqual(d["symmetry"], [{
            "entry_id": "4OPJ", "space_group_name_H-M": "P 1",
            "pdbx_full_space_group_name_H-M": "?", "cell_setting": "?",
            "Int_Tables_number": "?",
        }])
        self.assertEqual(d["atom_sites"], [{
            "entry_id": "4OPJ", "fract_transf_matrix[1][1]": "0.023592",
            "fract_transf_matrix[1][2]": "-0.000103",
            "fract_transf_matrix[1][3]": "0.004989",
            "fract_transf_matrix[2][1]": "0.000000",
            "fract_transf_matrix[2][2]": "0.021052",
            "fract_transf_matrix[2][3]": "0.004113",
            "fract_transf_matrix[3][1]": "0.000000",
            "fract_transf_matrix[3][2]": "0.000000",
            "fract_transf_matrix[3][3]": "0.018849",
            "fract_transf_vector[1]": "0.00000",
            "fract_transf_vector[2]": "0.00000",
            "fract_transf_vector[3]": "0.00000",
        }])
        self.assertEqual(d["database_PDB_matrix"], [{
            "entry_id": "4OPJ", "origx[1][1]": "1.000000",
            "origx[1][2]": "0.000000", "origx[1][3]": "0.000000",
            "origx[2][1]": "0.000000", "origx[2][2]": "1.000000",
            "origx[2][3]": "0.000000", "origx[3][1]": "0.000000",
            "origx[3][2]": "0.000000", "origx[3][3]": "1.000000",
            "origx_vector[1]": "0.00000", "origx_vector[2]": "0.00000",
            "origx_vector[3]": "0.00000",
        }])
        self.assertNotIn("struct_ncs_oper", d)

        # Citation categories
        self.assertEqual(d["audit_author"], [
            {"name": "Pallan, P.S", "pdbx_ordinal": "1"},
            {"name": "Egli, M", "pdbx_ordinal": "2"},
        ])
        self.assertEqual(d["citation_author"], [
            {"citation_id": "primary", "name": "Egli, M", "pdbx_ordinal": "1"},
            {"citation_id": "primary", "name": "Pallan, P.S", "pdbx_ordinal": "2"},
        ])
        self.assertNotIn("citation_editor", d)
        self.assertEqual(d["citation"], [{
            "id": "primary",
            "title": "GENERATING CRYSTALLOGRAPHIC MODELS OF DNA DODECAMERS FROM STRUCTURES OF RNASE H:DNA COMPLEXES.",
            "journal_abbrev": "Methods Mol.Biol.", "journal_volume": "1320",
            "page_first": "111", "page_last": "?", "year": "?",
            "journal_id_ASTM": "?", "country": "?", "journal_id_ISSN": "1064-3745",
            "journal_id_CSD": "?", "book_publisher": "?",
            "pdbx_database_id_PubMed": "26227040",
            "pdbx_database_id_DOI": "10.1007/978-1-4939-2763-0_8",
        }])

        # Experimental categories
        self.assertEqual(d["reflns"], [{
            "entry_id": "4OPJ", "observed_criterion_sigma_I": "?",
            "observed_criterion_sigma_F": "?", "d_resolution_low": "?",
            "d_resolution_high": "1.54", "number_obs": "?", "number_all": "?",
            "percent_possible_obs": "?", "pdbx_Rmerge_I_obs": "?",
            "pdbx_Rsym_value": "?", "pdbx_netI_over_sigmaI": "?",
            "B_iso_Wilson_estimate": "?", "pdbx_redundancy": "?",
            "R_free_details": "?", "limit_h_max": "?", "limit_h_min": "?",
            "limit_k_max": "?", "limit_k_min": "?", "limit_l_max": "?",
            "limit_l_min": "?", "observed_criterion_F_max": "?",
            "observed_criterion_F_min": "?", "pdbx_chi_squared": "?",
            "pdbx_scaling_rejects": "?", "pdbx_diffrn_id": "?", "pdbx_ordinal": "?",
        }])
        self.assertEqual(d["refine"], [{
            "entry_id": "4OPJ", "ls_number_reflns_obs": "52150",
            "ls_number_reflns_all": "?", "pdbx_ls_sigma_I": "?",
            "pdbx_ls_sigma_F": "?", "pdbx_data_cutoff_high_absF": "?",
            "pdbx_data_cutoff_low_absF": "?", "ls_d_res_low": "53.05",
            "ls_d_res_high": "1.54", "ls_percent_reflns_obs": "79.33",
            "ls_R_factor_obs": "0.162", "ls_R_factor_all": "0.162",
            "ls_R_factor_R_work": "0.162", "ls_R_factor_R_free": "0.220",
            "ls_R_factor_R_free_error": "?",
            "ls_R_factor_R_free_error_details": "?",
            "ls_percent_reflns_R_free": "7.800", "ls_number_reflns_R_free": "4392",
            "ls_number_parameters": "?", "ls_number_restraints": "?",
            "occupancy_min": "?", "occupancy_max": "?",
            "correlation_coeff_Fo_to_Fc": "?",
            "correlation_coeff_Fo_to_Fc_free": "?", "B_iso_mean": "32.70",
            "aniso_B[1][1]": "-2.06000", "aniso_B[2][2]": "3.82000",
            "aniso_B[3][3]": "-0.26000", "aniso_B[1][2]": "-0.60000",
            "aniso_B[1][3]": "-0.18000", "aniso_B[2][3]": "-3.57000",
            "solvent_model_details": "MASK", "solvent_model_param_ksol": "?",
            "solvent_model_param_bsol": "?", "pdbx_solvent_vdw_probe_radii": "?",
            "pdbx_solvent_ion_probe_radii": "?",
            "pdbx_solvent_shrinkage_radii": "?",
            "pdbx_ls_cross_valid_method": "THROUGHOUT", "details": "?",
            "pdbx_starting_model": "?", "pdbx_method_to_determine_struct": "?",
            "pdbx_isotropic_thermal_model": "?",
            "pdbx_stereochemistry_target_values": "MAXIMUM LIKELIHOOD",
            "pdbx_stereochem_target_val_spec_case": "?",
            "pdbx_R_Free_selection_details": "RANDOM",
            "pdbx_overall_ESU_R_Free": "?", "overall_SU_B": "?",
            "ls_redundancy_reflns_obs": "?", "B_iso_min": "?", "B_iso_max": "?",
            "overall_SU_R_Cruickshank_DPI": "?", "overall_SU_R_free": "?",
            "overall_SU_ML": "?", "pdbx_overall_ESU_R": "?",
            "pdbx_data_cutoff_high_rms_absF": "?", "pdbx_refine_id": "?",
            "pdbx_overall_phase_error": "?", "ls_wR_factor_R_free": "?",
            "ls_wR_factor_R_work": "?", "overall_FOM_free_R_set": "?",
            "overall_FOM_work_R_set": "?", "pdbx_diffrn_id": "?",
            "pdbx_TLS_residual_ADP_flag": "?",
            "pdbx_overall_SU_R_free_Cruickshank_DPI": "?",
            "pdbx_overall_SU_R_Blow_DPI": "?",
            "pdbx_overall_SU_R_free_Blow_DPI": "?",
        }])

        # Missing items categories
        self.assertEqual(len(d["pdbx_unobs_or_zero_occ_residues"]), 19)
        self.assertEqual(d["pdbx_unobs_or_zero_occ_residues"][0], {
            "id": "1", "PDB_model_num": "1", "polymer_flag": "Y",
            "occupancy_flag": "1", "auth_asym_id": "A", "auth_comp_id": "GLY",
            "auth_seq_id": "55", "PDB_ins_code": "?", "label_asym_id": "A",
            "label_comp_id": "GLY", "label_seq_id": "55",
        })
        self.assertEqual(d["pdbx_unobs_or_zero_occ_residues"][-1], {
            "id": "19", "PDB_model_num": "1", "polymer_flag": "Y",
            "occupancy_flag": "1", "auth_asym_id": "C", "auth_comp_id": "LYS",
            "auth_seq_id": "196", "PDB_ins_code": "?", "label_asym_id": "C",
            "label_comp_id": "LYS", "label_seq_id": "196",
        })

        # Assembly categories
        self.assertEqual(d["pdbx_struct_assembly"], [
            {"id": "1", "details": "?", "method_details": "PISA", "oligomeric_details": "?", "oligomeric_count": "?"},
            {"id": "2", "details": "?", "method_details": "PISA", "oligomeric_details": "?", "oligomeric_count": "?"},
        ])
        self.assertEqual(d["pdbx_struct_assembly_gen"], [
            {"assembly_id": "1", "oper_expression": "1", "asym_id_list": "D,H"},
            {"assembly_id": "1", "oper_expression": "2", "asym_id_list": "A,B,E,G,I"},
            {"assembly_id": "2", "oper_expression": "3", "asym_id_list": "B,E,G"},
            {"assembly_id": "2", "oper_expression": "2", "asym_id_list": "C,D,F,H,J"},
        ])
        self.assertEqual(d["pdbx_struct_assembly_prop"], [
            {"biol_id": "1", "type": "ABSA (A^2)", "value": "2660", "details": "?"},
            {"biol_id": "1", "type": "MORE", "value": "-11.0", "details": "?"},
            {"biol_id": "1", "type": "SSA (A^2)", "value": "10680", "details": "?"},
            {"biol_id": "2", "type": "ABSA (A^2)", "value": "2760", "details": "?"},
            {"biol_id": "2", "type": "MORE", "value": "-11.0", "details": "?"},
            {"biol_id": "2", "type": "SSA (A^2)", "value": "10770", "details": "?"},
        ])
        self.assertEqual(d["pdbx_struct_oper_list"], [{
            "id": "1", "type": "?", "name": "?", "symmetry_operation": "x,y,z",
            "matrix[1][1]": "1.0", "matrix[1][2]": "0.0", "matrix[1][3]": "0.0",
            "vector[1]": "42.387", "matrix[2][1]": "0.0", "matrix[2][2]": "1.0",
            "matrix[2][3]": "0.0", "vector[2]": "0.0", "matrix[3][1]": "0.0",
            "matrix[3][2]": "0.0", "matrix[3][3]": "1.0", "vector[3]": "0.0",
        }, {
            "id": "2", "type": "?", "name": "?", "symmetry_operation": "x,y,z",
            "matrix[1][1]": "1.0", "matrix[1][2]": "0.0", "matrix[1][3]": "0.0",
            "vector[1]": "0.0", "matrix[2][1]": "0.0", "matrix[2][2]": "1.0",
            "matrix[2][3]": "0.0", "vector[2]": "0.0", "matrix[3][1]": "0.0",
            "matrix[3][2]": "0.0", "matrix[3][3]": "1.0", "vector[3]": "0.0",
        }, {
            "id": "3", "type": "?", "name": "?", "symmetry_operation": "x,y,z",
            "matrix[1][1]": "1.0", "matrix[1][2]": "0.0", "matrix[1][3]": "0.0",
            "vector[1]": "-42.387", "matrix[2][1]": "0.0", "matrix[2][2]": "1.0",
            "matrix[2][3]": "0.0", "vector[2]": "0.0", "matrix[3][1]": "0.0",
            "matrix[3][2]": "0.0", "matrix[3][3]": "1.0", "vector[3]": "0.0",
        }])


    def test_3oaq(self):
        d = atomium.open("tests/integration/files/3oaq.pdb", dictionary=True)
        self.assertEqual(len(d.keys()), 40)

        # Metadata categories
        self.assertEqual(d["entry"], [{"id": "3OAQ"}])
        self.assertEqual(d["pdbx_database_status"], [
            {"status_code": "REL", "entry_id": "3OAQ", "recvd_initial_deposition_date": "2010-08-05"},
        ])
        self.assertEqual(d["struct_keywords"], [{
            "entry_id": "3OAQ", "pdbx_keywords": "RIBOSOME/ANTIBIOTIC",
            "text": "PROTEIN BIOSYNTHESIS, RIBOSOMES, RNA, TRNA, TRANSFER, TELITHROMYCIN, KETOLIDE, MACROLIDE, ANTIBIOTIC, EXIT, PEPTIDYL, 30S, 70S, 16S, RIBOSOMAL SUBUNIT, SMALL, RIBOSOME, RIBOSOME-ANTIBIOTIC COMPLEX",
        }])
        self.assertEqual(d["pdbx_database_PDB_obs_spr"], [
            {"id": "OBSLTE", "details": "?", "date": "2014-12-10", "pdb_id": "4V7S", "replace_pdb_id": "3OAQ"},
        ])
        self.assertEqual(d["struct"], [{
            "entry_id": "3OAQ",
            "title": "CRYSTAL STRUCTURE OF THE E. COLI RIBOSOME BOUND TO TELITHROMYCIN. THIS FILE CONTAINS THE 30S SUBUNIT OF THE FIRST 70S RIBOSOME.",
            "pdbx_descriptor": "?", "pdbx_model_details": "?",
            "pdbx_CASP_flag": "?", "pdbx_model_type_details": "?",
        }])
        self.assertEqual(d["pdbx_database_related"], [
            {"db_name": "PDB", "db_id": "3OAQ", "content_type": "split", "details": "Split 1"},
            {"db_name": "PDB", "db_id": "3OAR", "content_type": "split", "details": "Split 2"},
            {"db_name": "PDB", "db_id": "3OAS", "content_type": "split", "details": "Split 3"},
            {"db_name": "PDB", "db_id": "3OAT", "content_type": "split", "details": "Split 4"},
        ])
        self.assertNotIn("database_PDB_caveat", d)
        self.assertEqual(d["exptl"], [{"entry_id": "3OAQ", "method": "X-RAY DIFFRACTION"}])
        self.assertNotIn("pdbx_coordinate_model", d)
        self.assertEqual(d["pdbx_audit_revision_history"], [{
            "ordinal": "1", "data_content_type": "Structure model",
            "major_revision": "1", "minor_revision": "0",
            "revision_date": "2010-12-15",
        }, {
            "ordinal": "2", "data_content_type": "Structure model",
            "major_revision": "2", "minor_revision": "0",
            "revision_date": "2013-03-06",
        }, {
            "ordinal": "3", "data_content_type": "Structure model",
            "major_revision": "3", "minor_revision": "0",
            "revision_date": "2014-12-10",
        }])

        # Entity categories
        self.assertEqual(len(d["entity"]), 23)
        self.assertEqual(d["entity"][0], {
            "id": "1", "type": "polymer", "src_method": "nat",
            "pdbx_description": "16S RRNA", "formula_weight": "?",
            "pdbx_number_of_molecules": "1", "pdbx_ec": "?", "pdbx_mutation": "?",
            "pdbx_fragment": "?", "details": "?",
        })
        self.assertEqual(d["entity"][-1], {
            "id": "23", "type": "water", "src_method": "nat",
            "pdbx_description": "water", "formula_weight": "?",
            "pdbx_number_of_molecules": "208", "pdbx_ec": "?", "pdbx_mutation": "?",
            "pdbx_fragment": "?", "details": "?",
        })
        self.assertNotIn("entity_name_com", d)
        self.assertEqual(len(d["entity_poly"]), 21)
        self.assertEqual(d["entity_poly"][0], {
            "entity_id": "1", "type": "polyribonucleotide", "nstd_linkage": "no",
            "nstd_monomer": "no",
            "pdbx_seq_one_letter_code": "AAUUGAAGAGUUUGAUCAUGGCUCAGAUUGAACGCUGGCGGCAGGCCUAACACAUGCAAGUCGAACGGUAACAGGAAGAAGCUUGCUUCUUUGCUGACGAGUGGCGGACGGGUGAGUAAUGUCUGGGAAACUGCCUGAUGGAGGGGGAUAACUACUGGAAACGGUAGCUAAUACCGCAUAACGUCGCAAGACCAAAGAGGGGGACCUUCGGGCCUCUUGCCAUCGGAUGUGCCCAGAUGGGAUUAGCUAGUAGGUGGGGUAACGGCUCACCUAGGCGACGAUCCCUAGCUGGUCUGAGAGGAUGACCAGCCACACUGGAACUGAGACACGGUCCAGACUCCUACGGGAGGCAGCAGUGGGGAAUAUUGCACAAUGGGCGCAAGCCUGAUGCAGCCAUGCCGCGUGUAUGAAGAAGGCCUUCGGGUUGUAAAGUACUUUCAGCGGGGAGGAAGGGAGUAAAGUUAAUACCUUUGCUCAUUGACGUUACCCGCAGAAGAAGCACCGGCUAACUCCGUGCCAGCAGCCGCGGUAAUACGGAGGGUGCAAGCGUUAAUCGGAAUUACUGGGCGUAAAGCGCACGCAGGCGGUUUGUUAAGUCAGAUGUGAAAUCCCCGGGCUCAACCUGGGAACUGCAUCUGAUACUGGCAAGCUUGAGUCUCGUAGAGGGGGGUAGAAUUCCAGGUGUAGCGGUGAAAUGCGUAGAGAUCUGGAGGAAUACCGGUGGCGAAGGCGGCCCCCUGGACGAAGACUGACGCUCAGGUGCGAAAGCGUGGGGAGCAAACAGGAUUAGAUACCCUGGUAGUCCACGCCGUAAACGAUGUCGACUUGGAGGUUGUGCCCUUGAGGCGUGGCUUCCGGAGCUAACGCGUUAAGUCGACCGCCUGGGGAGUACGGCCGCAAGGUUAAAACUCAAAUGAAUUGACGGGGGCCCGCACAAGCGGUGGAGCAUGUGGUUUAAUUCGAUGCAACGCGAAGAACCUUACCUGGUCUUGACAUCCACGGAAGUUUUCAGAGAUGAGAAUGUGCCUUCGGGAACCGUGAGACAGGUGCUGCAUGGCUGUCGUCAGCUCGUGUUGUGAAAUGUUGGGUUAAGUCCCGCAACGAGCGCAACCCUUAUCCUUUGUUGCCAGCGGUCCGGCCGGGAACUCAAAGGAGACUGCCAGUGAUAAACUGGAGGAAGGUGGGGAUGACGUCAAGUCAUCAUGGCCCUUACGACCAGGGCUACACACGUGCUACAAUGGCGCAUACAAAGAGAAGCGACCUCGCGAGAGCAAGCGGACCUCAUAAAGUGCGUCGUAGUCCGGAUUGGAGUCUGCAACUCGACUCCAUGAAGUCGGAAUCGCUAGUAAUCGUGGAUCAGAAUGCCACGGUGAAUACGUUCCCGGGCCUUGUACACACCGCCCGUCACACCAUGGGAGUGGGUUGCAAAAGAAGUAGGUAGCUUAACCUUCGGGAGGGCGCUUACCACUUUGUGAUUCAUGACUGGGGUGAAGUCGUAACAAGGUAACCGUAGGGGAACCUGCGGUUGGAUCA",
            "pdbx_seq_one_letter_code_can": "AAUUGAAGAGUUUGAUCAUGGCUCAGAUUGAACGCUGGCGGCAGGCCUAACACAUGCAAGUCGAACGGUAACAGGAAGAAGCUUGCUUCUUUGCUGACGAGUGGCGGACGGGUGAGUAAUGUCUGGGAAACUGCCUGAUGGAGGGGGAUAACUACUGGAAACGGUAGCUAAUACCGCAUAACGUCGCAAGACCAAAGAGGGGGACCUUCGGGCCUCUUGCCAUCGGAUGUGCCCAGAUGGGAUUAGCUAGUAGGUGGGGUAACGGCUCACCUAGGCGACGAUCCCUAGCUGGUCUGAGAGGAUGACCAGCCACACUGGAACUGAGACACGGUCCAGACUCCUACGGGAGGCAGCAGUGGGGAAUAUUGCACAAUGGGCGCAAGCCUGAUGCAGCCAUGCCGCGUGUAUGAAGAAGGCCUUCGGGUUGUAAAGUACUUUCAGCGGGGAGGAAGGGAGUAAAGUUAAUACCUUUGCUCAUUGACGUUACCCGCAGAAGAAGCACCGGCUAACUCCGUGCCAGCAGCCGCGGUAAUACGGAGGGUGCAAGCGUUAAUCGGAAUUACUGGGCGUAAAGCGCACGCAGGCGGUUUGUUAAGUCAGAUGUGAAAUCCCCGGGCUCAACCUGGGAACUGCAUCUGAUACUGGCAAGCUUGAGUCUCGUAGAGGGGGGUAGAAUUCCAGGUGUAGCGGUGAAAUGCGUAGAGAUCUGGAGGAAUACCGGUGGCGAAGGCGGCCCCCUGGACGAAGACUGACGCUCAGGUGCGAAAGCGUGGGGAGCAAACAGGAUUAGAUACCCUGGUAGUCCACGCCGUAAACGAUGUCGACUUGGAGGUUGUGCCCUUGAGGCGUGGCUUCCGGAGCUAACGCGUUAAGUCGACCGCCUGGGGAGUACGGCCGCAAGGUUAAAACUCAAAUGAAUUGACGGGGGCCCGCACAAGCGGUGGAGCAUGUGGUUUAAUUCGAUGCAACGCGAAGAACCUUACCUGGUCUUGACAUCCACGGAAGUUUUCAGAGAUGAGAAUGUGCCUUCGGGAACCGUGAGACAGGUGCUGCAUGGCUGUCGUCAGCUCGUGUUGUGAAAUGUUGGGUUAAGUCCCGCAACGAGCGCAACCCUUAUCCUUUGUUGCCAGCGGUCCGGCCGGGAACUCAAAGGAGACUGCCAGUGAUAAACUGGAGGAAGGUGGGGAUGACGUCAAGUCAUCAUGGCCCUUACGACCAGGGCUACACACGUGCUACAAUGGCGCAUACAAAGAGAAGCGACCUCGCGAGAGCAAGCGGACCUCAUAAAGUGCGUCGUAGUCCGGAUUGGAGUCUGCAACUCGACUCCAUGAAGUCGGAAUCGCUAGUAAUCGUGGAUCAGAAUGCCACGGUGAAUACGUUCCCGGGCCUUGUACACACCGCCCGUCACACCAUGGGAGUGGGUUGCAAAAGAAGUAGGUAGCUUAACCUUCGGGAGGGCGCUUACCACUUUGUGAUUCAUGACUGGGGUGAAGUCGUAACAAGGUAACCGUAGGGGAACCUGCGGUUGGAUCA",
            "pdbx_strand_id": "A", "pdbx_target_identifier": "?",
        })
        self.assertEqual(d["entity_poly"][-1], {
            "entity_id": "21", "type": "polypeptide(L)", "nstd_linkage": "no",
            "nstd_monomer": "no",
            "pdbx_seq_one_letter_code": "IKVRENEPFDVALRRFKRSCEKAGVLAEVRRREFYEKPTTERKRAKASAVK",
            "pdbx_seq_one_letter_code_can": "IKVRENEPFDVALRRFKRSCEKAGVLAEVRRREFYEKPTTERKRAKASAVK",
            "pdbx_strand_id": "U", "pdbx_target_identifier": "?",
        })
        self.assertEqual(len(d["entity_poly_seq"]), 3891)
        self.assertEqual(d["entity_poly_seq"][0], {
            "entity_id": "1", "num": "1", "mon_id": "A", "hetero": "n",
        })
        self.assertEqual(d["entity_poly_seq"][-1], {
            "entity_id": "21", "num": "51", "mon_id": "LYS", "hetero": "n",
        })
        self.assertEqual(len(d["struct_ref"]), 21)
        self.assertEqual(d["struct_ref"][0], {
            "id": "1", "db_name": "GB", "db_code": "AP009048",
            "pdbx_db_accession": "AP009048.1", "pdbx_db_isoform": "?",
            "entity_id": "1", "pdbx_seq_one_letter_code": "?",
            "pdbx_align_begin": "?",
        })
        self.assertEqual(d["struct_ref"][-1], {
            "id": "21", "db_name": "UNP", "db_code": "RS21_ECOLI",
            "pdbx_db_accession": "P68679", "pdbx_db_isoform": "?",
            "entity_id": "21", "pdbx_seq_one_letter_code": "?",
            "pdbx_align_begin": "?",
        })
        self.assertEqual(len(d["struct_ref_seq"]), 21)
        self.assertEqual(d["struct_ref_seq"][0], {
            "align_id": "1", "ref_id": "1", "pdbx_PDB_id_code": "3OAQ",
            "pdbx_strand_id": "A", "seq_align_beg": "?",
            "pdbx_seq_align_beg_ins_code": "?", "seq_align_end": "?",
            "pdbx_seq_align_end_ins_code": "?", "pdbx_db_accession": "AP009048.1",
            "db_align_beg": "3427001", "pdbx_db_align_beg_ins_code": "?",
            "db_align_end": "3428533", "pdbx_db_align_end_ins_code": "?",
            "pdbx_auth_seq_align_beg": "2", "pdbx_auth_seq_align_end": "1534",
        })
        self.assertEqual(d["struct_ref_seq"][-1], {
            "align_id": "21", "ref_id": "21", "pdbx_PDB_id_code": "3OAQ",
            "pdbx_strand_id": "U", "seq_align_beg": "?",
            "pdbx_seq_align_beg_ins_code": "?", "seq_align_end": "?",
            "pdbx_seq_align_end_ins_code": "?", "pdbx_db_accession": "P68679",
            "db_align_beg": "4", "pdbx_db_align_beg_ins_code": "?",
            "db_align_end": "54", "pdbx_db_align_end_ins_code": "?",
            "pdbx_auth_seq_align_beg": "3", "pdbx_auth_seq_align_end": "53",
        })
        self.assertNotIn("struct_ref_seq_dif", d)
        self.assertNotIn("pdbx_struct_mod_residue", d)
        self.assertEqual(len(d["chem_comp"]), 26)
        self.assertEqual(d["chem_comp"][0], {
            "id": "A", "type": "L-peptide linking", "mon_nstd_flag": "y",
            "name": "ADENOSINE-5'-MONOPHOSPHATE", "pdbx_synonyms": "?",
            "formula": "C10 H14 N5 O7 P", "formula_weight": "347.220",
        })
        self.assertEqual(d["chem_comp"][12], {
            "id": "HOH", "type": "non-polymer", "mon_nstd_flag": ".",
            "name": "WATER", "pdbx_synonyms": "?", "formula": "H2 O",
            "formula_weight": "18.015",
        })
        self.assertEqual(d["chem_comp"][17], {
            "id": "MG", "type": "non-polymer", "mon_nstd_flag": ".",
            "name": "MAGNESIUM ION", "pdbx_synonyms": "?", "formula": "MG 2+",
            "formula_weight": "24.305",
        })
        self.assertEqual(d["chem_comp"][-1], {
            "id": "VAL", "type": "L-peptide linking", "mon_nstd_flag": "y",
            "name": "VALINE", "pdbx_synonyms": "?", "formula": "C5 H11 N O2",
            "formula_weight": "117.146",
        })
        self.assertEqual(d["pdbx_entity_nonpoly"], [
            {"entity_id": "22", "name": "MAGNESIUM ION", "comp_id": "MG"},
            {"entity_id": "23", "name": "water", "comp_id": "HOH"},
        ])
        self.assertEqual(len(d["struct_asym"]), 68)
        self.assertEqual(d["struct_asym"][0], {
            "id": "A", "pdbx_blank_PDB_chainid_flag": "N", "pdbx_modified": "N",
            "entity_id": "1", "details": "?",
        })
        self.assertEqual(d["struct_asym"][-1], {
            "id": "PB", "pdbx_blank_PDB_chainid_flag": "N", "pdbx_modified": "N",
            "entity_id": "23", "details": "?",
        })

        # Atom categories
        self.assertEqual(d["atom_type"], [
            {"symbol": "C"},
            {"symbol": "MG"},
            {"symbol": "N"},
            {"symbol": "O"},
            {"symbol": "P"},
            {"symbol": "S"},
        ])
        self.assertEqual(len(d["atom_site"]), 51700)
        self.assertEqual(d["atom_site"][0], {
            "group_PDB": "ATOM", "id": "1", "type_symbol": "P",
            "label_atom_id": "P", "label_alt_id": ".", "label_comp_id": "A",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-61.963", "Cartn_y": "54.361",
            "Cartn_z": "-72.804", "occupancy": "1.00", "B_iso_or_equiv": "121.71",
            "pdbx_formal_charge": "?", "auth_seq_id": "2", "auth_comp_id": "A",
            "auth_asym_id": "A", "auth_atom_id": "P", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][32895], {
            "group_PDB": "ATOM", "id": "32896", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "MET",
            "label_asym_id": "B", "label_entity_id": "2", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-69.701", "Cartn_y": "100.582",
            "Cartn_z": "-17.507", "occupancy": "1.00", "B_iso_or_equiv": "191.28",
            "pdbx_formal_charge": "?", "auth_seq_id": "8", "auth_comp_id": "MET",
            "auth_asym_id": "B", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][34600], {
            "group_PDB": "ATOM", "id": "34601", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "GLY",
            "label_asym_id": "C", "label_entity_id": "3", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-129.874", "Cartn_y": "42.130",
            "Cartn_z": "-20.099", "occupancy": "1.00", "B_iso_or_equiv": "99.55",
            "pdbx_formal_charge": "?", "auth_seq_id": "1", "auth_comp_id": "GLY",
            "auth_asym_id": "C", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][36225], {
            "group_PDB": "ATOM", "id": "36226", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "ALA",
            "label_asym_id": "D", "label_entity_id": "4", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-94.026", "Cartn_y": "28.033",
            "Cartn_z": "-77.735", "occupancy": "1.00", "B_iso_or_equiv": "52.55",
            "pdbx_formal_charge": "?", "auth_seq_id": "1", "auth_comp_id": "ALA",
            "auth_asym_id": "D", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][37868], {
            "group_PDB": "ATOM", "id": "37869", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "GLU",
            "label_asym_id": "E", "label_entity_id": "5", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-95.727", "Cartn_y": "67.753",
            "Cartn_z": "-56.981", "occupancy": "1.00", "B_iso_or_equiv": "154.87",
            "pdbx_formal_charge": "?", "auth_seq_id": "9", "auth_comp_id": "GLU",
            "auth_asym_id": "E", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][38974], {
            "group_PDB": "ATOM", "id": "38975", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "MET",
            "label_asym_id": "F", "label_entity_id": "6", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-35.935", "Cartn_y": "46.482",
            "Cartn_z": "19.922", "occupancy": "1.00", "B_iso_or_equiv": "103.37",
            "pdbx_formal_charge": "?", "auth_seq_id": "1", "auth_comp_id": "MET",
            "auth_asym_id": "F", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][39792], {
            "group_PDB": "ATOM", "id": "39793", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "PRO",
            "label_asym_id": "G", "label_entity_id": "7", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-122.467", "Cartn_y": "54.809",
            "Cartn_z": "9.412", "occupancy": "1.00", "B_iso_or_equiv": "79.06",
            "pdbx_formal_charge": "?", "auth_seq_id": "1", "auth_comp_id": "PRO",
            "auth_asym_id": "G", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][40974], {
            "group_PDB": "ATOM", "id": "40975", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "SER",
            "label_asym_id": "H", "label_entity_id": "8", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-53.354", "Cartn_y": "50.661",
            "Cartn_z": "-24.385", "occupancy": "1.00", "B_iso_or_equiv": "197.70",
            "pdbx_formal_charge": "?", "auth_seq_id": "1", "auth_comp_id": "SER",
            "auth_asym_id": "H", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][41953], {
            "group_PDB": "ATOM", "id": "41954", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "ASN",
            "label_asym_id": "I", "label_entity_id": "9", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-166.066", "Cartn_y": "87.462",
            "Cartn_z": "14.008", "occupancy": "1.00", "B_iso_or_equiv": "183.81",
            "pdbx_formal_charge": "?", "auth_seq_id": "3", "auth_comp_id": "ASN",
            "auth_asym_id": "I", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][42975], {
            "group_PDB": "ATOM", "id": "42976", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "ARG",
            "label_asym_id": "J", "label_entity_id": "10", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-178.616", "Cartn_y": "71.346",
            "Cartn_z": "-15.731", "occupancy": "1.00", "B_iso_or_equiv": "84.95",
            "pdbx_formal_charge": "?", "auth_seq_id": "5", "auth_comp_id": "ARG",
            "auth_asym_id": "J", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][43762], {
            "group_PDB": "ATOM", "id": "43763", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "ARG",
            "label_asym_id": "K", "label_entity_id": "11", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-94.640", "Cartn_y": "30.813",
            "Cartn_z": "60.611", "occupancy": "1.00", "B_iso_or_equiv": "174.13",
            "pdbx_formal_charge": "?", "auth_seq_id": "12", "auth_comp_id": "ARG",
            "auth_asym_id": "K", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][44639], {
            "group_PDB": "ATOM", "id": "44640", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "ALA",
            "label_asym_id": "L", "label_entity_id": "12", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-63.459", "Cartn_y": "37.369",
            "Cartn_z": "-31.396", "occupancy": "1.00", "B_iso_or_equiv": "69.28",
            "pdbx_formal_charge": "?", "auth_seq_id": "1", "auth_comp_id": "ALA",
            "auth_asym_id": "L", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][45594], {
            "group_PDB": "ATOM", "id": "45595", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "ALA",
            "label_asym_id": "M", "label_entity_id": "13", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-169.014", "Cartn_y": "13.486",
            "Cartn_z": "40.199", "occupancy": "1.00", "B_iso_or_equiv": "114.29",
            "pdbx_formal_charge": "?", "auth_seq_id": "1", "auth_comp_id": "ALA",
            "auth_asym_id": "M", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][46478], {
            "group_PDB": "ATOM", "id": "46479", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "ALA",
            "label_asym_id": "N", "label_entity_id": "14", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-148.648", "Cartn_y": "29.717",
            "Cartn_z": "-30.517", "occupancy": "1.00", "B_iso_or_equiv": "121.86",
            "pdbx_formal_charge": "?", "auth_seq_id": "1", "auth_comp_id": "ALA",
            "auth_asym_id": "N", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][47252], {
            "group_PDB": "ATOM", "id": "47253", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "SER",
            "label_asym_id": "O", "label_entity_id": "15", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-40.258", "Cartn_y": "44.183",
            "Cartn_z": "7.417", "occupancy": "1.00", "B_iso_or_equiv": "113.34",
            "pdbx_formal_charge": "?", "auth_seq_id": "1", "auth_comp_id": "SER",
            "auth_asym_id": "O", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][47966], {
            "group_PDB": "ATOM", "id": "47967", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "MET",
            "label_asym_id": "P", "label_entity_id": "16", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-39.701", "Cartn_y": "14.521",
            "Cartn_z": "-82.507", "occupancy": "1.00", "B_iso_or_equiv": "89.08",
            "pdbx_formal_charge": "?", "auth_seq_id": "1", "auth_comp_id": "MET",
            "auth_asym_id": "P", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][48615], {
            "group_PDB": "ATOM", "id": "48616", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "LYS",
            "label_asym_id": "Q", "label_entity_id": "17", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-22.841", "Cartn_y": "36.605",
            "Cartn_z": "-56.633", "occupancy": "1.00", "B_iso_or_equiv": "145.60",
            "pdbx_formal_charge": "?", "auth_seq_id": "3", "auth_comp_id": "LYS",
            "auth_asym_id": "Q", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][49264], {
            "group_PDB": "ATOM", "id": "49265", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "GLU",
            "label_asym_id": "R", "label_entity_id": "18", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-60.152", "Cartn_y": "56.192",
            "Cartn_z": "25.587", "occupancy": "1.00", "B_iso_or_equiv": "123.35",
            "pdbx_formal_charge": "?", "auth_seq_id": "19", "auth_comp_id": "GLU",
            "auth_asym_id": "R", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][49720], {
            "group_PDB": "ATOM", "id": "49721", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "ARG",
            "label_asym_id": "S", "label_entity_id": "19", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-172.218", "Cartn_y": "16.502",
            "Cartn_z": "0.639", "occupancy": "1.00", "B_iso_or_equiv": "203.12",
            "pdbx_formal_charge": "?", "auth_seq_id": "2", "auth_comp_id": "ARG",
            "auth_asym_id": "S", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][50358], {
            "group_PDB": "ATOM", "id": "50359", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "ASN",
            "label_asym_id": "T", "label_entity_id": "20", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-49.104", "Cartn_y": "-12.737",
            "Cartn_z": "-76.344", "occupancy": "1.00", "B_iso_or_equiv": "40.66",
            "pdbx_formal_charge": "?", "auth_seq_id": "2", "auth_comp_id": "ASN",
            "auth_asym_id": "T", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51023], {
            "group_PDB": "ATOM", "id": "51024", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "ILE",
            "label_asym_id": "U", "label_entity_id": "21", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-81.688", "Cartn_y": "48.668",
            "Cartn_z": "32.925", "occupancy": "1.00", "B_iso_or_equiv": "182.22",
            "pdbx_formal_charge": "?", "auth_seq_id": "3", "auth_comp_id": "ILE",
            "auth_asym_id": "U", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51449], {
            "group_PDB": "HETATM", "id": "51450", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "V", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-16.067", "Cartn_y": "7.183",
            "Cartn_z": "-56.354", "occupancy": "1.00", "B_iso_or_equiv": "77.18",
            "pdbx_formal_charge": "?", "auth_seq_id": "1", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51450], {
            "group_PDB": "HETATM", "id": "51451", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "W", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-47.877", "Cartn_y": "-5.110",
            "Cartn_z": "-61.679", "occupancy": "1.00", "B_iso_or_equiv": "152.08",
            "pdbx_formal_charge": "?", "auth_seq_id": "1535", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51451], {
            "group_PDB": "HETATM", "id": "51452", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "X", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-82.099", "Cartn_y": "-2.667",
            "Cartn_z": "-61.270", "occupancy": "1.00", "B_iso_or_equiv": "131.38",
            "pdbx_formal_charge": "?", "auth_seq_id": "1536", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51452], {
            "group_PDB": "HETATM", "id": "51453", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "Y", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-126.525", "Cartn_y": "30.947",
            "Cartn_z": "-71.222", "occupancy": "1.00", "B_iso_or_equiv": "121.27",
            "pdbx_formal_charge": "?", "auth_seq_id": "1537", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51453], {
            "group_PDB": "HETATM", "id": "51454", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "Z", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-102.787", "Cartn_y": "37.608",
            "Cartn_z": "-53.867", "occupancy": "1.00", "B_iso_or_equiv": "52.62",
            "pdbx_formal_charge": "?", "auth_seq_id": "1538", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51454], {
            "group_PDB": "HETATM", "id": "51455", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "AA", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-87.717", "Cartn_y": "26.895",
            "Cartn_z": "-75.864", "occupancy": "1.00", "B_iso_or_equiv": "72.18",
            "pdbx_formal_charge": "?", "auth_seq_id": "1539", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51455], {
            "group_PDB": "HETATM", "id": "51456", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "BA", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-71.593", "Cartn_y": "34.384",
            "Cartn_z": "-44.936", "occupancy": "1.00", "B_iso_or_equiv": "105.85",
            "pdbx_formal_charge": "?", "auth_seq_id": "1540", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51456], {
            "group_PDB": "HETATM", "id": "51457", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "CA", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-76.762", "Cartn_y": "34.006",
            "Cartn_z": "-19.550", "occupancy": "1.00", "B_iso_or_equiv": "47.06",
            "pdbx_formal_charge": "?", "auth_seq_id": "1541", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51457], {
            "group_PDB": "HETATM", "id": "51458", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "DA", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-65.215", "Cartn_y": "36.329",
            "Cartn_z": "-12.657", "occupancy": "1.00", "B_iso_or_equiv": "54.95",
            "pdbx_formal_charge": "?", "auth_seq_id": "1542", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51458], {
            "group_PDB": "HETATM", "id": "51459", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "EA", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-36.019", "Cartn_y": "51.486",
            "Cartn_z": "-35.109", "occupancy": "1.00", "B_iso_or_equiv": "201.45",
            "pdbx_formal_charge": "?", "auth_seq_id": "1543", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51459], {
            "group_PDB": "HETATM", "id": "51460", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "FA", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-64.338", "Cartn_y": "13.542",
            "Cartn_z": "9.781", "occupancy": "1.00", "B_iso_or_equiv": "60.92",
            "pdbx_formal_charge": "?", "auth_seq_id": "1544", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51460], {
            "group_PDB": "HETATM", "id": "51461", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "GA", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-68.639", "Cartn_y": "25.828",
            "Cartn_z": "-12.920", "occupancy": "1.00", "B_iso_or_equiv": "85.85",
            "pdbx_formal_charge": "?", "auth_seq_id": "1545", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51461], {
            "group_PDB": "HETATM", "id": "51462", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "HA", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-69.435", "Cartn_y": "14.051",
            "Cartn_z": "-16.923", "occupancy": "1.00", "B_iso_or_equiv": "50.36",
            "pdbx_formal_charge": "?", "auth_seq_id": "1546", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51462], {
            "group_PDB": "HETATM", "id": "51463", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "IA", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-128.900", "Cartn_y": "52.525",
            "Cartn_z": "3.230", "occupancy": "1.00", "B_iso_or_equiv": "156.21",
            "pdbx_formal_charge": "?", "auth_seq_id": "1547", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51463], {
            "group_PDB": "HETATM", "id": "51464", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "JA", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-124.828", "Cartn_y": "42.224",
            "Cartn_z": "13.312", "occupancy": "1.00", "B_iso_or_equiv": "128.95",
            "pdbx_formal_charge": "?", "auth_seq_id": "1548", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51464], {
            "group_PDB": "HETATM", "id": "51465", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "KA", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-159.987", "Cartn_y": "28.903",
            "Cartn_z": "-17.060", "occupancy": "1.00", "B_iso_or_equiv": "106.00",
            "pdbx_formal_charge": "?", "auth_seq_id": "1549", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51465], {
            "group_PDB": "HETATM", "id": "51466", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "LA", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-146.690", "Cartn_y": "28.556",
            "Cartn_z": "-33.937", "occupancy": "1.00", "B_iso_or_equiv": "147.63",
            "pdbx_formal_charge": "?", "auth_seq_id": "1550", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51466], {
            "group_PDB": "HETATM", "id": "51467", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "MA", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-125.619", "Cartn_y": "33.138",
            "Cartn_z": "-23.821", "occupancy": "1.00", "B_iso_or_equiv": "87.19",
            "pdbx_formal_charge": "?", "auth_seq_id": "1551", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51467], {
            "group_PDB": "HETATM", "id": "51468", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "NA", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-93.398", "Cartn_y": "59.701",
            "Cartn_z": "-31.516", "occupancy": "1.00", "B_iso_or_equiv": "209.26",
            "pdbx_formal_charge": "?", "auth_seq_id": "1552", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51468], {
            "group_PDB": "HETATM", "id": "51469", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "OA", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-158.866", "Cartn_y": "29.411",
            "Cartn_z": "20.838", "occupancy": "1.00", "B_iso_or_equiv": "126.40",
            "pdbx_formal_charge": "?", "auth_seq_id": "1553", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51469], {
            "group_PDB": "HETATM", "id": "51470", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "PA", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-76.073", "Cartn_y": "-2.011",
            "Cartn_z": "-29.290", "occupancy": "1.00", "B_iso_or_equiv": "143.19",
            "pdbx_formal_charge": "?", "auth_seq_id": "1554", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51470], {
            "group_PDB": "HETATM", "id": "51471", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "QA", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-93.584", "Cartn_y": "26.711",
            "Cartn_z": "-1.263", "occupancy": "1.00", "B_iso_or_equiv": "40.31",
            "pdbx_formal_charge": "?", "auth_seq_id": "1555", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51471], {
            "group_PDB": "HETATM", "id": "51472", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "RA", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-19.194", "Cartn_y": "8.882",
            "Cartn_z": "-51.653", "occupancy": "1.00", "B_iso_or_equiv": "89.09",
            "pdbx_formal_charge": "?", "auth_seq_id": "1556", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51472], {
            "group_PDB": "HETATM", "id": "51473", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "SA", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-86.737", "Cartn_y": "18.402",
            "Cartn_z": "15.009", "occupancy": "1.00", "B_iso_or_equiv": "102.28",
            "pdbx_formal_charge": "?", "auth_seq_id": "1557", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51473], {
            "group_PDB": "HETATM", "id": "51474", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "TA", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-74.853", "Cartn_y": "38.106",
            "Cartn_z": "-31.971", "occupancy": "1.00", "B_iso_or_equiv": "64.86",
            "pdbx_formal_charge": "?", "auth_seq_id": "1558", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51474], {
            "group_PDB": "HETATM", "id": "51475", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "UA", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-76.619", "Cartn_y": "56.385",
            "Cartn_z": "-11.903", "occupancy": "1.00", "B_iso_or_equiv": "52.72",
            "pdbx_formal_charge": "?", "auth_seq_id": "1559", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51475], {
            "group_PDB": "HETATM", "id": "51476", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "VA", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-112.074", "Cartn_y": "42.818",
            "Cartn_z": "0.439", "occupancy": "1.00", "B_iso_or_equiv": "137.27",
            "pdbx_formal_charge": "?", "auth_seq_id": "1560", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51476], {
            "group_PDB": "HETATM", "id": "51477", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "WA", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-130.327", "Cartn_y": "24.447",
            "Cartn_z": "-21.111", "occupancy": "1.00", "B_iso_or_equiv": "136.09",
            "pdbx_formal_charge": "?", "auth_seq_id": "1561", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51477], {
            "group_PDB": "HETATM", "id": "51478", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "XA", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-154.933", "Cartn_y": "18.071",
            "Cartn_z": "-18.638", "occupancy": "1.00", "B_iso_or_equiv": "84.22",
            "pdbx_formal_charge": "?", "auth_seq_id": "1562", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51478], {
            "group_PDB": "HETATM", "id": "51479", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "YA", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-33.548", "Cartn_y": "-18.980",
            "Cartn_z": "-55.205", "occupancy": "1.00", "B_iso_or_equiv": "189.36",
            "pdbx_formal_charge": "?", "auth_seq_id": "1563", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51479], {
            "group_PDB": "HETATM", "id": "51480", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "ZA", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-73.209", "Cartn_y": "38.963",
            "Cartn_z": "-48.672", "occupancy": "1.00", "B_iso_or_equiv": "168.77",
            "pdbx_formal_charge": "?", "auth_seq_id": "1564", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51480], {
            "group_PDB": "HETATM", "id": "51481", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "AB", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-31.094", "Cartn_y": "-1.463",
            "Cartn_z": "-70.789", "occupancy": "1.00", "B_iso_or_equiv": "90.82",
            "pdbx_formal_charge": "?", "auth_seq_id": "1565", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51481], {
            "group_PDB": "HETATM", "id": "51482", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "BB", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-118.919", "Cartn_y": "24.136",
            "Cartn_z": "-48.192", "occupancy": "1.00", "B_iso_or_equiv": "62.12",
            "pdbx_formal_charge": "?", "auth_seq_id": "1566", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51482], {
            "group_PDB": "HETATM", "id": "51483", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "CB", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-54.636", "Cartn_y": "32.146",
            "Cartn_z": "-70.174", "occupancy": "1.00", "B_iso_or_equiv": "62.19",
            "pdbx_formal_charge": "?", "auth_seq_id": "1567", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51483], {
            "group_PDB": "HETATM", "id": "51484", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "DB", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-129.596", "Cartn_y": "60.414",
            "Cartn_z": "-13.587", "occupancy": "1.00", "B_iso_or_equiv": "73.60",
            "pdbx_formal_charge": "?", "auth_seq_id": "1568", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51484], {
            "group_PDB": "HETATM", "id": "51485", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "EB", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-117.875", "Cartn_y": "64.747",
            "Cartn_z": "-14.844", "occupancy": "1.00", "B_iso_or_equiv": "199.41",
            "pdbx_formal_charge": "?", "auth_seq_id": "1569", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51485], {
            "group_PDB": "HETATM", "id": "51486", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "FB", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-81.840", "Cartn_y": "26.356",
            "Cartn_z": "2.744", "occupancy": "1.00", "B_iso_or_equiv": "96.47",
            "pdbx_formal_charge": "?", "auth_seq_id": "1570", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51486], {
            "group_PDB": "HETATM", "id": "51487", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "GB", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-97.160", "Cartn_y": "23.665",
            "Cartn_z": "-5.331", "occupancy": "1.00", "B_iso_or_equiv": "47.74",
            "pdbx_formal_charge": "?", "auth_seq_id": "1571", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51487], {
            "group_PDB": "HETATM", "id": "51488", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "HB", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-45.588", "Cartn_y": "-10.950",
            "Cartn_z": "-96.512", "occupancy": "1.00", "B_iso_or_equiv": "109.39",
            "pdbx_formal_charge": "?", "auth_seq_id": "1572", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51488], {
            "group_PDB": "HETATM", "id": "51489", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "IB", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-10.572", "Cartn_y": "-1.572",
            "Cartn_z": "-83.194", "occupancy": "1.00", "B_iso_or_equiv": "84.55",
            "pdbx_formal_charge": "?", "auth_seq_id": "1573", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51489], {
            "group_PDB": "HETATM", "id": "51490", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "JB", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-102.514", "Cartn_y": "21.401",
            "Cartn_z": "-53.132", "occupancy": "1.00", "B_iso_or_equiv": "148.54",
            "pdbx_formal_charge": "?", "auth_seq_id": "1574", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51490], {
            "group_PDB": "HETATM", "id": "51491", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "KB", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-60.260", "Cartn_y": "7.255",
            "Cartn_z": "-53.415", "occupancy": "1.00", "B_iso_or_equiv": "36.82",
            "pdbx_formal_charge": "?", "auth_seq_id": "1575", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51491], {
            "group_PDB": "HETATM", "id": "51492", "type_symbol": "MG",
            "label_atom_id": "MG", "label_alt_id": ".", "label_comp_id": "MG",
            "label_asym_id": "LB", "label_entity_id": "22", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-57.773", "Cartn_y": "-5.293",
            "Cartn_z": "-72.162", "occupancy": "1.00", "B_iso_or_equiv": "43.07",
            "pdbx_formal_charge": "?", "auth_seq_id": "1576", "auth_comp_id": "MG",
            "auth_asym_id": "A", "auth_atom_id": "MG", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51492], {
            "group_PDB": "HETATM", "id": "51493", "type_symbol": "O",
            "label_atom_id": "O", "label_alt_id": ".", "label_comp_id": "HOH",
            "label_asym_id": "MB", "label_entity_id": "23", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-15.156", "Cartn_y": "5.377",
            "Cartn_z": "-55.543", "occupancy": "1.00", "B_iso_or_equiv": "125.83",
            "pdbx_formal_charge": "?", "auth_seq_id": "1577", "auth_comp_id": "HOH",
            "auth_asym_id": "A", "auth_atom_id": "O", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51693], {
            "group_PDB": "HETATM", "id": "51694", "type_symbol": "O",
            "label_atom_id": "O", "label_alt_id": ".", "label_comp_id": "HOH",
            "label_asym_id": "NB", "label_entity_id": "23", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-83.359", "Cartn_y": "-1.796",
            "Cartn_z": "-59.860", "occupancy": "1.00", "B_iso_or_equiv": "135.92",
            "pdbx_formal_charge": "?", "auth_seq_id": "124", "auth_comp_id": "HOH",
            "auth_asym_id": "L", "auth_atom_id": "O", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51694], {
            "group_PDB": "HETATM", "id": "51695", "type_symbol": "O",
            "label_atom_id": "O", "label_alt_id": ".", "label_comp_id": "HOH",
            "label_asym_id": "OB", "label_entity_id": "23", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-159.202", "Cartn_y": "28.186",
            "Cartn_z": "-18.959", "occupancy": "1.00", "B_iso_or_equiv": "182.56",
            "pdbx_formal_charge": "?", "auth_seq_id": "101", "auth_comp_id": "HOH",
            "auth_asym_id": "N", "auth_atom_id": "O", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][51699], {
            "group_PDB": "HETATM", "id": "51700", "type_symbol": "O",
            "label_atom_id": "O", "label_alt_id": ".", "label_comp_id": "HOH",
            "label_asym_id": "PB", "label_entity_id": "23", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-82.874", "Cartn_y": "26.933",
            "Cartn_z": "4.575", "occupancy": "1.00", "B_iso_or_equiv": "60.23",
            "pdbx_formal_charge": "?", "auth_seq_id": "219", "auth_comp_id": "HOH",
            "auth_asym_id": "U", "auth_atom_id": "O", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][-1], {
            "group_PDB": "HETATM", "id": "51700", "type_symbol": "O",
            "label_atom_id": "O", "label_alt_id": ".", "label_comp_id": "HOH",
            "label_asym_id": "PB", "label_entity_id": "23", "label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-82.874", "Cartn_y": "26.933",
            "Cartn_z": "4.575", "occupancy": "1.00", "B_iso_or_equiv": "60.23",
            "pdbx_formal_charge": "?", "auth_seq_id": "219", "auth_comp_id": "HOH",
            "auth_asym_id": "U", "auth_atom_id": "O", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(len(d["atom_site_anisotrop"]), 50893)
        self.assertEqual(d["atom_site_anisotrop"][0], {
            "id": "1", "type_symbol": "P", "pdbx_label_atom_id": "P",
            "pdbx_label_alt_id": ".", "pdbx_label_comp_id": "A",
            "pdbx_label_asym_id": "A", "pdbx_label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "U[1][1]": "1.7681", "U[2][2]": "1.4015",
            "U[3][3]": "1.4547", "U[1][2]": "0.1896", "U[1][3]": "0.0851",
            "U[2][3]": "0.2584", "pdbx_auth_seq_id": "2", "pdbx_auth_comp_id": "A",
            "pdbx_auth_asym_id": "A", "pdbx_auth_atom_id ": "P",
        })
        self.assertEqual(d["atom_site_anisotrop"][-1], {
            "id": "51700", "type_symbol": "O", "pdbx_label_atom_id": "O",
            "pdbx_label_alt_id": ".", "pdbx_label_comp_id": "HOH",
            "pdbx_label_asym_id": "PB", "pdbx_label_seq_id": ".",
            "pdbx_PDB_ins_code": "?", "U[1][1]": "0.8884", "U[2][2]": "0.6037",
            "U[3][3]": "0.7963", "U[1][2]": "0.1388", "U[1][3]": "0.0125",
            "U[2][3]": "0.0596", "pdbx_auth_seq_id": "219",
            "pdbx_auth_comp_id": "HOH", "pdbx_auth_asym_id": "U",
            "pdbx_auth_atom_id ": "O",
        })

        # Annotation categories
        self.assertEqual(len(d["struct_conf"]), 80)
        self.assertEqual(d["struct_conf"][0], {
            "conf_type_id": "HELX_P", "id": "HELX_P1", "pdbx_PDB_helix_id": "1",
            "beg_label_comp_id": "ASN", "beg_label_asym_id": "B",
            "beg_label_seq_id": "16", "pdbx_beg_PDB_ins_code": "?",
            "end_label_comp_id": "LYS", "end_label_asym_id": "B",
            "end_label_seq_id": "20", "pdbx_end_PDB_ins_code": "?",
            "beg_auth_comp_id": "ASN", "beg_auth_asym_id": "B",
            "beg_auth_seq_id": "23", "end_auth_comp_id": "LYS",
            "end_auth_asym_id": "B", "end_auth_seq_id": "27",
            "pdbx_PDB_helix_class": "5", "details": "?",
            "pdbx_PDB_helix_length": "5",
        })
        self.assertEqual(d["struct_conf"][-1], {
            "conf_type_id": "HELX_P", "id": "HELX_P80", "pdbx_PDB_helix_id": "80",
            "beg_label_comp_id": "LYS", "beg_label_asym_id": "U",
            "beg_label_seq_id": "37", "pdbx_beg_PDB_ins_code": "?",
            "end_label_comp_id": "SER", "end_label_asym_id": "U",
            "end_label_seq_id": "48", "pdbx_end_PDB_ins_code": "?",
            "beg_auth_comp_id": "LYS", "beg_auth_asym_id": "U",
            "beg_auth_seq_id": "39", "end_auth_comp_id": "SER",
            "end_auth_asym_id": "U", "end_auth_seq_id": "50",
            "pdbx_PDB_helix_class": "1", "details": "?",
            "pdbx_PDB_helix_length": "12",
        })
        self.assertEqual(len(d["struct_sheet"]), 24)
        self.assertEqual(d["struct_sheet"][0], {
            "id": "A", "type": "?", "number_strands": "5", "details": "?",
        })
        self.assertEqual(d["struct_sheet"][-1], {
            "id": "X", "type": "?", "number_strands": "3", "details": "?",
        })
        self.assertEqual(len(d["struct_sheet_range"]), 76)
        self.assertEqual(d["struct_sheet_range"][0], {
            "sheet_id": "A", "id": "1", "beg_label_comp_id": "PHE",
            "beg_label_asym_id": "B", "beg_label_seq_id": "83",
            "pdbx_beg_PDB_ins_code": "?", "end_label_comp_id": "VAL",
            "end_label_asym_id": "B", "end_label_seq_id": "84",
            "pdbx_end_PDB_ins_code": "?", "beg_auth_comp_id": "PHE",
            "beg_auth_asym_id": "B", "beg_auth_seq_id": "90",
            "end_auth_comp_id": "VAL", "end_auth_asym_id": "B",
            "end_auth_seq_id": "91",
        })
        self.assertEqual(d["struct_sheet_range"][-1], {
            "sheet_id": "X", "id": "3", "beg_label_comp_id": "HIS",
            "beg_label_asym_id": "S", "beg_label_seq_id": "55",
            "pdbx_beg_PDB_ins_code": "?", "end_label_comp_id": "PRO",
            "end_label_asym_id": "S", "end_label_seq_id": "57",
            "pdbx_end_PDB_ins_code": "?", "beg_auth_comp_id": "HIS",
            "beg_auth_asym_id": "S", "beg_auth_seq_id": "56",
            "end_auth_comp_id": "PRO", "end_auth_asym_id": "S",
            "end_auth_seq_id": "58",
        })
        self.assertEqual(len(d["struct_sheet_order"]), 52)
        self.assertEqual(d["struct_sheet_order"][0], {
            "sheet_id": "A", "range_id_1": "1", "range_id_2": "2", "offset": "?",
            "sense": "parallel",
        })
        self.assertEqual(d["struct_sheet_order"][-1], {
            "sheet_id": "X", "range_id_1": "2", "range_id_2": "3", "offset": "?",
            "sense": "anti-parallel",
        })
        self.assertEqual(len(d["pdbx_struct_sheet_hbond"]), 52)
        self.assertEqual(d["pdbx_struct_sheet_hbond"][0], {
            "sheet_id": "A", "range_id_1": "1", "range_id_2": "2",
            "range_1_label_atom_id": "O", "range_1_label_comp_id": "VAL",
            "range_1_label_asym_id": "B", "range_1_label_seq_id": "84",
            "range_1_PDB_ins_code": "?", "range_1_auth_atom_id": "O",
            "range_1_auth_comp_id": "VAL", "range_1_auth_asym_id": "B",
            "range_1_auth_seq_id": "91", "range_2_label_atom_id": "N",
            "range_2_label_comp_id": "PHE", "range_2_label_asym_id": "B",
            "range_2_label_seq_id": "62", "range_2_PDB_ins_code": "?",
            "range_2_auth_atom_id": "N", "range_2_auth_comp_id": "PHE",
            "range_2_auth_asym_id": "B", "range_2_auth_seq_id": "69",
        })
        self.assertEqual(d["pdbx_struct_sheet_hbond"][-1], {
            "sheet_id": "X", "range_id_1": "2", "range_id_2": "3",
            "range_1_label_atom_id": "N", "range_1_label_comp_id": "VAL",
            "range_1_label_asym_id": "S", "range_1_label_seq_id": "49",
            "range_1_PDB_ins_code": "?", "range_1_auth_atom_id": "N",
            "range_1_auth_comp_id": "VAL", "range_1_auth_asym_id": "S",
            "range_1_auth_seq_id": "50", "range_2_label_atom_id": "O",
            "range_2_label_comp_id": "VAL", "range_2_label_asym_id": "S",
            "range_2_label_seq_id": "57", "range_2_PDB_ins_code": "?",
            "range_2_auth_atom_id": "O", "range_2_auth_comp_id": "VAL",
            "range_2_auth_asym_id": "S", "range_2_auth_seq_id": "58",
        })
        self.assertEqual(len(d["struct_conn"]), 250)
        self.assertEqual(d["struct_conn"][0], {
            "id": "covale1", "conn_type_id": "covale",
            "pdbx_leaving_atom_flag": "?", "pdbx_PDB_id": "?",
            "ptnr1_label_asym_id": "Z", "ptnr1_label_comp_id": "MG",
            "ptnr1_label_seq_id": ".", "ptnr1_label_atom_id": "MG",
            "pdbx_ptnr1_label_alt_id": "?", "pdbx_ptnr1_PDB_ins_code": "?",
            "pdbx_ptnr1_standard_comp_id": "?", "ptnr1_symmetry": "1555",
            "ptnr2_label_asym_id": "MB", "ptnr2_label_comp_id": "HOH",
            "ptnr2_label_seq_id": ".", "ptnr2_label_atom_id": "O",
            "pdbx_ptnr2_label_alt_id": "?", "pdbx_ptnr2_PDB_ins_code": "?",
            "ptnr1_auth_asym_id": "A", "ptnr1_auth_comp_id": "MG",
            "ptnr1_auth_seq_id": "1538", "ptnr2_auth_asym_id": "A",
            "ptnr2_auth_comp_id": "HOH", "ptnr2_auth_seq_id": "1599",
            "ptnr2_symmetry": "1555", "pdbx_ptnr3_label_atom_id": "?",
            "pdbx_ptnr3_label_seq_id": "?", "pdbx_ptnr3_label_comp_id": "?",
            "pdbx_ptnr3_label_asym_id": "?", "pdbx_ptnr3_label_alt_id": "?",
            "pdbx_ptnr3_PDB_ins_code": "?", "details": "?",
            "pdbx_dist_value": "2.06", "pdbx_value_order": "?",
        })
        self.assertEqual(d["struct_conn"][-1], {
            "id": "covale250", "conn_type_id": "covale",
            "pdbx_leaving_atom_flag": "?", "pdbx_PDB_id": "?",
            "ptnr1_label_asym_id": "A", "ptnr1_label_comp_id": "C",
            "ptnr1_label_seq_id": "351", "ptnr1_label_atom_id": "OP2",
            "pdbx_ptnr1_label_alt_id": "?", "pdbx_ptnr1_PDB_ins_code": "?",
            "pdbx_ptnr1_standard_comp_id": "?", "ptnr1_symmetry": "1555",
            "ptnr2_label_asym_id": "LB", "ptnr2_label_comp_id": "MG",
            "ptnr2_label_seq_id": ".", "ptnr2_label_atom_id": "MG",
            "pdbx_ptnr2_label_alt_id": "?", "pdbx_ptnr2_PDB_ins_code": "?",
            "ptnr1_auth_asym_id": "A", "ptnr1_auth_comp_id": "C",
            "ptnr1_auth_seq_id": "352", "ptnr2_auth_asym_id": "A",
            "ptnr2_auth_comp_id": "MG", "ptnr2_auth_seq_id": "1576",
            "ptnr2_symmetry": "1555", "pdbx_ptnr3_label_atom_id": "?",
            "pdbx_ptnr3_label_seq_id": "?", "pdbx_ptnr3_label_comp_id": "?",
            "pdbx_ptnr3_label_asym_id": "?", "pdbx_ptnr3_label_alt_id": "?",
            "pdbx_ptnr3_PDB_ins_code": "?", "details": "?",
            "pdbx_dist_value": "2.78", "pdbx_value_order": "?",
        })
        self.assertNotIn("struct_mon_prot_cis", d)
        self.assertEqual(len(d["struct_site"]), 43)
        self.assertEqual(d["struct_site"][0], {
            "id": "AC1", "pdbx_evidence_code": "Software", "pdbx_auth_asym_id": "?",
            "pdbx_auth_comp_id": "?", "pdbx_auth_seq_id": "?",
            "pdbx_auth_ins_code": "?", "pdbx_num_residues": "?",
            "details": "BINDING SITE FOR RESIDUE MG A 1",
        })
        self.assertEqual(d["struct_site"][-1], {
            "id": "EC7", "pdbx_evidence_code": "Software", "pdbx_auth_asym_id": "?",
            "pdbx_auth_comp_id": "?", "pdbx_auth_seq_id": "?",
            "pdbx_auth_ins_code": "?", "pdbx_num_residues": "?",
            "details": "BINDING SITE FOR RESIDUE MG A 1576",
        })
        self.assertEqual(len(d["struct_site_gen"]), 271)
        self.assertEqual(d["struct_site_gen"][0], {
            "id": "1", "site_id": "AC1", "pdbx_num_res": "6",
            "label_comp_id": "HOH", "label_asym_id": "MB", "label_seq_id": ".",
            "pdbx_auth_ins_code": "?", "auth_comp_id": "HOH", "auth_asym_id": "A",
            "auth_seq_id": "1577", "label_atom_id": ".", "label_alt_id": "?",
            "symmetry": "1_555", "details": "?",
        })
        self.assertEqual(d["struct_site_gen"][-1], {
            "id": "271", "site_id": "EC7", "pdbx_num_res": "7",
            "label_comp_id": "HOH", "label_asym_id": "MB", "label_seq_id": ".",
            "pdbx_auth_ins_code": "?", "auth_comp_id": "HOH", "auth_asym_id": "A",
            "auth_seq_id": "1777", "label_atom_id": ".", "label_alt_id": "?",
            "symmetry": "1_555", "details": "?",
        })

        # Crystal categories
        self.assertEqual(d["cell"], [{
            "entry_id": "3OAQ", "length_a": "210.759", "length_b": "433.272",
            "length_c": "618.863", "angle_alpha": "90.00", "angle_beta": "90.00",
            "angle_gamma": "90.00", "Z_pdb": "4", "pdbx_unique_axis": "?",
        }])
        self.assertEqual(d["symmetry"], [{
            "entry_id": "3OAQ", "space_group_name_H-M": "P 21 21 21",
            "pdbx_full_space_group_name_H-M": "?", "cell_setting": "?",
            "Int_Tables_number": "?",
        }])
        self.assertEqual(d["atom_sites"], [{
            "entry_id": "3OAQ", "fract_transf_matrix[1][1]": "0.004745",
            "fract_transf_matrix[1][2]": "0.000000",
            "fract_transf_matrix[1][3]": "0.000000",
            "fract_transf_matrix[2][1]": "0.000000",
            "fract_transf_matrix[2][2]": "0.002308",
            "fract_transf_matrix[2][3]": "0.000000",
            "fract_transf_matrix[3][1]": "0.000000",
            "fract_transf_matrix[3][2]": "0.000000",
            "fract_transf_matrix[3][3]": "0.001616",
            "fract_transf_vector[1]": "0.00000",
            "fract_transf_vector[2]": "0.00000",
            "fract_transf_vector[3]": "0.00000",
        }])
        self.assertEqual(d["database_PDB_matrix"], [{
            "entry_id": "3OAQ", "origx[1][1]": "1.000000",
            "origx[1][2]": "0.000000", "origx[1][3]": "0.000000",
            "origx[2][1]": "0.000000", "origx[2][2]": "1.000000",
            "origx[2][3]": "0.000000", "origx[3][1]": "0.000000",
            "origx[3][2]": "0.000000", "origx[3][3]": "1.000000",
            "origx_vector[1]": "0.00000", "origx_vector[2]": "0.00000",
            "origx_vector[3]": "0.00000",
        }])
        self.assertNotIn("struct_ncs_oper", d)

        # Citation categories
        self.assertEqual(d["audit_author"], [
            {"name": "Dunkle, J.A", "pdbx_ordinal": "1"},
            {"name": "Xiong, L", "pdbx_ordinal": "2"},
            {"name": "Mankin, A.S", "pdbx_ordinal": "3"},
            {"name": "Cate, J.H.D", "pdbx_ordinal": "4"},
        ])
        self.assertEqual(d["citation_author"], [
            {"citation_id": "primary", "name": "Dunkle, J.A", "pdbx_ordinal": "1"},
            {"citation_id": "primary", "name": "Xiong, L", "pdbx_ordinal": "2"},
            {"citation_id": "primary", "name": "Mankin, A.S", "pdbx_ordinal": "3"},
            {"citation_id": "primary", "name": "Cate, J.H", "pdbx_ordinal": "4"},
        ])
        self.assertNotIn("citation_editor", d)
        self.assertEqual(d["citation"], [{
            "id": "primary",
            "title": "STRUCTURES OF THE ESCHERICHIA COLI RIBOSOME WITH ANTIBIOTICS BOUND NEAR THE PEPTIDYL TRANSFERASE CENTER EXPLAIN SPECTRA OF DRUG ACTION.",
            "journal_abbrev": "Proc.Natl.Acad.Sci.Usa", "journal_volume": "107",
            "page_first": "17152", "page_last": "?", "year": "2010",
            "journal_id_ASTM": "?", "country": "?", "journal_id_ISSN": "0027-8424",
            "journal_id_CSD": "?", "book_publisher": "?",
            "pdbx_database_id_PubMed": "20876128",
            "pdbx_database_id_DOI": "10.1073/pnas.1007988107",
        }])

        # Experimental categories
        self.assertEqual(d["reflns"], [{
            "entry_id": "3OAQ", "observed_criterion_sigma_I": "?",
            "observed_criterion_sigma_F": "?", "d_resolution_low": "?",
            "d_resolution_high": "3.25", "number_obs": "?", "number_all": "?",
            "percent_possible_obs": "?", "pdbx_Rmerge_I_obs": "?",
            "pdbx_Rsym_value": "?", "pdbx_netI_over_sigmaI": "?",
            "B_iso_Wilson_estimate": "?", "pdbx_redundancy": "?",
            "R_free_details": "?", "limit_h_max": "?", "limit_h_min": "?",
            "limit_k_max": "?", "limit_k_min": "?", "limit_l_max": "?",
            "limit_l_min": "?", "observed_criterion_F_max": "?",
            "observed_criterion_F_min": "?", "pdbx_chi_squared": "?",
            "pdbx_scaling_rejects": "?", "pdbx_diffrn_id": "?", "pdbx_ordinal": "?",
        }])
        self.assertEqual(d["refine"], [{
            "entry_id": "3OAQ", "ls_number_reflns_obs": "754750",
            "ls_number_reflns_all": "?", "pdbx_ls_sigma_I": "?",
            "pdbx_ls_sigma_F": "?", "pdbx_data_cutoff_high_absF": "?",
            "pdbx_data_cutoff_low_absF": "?", "ls_d_res_low": "85.22",
            "ls_d_res_high": "3.25", "ls_percent_reflns_obs": "?",
            "ls_R_factor_obs": "0.193", "ls_R_factor_all": "0.193",
            "ls_R_factor_R_work": "0.193", "ls_R_factor_R_free": "0.245",
            "ls_R_factor_R_free_error": "?",
            "ls_R_factor_R_free_error_details": "?",
            "ls_percent_reflns_R_free": "2.020", "ls_number_reflns_R_free": "15234",
            "ls_number_parameters": "?", "ls_number_restraints": "?",
            "occupancy_min": "?", "occupancy_max": "?",
            "correlation_coeff_Fo_to_Fc": "?",
            "correlation_coeff_Fo_to_Fc_free": "?", "B_iso_mean": "?",
            "aniso_B[1][1]": "-2.25550", "aniso_B[2][2]": "-12.54200",
            "aniso_B[3][3]": "9.45610", "aniso_B[1][2]": "-0.00000",
            "aniso_B[1][3]": "0.00000", "aniso_B[2][3]": "0.00000",
            "solvent_model_details": "FLAT BULK SOLVENT MODEL",
            "solvent_model_param_ksol": "?", "solvent_model_param_bsol": "?",
            "pdbx_solvent_vdw_probe_radii": "?",
            "pdbx_solvent_ion_probe_radii": "?",
            "pdbx_solvent_shrinkage_radii": "?", "pdbx_ls_cross_valid_method": "?",
            "details": "?", "pdbx_starting_model": "?",
            "pdbx_method_to_determine_struct": "?",
            "pdbx_isotropic_thermal_model": "?",
            "pdbx_stereochemistry_target_values": "ML",
            "pdbx_stereochem_target_val_spec_case": "?",
            "pdbx_R_Free_selection_details": "?", "pdbx_overall_ESU_R_Free": "?",
            "overall_SU_B": "?", "ls_redundancy_reflns_obs": "?", "B_iso_min": "?",
            "B_iso_max": "?", "overall_SU_R_Cruickshank_DPI": "?",
            "overall_SU_R_free": "?", "overall_SU_ML": "?",
            "pdbx_overall_ESU_R": "?", "pdbx_data_cutoff_high_rms_absF": "?",
            "pdbx_refine_id": "?", "pdbx_overall_phase_error": "?",
            "ls_wR_factor_R_free": "?", "ls_wR_factor_R_work": "?",
            "overall_FOM_free_R_set": "?", "overall_FOM_work_R_set": "?",
            "pdbx_diffrn_id": "?", "pdbx_TLS_residual_ADP_flag": "?",
            "pdbx_overall_SU_R_free_Cruickshank_DPI": "?",
            "pdbx_overall_SU_R_Blow_DPI": "?",
            "pdbx_overall_SU_R_free_Blow_DPI": "?",
        }])

        # Missing items categories
        self.assertEqual(d["pdbx_unobs_or_zero_occ_residues"], [{
            "id": "1", "PDB_model_num": "1", "polymer_flag": "Y",
            "occupancy_flag": "1", "auth_asym_id": "N", "auth_comp_id": "SER",
            "auth_seq_id": "36", "PDB_ins_code": "?", "label_asym_id": "N",
            "label_comp_id": "SER", "label_seq_id": "36",
        }, {
            "id": "2", "PDB_model_num": "1", "polymer_flag": "Y",
            "occupancy_flag": "1", "auth_asym_id": "N", "auth_comp_id": "ASP",
            "auth_seq_id": "37", "PDB_ins_code": "?", "label_asym_id": "N",
            "label_comp_id": "ASP", "label_seq_id": "37",
        }, {
            "id": "3", "PDB_model_num": "1", "polymer_flag": "Y",
            "occupancy_flag": "1", "auth_asym_id": "N", "auth_comp_id": "GLU",
            "auth_seq_id": "38", "PDB_ins_code": "?", "label_asym_id": "N",
            "label_comp_id": "GLU", "label_seq_id": "38",
        }, {
            "id": "4", "PDB_model_num": "1", "polymer_flag": "Y",
            "occupancy_flag": "1", "auth_asym_id": "N", "auth_comp_id": "ASP",
            "auth_seq_id": "39", "PDB_ins_code": "?", "label_asym_id": "N",
            "label_comp_id": "ASP", "label_seq_id": "39",
        }])

        # Assembly categories
        self.assertEqual(d["pdbx_struct_assembly"], [
            {"id": "1", "details": "?", "method_details": "?", "oligomeric_details": "?", "oligomeric_count": "?"},
        ])
        self.assertEqual(d["pdbx_struct_assembly_gen"], [{
            "assembly_id": "1", "oper_expression": "1",
            "asym_id_list": "A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,AA,BA,CA,DA,EA,FA,GA,HA,IA,JA,KA,LA,MA,NA,OA,PA,QA,RA,SA,TA,UA,VA,WA,XA,YA,ZA,AB,BB,CB,DB,EB,FB,GB,HB,IB,JB,KB,LB,MB,NB,OB,PB",
        }])
        self.assertNotIn("pdbx_struct_assembly_prop", d)
        self.assertEqual(d["pdbx_struct_oper_list"], [{
            "id": "1", "type": "?", "name": "?", "symmetry_operation": "x,y,z",
            "matrix[1][1]": "1.0", "matrix[1][2]": "0.0", "matrix[1][3]": "0.0",
            "vector[1]": "0.0", "matrix[2][1]": "0.0", "matrix[2][2]": "1.0",
            "matrix[2][3]": "0.0", "vector[2]": "0.0", "matrix[3][1]": "0.0",
            "matrix[3][2]": "0.0", "matrix[3][3]": "1.0", "vector[3]": "0.0",
        }])


    def test_2ldc(self):
        d = atomium.open("tests/integration/files/2ldc.pdb", dictionary=True)
        self.assertEqual(len(d.keys()), 28)

        # Metadata categories
        self.assertEqual(d["entry"], [{"id": "2LDC"}])
        self.assertEqual(d["pdbx_database_status"], [
            {"status_code": "REL", "entry_id": "2LDC", "recvd_initial_deposition_date": "2011-05-20"},
        ])
        self.assertEqual(d["struct_keywords"], [
            {"entry_id": "2LDC", "pdbx_keywords": "DE NOVO PROTEIN", "text": "DE NOVO PROTEIN"},
        ])
        self.assertNotIn("pdbx_database_PDB_obs_spr", d)
        self.assertEqual(d["struct"], [{
            "entry_id": "2LDC",
            "title": "SOLUTION STRUCTURE OF THE ESTROGEN RECEPTOR-BINDING STAPLED PEPTIDE SP1 (AC-HXILHXLLQDS-NH2)",
            "pdbx_descriptor": "?", "pdbx_model_details": "?",
            "pdbx_CASP_flag": "?", "pdbx_model_type_details": "?",
        }])
        self.assertNotIn("pdbx_database_related", d)
        self.assertNotIn("database_PDB_caveat", d)
        self.assertEqual(d["exptl"], [{"entry_id": "2LDC", "method": "SOLUTION NMR"}])
        self.assertNotIn("pdbx_coordinate_model", d)
        self.assertEqual(d["pdbx_audit_revision_history"], [{
            "ordinal": "1", "data_content_type": "Structure model",
            "major_revision": "1", "minor_revision": "0",
            "revision_date": "2011-07-06",
        }])

        # Entity categories
        self.assertEqual(d["entity"], [{
            "id": "1", "type": "polymer", "src_method": "syn",
            "pdbx_description": "ESTROGEN RECEPTOR-BINDING STAPLED PEPTIDE SP1",
            "formula_weight": "?", "pdbx_number_of_molecules": "1", "pdbx_ec": "?",
            "pdbx_mutation": "?", "pdbx_fragment": "?", "details": "?",
        }])
        self.assertNotIn("entity_name_com", d)
        self.assertEqual(d["entity_poly"], [{
            "entity_id": "1", "type": "polypeptide(L)", "nstd_linkage": "no",
            "nstd_monomer": "yes",
            "pdbx_seq_one_letter_code": "(ACE)H(MK8)ILH(MK8)LLQDS(NH2)",
            "pdbx_seq_one_letter_code_can": "XHLILHLLLQDSX", "pdbx_strand_id": "A",
            "pdbx_target_identifier": "?",
        }])
        self.assertEqual(d["entity_poly_seq"], [
            {"entity_id": "1", "num": "1", "mon_id": "ACE", "hetero": "n"},
            {"entity_id": "1", "num": "2", "mon_id": "HIS", "hetero": "n"},
            {"entity_id": "1", "num": "3", "mon_id": "MK8", "hetero": "n"},
            {"entity_id": "1", "num": "4", "mon_id": "ILE", "hetero": "n"},
            {"entity_id": "1", "num": "5", "mon_id": "LEU", "hetero": "n"},
            {"entity_id": "1", "num": "6", "mon_id": "HIS", "hetero": "n"},
            {"entity_id": "1", "num": "7", "mon_id": "MK8", "hetero": "n"},
            {"entity_id": "1", "num": "8", "mon_id": "LEU", "hetero": "n"},
            {"entity_id": "1", "num": "9", "mon_id": "LEU", "hetero": "n"},
            {"entity_id": "1", "num": "10", "mon_id": "GLN", "hetero": "n"},
            {"entity_id": "1", "num": "11", "mon_id": "ASP", "hetero": "n"},
            {"entity_id": "1", "num": "12", "mon_id": "SER", "hetero": "n"},
            {"entity_id": "1", "num": "13", "mon_id": "NH2", "hetero": "n"},
        ])
        self.assertEqual(d["struct_ref"], [{
            "id": "1", "db_name": "PDB", "db_code": "2LDC",
            "pdbx_db_accession": "2LDC", "pdbx_db_isoform": "?", "entity_id": "1",
            "pdbx_seq_one_letter_code": "?", "pdbx_align_begin": "?",
        }])
        self.assertEqual(d["struct_ref_seq"], [{
            "align_id": "1", "ref_id": "1", "pdbx_PDB_id_code": "2LDC",
            "pdbx_strand_id": "A", "seq_align_beg": "?",
            "pdbx_seq_align_beg_ins_code": "?", "seq_align_end": "?",
            "pdbx_seq_align_end_ins_code": "?", "pdbx_db_accession": "2LDC",
            "db_align_beg": "0", "pdbx_db_align_beg_ins_code": "?",
            "db_align_end": "12", "pdbx_db_align_end_ins_code": "?",
            "pdbx_auth_seq_align_beg": "0", "pdbx_auth_seq_align_end": "12",
        }])
        self.assertNotIn("struct_ref_seq_dif", d)
        self.assertEqual(d["pdbx_struct_mod_residue"], [{
            "id": "1", "label_asym_id": "A", "label_comp_id": "MK8",
            "label_seq_id": "3", "auth_asym_id": "A", "auth_comp_id": "MK8",
            "auth_seq_id": "2", "PDB_ins_code": "?", "parent_comp_id": "LEU",
            "details": "2-METHYL-L-NORLEUCINE",
        }, {
            "id": "2", "label_asym_id": "A", "label_comp_id": "MK8",
            "label_seq_id": "7", "auth_asym_id": "A", "auth_comp_id": "MK8",
            "auth_seq_id": "6", "PDB_ins_code": "?", "parent_comp_id": "LEU",
            "details": "2-METHYL-L-NORLEUCINE",
        }])
        self.assertEqual(len(d["chem_comp"]), 9)
        self.assertEqual(d["chem_comp"][0], {
            "id": "ACE", "type": "L-peptide linking", "mon_nstd_flag": "n",
            "name": "ACETYL GROUP", "pdbx_synonyms": "?", "formula": "C2 H4 O",
            "formula_weight": "44.052",
        })
        self.assertEqual(d["chem_comp"][6], {
            "id": "MK8", "type": "L-peptide linking", "mon_nstd_flag": "n",
            "name": "2-METHYL-L-NORLEUCINE", "pdbx_synonyms": "?",
            "formula": "C7 H15 N O2", "formula_weight": "145.199",
        })
        self.assertEqual(d["chem_comp"][7], {
            "id": "NH2", "type": "L-peptide linking", "mon_nstd_flag": "n",
            "name": "AMINO GROUP", "pdbx_synonyms": "?", "formula": "H2 N",
            "formula_weight": "16.023",
        })
        self.assertEqual(d["chem_comp"][-1], {
            "id": "SER", "type": "L-peptide linking", "mon_nstd_flag": "y",
            "name": "SERINE", "pdbx_synonyms": "?", "formula": "C3 H7 N O3",
            "formula_weight": "115.130",
        })
        self.assertNotIn("pdbx_entity_nonpoly", d)
        self.assertEqual(d["struct_asym"], [
            {"id": "A", "pdbx_blank_PDB_chainid_flag": "N", "pdbx_modified": "N", "entity_id": "1", "details": "?"},
        ])

        # Atom categories
        self.assertEqual(d["atom_type"], [{"symbol": "C"}, {"symbol": "H"}, {"symbol": "N"}, {"symbol": "O"}])
        self.assertEqual(len(d["atom_site"]), 2010)
        self.assertEqual(d["atom_site"][0], {
            "group_PDB": "HETATM", "id": "1", "type_symbol": "C",
            "label_atom_id": "C", "label_alt_id": ".", "label_comp_id": "ACE",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-3.378", "Cartn_y": "-3.101",
            "Cartn_z": "7.021", "occupancy": "1.00", "B_iso_or_equiv": "0.00",
            "pdbx_formal_charge": "?", "auth_seq_id": "0", "auth_comp_id": "ACE",
            "auth_asym_id": "A", "auth_atom_id": "C", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][201], {
            "group_PDB": "HETATM", "id": "202", "type_symbol": "C",
            "label_atom_id": "C", "label_alt_id": ".", "label_comp_id": "ACE",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-3.166", "Cartn_y": "-2.498",
            "Cartn_z": "6.924", "occupancy": "1.00", "B_iso_or_equiv": "0.00",
            "pdbx_formal_charge": "?", "auth_seq_id": "0", "auth_comp_id": "ACE",
            "auth_asym_id": "A", "auth_atom_id": "C", "pdbx_PDB_model_num": "2",
        })
        self.assertEqual(d["atom_site"][402], {
            "group_PDB": "HETATM", "id": "403", "type_symbol": "C",
            "label_atom_id": "C", "label_alt_id": ".", "label_comp_id": "ACE",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-3.269", "Cartn_y": "-2.772",
            "Cartn_z": "6.874", "occupancy": "1.00", "B_iso_or_equiv": "0.00",
            "pdbx_formal_charge": "?", "auth_seq_id": "0", "auth_comp_id": "ACE",
            "auth_asym_id": "A", "auth_atom_id": "C", "pdbx_PDB_model_num": "3",
        })
        self.assertEqual(d["atom_site"][603], {
            "group_PDB": "HETATM", "id": "604", "type_symbol": "C",
            "label_atom_id": "C", "label_alt_id": ".", "label_comp_id": "ACE",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-3.290", "Cartn_y": "-2.872",
            "Cartn_z": "7.516", "occupancy": "1.00", "B_iso_or_equiv": "0.00",
            "pdbx_formal_charge": "?", "auth_seq_id": "0", "auth_comp_id": "ACE",
            "auth_asym_id": "A", "auth_atom_id": "C", "pdbx_PDB_model_num": "4",
        })
        self.assertEqual(d["atom_site"][804], {
            "group_PDB": "HETATM", "id": "805", "type_symbol": "C",
            "label_atom_id": "C", "label_alt_id": ".", "label_comp_id": "ACE",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-2.757", "Cartn_y": "-2.786",
            "Cartn_z": "7.770", "occupancy": "1.00", "B_iso_or_equiv": "0.00",
            "pdbx_formal_charge": "?", "auth_seq_id": "0", "auth_comp_id": "ACE",
            "auth_asym_id": "A", "auth_atom_id": "C", "pdbx_PDB_model_num": "5",
        })
        self.assertEqual(d["atom_site"][1005], {
            "group_PDB": "HETATM", "id": "1006", "type_symbol": "C",
            "label_atom_id": "C", "label_alt_id": ".", "label_comp_id": "ACE",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-2.885", "Cartn_y": "-2.441",
            "Cartn_z": "7.465", "occupancy": "1.00", "B_iso_or_equiv": "0.00",
            "pdbx_formal_charge": "?", "auth_seq_id": "0", "auth_comp_id": "ACE",
            "auth_asym_id": "A", "auth_atom_id": "C", "pdbx_PDB_model_num": "6",
        })
        self.assertEqual(d["atom_site"][1206], {
            "group_PDB": "HETATM", "id": "1207", "type_symbol": "C",
            "label_atom_id": "C", "label_alt_id": ".", "label_comp_id": "ACE",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-3.097", "Cartn_y": "-2.892",
            "Cartn_z": "7.467", "occupancy": "1.00", "B_iso_or_equiv": "0.00",
            "pdbx_formal_charge": "?", "auth_seq_id": "0", "auth_comp_id": "ACE",
            "auth_asym_id": "A", "auth_atom_id": "C", "pdbx_PDB_model_num": "7",
        })
        self.assertEqual(d["atom_site"][1407], {
            "group_PDB": "HETATM", "id": "1408", "type_symbol": "C",
            "label_atom_id": "C", "label_alt_id": ".", "label_comp_id": "ACE",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-3.373", "Cartn_y": "-2.603",
            "Cartn_z": "6.700", "occupancy": "1.00", "B_iso_or_equiv": "0.00",
            "pdbx_formal_charge": "?", "auth_seq_id": "0", "auth_comp_id": "ACE",
            "auth_asym_id": "A", "auth_atom_id": "C", "pdbx_PDB_model_num": "8",
        })
        self.assertEqual(d["atom_site"][1608], {
            "group_PDB": "HETATM", "id": "1609", "type_symbol": "C",
            "label_atom_id": "C", "label_alt_id": ".", "label_comp_id": "ACE",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-3.230", "Cartn_y": "-2.552",
            "Cartn_z": "6.751", "occupancy": "1.00", "B_iso_or_equiv": "0.00",
            "pdbx_formal_charge": "?", "auth_seq_id": "0", "auth_comp_id": "ACE",
            "auth_asym_id": "A", "auth_atom_id": "C", "pdbx_PDB_model_num": "9",
        })
        self.assertEqual(d["atom_site"][1809], {
            "group_PDB": "HETATM", "id": "1810", "type_symbol": "C",
            "label_atom_id": "C", "label_alt_id": ".", "label_comp_id": "ACE",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-3.371", "Cartn_y": "-3.096",
            "Cartn_z": "6.994", "occupancy": "1.00", "B_iso_or_equiv": "0.00",
            "pdbx_formal_charge": "?", "auth_seq_id": "0", "auth_comp_id": "ACE",
            "auth_asym_id": "A", "auth_atom_id": "C", "pdbx_PDB_model_num": "10",
        })
        self.assertEqual(d["atom_site"][-1], {
            "group_PDB": "HETATM", "id": "2010", "type_symbol": "H",
            "label_atom_id": "HN2", "label_alt_id": ".", "label_comp_id": "NH2",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "13",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "2.870", "Cartn_y": "3.746",
            "Cartn_z": "-9.328", "occupancy": "1.00", "B_iso_or_equiv": "0.00",
            "pdbx_formal_charge": "?", "auth_seq_id": "12", "auth_comp_id": "NH2",
            "auth_asym_id": "A", "auth_atom_id": "HN2", "pdbx_PDB_model_num": "10",
        })
        self.assertNotIn("atom_site_anisotrop", d)

        # Annotation categories
        self.assertEqual(d["struct_conf"], [{
            "conf_type_id": "HELX_P", "id": "HELX_P1", "pdbx_PDB_helix_id": "1",
            "beg_label_comp_id": "HIS", "beg_label_asym_id": "A",
            "beg_label_seq_id": "2", "pdbx_beg_PDB_ins_code": "?",
            "end_label_comp_id": "GLN", "end_label_asym_id": "A",
            "end_label_seq_id": "10", "pdbx_end_PDB_ins_code": "?",
            "beg_auth_comp_id": "HIS", "beg_auth_asym_id": "A",
            "beg_auth_seq_id": "1", "end_auth_comp_id": "GLN",
            "end_auth_asym_id": "A", "end_auth_seq_id": "9",
            "pdbx_PDB_helix_class": "1", "details": "?",
            "pdbx_PDB_helix_length": "9",
        }])
        self.assertNotIn("struct_sheet", d)
        self.assertNotIn("struct_sheet_range", d)
        self.assertNotIn("struct_sheet_order", d)
        self.assertNotIn("pdbx_struct_sheet_hbond", d)
        self.assertEqual(len(d["struct_conn"]), 7)
        self.assertEqual(d["struct_conn"][0], {
            "id": "covale1", "conn_type_id": "covale",
            "pdbx_leaving_atom_flag": "?", "pdbx_PDB_id": "?",
            "ptnr1_label_asym_id": "A", "ptnr1_label_comp_id": "ACE",
            "ptnr1_label_seq_id": "1", "ptnr1_label_atom_id": "C",
            "pdbx_ptnr1_label_alt_id": "?", "pdbx_ptnr1_PDB_ins_code": "?",
            "pdbx_ptnr1_standard_comp_id": "?", "ptnr1_symmetry": "1555",
            "ptnr2_label_asym_id": "A", "ptnr2_label_comp_id": "HIS",
            "ptnr2_label_seq_id": "2", "ptnr2_label_atom_id": "N",
            "pdbx_ptnr2_label_alt_id": "?", "pdbx_ptnr2_PDB_ins_code": "?",
            "ptnr1_auth_asym_id": "A", "ptnr1_auth_comp_id": "ACE",
            "ptnr1_auth_seq_id": "0", "ptnr2_auth_asym_id": "A",
            "ptnr2_auth_comp_id": "HIS", "ptnr2_auth_seq_id": "1",
            "ptnr2_symmetry": "1555", "pdbx_ptnr3_label_atom_id": "?",
            "pdbx_ptnr3_label_seq_id": "?", "pdbx_ptnr3_label_comp_id": "?",
            "pdbx_ptnr3_label_asym_id": "?", "pdbx_ptnr3_label_alt_id": "?",
            "pdbx_ptnr3_PDB_ins_code": "?", "details": "?",
            "pdbx_dist_value": "1.38", "pdbx_value_order": "?",
        })
        self.assertEqual(d["struct_conn"][-1], {
            "id": "covale7", "conn_type_id": "covale",
            "pdbx_leaving_atom_flag": "?", "pdbx_PDB_id": "?",
            "ptnr1_label_asym_id": "A", "ptnr1_label_comp_id": "MK8",
            "ptnr1_label_seq_id": "3", "ptnr1_label_atom_id": "CE",
            "pdbx_ptnr1_label_alt_id": "?", "pdbx_ptnr1_PDB_ins_code": "?",
            "pdbx_ptnr1_standard_comp_id": "?", "ptnr1_symmetry": "1555",
            "ptnr2_label_asym_id": "A", "ptnr2_label_comp_id": "MK8",
            "ptnr2_label_seq_id": "7", "ptnr2_label_atom_id": "CE",
            "pdbx_ptnr2_label_alt_id": "?", "pdbx_ptnr2_PDB_ins_code": "?",
            "ptnr1_auth_asym_id": "A", "ptnr1_auth_comp_id": "MK8",
            "ptnr1_auth_seq_id": "2", "ptnr2_auth_asym_id": "A",
            "ptnr2_auth_comp_id": "MK8", "ptnr2_auth_seq_id": "6",
            "ptnr2_symmetry": "1555", "pdbx_ptnr3_label_atom_id": "?",
            "pdbx_ptnr3_label_seq_id": "?", "pdbx_ptnr3_label_comp_id": "?",
            "pdbx_ptnr3_label_asym_id": "?", "pdbx_ptnr3_label_alt_id": "?",
            "pdbx_ptnr3_PDB_ins_code": "?", "details": "?",
            "pdbx_dist_value": "1.39", "pdbx_value_order": "?",
        })
        self.assertNotIn("struct_mon_prot_cis", d)
        self.assertNotIn("struct_site", d)
        self.assertNotIn("struct_site_gen", d)

        # Crystal categories
        self.assertEqual(d["cell"], [{
            "entry_id": "2LDC", "length_a": "1.000", "length_b": "1.000",
            "length_c": "1.000", "angle_alpha": "90.00", "angle_beta": "90.00",
            "angle_gamma": "90.00", "Z_pdb": "1", "pdbx_unique_axis": "?",
        }])
        self.assertEqual(d["symmetry"], [{
            "entry_id": "2LDC", "space_group_name_H-M": "P 1",
            "pdbx_full_space_group_name_H-M": "?", "cell_setting": "?",
            "Int_Tables_number": "?",
        }])
        self.assertEqual(d["atom_sites"], [{
            "entry_id": "2LDC", "fract_transf_matrix[1][1]": "1.000000",
            "fract_transf_matrix[1][2]": "0.000000",
            "fract_transf_matrix[1][3]": "0.000000",
            "fract_transf_matrix[2][1]": "0.000000",
            "fract_transf_matrix[2][2]": "1.000000",
            "fract_transf_matrix[2][3]": "0.000000",
            "fract_transf_matrix[3][1]": "0.000000",
            "fract_transf_matrix[3][2]": "0.000000",
            "fract_transf_matrix[3][3]": "1.000000",
            "fract_transf_vector[1]": "0.00000",
            "fract_transf_vector[2]": "0.00000",
            "fract_transf_vector[3]": "0.00000",
        }])
        self.assertEqual(d["database_PDB_matrix"], [{
            "entry_id": "2LDC", "origx[1][1]": "1.000000",
            "origx[1][2]": "0.000000", "origx[1][3]": "0.000000",
            "origx[2][1]": "0.000000", "origx[2][2]": "1.000000",
            "origx[2][3]": "0.000000", "origx[3][1]": "0.000000",
            "origx[3][2]": "0.000000", "origx[3][3]": "1.000000",
            "origx_vector[1]": "0.00000", "origx_vector[2]": "0.00000",
            "origx_vector[3]": "0.00000",
        }])
        self.assertNotIn("struct_ncs_oper", d)

        # Citation categories
        self.assertEqual(d["audit_author"], [
            {"name": "Phillips, C", "pdbx_ordinal": "1"},
            {"name": "Bazin, R", "pdbx_ordinal": "2"},
            {"name": "Bent, A", "pdbx_ordinal": "3"},
            {"name": "Davies, N", "pdbx_ordinal": "4"},
            {"name": "Moore, R", "pdbx_ordinal": "5"},
            {"name": "Pannifer, A", "pdbx_ordinal": "6"},
            {"name": "Pickford, A", "pdbx_ordinal": "7"},
            {"name": "Prior, S", "pdbx_ordinal": "8"},
            {"name": "Read, C", "pdbx_ordinal": "9"},
            {"name": "Roberts, L", "pdbx_ordinal": "10"},
            {"name": "Schade, M", "pdbx_ordinal": "11"},
            {"name": "Scott, A", "pdbx_ordinal": "12"},
            {"name": "Brown, D", "pdbx_ordinal": "13"},
            {"name": "Xu, B", "pdbx_ordinal": "14"},
            {"name": "Irving, S", "pdbx_ordinal": "15"},
        ])
        self.assertEqual(len(d["citation_author"]), 15)
        self.assertEqual(d["citation_author"][0], {
            "citation_id": "primary", "name": "Phillips, C", "pdbx_ordinal": "1",
        })
        self.assertEqual(d["citation_author"][-1], {
            "citation_id": "primary", "name": "Irving, S.L", "pdbx_ordinal": "15",
        })
        self.assertNotIn("citation_editor", d)
        self.assertEqual(d["citation"], [{
            "id": "primary",
            "title": "DESIGN AND STRUCTURE OF STAPLED PEPTIDES BINDING TO ESTROGEN RECEPTORS.",
            "journal_abbrev": "J.Am.Chem.Soc.", "journal_volume": "133",
            "page_first": "9696", "page_last": "?", "year": "2011",
            "journal_id_ASTM": "?", "country": "?", "journal_id_ISSN": "0002-7863",
            "journal_id_CSD": "?", "book_publisher": "?",
            "pdbx_database_id_PubMed": "21612236",
            "pdbx_database_id_DOI": "10.1021/ja202946k",
        }])

        # Experimental categories
        self.assertNotIn("reflns", d)
        self.assertNotIn("refine", d)

        # Missing items categories
        self.assertNotIn("pdbx_unobs_or_zero_occ_residues", d)

        # Assembly categories
        self.assertEqual(d["pdbx_struct_assembly"], [
            {"id": "1", "details": "?", "method_details": "?", "oligomeric_details": "?", "oligomeric_count": "?"},
        ])
        self.assertEqual(d["pdbx_struct_assembly_gen"], [{"assembly_id": "1", "oper_expression": "1", "asym_id_list": "A"}])
        self.assertNotIn("pdbx_struct_assembly_prop", d)
        self.assertEqual(d["pdbx_struct_oper_list"], [{
            "id": "1", "type": "?", "name": "?", "symmetry_operation": "x,y,z",
            "matrix[1][1]": "1.0", "matrix[1][2]": "0.0", "matrix[1][3]": "0.0",
            "vector[1]": "0.0", "matrix[2][1]": "0.0", "matrix[2][2]": "1.0",
            "matrix[2][3]": "0.0", "vector[2]": "0.0", "matrix[3][1]": "0.0",
            "matrix[3][2]": "0.0", "matrix[3][3]": "1.0", "vector[3]": "0.0",
        }])


    '''

    def test_1CK8_pdb(self):
        # Tests OBSLTE, CAVEAT, AUTHOR
        d = atomium.open("tests/integration/files/1ck8.pdb", dictionary=True)
        self.assertEqual(d["pdbx_database_PDB_obs_spr"], [{
            "id": "OBSLTE", "date": "2006-07-25", "pdb_id": "1T0K",
            "replace_pdb_id": "1CK8", "details": "?"
        }])
        self.assertEqual(d["database_PDB_caveat"], [{
            "id": "1", "text": "THERE ARE CHIRALITY ERRORS IN C-ALPHA CENTERS"
        }])
        self.assertEqual(d["audit_author"], [
            {"name": "Mao, H", "pdbx_ordinal": "1"},
            {"name": "White, S.A", "pdbx_ordinal": "2"},
            {"name": "Willamson, J.R", "pdbx_ordinal": "3"},
        ])
    

    def test_4HHB_pdb(self):
        # Tests SPRSDE, CAVEAT, REVDAT
        d = atomium.open("tests/integration/files/4hhb.pdb", dictionary=True)
        self.assertEqual(d["pdbx_database_PDB_obs_spr"], [{
            "id": "SPRSDE", "date": "1984-07-17", "pdb_id": "4HHB",
            "replace_pdb_id": "1HHB", "details": "?"
        }])
        self.assertEqual(d["database_PDB_caveat"], [{
            "id": "1",
            "text": "THR A 137 HAS WRONG CHIRALITY AT ATOM CB "
            "THR B 12 HAS WRONG CHIRALITY AT ATOM CB "
            "THR B 50 HAS WRONG CHIRALITY AT ATOM CB "
            "ASN C 78 HAS WRONG CHIRALITY AT ATOM CA "
            "THR C 118 HAS WRONG CHIRALITY AT ATOM CB "
            "HIS D 2 HAS WRONG CHIRALITY AT ATOM CA "
            "SER D 72 HAS WRONG CHIRALITY AT ATOM CA "
            "ASP D 73 HAS WRONG CHIRALITY AT ATOM CA "
            "LEU D 78 HAS WRONG CHIRALITY AT ATOM CA "
            "LYS D 144 HAS WRONG CHIRALITY AT ATOM CA"
        }])
        self.assertEqual(d["pdbx_audit_revision_history"], [{
           "ordinal": "1", "data_content_type": "Structure model",
           "major_revision": "1", "minor_revision": "0", "revision_date": "1984-07-17"
        }, {
            "ordinal": "2", "data_content_type": "Structure model",
            "major_revision": "2", "minor_revision": "0", "revision_date": "1989-10-15"
        }, {
            "ordinal": "3", "data_content_type": "Structure model",
            "major_revision": "3", "minor_revision": "0", "revision_date": "2003-04-01"
        }, {
            "ordinal": "4", "data_content_type": "Structure model",
            "major_revision": "4", "minor_revision": "0", "revision_date": "2009-02-24"
        }, {
            "ordinal": "5", "data_content_type": "Structure model",
            "major_revision": "5", "minor_revision": "0", "revision_date": "2011-07-13"
        }, {
            "ordinal": "6", "data_content_type": "Structure model",
            "major_revision": "6", "minor_revision": "0", "revision_date": "2020-06-17"
        }, {
            "ordinal": "7", "data_content_type": "Structure model",
            "major_revision": "7", "minor_revision": "0", "revision_date": "2021-03-31"
        }])
    

    def test_1GDJ_pdb(self):
        # Tests SPRSDE
        d = atomium.open("tests/integration/files/1gdj.pdb", dictionary=True)
        self.assertEqual(d["pdbx_database_PDB_obs_spr"], [{
            "id": "SPRSDE", "date": "1995-02-27", "pdb_id": "1GDJ",
            "replace_pdb_id": "1LH4 2LH4", "details": "?"
        }])
    

    def test_4CWU_pdb(self):
        # Tests OBSLTE, SPRSDE
        d = atomium.open("tests/integration/files/4cwu.pdb", dictionary=True)
        self.assertEqual(d["pdbx_database_PDB_obs_spr"], [{
            "id": "SPRSDE", "date": "2014-08-06", "pdb_id": "4CWU",
            "replace_pdb_id": "1VSZ", "details": "?"
        }, {
            "id": "OBSLTE", "date": "2018-05-02", "pdb_id": "6CGV",
            "replace_pdb_id": "4CWU", "details": "?"
        }])
    

    def test_1VUA_pdb(self):
        # Tests SPLIT
        d = atomium.open("tests/integration/files/1vua.pdb", dictionary=True)
        self.assertEqual(d["pdbx_database_related"], [{
            "db_name": "PDB", "db_id": code, "content_type": "split",
            "details": f"Split {n}"
        } for n, code in enumerate([
            "1VU4", "1VU5", "1VU6", "1VU7", "1VU8", "1VU9", "1VUA", "1VUC",
            "1VUD", "1VUE", "1VUF", "1VUG", "1VUH", "1VUI", "1VUJ", "1VUK",
            "1VUL", "1VUM", "1VUN", "1VUO", "1VUP", "1VUQ", "1VUR", "1VUS", "1VUT"
        ], start=1)])
    

    def test_2CSE_pdb(self):
        # Tests MDLTYP
        d = atomium.open("tests/integration/files/2cse.pdb", dictionary=True)
        self.assertEqual(d["struct"], [{
            "entry_id": "2CSE", "pdbx_CASP_flag": "?", "pdbx_descriptor": "?",
            "pdbx_model_details": "?",
            "pdbx_model_type_details": "?",
            "title": "FEATURES OF REOVIRUS OUTER-CAPSID PROTEIN MU1 REVEALED BY "
            "ELECTRON AND IMAGE RECONSTRUCTION OF THE VIRION AT 7.0-A RESOLUTION"
        }])
        self.assertEqual(d["pdbx_coordinate_model"], [{
            "asym_id": char, "type": "CA ATOMS ONLY"
        } for char in "ABCSDEFMNOGHIPQRJKLTUVWXYZ1"
        ])
    

    def test_1GIY_pdb(self):
        # Tests MDLTYP
        d = atomium.open("tests/integration/files/1giy.pdb", dictionary=True)
        self.assertEqual(d["struct"][0]["pdbx_model_type_details"], "?")
        self.assertEqual(d["pdbx_coordinate_model"], [{
            "asym_id": char, "type": "CA ATOMS ONLY"
        } for char in "CDEFGHIJKLMNOPQRSTUVWX"] + [{
            "asym_id": char, "type": "P ATOMS ONLY"
        } for char in "AB"])
    

    def test_2JRQ_pdb(self):
        # Tests MDLTYP
        d = atomium.open("tests/integration/files/2jrq.pdb", dictionary=True)
        self.assertEqual(d["struct"], [{
            "entry_id": "2JRQ", "pdbx_CASP_flag": "?", "pdbx_descriptor": "?",
            "pdbx_model_details": "?",
            "pdbx_model_type_details": "MINIMIZED AVERAGE",
            "title": "NMR SOLUTION STRUCTURE OF THE ANTICODON OF E. COLI TRNA-VAL3 WITH 1 MODIFICATION (CMO5U34)"
        }])
        self.assertNotIn("pdbx_coordinate_model", d)
    

    def test_1YGT_pdb(self):
        # Tests REMARK 200
        d = atomium.open("tests/integration/files/1ygt.pdb", dictionary=True)
        self.assertEqual(d["diffrn"], [{
            "id": "1", "ambient_temp": "100", "ambient_temp_details": "?", "crystal_id": "1"
        }])
        self.assertEqual(d["diffrn_detector"], [{
            "diffrn_id": "1", "detector": "CCD", "type": "ADSC QUANTUM 4",
            "pdbx_collection_date": "2000-09-04", "details": "NULL"
        }])
    

    def test_1COJ_pdb(self):
        # Tests JRNL
        d = atomium.open("tests/integration/files/1coj.pdb", dictionary=True)
        self.assertEqual(d["citation"], [{
            "id": "primary",
            "title": "THE CRYSTAL STRUCTURE OF AN FE-SUPEROXIDE DISMUTASE FROM THE HYPERTHERMOPHILE AQUIFEX PYROPHILUS AT 1.9 A RESOLUTION: STRUCTURAL BASIS FOR THERMOSTABILITY.",
            "journal_abbrev": "J.Mol.Biol.", "journal_volume": "270",
            "page_first": "259", "page_last": "?", "year": "1997",
            "journal_id_ASTM": "?", "country": "?", "journal_id_ISSN": "0022-2836",
            "journal_id_CSD": "?", "book_publisher": "U.K. : ACADEMIC PRESS",
            "pdbx_database_id_PubMed": "9236127",
            "pdbx_database_id_DOI": "10.1006/jmbi.1997.1105"
        }])
        self.assertEqual(d["citation_author"], [
            {"citation_id": "primary", "name": "Lim, J.H", "pdbx_ordinal": "1"},
            {"citation_id": "primary", "name": "Yu, Y.G", "pdbx_ordinal": "2"},
            {"citation_id": "primary", "name": "Han, Y.S", "pdbx_ordinal": "3"},
            {"citation_id": "primary", "name": "Cho, S", "pdbx_ordinal": "4"},
            {"citation_id": "primary", "name": "Ahn, B.Y", "pdbx_ordinal": "5"},
            {"citation_id": "primary", "name": "Kim, S.H", "pdbx_ordinal": "6"},
            {"citation_id": "primary", "name": "Cho, Y", "pdbx_ordinal": "7"}
        ])
        self.assertEqual(d["citation_editor"], [
            {"citation_id": "primary", "name": "Peter Wright", "pdbx_ordinal": "1"},
        ])
    

    def test_2NQ2_pdb(self):
        # Tests REMARK 470,480
        d = atomium.open("tests/integration/files/2nq2.pdb", dictionary=True)
        self.assertEqual(len(d["pdbx_unobs_or_zero_occ_atoms"]), 129)
        self.assertEqual(d["pdbx_unobs_or_zero_occ_atoms"][0], {
            "id": "1", "PDB_model_num": "1", "polymer_flag": "Y", "occupancy_flag": "0",
            "auth_asym_id": "B", "auth_comp_id": "ILE", "auth_seq_id": "9", "PDB_ins_code": "?",
            "auth_atom_id": "CD1", "label_alt_id": "?", "label_asym_id": "B",
            "label_comp_id": "ILE", "label_seq_id": "9", "label_atom_id": "CD1"
        })
        self.assertEqual(d["pdbx_unobs_or_zero_occ_atoms"][-1], {
            "id": "129", "PDB_model_num": "1", "polymer_flag": "Y", "occupancy_flag": "0",
            "auth_asym_id": "B", "auth_comp_id": "LEU", "auth_seq_id": "330", "PDB_ins_code": "?",
            "auth_atom_id": "O", "label_alt_id": "?", "label_asym_id": "B",
            "label_comp_id": "LEU", "label_seq_id": "330", "label_atom_id": "O"
        })
    

    def test_1NCF_pdb(self):
        # Tests MTRIXn
        d = atomium.open("tests/integration/files/1ncf.pdb", dictionary=True)
        self.assertEqual(d["struct_ncs_oper"], [{
            "id": "1", "code": "given", "details": "?",
            "matrix[1][1]": "-0.745000", "matrix[1][2]": "0.008000", "matrix[1][3]": "0.667000",
            "matrix[2][1]": "0.020000", "matrix[2][2]": "-0.999000", "matrix[2][3]": "0.035000",
            "matrix[3][1]": "0.667000", "matrix[3][2]": "0.040000", "matrix[3][3]": "0.744000",
            "vector[1]": "14.56700", "vector[2]": "27.10700", "vector[3]": "-6.04300"
        }, {
            "id": "2", "code": "given", "details": "?",
            "matrix[1][1]": "-0.746000", "matrix[1][2]": "-0.029000", "matrix[1][3]": "0.665000",
            "matrix[2][1]": "0.045000", "matrix[2][2]": "-0.999000", "matrix[2][3]": "0.007000",
            "matrix[3][1]": "0.664000", "matrix[3][2]": "0.035000", "matrix[3][3]": "0.747000",
            "vector[1]": "14.68300", "vector[2]": "27.24500", "vector[3]": "-5.93200"
        }, {
            "id": "3", "code": "given", "details": "?",
            "matrix[1][1]": "-0.814000", "matrix[1][2]": "-0.073000", "matrix[1][3]": "0.576000",
            "matrix[2][1]": "0.147000", "matrix[2][2]": "-0.986000", "matrix[2][3]": "0.083000",
            "matrix[3][1]": "0.562000", "matrix[3][2]": "0.153000", "matrix[3][3]": "0.813000",
            "vector[1]": "19.59000", "vector[2]": "22.20100", "vector[3]": "-7.71900"
        }, {
            "id": "4", "code": "given", "details": "?",
            "matrix[1][1]": "-0.862000", "matrix[1][2]": "0.043000", "matrix[1][3]": "0.505000",
            "matrix[2][1]": "0.102000", "matrix[2][2]": "-0.961000", "matrix[2][3]": "0.255000",
            "matrix[3][1]": "0.496000", "matrix[3][2]": "0.272000", "matrix[3][3]": "0.825000",
            "vector[1]": "23.91900", "vector[2]": "11.13100", "vector[3]": "-7.60600"
        }])
    

    def test_1GSG_pdb(self):
        # Tests MODRES
        d = atomium.open("tests/integration/files/1gsg.pdb", dictionary=True)
        self.assertEqual(d["entry"], [{"id": "1GSG"}])
        self.assertEqual(d["pdbx_struct_mod_residue"][0], {
            "id": "1", "label_asym_id": "?", "label_comp_id": "4SU", "label_seq_id": "?",
            "auth_asym_id": "T", "auth_comp_id": "4SU", "auth_seq_id": "8", "PDB_ins_code": "?",
            "parent_comp_id": "U", "details": "4-THIOURIDINE-5'-MONOPHOSPHATE"
        })
    

    def test_4y60_pdb(self):
        # Test ANISOU
        d = atomium.open("tests/integration/files/4y60.pdb", dictionary=True)
        
        self.assertEqual(d["entry"], [{"id": "4Y60"}])

        self.assertEqual(d["atom_site"][0], {
            "group_PDB": "ATOM", "id": "1", "type_symbol": "N", "label_atom_id": "N",
            "label_alt_id": ".", "label_comp_id": "GLY", "label_asym_id": "A",
            "label_entity_id": "1", "label_seq_id": "1", "pdbx_PDB_ins_code": "?",
            "Cartn_x": "43.447", "Cartn_y": "-56.622", "Cartn_z": "-20.561",
            "occupancy": "1.00", "B_iso_or_equiv": "56.53", "pdbx_formal_charge": "?",
            "auth_seq_id": "0", "auth_comp_id": "GLY", "auth_asym_id": "C",
            "auth_atom_id": "N", "pdbx_PDB_model_num": "1"
        })
        self.assertEqual(d["atom_site_anisotrop"][0], {
            "id": "1", "type_symbol": "N", "pdbx_label_atom_id": "N",
            "pdbx_label_alt_id": ".", "pdbx_label_comp_id": "GLY",
            "pdbx_label_asym_id": "A", "pdbx_label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "U[1][1]": "0.9838", "U[2][2]": "0.7489",
            "U[3][3]": "0.4152", "U[1][2]": "-0.1159", "U[1][3]": "-0.0115",
            "U[2][3]": "-0.2655", "pdbx_auth_seq_id": "0", "pdbx_auth_comp_id": "GLY",
            "pdbx_auth_asym_id": "C", "pdbx_auth_atom_id ": "N"
        })
        self.assertEqual(d["atom_site_anisotrop"][-1], {
            "id": "1319", "type_symbol": "C", "pdbx_label_atom_id": "C4",
            "pdbx_label_alt_id": ".", "pdbx_label_comp_id": "DG",
            "pdbx_label_asym_id": "C", "pdbx_label_seq_id": "16",
            "pdbx_PDB_ins_code": "?", "U[1][1]": "0.4728", "U[2][2]": "0.6687",
            "U[3][3]": "0.5867", "U[1][2]": "0.2574", "U[1][3]": "-0.0037",
            "U[2][3]": "-0.1117", "pdbx_auth_seq_id": "15", "pdbx_auth_comp_id": "DG",
            "pdbx_auth_asym_id": "B", "pdbx_auth_atom_id ": "C4"
        })
    

    def test_5XME_pdb(self):
        # Test multiple models
        d = atomium.open("tests/integration/files/5xme.pdb", dictionary=True)
        
        self.assertEqual(d["entry"], [{"id": "5XME"}])

        
        self.assertEqual(d["atom_site"][0], {
            "group_PDB": "ATOM", "id": "1", "type_symbol": "N", "label_atom_id": "N",
            "label_alt_id": ".", "label_comp_id": "ALA", "label_asym_id": "A",
            "label_entity_id": "1", "label_seq_id": "14", "pdbx_PDB_ins_code": "?",
            "Cartn_x": "33.969", "Cartn_y": "-8.430", "Cartn_z": "-0.271",
            "occupancy": "1.00", "B_iso_or_equiv": "0.00", "pdbx_formal_charge": "?",
            "auth_seq_id": "199", "auth_comp_id": "ALA", "auth_asym_id": "A",
            "auth_atom_id": "N", "pdbx_PDB_model_num": "1"
        })

        self.assertEqual(d["atom_site"][1827], {
            "group_PDB": "ATOM", "id": "1828", "type_symbol": "N", "label_atom_id": "N",
            "label_alt_id": ".", "label_comp_id": "ALA", "label_asym_id": "A",
            "label_entity_id": "1", "label_seq_id": "14", "pdbx_PDB_ins_code": "?",
            "Cartn_x": "34.064", "Cartn_y": "-8.092", "Cartn_z": "-0.062",
            "occupancy": "1.00", "B_iso_or_equiv": "0.00", "pdbx_formal_charge": "?",
            "auth_seq_id": "199", "auth_comp_id": "ALA", "auth_asym_id": "A",
            "auth_atom_id": "N", "pdbx_PDB_model_num": "2"
        })



class FullParsingTests(TestCase):
    
    def test_1LOL_pdb(self):
        pdb = atomium.open("tests/integration/files/1lol.pdb")
        self.assertEqual(pdb.name, "1LOL")
        self.assertEqual(pdb.title, "CRYSTAL STRUCTURE OF OROTIDINE MONOPHOSPHATE DECARBOXYLASE COMPLEX WITH XMP")
        self.assertEqual(pdb.keywords, ["TIM BARREL", "LYASE"])
        self.assertEqual(str(pdb.model), "<Model (2 polymers, 4 non-polymers)>")'''



class DictToFileTests(TestCase):

    def setUp(self):
        if os.path.exists("tests/integration/files/output/"):
            shutil.rmtree("tests/integration/files/output/")
        os.mkdir("tests/integration/files/output")


    def tearDown(self):
        if os.path.exists("tests/integration/files/output/"):
            shutil.rmtree("tests/integration/files/output/")
    

    def save(self, code):
        original = atomium.open(f"tests/integration/files/{code}.pdb", dictionary=True)
        atomium.save_dictionary(original, f"tests/integration/files/output/{code}.pdb")
        saved = atomium.open(f"tests/integration/files/output/{code}.pdb", dictionary=True)
        self.assertEqual(original.keys(), saved.keys())
        for key in original:
            self.assertEqual(original[key], saved[key])
        self.assertEqual(original, saved)
    

    def test_1lol(self):
        self.save("1lol")
    

    def test_5xme(self):
        self.save("5xme")
    

    def test_1gsg(self):
        self.save("1gsg")