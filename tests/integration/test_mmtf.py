import os
import shutil
from unittest import TestCase
import atomium

class FileToDictTests(TestCase):

    def test_1lol(self):
        d = atomium.open("tests/integration/files/1lol.mmtf", dictionary=True)
        self.assertEqual(len(d.keys()), 18)

        # Header categories
        self.assertEqual(d["entry"], [{"id": "1LOL"}])
        self.assertEqual(d["struct"], [
            {"entry_id": "1LOL", "title": "Crystal structure of orotidine monophosphate decarboxylase complex with XMP"},
        ])
        self.assertEqual(d["pdbx_database_status"], [{"entry_id": "1LOL", "recvd_initial_deposition_date": "2002-05-06"}])
        self.assertEqual(d["exptl"], [{"entry_id": "1LOL", "method": "X-RAY DIFFRACTION"}])

        # Quality categories
        self.assertEqual(d["refine"], [{
            "entry_id": "1LOL", "ls_d_res_high": "1.9", "ls_R_factor_obs": "0.193",
            "ls_R_factor_all": "0.193", "ls_R_factor_R_work": "0.193",
            "ls_R_factor_R_free": "0.229",
        }])

        # Crystal categories
        self.assertEqual(d["symmetry"], [{"entry_id": "1LOL", "space_group_name_H-M": "P 1 21 1"}])
        self.assertEqual(d["cell"], [{
            "entry_id": "1LOL", "length_a": "57.57", "length_b": "55.482",
            "length_c": "66.129", "length_alpha": "90.0", "length_beta": "94.28",
            "length_gamma": "90.0",
        }])

        # Assembly categories
        self.assertEqual(d["pdbx_struct_assembly"], [{"id": "1"}])
        self.assertEqual(d["pdbx_struct_assembly_gen"], [
            {"assembly_id": "1", "oper_expression": "1", "asym_id_list": "A,B,C,D,E,F,G,H"},
        ])
        self.assertEqual(d["pdbx_struct_oper_list"], [{
            "id": "1", "matrix[1][1]": "1.0000000000",
            "matrix[1][2]": "0.0000000000", "matrix[1][3]": "0.0000000000",
            "vector[1]": "0.0000000000", "matrix[2][1]": "0.0000000000",
            "matrix[2][2]": "1.0000000000", "matrix[2][3]": "0.0000000000",
            "vector[2]": "0.0000000000", "matrix[3][1]": "0.0000000000",
            "matrix[3][2]": "0.0000000000", "matrix[3][3]": "1.0000000000",
            "vector[3]": "0.0000000000",
        }])

        # Entity categories
        self.assertEqual(d["entity"], [{
            "id": "1", "type": "polymer",
            "pdbx_description": "orotidine 5'-monophosphate decarboxylase",
            "pdbx_number_of_molecules": "2",
        }, {
            "id": "2", "type": "non-polymer", "pdbx_description": "1,3-BUTANEDIOL",
            "pdbx_number_of_molecules": "2",
        }, {
            "id": "3", "type": "non-polymer",
            "pdbx_description": "XANTHOSINE-5'-MONOPHOSPHATE",
            "pdbx_number_of_molecules": "2",
        }, {
            "id": "4", "type": "water", "pdbx_description": "water",
            "pdbx_number_of_molecules": "2",
        }])
        self.assertEqual(d["entity_poly"], [{
            "entity_id": "1",
            "pdbx_seq_one_letter_code": "LRSRRVDVMDVMNRLILAMDLMNRDDALRVTGEVREYIDTVKIGYPLVLSEGMDIIAEFRKRFGCRIIADFKVADIPETNEKICRATFKAGADAIIVHGFPGADSVRACLNVAEEMGREVFLLTEMSHPGAEMFIQGAADEIARMGVDLGVKNYVGPSTRPERLSRLREIIGQDSFLISPGVGAQGGDPGETLRFADAIIVGRSIYLADNPAAAAAGIIESIKDLLIPE",
            "pdbx_seq_one_letter_code_can": "LRSRRVDVMDVMNRLILAMDLMNRDDALRVTGEVREYIDTVKIGYPLVLSEGMDIIAEFRKRFGCRIIADFKVADIPETNEKICRATFKAGADAIIVHGFPGADSVRACLNVAEEMGREVFLLTEMSHPGAEMFIQGAADEIARMGVDLGVKNYVGPSTRPERLSRLREIIGQDSFLISPGVGAQGGDPGETLRFADAIIVGRSIYLADNPAAAAAGIIESIKDLLIPE",
        }])
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

        # Compound categories
        self.assertEqual(len(d["chem_comp"]), 22)
        self.assertEqual(d["chem_comp"][0], {
            "id": "ALA", "type": "L-PEPTIDE LINKING", "mon_nstd_flag": "y",
            "name": "alanine", "pdbx_synonyms": "?", "formula": "?",
            "formula_weight": "?",
        })
        self.assertEqual(d["chem_comp"][4], {
            "id": "BU2", "type": "NON-POLYMER", "mon_nstd_flag": ".", "name": "?",
            "pdbx_synonyms": "?", "formula": "?", "formula_weight": "?",
        })
        self.assertEqual(d["chem_comp"][10], {
            "id": "HOH", "type": "NON-POLYMER", "mon_nstd_flag": ".", "name": "?",
            "pdbx_synonyms": "?", "formula": "?", "formula_weight": "?",
        })
        self.assertEqual(d["chem_comp"][-1], {
            "id": "XMP", "type": "NON-POLYMER", "mon_nstd_flag": ".", "name": "?",
            "pdbx_synonyms": "?", "formula": "?", "formula_weight": "?",
        })

        # Model categories
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
            "Cartn_z": "63.219", "occupancy": "1.0", "B_iso_or_equiv": "21.5",
            "pdbx_formal_charge": "0", "auth_seq_id": "11", "auth_comp_id": "VAL",
            "auth_asym_id": "A", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1557], {
            "group_PDB": "ATOM", "id": "1558", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "VAL",
            "label_asym_id": "B", "label_entity_id": "1", "label_seq_id": "1011",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-26.384", "Cartn_y": "61.433",
            "Cartn_z": "36.898", "occupancy": "1.0", "B_iso_or_equiv": "39.3",
            "pdbx_formal_charge": "0", "auth_seq_id": "1011", "auth_comp_id": "VAL",
            "auth_asym_id": "B", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][3191], {
            "group_PDB": "HETATM", "id": "3192", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "BU2",
            "label_asym_id": "C", "label_entity_id": "2", "label_seq_id": "5001",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "2.646", "Cartn_y": "45.112",
            "Cartn_z": "48.995", "occupancy": "1.0", "B_iso_or_equiv": "43.24",
            "pdbx_formal_charge": "0", "auth_seq_id": "5001", "auth_comp_id": "BU2",
            "auth_asym_id": "A", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][3197], {
            "group_PDB": "HETATM", "id": "3198", "type_symbol": "P",
            "label_atom_id": "P", "label_alt_id": ".", "label_comp_id": "XMP",
            "label_asym_id": "D", "label_entity_id": "3", "label_seq_id": "2001",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "3.293", "Cartn_y": "36.948",
            "Cartn_z": "44.605", "occupancy": "1.0", "B_iso_or_equiv": "20.47",
            "pdbx_formal_charge": "0", "auth_seq_id": "2001", "auth_comp_id": "XMP",
            "auth_asym_id": "A", "auth_atom_id": "P", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][3221], {
            "group_PDB": "HETATM", "id": "3222", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "BU2",
            "label_asym_id": "E", "label_entity_id": "2", "label_seq_id": "5002",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-14.563", "Cartn_y": "61.208",
            "Cartn_z": "49.005", "occupancy": "1.0", "B_iso_or_equiv": "45.5",
            "pdbx_formal_charge": "0", "auth_seq_id": "5002", "auth_comp_id": "BU2",
            "auth_asym_id": "B", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][3227], {
            "group_PDB": "HETATM", "id": "3228", "type_symbol": "P",
            "label_atom_id": "P", "label_alt_id": ".", "label_comp_id": "XMP",
            "label_asym_id": "F", "label_entity_id": "3", "label_seq_id": "2002",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-23.846", "Cartn_y": "63.007",
            "Cartn_z": "53.02", "occupancy": "1.0", "B_iso_or_equiv": "33.88",
            "pdbx_formal_charge": "0", "auth_seq_id": "2002", "auth_comp_id": "XMP",
            "auth_asym_id": "B", "auth_atom_id": "P", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][3251], {
            "group_PDB": "HETATM", "id": "3252", "type_symbol": "O",
            "label_atom_id": "O", "label_alt_id": ".", "label_comp_id": "HOH",
            "label_asym_id": "G", "label_entity_id": "4", "label_seq_id": "3005",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-7.624", "Cartn_y": "43.71",
            "Cartn_z": "47.691", "occupancy": "1.0", "B_iso_or_equiv": "15.72",
            "pdbx_formal_charge": "0", "auth_seq_id": "3005", "auth_comp_id": "HOH",
            "auth_asym_id": "A", "auth_atom_id": "O", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][3347], {
            "group_PDB": "HETATM", "id": "3348", "type_symbol": "O",
            "label_atom_id": "O", "label_alt_id": ".", "label_comp_id": "HOH",
            "label_asym_id": "H", "label_entity_id": "4", "label_seq_id": "3001",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-9.22", "Cartn_y": "50.8",
            "Cartn_z": "49.092", "occupancy": "1.0", "B_iso_or_equiv": "17.69",
            "pdbx_formal_charge": "0", "auth_seq_id": "3001", "auth_comp_id": "HOH",
            "auth_asym_id": "B", "auth_atom_id": "O", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][-1], {
            "group_PDB": "HETATM", "id": "3431", "type_symbol": "O",
            "label_atom_id": "O", "label_alt_id": ".", "label_comp_id": "HOH",
            "label_asym_id": "H", "label_entity_id": "4", "label_seq_id": "3180",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-39.239", "Cartn_y": "51.357",
            "Cartn_z": "40.064", "occupancy": "1.0", "B_iso_or_equiv": "37.93",
            "pdbx_formal_charge": "0", "auth_seq_id": "3180", "auth_comp_id": "HOH",
            "auth_asym_id": "B", "auth_atom_id": "O", "pdbx_PDB_model_num": "1",
        })
    
    
    def test_1xda(self):
        d = atomium.open("tests/integration/files/1xda.mmtf", dictionary=True)
        self.assertEqual(len(d.keys()), 18)

        # Header categories
        self.assertEqual(d["entry"], [{"id": "1XDA"}])
        self.assertEqual(d["struct"], [{"entry_id": "1XDA", "title": "STRUCTURE OF INSULIN"}])
        self.assertEqual(d["pdbx_database_status"], [{"entry_id": "1XDA", "recvd_initial_deposition_date": "1996-12-18"}])
        self.assertEqual(d["exptl"], [{"entry_id": "1XDA", "method": "X-RAY DIFFRACTION"}])

        # Quality categories
        self.assertEqual(d["refine"], [{
            "entry_id": "1XDA", "ls_d_res_high": "1.8", "ls_R_factor_obs": "0.174",
            "ls_R_factor_all": "0.174", "ls_R_factor_R_work": "0.174",
            "ls_R_factor_R_free": "?",
        }])

        # Crystal categories
        self.assertEqual(d["symmetry"], [{"entry_id": "1XDA", "space_group_name_H-M": "H 3"}])
        self.assertEqual(d["cell"], [{
            "entry_id": "1XDA", "length_a": "78.752", "length_b": "78.752",
            "length_c": "79.199", "length_alpha": "90.0", "length_beta": "90.0",
            "length_gamma": "120.0",
        }])

        # Assembly categories
        self.assertEqual(d["pdbx_struct_assembly"], [
            {"id": "1"},
            {"id": "2"},
            {"id": "3"},
            {"id": "4"},
            {"id": "5"},
            {"id": "6"},
            {"id": "7"},
            {"id": "8"},
            {"id": "9"},
            {"id": "10"},
            {"id": "11"},
            {"id": "12"},
        ])
        self.assertEqual(len(d["pdbx_struct_assembly_gen"]), 24)
        self.assertEqual(d["pdbx_struct_assembly_gen"][0], {
            "assembly_id": "1", "oper_expression": "1",
            "asym_id_list": "A,B,I,J,K,L,Y,Z",
        })
        self.assertEqual(d["pdbx_struct_assembly_gen"][-1], {
            "assembly_id": "12", "oper_expression": "24",
            "asym_id_list": "A,B,C,D,I,J,K,L,M,N,O,P,Y,Z,AA,BA",
        })
        self.assertEqual(len(d["pdbx_struct_oper_list"]), 24)
        self.assertEqual(d["pdbx_struct_oper_list"][0], {
            "id": "1", "matrix[1][1]": "1.0000000000",
            "matrix[1][2]": "0.0000000000", "matrix[1][3]": "0.0000000000",
            "vector[1]": "0.0000000000", "matrix[2][1]": "0.0000000000",
            "matrix[2][2]": "1.0000000000", "matrix[2][3]": "0.0000000000",
            "vector[2]": "0.0000000000", "matrix[3][1]": "0.0000000000",
            "matrix[3][2]": "0.0000000000", "matrix[3][3]": "1.0000000000",
            "vector[3]": "0.0000000000",
        })
        self.assertEqual(d["pdbx_struct_oper_list"][-1], {
            "id": "24", "matrix[1][1]": "1.0000000000",
            "matrix[1][2]": "0.0000000000", "matrix[1][3]": "0.0000000000",
            "vector[1]": "0.0000000000", "matrix[2][1]": "0.0000000000",
            "matrix[2][2]": "1.0000000000", "matrix[2][3]": "0.0000000000",
            "vector[2]": "0.0000000000", "matrix[3][1]": "0.0000000000",
            "matrix[3][2]": "0.0000000000", "matrix[3][3]": "1.0000000000",
            "vector[3]": "0.0000000000",
        })

        # Entity categories
        self.assertEqual(d["entity"], [
            {"id": "1", "type": "polymer", "pdbx_description": "FATTY ACID ACYLATED INSULIN", "pdbx_number_of_molecules": "4"},
            {"id": "2", "type": "polymer", "pdbx_description": "FATTY ACID ACYLATED INSULIN", "pdbx_number_of_molecules": "4"},
            {"id": "3", "type": "non-polymer", "pdbx_description": "PHENOL", "pdbx_number_of_molecules": "4"},
            {"id": "4", "type": "non-polymer", "pdbx_description": "ZINC ION", "pdbx_number_of_molecules": "4"},
            {"id": "5", "type": "non-polymer", "pdbx_description": "CHLORIDE ION", "pdbx_number_of_molecules": "4"},
            {"id": "6", "type": "non-polymer", "pdbx_description": "MYRISTIC ACID", "pdbx_number_of_molecules": "4"},
            {"id": "7", "type": "water", "pdbx_description": "water", "pdbx_number_of_molecules": "8"},
        ])
        self.assertEqual(d["entity_poly"], [{
            "entity_id": "1", "pdbx_seq_one_letter_code": "GIVEQCCTSICSLYQLENYCN",
            "pdbx_seq_one_letter_code_can": "GIVEQCCTSICSLYQLENYCN",
        }, {
            "entity_id": "2",
            "pdbx_seq_one_letter_code": "FVNQHLCGSHLVEALYLVCGERGFFYTPK",
            "pdbx_seq_one_letter_code_can": "FVNQHLCGSHLVEALYLVCGERGFFYTPK",
        }])
        self.assertEqual(len(d["struct_asym"]), 32)
        self.assertEqual(d["struct_asym"][0], {
            "id": "A", "pdbx_blank_PDB_chainid_flag": "N", "pdbx_modified": "N",
            "entity_id": "1", "details": "?",
        })
        self.assertEqual(d["struct_asym"][-1], {
            "id": "FA", "pdbx_blank_PDB_chainid_flag": "N", "pdbx_modified": "N",
            "entity_id": "7", "details": "?",
        })

        # Compound categories
        self.assertEqual(len(d["chem_comp"]), 22)
        self.assertEqual(d["chem_comp"][0], {
            "id": "ALA", "type": "L-PEPTIDE LINKING", "mon_nstd_flag": "y",
            "name": "alanine", "pdbx_synonyms": "?", "formula": "?",
            "formula_weight": "?",
        })
        self.assertEqual(d["chem_comp"][3], {
            "id": "CL", "type": "NON-POLYMER", "mon_nstd_flag": ".", "name": "?",
            "pdbx_synonyms": "?", "formula": "?", "formula_weight": "?",
        })
        self.assertEqual(d["chem_comp"][9], {
            "id": "HOH", "type": "NON-POLYMER", "mon_nstd_flag": ".", "name": "?",
            "pdbx_synonyms": "?", "formula": "?", "formula_weight": "?",
        })
        self.assertEqual(d["chem_comp"][11], {
            "id": "IPH", "type": "NON-POLYMER", "mon_nstd_flag": ".", "name": "?",
            "pdbx_synonyms": "?", "formula": "?", "formula_weight": "?",
        })
        self.assertEqual(d["chem_comp"][14], {
            "id": "MYR", "type": "NON-POLYMER", "mon_nstd_flag": ".", "name": "?",
            "pdbx_synonyms": "?", "formula": "?", "formula_weight": "?",
        })
        self.assertEqual(d["chem_comp"][-1], {
            "id": "ZN", "type": "NON-POLYMER", "mon_nstd_flag": ".", "name": "?",
            "pdbx_synonyms": "?", "formula": "?", "formula_weight": "?",
        })

        # Model categories
        self.assertEqual(d["atom_type"], [
            {"symbol": "C"},
            {"symbol": "Cl"},
            {"symbol": "N"},
            {"symbol": "O"},
            {"symbol": "S"},
            {"symbol": "Zn"},
        ])
        self.assertEqual(len(d["atom_site"]), 1860)
        self.assertEqual(d["atom_site"][0], {
            "group_PDB": "ATOM", "id": "1", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "GLY",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "16.136", "Cartn_y": "11.725",
            "Cartn_z": "12.834", "occupancy": "1.0", "B_iso_or_equiv": "17.65",
            "pdbx_formal_charge": "0", "auth_seq_id": "1", "auth_comp_id": "GLY",
            "auth_asym_id": "A", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][163], {
            "group_PDB": "ATOM", "id": "164", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "PHE",
            "label_asym_id": "B", "label_entity_id": "2", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-0.033", "Cartn_y": "11.515",
            "Cartn_z": "9.458", "occupancy": "1.0", "B_iso_or_equiv": "19.34",
            "pdbx_formal_charge": "0", "auth_seq_id": "1", "auth_comp_id": "PHE",
            "auth_asym_id": "B", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][405], {
            "group_PDB": "ATOM", "id": "406", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "GLY",
            "label_asym_id": "C", "label_entity_id": "1", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "9.561", "Cartn_y": "16.902",
            "Cartn_z": "40.4", "occupancy": "1.0", "B_iso_or_equiv": "25.15",
            "pdbx_formal_charge": "0", "auth_seq_id": "1", "auth_comp_id": "GLY",
            "auth_asym_id": "C", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][568], {
            "group_PDB": "ATOM", "id": "569", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "PHE",
            "label_asym_id": "D", "label_entity_id": "2", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "11.9", "Cartn_y": "0.573",
            "Cartn_z": "42.957", "occupancy": "1.0", "B_iso_or_equiv": "26.7",
            "pdbx_formal_charge": "0", "auth_seq_id": "1", "auth_comp_id": "PHE",
            "auth_asym_id": "D", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][806], {
            "group_PDB": "ATOM", "id": "807", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "GLY",
            "label_asym_id": "E", "label_entity_id": "1", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "17.335", "Cartn_y": "9.341",
            "Cartn_z": "52.401", "occupancy": "1.0", "B_iso_or_equiv": "19.16",
            "pdbx_formal_charge": "0", "auth_seq_id": "1", "auth_comp_id": "GLY",
            "auth_asym_id": "E", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][969], {
            "group_PDB": "ATOM", "id": "970", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "PHE",
            "label_asym_id": "F", "label_entity_id": "2", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "1.175", "Cartn_y": "11.882",
            "Cartn_z": "49.858", "occupancy": "1.0", "B_iso_or_equiv": "28.7",
            "pdbx_formal_charge": "0", "auth_seq_id": "1", "auth_comp_id": "PHE",
            "auth_asym_id": "F", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1207], {
            "group_PDB": "ATOM", "id": "1208", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "GLY",
            "label_asym_id": "G", "label_entity_id": "1", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "12.506", "Cartn_y": "15.301",
            "Cartn_z": "80.111", "occupancy": "1.0", "B_iso_or_equiv": "22.75",
            "pdbx_formal_charge": "0", "auth_seq_id": "1", "auth_comp_id": "GLY",
            "auth_asym_id": "G", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1373], {
            "group_PDB": "ATOM", "id": "1374", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "PHE",
            "label_asym_id": "H", "label_entity_id": "2", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "11.577", "Cartn_y": "-0.737",
            "Cartn_z": "83.108", "occupancy": "1.0", "B_iso_or_equiv": "17.85",
            "pdbx_formal_charge": "0", "auth_seq_id": "1", "auth_comp_id": "PHE",
            "auth_asym_id": "H", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1610], {
            "group_PDB": "HETATM", "id": "1611", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "IPH",
            "label_asym_id": "I", "label_entity_id": "3", "label_seq_id": "22",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-5.661", "Cartn_y": "9.282",
            "Cartn_z": "16.887", "occupancy": "1.0", "B_iso_or_equiv": "13.71",
            "pdbx_formal_charge": "0", "auth_seq_id": "22", "auth_comp_id": "IPH",
            "auth_asym_id": "A", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1617], {
            "group_PDB": "HETATM", "id": "1618", "type_symbol": "Zn",
            "label_atom_id": "ZN", "label_alt_id": ".", "label_comp_id": "ZN",
            "label_asym_id": "J", "label_entity_id": "4", "label_seq_id": "30",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "0.0", "Cartn_y": "0.0",
            "Cartn_z": "18.779", "occupancy": "0.33", "B_iso_or_equiv": "10.69",
            "pdbx_formal_charge": "2", "auth_seq_id": "30", "auth_comp_id": "ZN",
            "auth_asym_id": "B", "auth_atom_id": "ZN", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1618], {
            "group_PDB": "HETATM", "id": "1619", "type_symbol": "Cl",
            "label_atom_id": "CL", "label_alt_id": ".", "label_comp_id": "CL",
            "label_asym_id": "K", "label_entity_id": "5", "label_seq_id": "31",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "0.0", "Cartn_y": "0.0",
            "Cartn_z": "16.611", "occupancy": "0.33", "B_iso_or_equiv": "12.22",
            "pdbx_formal_charge": "-1", "auth_seq_id": "31", "auth_comp_id": "CL",
            "auth_asym_id": "B", "auth_atom_id": "CL", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1619], {
            "group_PDB": "HETATM", "id": "1620", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "MYR",
            "label_asym_id": "L", "label_entity_id": "6", "label_seq_id": "39",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "7.835", "Cartn_y": "13.475",
            "Cartn_z": "11.645", "occupancy": "0.0", "B_iso_or_equiv": "20.0",
            "pdbx_formal_charge": "0", "auth_seq_id": "39", "auth_comp_id": "MYR",
            "auth_asym_id": "B", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1634], {
            "group_PDB": "HETATM", "id": "1635", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "IPH",
            "label_asym_id": "M", "label_entity_id": "3", "label_seq_id": "22",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "9.806", "Cartn_y": "-4.594",
            "Cartn_z": "35.946", "occupancy": "1.0", "B_iso_or_equiv": "18.13",
            "pdbx_formal_charge": "0", "auth_seq_id": "22", "auth_comp_id": "IPH",
            "auth_asym_id": "C", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1641], {
            "group_PDB": "HETATM", "id": "1642", "type_symbol": "Zn",
            "label_atom_id": "ZN", "label_alt_id": ".", "label_comp_id": "ZN",
            "label_asym_id": "N", "label_entity_id": "4", "label_seq_id": "30",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "0.0", "Cartn_y": "0.0",
            "Cartn_z": "34.196", "occupancy": "0.33", "B_iso_or_equiv": "12.34",
            "pdbx_formal_charge": "2", "auth_seq_id": "30", "auth_comp_id": "ZN",
            "auth_asym_id": "D", "auth_atom_id": "ZN", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1642], {
            "group_PDB": "HETATM", "id": "1643", "type_symbol": "Cl",
            "label_atom_id": "CL", "label_alt_id": ".", "label_comp_id": "CL",
            "label_asym_id": "O", "label_entity_id": "5", "label_seq_id": "31",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "0.0", "Cartn_y": "0.0",
            "Cartn_z": "36.366", "occupancy": "0.33", "B_iso_or_equiv": "11.37",
            "pdbx_formal_charge": "-1", "auth_seq_id": "31", "auth_comp_id": "CL",
            "auth_asym_id": "D", "auth_atom_id": "CL", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1643], {
            "group_PDB": "HETATM", "id": "1644", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "MYR",
            "label_asym_id": "P", "label_entity_id": "6", "label_seq_id": "39",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "13.145", "Cartn_y": "5.879",
            "Cartn_z": "44.308", "occupancy": "1.0", "B_iso_or_equiv": "43.86",
            "pdbx_formal_charge": "0", "auth_seq_id": "39", "auth_comp_id": "MYR",
            "auth_asym_id": "D", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1658], {
            "group_PDB": "HETATM", "id": "1659", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "IPH",
            "label_asym_id": "Q", "label_entity_id": "3", "label_seq_id": "22",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-4.433", "Cartn_y": "9.971",
            "Cartn_z": "56.719", "occupancy": "1.0", "B_iso_or_equiv": "15.57",
            "pdbx_formal_charge": "0", "auth_seq_id": "22", "auth_comp_id": "IPH",
            "auth_asym_id": "E", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1665], {
            "group_PDB": "HETATM", "id": "1666", "type_symbol": "Zn",
            "label_atom_id": "ZN", "label_alt_id": ".", "label_comp_id": "ZN",
            "label_asym_id": "R", "label_entity_id": "4", "label_seq_id": "30",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "0.0", "Cartn_y": "0.0",
            "Cartn_z": "58.639", "occupancy": "0.33", "B_iso_or_equiv": "11.99",
            "pdbx_formal_charge": "2", "auth_seq_id": "30", "auth_comp_id": "ZN",
            "auth_asym_id": "F", "auth_atom_id": "ZN", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1666], {
            "group_PDB": "HETATM", "id": "1667", "type_symbol": "Cl",
            "label_atom_id": "CL", "label_alt_id": ".", "label_comp_id": "CL",
            "label_asym_id": "S", "label_entity_id": "5", "label_seq_id": "31",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "0.0", "Cartn_y": "0.0",
            "Cartn_z": "56.498", "occupancy": "0.33", "B_iso_or_equiv": "13.93",
            "pdbx_formal_charge": "-1", "auth_seq_id": "31", "auth_comp_id": "CL",
            "auth_asym_id": "F", "auth_atom_id": "CL", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1667], {
            "group_PDB": "HETATM", "id": "1668", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "MYR",
            "label_asym_id": "T", "label_entity_id": "6", "label_seq_id": "39",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "12.146", "Cartn_y": "13.183",
            "Cartn_z": "45.515", "occupancy": "0.0", "B_iso_or_equiv": "20.0",
            "pdbx_formal_charge": "0", "auth_seq_id": "39", "auth_comp_id": "MYR",
            "auth_asym_id": "F", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1682], {
            "group_PDB": "HETATM", "id": "1683", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "IPH",
            "label_asym_id": "U", "label_entity_id": "3", "label_seq_id": "22",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "8.823", "Cartn_y": "-6.017",
            "Cartn_z": "75.855", "occupancy": "1.0", "B_iso_or_equiv": "15.74",
            "pdbx_formal_charge": "0", "auth_seq_id": "22", "auth_comp_id": "IPH",
            "auth_asym_id": "G", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1689], {
            "group_PDB": "HETATM", "id": "1690", "type_symbol": "Zn",
            "label_atom_id": "ZN", "label_alt_id": ".", "label_comp_id": "ZN",
            "label_asym_id": "V", "label_entity_id": "4", "label_seq_id": "30",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "0.0", "Cartn_y": "0.0",
            "Cartn_z": "74.024", "occupancy": "0.33", "B_iso_or_equiv": "14.68",
            "pdbx_formal_charge": "2", "auth_seq_id": "30", "auth_comp_id": "ZN",
            "auth_asym_id": "H", "auth_atom_id": "ZN", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1690], {
            "group_PDB": "HETATM", "id": "1691", "type_symbol": "Cl",
            "label_atom_id": "CL", "label_alt_id": ".", "label_comp_id": "CL",
            "label_asym_id": "W", "label_entity_id": "5", "label_seq_id": "31",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "0.0", "Cartn_y": "0.0",
            "Cartn_z": "76.191", "occupancy": "0.33", "B_iso_or_equiv": "15.86",
            "pdbx_formal_charge": "-1", "auth_seq_id": "31", "auth_comp_id": "CL",
            "auth_asym_id": "H", "auth_atom_id": "CL", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1691], {
            "group_PDB": "HETATM", "id": "1692", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "MYR",
            "label_asym_id": "X", "label_entity_id": "6", "label_seq_id": "39",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "14.782", "Cartn_y": "6.834",
            "Cartn_z": "84.17", "occupancy": "1.0", "B_iso_or_equiv": "43.17",
            "pdbx_formal_charge": "0", "auth_seq_id": "39", "auth_comp_id": "MYR",
            "auth_asym_id": "H", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1706], {
            "group_PDB": "HETATM", "id": "1707", "type_symbol": "O",
            "label_atom_id": "O", "label_alt_id": ".", "label_comp_id": "HOH",
            "label_asym_id": "Y", "label_entity_id": "7", "label_seq_id": "23",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "17.639", "Cartn_y": "-5.979",
            "Cartn_z": "20.33", "occupancy": "1.0", "B_iso_or_equiv": "27.85",
            "pdbx_formal_charge": "0", "auth_seq_id": "23", "auth_comp_id": "HOH",
            "auth_asym_id": "A", "auth_atom_id": "O", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1721], {
            "group_PDB": "HETATM", "id": "1722", "type_symbol": "O",
            "label_atom_id": "O", "label_alt_id": ".", "label_comp_id": "HOH",
            "label_asym_id": "Z", "label_entity_id": "7", "label_seq_id": "40",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "13.817", "Cartn_y": "11.255",
            "Cartn_z": "20.074", "occupancy": "1.0", "B_iso_or_equiv": "18.8",
            "pdbx_formal_charge": "0", "auth_seq_id": "40", "auth_comp_id": "HOH",
            "auth_asym_id": "B", "auth_atom_id": "O", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1746], {
            "group_PDB": "HETATM", "id": "1747", "type_symbol": "O",
            "label_atom_id": "O", "label_alt_id": ".", "label_comp_id": "HOH",
            "label_asym_id": "AA", "label_entity_id": "7", "label_seq_id": "42",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "10.432", "Cartn_y": "18.286",
            "Cartn_z": "28.169", "occupancy": "1.0", "B_iso_or_equiv": "17.65",
            "pdbx_formal_charge": "0", "auth_seq_id": "42", "auth_comp_id": "HOH",
            "auth_asym_id": "C", "auth_atom_id": "O", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1762], {
            "group_PDB": "HETATM", "id": "1763", "type_symbol": "O",
            "label_atom_id": "O", "label_alt_id": ".", "label_comp_id": "HOH",
            "label_asym_id": "BA", "label_entity_id": "7", "label_seq_id": "40",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "-1.686", "Cartn_y": "4.748",
            "Cartn_z": "29.405", "occupancy": "1.0", "B_iso_or_equiv": "28.76",
            "pdbx_formal_charge": "0", "auth_seq_id": "40", "auth_comp_id": "HOH",
            "auth_asym_id": "D", "auth_atom_id": "O", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1779], {
            "group_PDB": "HETATM", "id": "1780", "type_symbol": "O",
            "label_atom_id": "O", "label_alt_id": ".", "label_comp_id": "HOH",
            "label_asym_id": "CA", "label_entity_id": "7", "label_seq_id": "23",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "19.132", "Cartn_y": "10.418",
            "Cartn_z": "56.864", "occupancy": "1.0", "B_iso_or_equiv": "17.46",
            "pdbx_formal_charge": "0", "auth_seq_id": "23", "auth_comp_id": "HOH",
            "auth_asym_id": "E", "auth_atom_id": "O", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1794], {
            "group_PDB": "HETATM", "id": "1795", "type_symbol": "O",
            "label_atom_id": "O", "label_alt_id": ".", "label_comp_id": "HOH",
            "label_asym_id": "DA", "label_entity_id": "7", "label_seq_id": "40",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "14.151", "Cartn_y": "4.653",
            "Cartn_z": "70.826", "occupancy": "1.0", "B_iso_or_equiv": "18.08",
            "pdbx_formal_charge": "0", "auth_seq_id": "40", "auth_comp_id": "HOH",
            "auth_asym_id": "F", "auth_atom_id": "O", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1810], {
            "group_PDB": "HETATM", "id": "1811", "type_symbol": "O",
            "label_atom_id": "O", "label_alt_id": ".", "label_comp_id": "HOH",
            "label_asym_id": "EA", "label_entity_id": "7", "label_seq_id": "23",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "5.232", "Cartn_y": "20.03",
            "Cartn_z": "75.975", "occupancy": "1.0", "B_iso_or_equiv": "25.28",
            "pdbx_formal_charge": "0", "auth_seq_id": "23", "auth_comp_id": "HOH",
            "auth_asym_id": "G", "auth_atom_id": "O", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1840], {
            "group_PDB": "HETATM", "id": "1841", "type_symbol": "O",
            "label_atom_id": "O", "label_alt_id": ".", "label_comp_id": "HOH",
            "label_asym_id": "FA", "label_entity_id": "7", "label_seq_id": "40",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "12.513", "Cartn_y": "14.38",
            "Cartn_z": "75.003", "occupancy": "1.0", "B_iso_or_equiv": "15.41",
            "pdbx_formal_charge": "0", "auth_seq_id": "40", "auth_comp_id": "HOH",
            "auth_asym_id": "H", "auth_atom_id": "O", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][-1], {
            "group_PDB": "HETATM", "id": "1860", "type_symbol": "O",
            "label_atom_id": "O", "label_alt_id": ".", "label_comp_id": "HOH",
            "label_asym_id": "FA", "label_entity_id": "7", "label_seq_id": "59",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "11.353", "Cartn_y": "5.616",
            "Cartn_z": "75.248", "occupancy": "1.0", "B_iso_or_equiv": "33.51",
            "pdbx_formal_charge": "0", "auth_seq_id": "59", "auth_comp_id": "HOH",
            "auth_asym_id": "H", "auth_atom_id": "O", "pdbx_PDB_model_num": "1",
        })


    def test_5xme(self):
        d = atomium.open("tests/integration/files/5xme.mmtf", dictionary=True)
        self.assertEqual(len(d.keys()), 16)

        # Header categories
        self.assertEqual(d["entry"], [{"id": "5XME"}])
        self.assertEqual(d["struct"], [
            {"entry_id": "5XME", "title": "Solution structure of C-terminal domain of TRADD"},
        ])
        self.assertEqual(d["pdbx_database_status"], [{"entry_id": "5XME", "recvd_initial_deposition_date": "2017-05-15"}])
        self.assertEqual(d["exptl"], [{"entry_id": "5XME", "method": "SOLUTION NMR"}])

        # Quality categories
        self.assertNotIn("refine", d)

        # Crystal categories
        self.assertEqual(d["symmetry"], [{"entry_id": "5XME", "space_group_name_H-M": "NA"}])
        self.assertNotIn("cell", d)

        # Assembly categories
        self.assertEqual(d["pdbx_struct_assembly"], [{"id": "1"}])
        self.assertEqual(d["pdbx_struct_assembly_gen"], [{"assembly_id": "1", "oper_expression": "1", "asym_id_list": "A"}])
        self.assertEqual(d["pdbx_struct_oper_list"], [{
            "id": "1", "matrix[1][1]": "1.0000000000",
            "matrix[1][2]": "0.0000000000", "matrix[1][3]": "0.0000000000",
            "vector[1]": "0.0000000000", "matrix[2][1]": "0.0000000000",
            "matrix[2][2]": "1.0000000000", "matrix[2][3]": "0.0000000000",
            "vector[2]": "0.0000000000", "matrix[3][1]": "0.0000000000",
            "matrix[3][2]": "0.0000000000", "matrix[3][3]": "1.0000000000",
            "vector[3]": "0.0000000000",
        }])

        # Entity categories
        self.assertEqual(d["entity"], [{
            "id": "1", "type": "polymer",
            "pdbx_description": "Tumor necrosis factor receptor type 1-associated DEATH domain protein",
            "pdbx_number_of_molecules": "1",
        }])
        self.assertEqual(d["entity_poly"], [{
            "entity_id": "1",
            "pdbx_seq_one_letter_code": "MHHHHHHSSGRGSAQTFLFQGQPVVNRPLSLKDQQTFARSVGLKWRKVGRSLQRGCRALRDPALDSLAYEYEREGLYEQAFQLLRRFVQAEGRRATLQRLVEALEENELTSLAEDLLGLTDPNGGLA",
            "pdbx_seq_one_letter_code_can": "MHHHHHHSSGRGSAQTFLFQGQPVVNRPLSLKDQQTFARSVGLKWRKVGRSLQRGCRALRDPALDSLAYEYEREGLYEQAFQLLRRFVQAEGRRATLQRLVEALEENELTSLAEDLLGLTDPNGGLA",
        }])
        self.assertEqual(len(d["struct_asym"]), 10)
        self.assertEqual(d["struct_asym"][0], {
            "id": "A", "pdbx_blank_PDB_chainid_flag": "N", "pdbx_modified": "N",
            "entity_id": "1", "details": "?",
        })
        self.assertEqual(d["struct_asym"][-1], {
            "id": "A", "pdbx_blank_PDB_chainid_flag": "N", "pdbx_modified": "N",
            "entity_id": "1", "details": "?",
        })

        # Compound categories
        self.assertEqual(len(d["chem_comp"]), 17)
        self.assertEqual(d["chem_comp"][0], {
            "id": "ALA", "type": "L-PEPTIDE LINKING", "mon_nstd_flag": "y",
            "name": "alanine", "pdbx_synonyms": "?", "formula": "?",
            "formula_weight": "?",
        })
        self.assertEqual(d["chem_comp"][-1], {
            "id": "VAL", "type": "L-PEPTIDE LINKING", "mon_nstd_flag": "y",
            "name": "valine", "pdbx_synonyms": "?", "formula": "?",
            "formula_weight": "?",
        })

        # Model categories
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
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "199",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "33.969", "Cartn_y": "-8.43",
            "Cartn_z": "-0.271", "occupancy": "1.0", "B_iso_or_equiv": "0.0",
            "pdbx_formal_charge": "0", "auth_seq_id": "199", "auth_comp_id": "ALA",
            "auth_asym_id": "A", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][1827], {
            "group_PDB": "ATOM", "id": "1828", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "ALA",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "199",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "34.064", "Cartn_y": "-8.092",
            "Cartn_z": "-0.062", "occupancy": "1.0", "B_iso_or_equiv": "0.0",
            "pdbx_formal_charge": "0", "auth_seq_id": "199", "auth_comp_id": "ALA",
            "auth_asym_id": "A", "auth_atom_id": "N", "pdbx_PDB_model_num": "2",
        })
        self.assertEqual(d["atom_site"][3654], {
            "group_PDB": "ATOM", "id": "3655", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "ALA",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "199",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "37.369", "Cartn_y": "-7.973",
            "Cartn_z": "-0.242", "occupancy": "1.0", "B_iso_or_equiv": "0.0",
            "pdbx_formal_charge": "0", "auth_seq_id": "199", "auth_comp_id": "ALA",
            "auth_asym_id": "A", "auth_atom_id": "N", "pdbx_PDB_model_num": "3",
        })
        self.assertEqual(d["atom_site"][5481], {
            "group_PDB": "ATOM", "id": "5482", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "ALA",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "199",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "36.023", "Cartn_y": "-7.429",
            "Cartn_z": "-0.637", "occupancy": "1.0", "B_iso_or_equiv": "0.0",
            "pdbx_formal_charge": "0", "auth_seq_id": "199", "auth_comp_id": "ALA",
            "auth_asym_id": "A", "auth_atom_id": "N", "pdbx_PDB_model_num": "4",
        })
        self.assertEqual(d["atom_site"][7308], {
            "group_PDB": "ATOM", "id": "7309", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "ALA",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "199",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "35.245", "Cartn_y": "-9.094",
            "Cartn_z": "0.245", "occupancy": "1.0", "B_iso_or_equiv": "0.0",
            "pdbx_formal_charge": "0", "auth_seq_id": "199", "auth_comp_id": "ALA",
            "auth_asym_id": "A", "auth_atom_id": "N", "pdbx_PDB_model_num": "5",
        })
        self.assertEqual(d["atom_site"][9135], {
            "group_PDB": "ATOM", "id": "9136", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "ALA",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "199",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "35.835", "Cartn_y": "-7.648",
            "Cartn_z": "-0.888", "occupancy": "1.0", "B_iso_or_equiv": "0.0",
            "pdbx_formal_charge": "0", "auth_seq_id": "199", "auth_comp_id": "ALA",
            "auth_asym_id": "A", "auth_atom_id": "N", "pdbx_PDB_model_num": "6",
        })
        self.assertEqual(d["atom_site"][10962], {
            "group_PDB": "ATOM", "id": "10963", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "ALA",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "199",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "37.525", "Cartn_y": "-7.759",
            "Cartn_z": "-0.297", "occupancy": "1.0", "B_iso_or_equiv": "0.0",
            "pdbx_formal_charge": "0", "auth_seq_id": "199", "auth_comp_id": "ALA",
            "auth_asym_id": "A", "auth_atom_id": "N", "pdbx_PDB_model_num": "7",
        })
        self.assertEqual(d["atom_site"][12789], {
            "group_PDB": "ATOM", "id": "12790", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "ALA",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "199",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "35.062", "Cartn_y": "-9.22",
            "Cartn_z": "0.031", "occupancy": "1.0", "B_iso_or_equiv": "0.0",
            "pdbx_formal_charge": "0", "auth_seq_id": "199", "auth_comp_id": "ALA",
            "auth_asym_id": "A", "auth_atom_id": "N", "pdbx_PDB_model_num": "8",
        })
        self.assertEqual(d["atom_site"][14616], {
            "group_PDB": "ATOM", "id": "14617", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "ALA",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "199",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "36.244", "Cartn_y": "-7.381",
            "Cartn_z": "-0.745", "occupancy": "1.0", "B_iso_or_equiv": "0.0",
            "pdbx_formal_charge": "0", "auth_seq_id": "199", "auth_comp_id": "ALA",
            "auth_asym_id": "A", "auth_atom_id": "N", "pdbx_PDB_model_num": "9",
        })
        self.assertEqual(d["atom_site"][16443], {
            "group_PDB": "ATOM", "id": "16444", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "ALA",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "199",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "37.677", "Cartn_y": "-7.651",
            "Cartn_z": "-0.218", "occupancy": "1.0", "B_iso_or_equiv": "0.0",
            "pdbx_formal_charge": "0", "auth_seq_id": "199", "auth_comp_id": "ALA",
            "auth_asym_id": "A", "auth_atom_id": "N", "pdbx_PDB_model_num": "10",
        })
        self.assertEqual(d["atom_site"][-1], {
            "group_PDB": "ATOM", "id": "18270", "type_symbol": "H",
            "label_atom_id": "HB3", "label_alt_id": ".", "label_comp_id": "ALA",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "312",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "22.641", "Cartn_y": "-15.379",
            "Cartn_z": "-10.884", "occupancy": "1.0", "B_iso_or_equiv": "0.0",
            "pdbx_formal_charge": "0", "auth_seq_id": "312", "auth_comp_id": "ALA",
            "auth_asym_id": "A", "auth_atom_id": "HB3", "pdbx_PDB_model_num": "10",
        })
    

    def test_1cbn(self):
        d = atomium.open("tests/integration/files/1cbn.mmtf", dictionary=True)
        self.assertEqual(len(d.keys()), 18)

        # Header categories
        self.assertEqual(d["entry"], [{"id": "1CBN"}])
        self.assertEqual(d["struct"], [{
            "entry_id": "1CBN",
            "title": "ATOMIC RESOLUTION (0.83 ANGSTROMS) CRYSTAL STRUCTURE OF THE HYDROPHOBIC PROTEIN CRAMBIN AT 130 K",
        }])
        self.assertEqual(d["pdbx_database_status"], [{"entry_id": "1CBN", "recvd_initial_deposition_date": "1991-10-11"}])
        self.assertEqual(d["exptl"], [{"entry_id": "1CBN", "method": "X-RAY DIFFRACTION"}])

        # Quality categories
        self.assertEqual(d["refine"], [{
            "entry_id": "1CBN", "ls_d_res_high": "0.83", "ls_R_factor_obs": "?",
            "ls_R_factor_all": "?", "ls_R_factor_R_work": "?",
            "ls_R_factor_R_free": "?",
        }])

        # Crystal categories
        self.assertEqual(d["symmetry"], [{"entry_id": "1CBN", "space_group_name_H-M": "P 1 21 1"}])
        self.assertEqual(d["cell"], [{
            "entry_id": "1CBN", "length_a": "40.763", "length_b": "18.492",
            "length_c": "22.333", "length_alpha": "90.0", "length_beta": "90.61",
            "length_gamma": "90.0",
        }])

        # Assembly categories
        self.assertEqual(d["pdbx_struct_assembly"], [{"id": "1"}])
        self.assertEqual(d["pdbx_struct_assembly_gen"], [{"assembly_id": "1", "oper_expression": "1", "asym_id_list": "A,B"}])
        self.assertEqual(d["pdbx_struct_oper_list"], [{
            "id": "1", "matrix[1][1]": "1.0000000000",
            "matrix[1][2]": "0.0000000000", "matrix[1][3]": "0.0000000000",
            "vector[1]": "0.0000000000", "matrix[2][1]": "0.0000000000",
            "matrix[2][2]": "1.0000000000", "matrix[2][3]": "0.0000000000",
            "vector[2]": "0.0000000000", "matrix[3][1]": "0.0000000000",
            "matrix[3][2]": "0.0000000000", "matrix[3][3]": "1.0000000000",
            "vector[3]": "0.0000000000",
        }])

        # Entity categories
        self.assertEqual(d["entity"], [
            {"id": "1", "type": "polymer", "pdbx_description": "CRAMBIN", "pdbx_number_of_molecules": "1"},
            {"id": "2", "type": "non-polymer", "pdbx_description": "ETHANOL", "pdbx_number_of_molecules": "1"},
        ])
        self.assertEqual(d["entity_poly"], [{
            "entity_id": "1",
            "pdbx_seq_one_letter_code": "TTCCPSIVARSNFNVCRLPGTSEAICATYTGCIIIPGATCPGDYAN",
            "pdbx_seq_one_letter_code_can": "TTCCPSIVARSNFNVCRLPGTSEAICATYTGCIIIPGATCPGDYAN",
        }])
        self.assertEqual(d["struct_asym"], [
            {"id": "A", "pdbx_blank_PDB_chainid_flag": "N", "pdbx_modified": "N", "entity_id": "1", "details": "?"},
            {"id": "B", "pdbx_blank_PDB_chainid_flag": "N", "pdbx_modified": "N", "entity_id": "2", "details": "?"},
        ])

        # Compound categories
        self.assertEqual(len(d["chem_comp"]), 16)
        self.assertEqual(d["chem_comp"][0], {
            "id": "ALA", "type": "L-PEPTIDE LINKING", "mon_nstd_flag": "y",
            "name": "alanine", "pdbx_synonyms": "?", "formula": "?",
            "formula_weight": "?",
        })
        self.assertEqual(d["chem_comp"][5], {
            "id": "EOH", "type": "NON-POLYMER", "mon_nstd_flag": ".", "name": "?",
            "pdbx_synonyms": "?", "formula": "?", "formula_weight": "?",
        })
        self.assertEqual(d["chem_comp"][-1], {
            "id": "VAL", "type": "L-PEPTIDE LINKING", "mon_nstd_flag": "y",
            "name": "valine", "pdbx_synonyms": "?", "formula": "?",
            "formula_weight": "?",
        })

        # Model categories
        self.assertEqual(d["atom_type"], [
            {"symbol": "C"},
            {"symbol": "H"},
            {"symbol": "N"},
            {"symbol": "O"},
            {"symbol": "S"},
        ])
        self.assertEqual(len(d["atom_site"]), 777)
        self.assertEqual(d["atom_site"][0], {
            "group_PDB": "ATOM", "id": "1", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": "A", "label_comp_id": "THR",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "16.864", "Cartn_y": "14.059",
            "Cartn_z": "3.442", "occupancy": "0.8", "B_iso_or_equiv": "6.22",
            "pdbx_formal_charge": "0", "auth_seq_id": "1", "auth_comp_id": "THR",
            "auth_asym_id": "A", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][772], {
            "group_PDB": "HETATM", "id": "773", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "EOH",
            "label_asym_id": "B", "label_entity_id": "2", "label_seq_id": "66",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "15.702", "Cartn_y": "0.904",
            "Cartn_z": "12.771", "occupancy": "1.0", "B_iso_or_equiv": "17.9",
            "pdbx_formal_charge": "0", "auth_seq_id": "66", "auth_comp_id": "EOH",
            "auth_asym_id": "A", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][-1], {
            "group_PDB": "HETATM", "id": "777", "type_symbol": "O",
            "label_atom_id": "O", "label_alt_id": "B", "label_comp_id": "EOH",
            "label_asym_id": "B", "label_entity_id": "2", "label_seq_id": "66",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "14.54", "Cartn_y": "1.708",
            "Cartn_z": "12.931", "occupancy": "0.3", "B_iso_or_equiv": "6.09",
            "pdbx_formal_charge": "0", "auth_seq_id": "66", "auth_comp_id": "EOH",
            "auth_asym_id": "A", "auth_atom_id": "O", "pdbx_PDB_model_num": "1",
        })


    def test_1m4x(self):
        d = atomium.open("tests/integration/files/1m4x.mmtf", dictionary=True)
        self.assertEqual(len(d.keys()), 17)

        # Header categories
        self.assertEqual(d["entry"], [{"id": "1M4X"}])
        self.assertEqual(d["struct"], [{"entry_id": "1M4X", "title": "PBCV-1 virus capsid, quasi-atomic model"}])
        self.assertEqual(d["pdbx_database_status"], [{"entry_id": "1M4X", "recvd_initial_deposition_date": "2002-07-05"}])
        self.assertEqual(d["exptl"], [{"entry_id": "1M4X", "method": "ELECTRON MICROSCOPY"}])

        # Quality categories
        self.assertNotIn("refine", d)

        # Crystal categories
        self.assertEqual(d["symmetry"], [{"entry_id": "1M4X", "space_group_name_H-M": "P 2 3"}])
        self.assertEqual(d["cell"], [{
            "entry_id": "1M4X", "length_a": "1927.0", "length_b": "1927.0",
            "length_c": "1927.0", "length_alpha": "90.0", "length_beta": "90.0",
            "length_gamma": "90.0",
        }])

        # Assembly categories
        self.assertEqual(d["pdbx_struct_assembly"], [
            {"id": "1"},
            {"id": "2"},
            {"id": "3"},
            {"id": "4"},
            {"id": "5"},
            {"id": "6"},
            {"id": "7"},
        ])
        self.assertEqual(len(d["pdbx_struct_assembly_gen"]), 2140)
        self.assertEqual(d["pdbx_struct_assembly_gen"][0], {
            "assembly_id": "1", "oper_expression": "1", "asym_id_list": "A,B,C",
        })
        self.assertEqual(d["pdbx_struct_assembly_gen"][-1], {
            "assembly_id": "7", "oper_expression": "2140", "asym_id_list": "A,B,C",
        })
        self.assertEqual(len(d["pdbx_struct_oper_list"]), 2140)
        self.assertEqual(d["pdbx_struct_oper_list"][0], {
            "id": "1", "matrix[1][1]": "1.0000000000",
            "matrix[1][2]": "0.0000000000", "matrix[1][3]": "0.0000000000",
            "vector[1]": "0.0000000000", "matrix[2][1]": "0.0000000000",
            "matrix[2][2]": "1.0000000000", "matrix[2][3]": "0.0000000000",
            "vector[2]": "0.0000000000", "matrix[3][1]": "0.0000000000",
            "matrix[3][2]": "0.0000000000", "matrix[3][3]": "1.0000000000",
            "vector[3]": "0.0000000000",
        })
        self.assertEqual(d["pdbx_struct_oper_list"][-1], {
            "id": "2140", "matrix[1][1]": "0.5313321011",
            "matrix[1][2]": "-0.3081559053", "matrix[1][3]": "0.7891296276",
            "vector[1]": "-89.3238998160", "matrix[2][1]": "0.7551163827",
            "matrix[2][2]": "-0.2499963986", "matrix[2][3]": "-0.6060533418",
            "vector[2]": "-134.0685044124", "matrix[3][1]": "0.3840387183",
            "matrix[3][2]": "0.9179014933", "matrix[3][3]": "0.0998619693",
            "vector[3]": "11.9906045963",
        })

        # Entity categories
        self.assertEqual(d["entity"], [
            {"id": "1", "type": "polymer", "pdbx_description": "PBCV-1 virus capsid", "pdbx_number_of_molecules": "3"},
        ])
        self.assertEqual(d["entity_poly"], [{
            "entity_id": "1",
            "pdbx_seq_one_letter_code": "TFFKTVYRRYTNFAIESIQQTINGSVGFGNKVSTQISRNGDLITDIVVEFVLTKGGNGGTTYYPAEELLQDVELEIGGQRIDKHYNDWFRTYDALFRMNDDRYNYRRMTDWVNNELVGAQKRFYVPLIFFFNQTPGLALPLIALQYHEVKLYFTLASQVQGVNYNGSSAIAGAAQPTMSVWVDYIFLDTQERTRFAQLPHEYLIEQLQFTGSETATPSATTQAAQNIRLNFNHPTKYLAWNFNNPTNYGQYTALANIPGACSGAGTAAATVTTPDYGNTGTYNEQLAVLDSAKIQLNGQDRFATRKGSYFNKVQPYQSIGGVTPAGVYLYSFALKPAGRQPSGTCNFSRIDNATLSLTYKTCSIDATSPAAVLGNTETVTANTATLLTALNIYAKNYNVLRIMSGMGGLAYAN",
            "pdbx_seq_one_letter_code_can": "TFFKTVYRRYTNFAIESIQQTINGSVGFGNKVSTQISRNGDLITDIVVEFVLTKGGNGGTTYYPAEELLQDVELEIGGQRIDKHYNDWFRTYDALFRMNDDRYNYRRMTDWVNNELVGAQKRFYVPLIFFFNQTPGLALPLIALQYHEVKLYFTLASQVQGVNYNGSSAIAGAAQPTMSVWVDYIFLDTQERTRFAQLPHEYLIEQLQFTGSETATPSATTQAAQNIRLNFNHPTKYLAWNFNNPTNYGQYTALANIPGACSGAGTAAATVTTPDYGNTGTYNEQLAVLDSAKIQLNGQDRFATRKGSYFNKVQPYQSIGGVTPAGVYLYSFALKPAGRQPSGTCNFSRIDNATLSLTYKTCSIDATSPAAVLGNTETVTANTATLLTALNIYAKNYNVLRIMSGMGGLAYAN",
        }])
        self.assertEqual(d["struct_asym"], [
            {"id": "A", "pdbx_blank_PDB_chainid_flag": "N", "pdbx_modified": "N", "entity_id": "1", "details": "?"},
            {"id": "B", "pdbx_blank_PDB_chainid_flag": "N", "pdbx_modified": "N", "entity_id": "1", "details": "?"},
            {"id": "C", "pdbx_blank_PDB_chainid_flag": "N", "pdbx_modified": "N", "entity_id": "1", "details": "?"},
        ])

        # Compound categories
        self.assertEqual(len(d["chem_comp"]), 20)
        self.assertEqual(d["chem_comp"][0], {
            "id": "ALA", "type": "L-PEPTIDE LINKING", "mon_nstd_flag": "y",
            "name": "alanine", "pdbx_synonyms": "?", "formula": "?",
            "formula_weight": "?",
        })
        self.assertEqual(d["chem_comp"][-1], {
            "id": "VAL", "type": "L-PEPTIDE LINKING", "mon_nstd_flag": "y",
            "name": "valine", "pdbx_synonyms": "?", "formula": "?",
            "formula_weight": "?",
        })

        # Model categories
        self.assertEqual(d["atom_type"], [{"symbol": "C"}, {"symbol": "N"}, {"symbol": "O"}, {"symbol": "S"}])
        self.assertEqual(len(d["atom_site"]), 9693)
        self.assertEqual(d["atom_site"][0], {
            "group_PDB": "ATOM", "id": "1", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "THR",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "25",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "529.72", "Cartn_y": "577.569",
            "Cartn_z": "66.498", "occupancy": "1.0", "B_iso_or_equiv": "32.89",
            "pdbx_formal_charge": "0", "auth_seq_id": "25", "auth_comp_id": "THR",
            "auth_asym_id": "A", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][3231], {
            "group_PDB": "ATOM", "id": "3232", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "THR",
            "label_asym_id": "B", "label_entity_id": "1", "label_seq_id": "25",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "528.258", "Cartn_y": "582.202",
            "Cartn_z": "44.416", "occupancy": "1.0", "B_iso_or_equiv": "32.89",
            "pdbx_formal_charge": "0", "auth_seq_id": "25", "auth_comp_id": "THR",
            "auth_asym_id": "B", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][6462], {
            "group_PDB": "ATOM", "id": "6463", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "THR",
            "label_asym_id": "C", "label_entity_id": "1", "label_seq_id": "25",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "546.577", "Cartn_y": "571.77",
            "Cartn_z": "52.59", "occupancy": "1.0", "B_iso_or_equiv": "32.89",
            "pdbx_formal_charge": "0", "auth_seq_id": "25", "auth_comp_id": "THR",
            "auth_asym_id": "C", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][-1], {
            "group_PDB": "ATOM", "id": "9693", "type_symbol": "O",
            "label_atom_id": "OXT", "label_alt_id": ".", "label_comp_id": "ASN",
            "label_asym_id": "C", "label_entity_id": "1", "label_seq_id": "437",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "511.742", "Cartn_y": "599.344",
            "Cartn_z": "66.077", "occupancy": "1.0", "B_iso_or_equiv": "52.28",
            "pdbx_formal_charge": "0", "auth_seq_id": "437", "auth_comp_id": "ASN",
            "auth_asym_id": "C", "auth_atom_id": "OXT", "pdbx_PDB_model_num": "1",
        })


    def test_6xlu(self):
        d = atomium.open("tests/integration/files/6xlu.mmtf", dictionary=True)
        self.assertEqual(len(d.keys()), 16)

        # Header categories
        self.assertEqual(d["entry"], [{"id": "6XLU"}])
        self.assertEqual(d["struct"], [{"entry_id": "6XLU", "title": "Structure of SARS-CoV-2 spike at pH 4.0"}])
        self.assertEqual(d["pdbx_database_status"], [{"entry_id": "6XLU", "recvd_initial_deposition_date": "2020-06-29"}])
        self.assertEqual(d["exptl"], [{"entry_id": "6XLU", "method": "ELECTRON MICROSCOPY"}])

        # Quality categories
        self.assertNotIn("refine", d)

        # Crystal categories
        self.assertEqual(d["symmetry"], [{"entry_id": "6XLU", "space_group_name_H-M": "P 1"}])
        self.assertNotIn("cell", d)

        # Assembly categories
        self.assertEqual(d["pdbx_struct_assembly"], [{"id": "1"}])
        self.assertEqual(d["pdbx_struct_assembly_gen"], [{
            "assembly_id": "1", "oper_expression": "1",
            "asym_id_list": "A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,AA,BA,CA,DA,EA,FA,GA,HA,IA,JA,KA,LA,MA,NA,OA,PA,QA,RA,SA,TA,UA,VA,WA,XA,YA,ZA,AB",
        }])
        self.assertEqual(d["pdbx_struct_oper_list"], [{
            "id": "1", "matrix[1][1]": "1.0000000000",
            "matrix[1][2]": "0.0000000000", "matrix[1][3]": "0.0000000000",
            "vector[1]": "0.0000000000", "matrix[2][1]": "0.0000000000",
            "matrix[2][2]": "1.0000000000", "matrix[2][3]": "0.0000000000",
            "vector[2]": "0.0000000000", "matrix[3][1]": "0.0000000000",
            "matrix[3][2]": "0.0000000000", "matrix[3][3]": "1.0000000000",
            "vector[3]": "0.0000000000",
        }])

        # Entity categories
        self.assertEqual(d["entity"], [{
            "id": "1", "type": "polymer", "pdbx_description": "Spike glycoprotein",
            "pdbx_number_of_molecules": "3",
        }, {
            "id": "2", "type": "branched",
            "pdbx_description": "2-acetamido-2-deoxy-beta-D-glucopyranose-(1-4)-2-acetamido-2-deoxy-beta-D-glucopyranose",
            "pdbx_number_of_molecules": "15",
        }, {
            "id": "3", "type": "non-polymer",
            "pdbx_description": "2-acetamido-2-deoxy-beta-D-glucopyranose",
            "pdbx_number_of_molecules": "32",
        }, {
            "id": "4", "type": "water", "pdbx_description": "water",
            "pdbx_number_of_molecules": "3",
        }])
        self.assertEqual(d["entity_poly"], [{
            "entity_id": "1",
            "pdbx_seq_one_letter_code": "QCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPGSASSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDPPEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQGSGYIPEAPRDGQAYVRKDGEWVLLSTFLGRSLEVLFQGPGHHHHHHHHSAWSHPQFEKGGGSGGGGSGGSAWSHPQFEK",
            "pdbx_seq_one_letter_code_can": "QCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPGSASSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDPPEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQGSGYIPEAPRDGQAYVRKDGEWVLLSTFLGRSLEVLFQGPGHHHHHHHHSAWSHPQFEKGGGSGGGGSGGSAWSHPQFEK",
        }])
        self.assertEqual(len(d["struct_asym"]), 53)
        self.assertEqual(d["struct_asym"][0], {
            "id": "A", "pdbx_blank_PDB_chainid_flag": "N", "pdbx_modified": "N",
            "entity_id": "1", "details": "?",
        })
        self.assertEqual(d["struct_asym"][-1], {
            "id": "AB", "pdbx_blank_PDB_chainid_flag": "N", "pdbx_modified": "N",
            "entity_id": "4", "details": "?",
        })

        # Compound categories
        self.assertEqual(len(d["chem_comp"]), 22)
        self.assertEqual(d["chem_comp"][0], {
            "id": "ALA", "type": "L-PEPTIDE LINKING", "mon_nstd_flag": "y",
            "name": "alanine", "pdbx_synonyms": "?", "formula": "?",
            "formula_weight": "?",
        })
        self.assertEqual(d["chem_comp"][9], {
            "id": "HOH", "type": "NON-POLYMER", "mon_nstd_flag": ".", "name": "?",
            "pdbx_synonyms": "?", "formula": "?", "formula_weight": "?",
        })
        self.assertEqual(d["chem_comp"][14], {
            "id": "NAG", "type": "D-SACCHARIDE", "mon_nstd_flag": ".", "name": "?",
            "pdbx_synonyms": "?", "formula": "?", "formula_weight": "?",
        })
        self.assertEqual(d["chem_comp"][-1], {
            "id": "VAL", "type": "L-PEPTIDE LINKING", "mon_nstd_flag": "y",
            "name": "valine", "pdbx_synonyms": "?", "formula": "?",
            "formula_weight": "?",
        })

        # Model categories
        self.assertEqual(d["atom_type"], [{"symbol": "C"}, {"symbol": "N"}, {"symbol": "O"}, {"symbol": "S"}])
        self.assertEqual(len(d["atom_site"]), 25990)
        self.assertEqual(d["atom_site"][0], {
            "group_PDB": "ATOM", "id": "1", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "ALA",
            "label_asym_id": "A", "label_entity_id": "1", "label_seq_id": "27",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "137.193", "Cartn_y": "160.456",
            "Cartn_z": "208.534", "occupancy": "1.0", "B_iso_or_equiv": "62.75",
            "pdbx_formal_charge": "0", "auth_seq_id": "27", "auth_comp_id": "ALA",
            "auth_asym_id": "A", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][8276], {
            "group_PDB": "ATOM", "id": "8277", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "ALA",
            "label_asym_id": "B", "label_entity_id": "1", "label_seq_id": "27",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "234.997", "Cartn_y": "157.197",
            "Cartn_z": "208.5", "occupancy": "1.0", "B_iso_or_equiv": "60.77",
            "pdbx_formal_charge": "0", "auth_seq_id": "27", "auth_comp_id": "ALA",
            "auth_asym_id": "B", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][16578], {
            "group_PDB": "ATOM", "id": "16579", "type_symbol": "N",
            "label_atom_id": "N", "label_alt_id": ".", "label_comp_id": "ALA",
            "label_asym_id": "C", "label_entity_id": "1", "label_seq_id": "27",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "188.766", "Cartn_y": "243.345",
            "Cartn_z": "208.464", "occupancy": "1.0", "B_iso_or_equiv": "60.93",
            "pdbx_formal_charge": "0", "auth_seq_id": "27", "auth_comp_id": "ALA",
            "auth_asym_id": "C", "auth_atom_id": "N", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][24872], {
            "group_PDB": "HETATM", "id": "24873", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "D", "label_entity_id": "2", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "158.292", "Cartn_y": "167.896",
            "Cartn_z": "226.28", "occupancy": "1.0", "B_iso_or_equiv": "68.52",
            "pdbx_formal_charge": "0", "auth_seq_id": "1", "auth_comp_id": "NAG",
            "auth_asym_id": "D", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][24900], {
            "group_PDB": "HETATM", "id": "24901", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "E", "label_entity_id": "2", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "159.383", "Cartn_y": "181.948",
            "Cartn_z": "133.968", "occupancy": "1.0", "B_iso_or_equiv": "39.16",
            "pdbx_formal_charge": "0", "auth_seq_id": "1", "auth_comp_id": "NAG",
            "auth_asym_id": "E", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][24928], {
            "group_PDB": "HETATM", "id": "24929", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "F", "label_entity_id": "2", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "155.569", "Cartn_y": "201.666",
            "Cartn_z": "144.621", "occupancy": "1.0", "B_iso_or_equiv": "37.3",
            "pdbx_formal_charge": "0", "auth_seq_id": "1", "auth_comp_id": "NAG",
            "auth_asym_id": "F", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][24956], {
            "group_PDB": "HETATM", "id": "24957", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "G", "label_entity_id": "2", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "166.347", "Cartn_y": "170.261",
            "Cartn_z": "117.882", "occupancy": "1.0", "B_iso_or_equiv": "55.39",
            "pdbx_formal_charge": "0", "auth_seq_id": "1", "auth_comp_id": "NAG",
            "auth_asym_id": "G", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][24984], {
            "group_PDB": "HETATM", "id": "24985", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "H", "label_entity_id": "2", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "183.938", "Cartn_y": "159.964",
            "Cartn_z": "113.723", "occupancy": "1.0", "B_iso_or_equiv": "67.06",
            "pdbx_formal_charge": "0", "auth_seq_id": "1", "auth_comp_id": "NAG",
            "auth_asym_id": "H", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25012], {
            "group_PDB": "HETATM", "id": "25013", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "I", "label_entity_id": "2", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "217.803", "Cartn_y": "171.699",
            "Cartn_z": "226.262", "occupancy": "1.0", "B_iso_or_equiv": "68.57",
            "pdbx_formal_charge": "0", "auth_seq_id": "1", "auth_comp_id": "NAG",
            "auth_asym_id": "I", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25040], {
            "group_PDB": "HETATM", "id": "25041", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "J", "label_entity_id": "2", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "205.173", "Cartn_y": "165.618",
            "Cartn_z": "133.979", "occupancy": "1.0", "B_iso_or_equiv": "40.85",
            "pdbx_formal_charge": "0", "auth_seq_id": "1", "auth_comp_id": "NAG",
            "auth_asym_id": "J", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25068], {
            "group_PDB": "HETATM", "id": "25069", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "K", "label_entity_id": "2", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "189.993", "Cartn_y": "152.448",
            "Cartn_z": "144.612", "occupancy": "1.0", "B_iso_or_equiv": "37.36",
            "pdbx_formal_charge": "0", "auth_seq_id": "1", "auth_comp_id": "NAG",
            "auth_asym_id": "K", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25096], {
            "group_PDB": "HETATM", "id": "25097", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "L", "label_entity_id": "2", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "211.852", "Cartn_y": "177.523",
            "Cartn_z": "117.872", "occupancy": "1.0", "B_iso_or_equiv": "52.88",
            "pdbx_formal_charge": "0", "auth_seq_id": "1", "auth_comp_id": "NAG",
            "auth_asym_id": "L", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25124], {
            "group_PDB": "HETATM", "id": "25125", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "M", "label_entity_id": "2", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "211.963", "Cartn_y": "197.857",
            "Cartn_z": "113.701", "occupancy": "1.0", "B_iso_or_equiv": "69.32",
            "pdbx_formal_charge": "0", "auth_seq_id": "1", "auth_comp_id": "NAG",
            "auth_asym_id": "M", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25152], {
            "group_PDB": "HETATM", "id": "25153", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "N", "label_entity_id": "2", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "184.802", "Cartn_y": "221.328",
            "Cartn_z": "226.28", "occupancy": "1.0", "B_iso_or_equiv": "67.25",
            "pdbx_formal_charge": "0", "auth_seq_id": "1", "auth_comp_id": "NAG",
            "auth_asym_id": "N", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25180], {
            "group_PDB": "HETATM", "id": "25181", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "O", "label_entity_id": "2", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "196.431", "Cartn_y": "213.401",
            "Cartn_z": "133.99", "occupancy": "1.0", "B_iso_or_equiv": "39.07",
            "pdbx_formal_charge": "0", "auth_seq_id": "1", "auth_comp_id": "NAG",
            "auth_asym_id": "O", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25208], {
            "group_PDB": "HETATM", "id": "25209", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "P", "label_entity_id": "2", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "215.378", "Cartn_y": "206.864",
            "Cartn_z": "144.659", "occupancy": "1.0", "B_iso_or_equiv": "36.31",
            "pdbx_formal_charge": "0", "auth_seq_id": "1", "auth_comp_id": "NAG",
            "auth_asym_id": "P", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25236], {
            "group_PDB": "HETATM", "id": "25237", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "Q", "label_entity_id": "2", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "182.774", "Cartn_y": "213.299",
            "Cartn_z": "117.869", "occupancy": "1.0", "B_iso_or_equiv": "57.43",
            "pdbx_formal_charge": "0", "auth_seq_id": "1", "auth_comp_id": "NAG",
            "auth_asym_id": "Q", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25264], {
            "group_PDB": "HETATM", "id": "25265", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "R", "label_entity_id": "2", "label_seq_id": "1",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "165.065", "Cartn_y": "203.25",
            "Cartn_z": "113.813", "occupancy": "1.0", "B_iso_or_equiv": "62.93",
            "pdbx_formal_charge": "0", "auth_seq_id": "1", "auth_comp_id": "NAG",
            "auth_asym_id": "R", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25292], {
            "group_PDB": "HETATM", "id": "25293", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "S", "label_entity_id": "3", "label_seq_id": "1301",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "145.422", "Cartn_y": "161.696",
            "Cartn_z": "201.354", "occupancy": "1.0", "B_iso_or_equiv": "81.13",
            "pdbx_formal_charge": "0", "auth_seq_id": "1301", "auth_comp_id": "NAG",
            "auth_asym_id": "A", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25306], {
            "group_PDB": "HETATM", "id": "25307", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "T", "label_entity_id": "3", "label_seq_id": "1302",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "134.38", "Cartn_y": "185.068",
            "Cartn_z": "231.558", "occupancy": "1.0", "B_iso_or_equiv": "76.02",
            "pdbx_formal_charge": "0", "auth_seq_id": "1302", "auth_comp_id": "NAG",
            "auth_asym_id": "A", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25320], {
            "group_PDB": "HETATM", "id": "25321", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "U", "label_entity_id": "3", "label_seq_id": "1303",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "143.911", "Cartn_y": "198.545",
            "Cartn_z": "195.564", "occupancy": "1.0", "B_iso_or_equiv": "75.72",
            "pdbx_formal_charge": "0", "auth_seq_id": "1303", "auth_comp_id": "NAG",
            "auth_asym_id": "A", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25334], {
            "group_PDB": "HETATM", "id": "25335", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "V", "label_entity_id": "3", "label_seq_id": "1304",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "183.3", "Cartn_y": "142.267",
            "Cartn_z": "220.784", "occupancy": "1.0", "B_iso_or_equiv": "67.59",
            "pdbx_formal_charge": "0", "auth_seq_id": "1304", "auth_comp_id": "NAG",
            "auth_asym_id": "A", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25348], {
            "group_PDB": "HETATM", "id": "25349", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "W", "label_entity_id": "3", "label_seq_id": "1305",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "184.641", "Cartn_y": "162.816",
            "Cartn_z": "244.044", "occupancy": "1.0", "B_iso_or_equiv": "69.85",
            "pdbx_formal_charge": "0", "auth_seq_id": "1305", "auth_comp_id": "NAG",
            "auth_asym_id": "A", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25362], {
            "group_PDB": "HETATM", "id": "25363", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "X", "label_entity_id": "3", "label_seq_id": "1306",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "150.119", "Cartn_y": "183.863",
            "Cartn_z": "177.16", "occupancy": "1.0", "B_iso_or_equiv": "69.9",
            "pdbx_formal_charge": "0", "auth_seq_id": "1306", "auth_comp_id": "NAG",
            "auth_asym_id": "A", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25376], {
            "group_PDB": "HETATM", "id": "25377", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "Y", "label_entity_id": "3", "label_seq_id": "1307",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "169.682", "Cartn_y": "149.237",
            "Cartn_z": "182.593", "occupancy": "1.0", "B_iso_or_equiv": "67.74",
            "pdbx_formal_charge": "0", "auth_seq_id": "1307", "auth_comp_id": "NAG",
            "auth_asym_id": "A", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25390], {
            "group_PDB": "HETATM", "id": "25391", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "Z", "label_entity_id": "3", "label_seq_id": "1308",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "158.1", "Cartn_y": "151.419",
            "Cartn_z": "166.545", "occupancy": "1.0", "B_iso_or_equiv": "59.37",
            "pdbx_formal_charge": "0", "auth_seq_id": "1308", "auth_comp_id": "NAG",
            "auth_asym_id": "A", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25404], {
            "group_PDB": "HETATM", "id": "25405", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "AA", "label_entity_id": "3", "label_seq_id": "1309",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "180.077", "Cartn_y": "156.831",
            "Cartn_z": "129.475", "occupancy": "1.0", "B_iso_or_equiv": "44.61",
            "pdbx_formal_charge": "0", "auth_seq_id": "1309", "auth_comp_id": "NAG",
            "auth_asym_id": "A", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25418], {
            "group_PDB": "HETATM", "id": "25419", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "BA", "label_entity_id": "3", "label_seq_id": "1310",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "166.689", "Cartn_y": "163.753",
            "Cartn_z": "134.108", "occupancy": "1.0", "B_iso_or_equiv": "33.59",
            "pdbx_formal_charge": "0", "auth_seq_id": "1310", "auth_comp_id": "NAG",
            "auth_asym_id": "A", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25432], {
            "group_PDB": "HETATM", "id": "25433", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "CA", "label_entity_id": "3", "label_seq_id": "1311",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "151.958", "Cartn_y": "175.259",
            "Cartn_z": "243.466", "occupancy": "1.0", "B_iso_or_equiv": "74.39",
            "pdbx_formal_charge": "0", "auth_seq_id": "1311", "auth_comp_id": "NAG",
            "auth_asym_id": "A", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25446], {
            "group_PDB": "HETATM", "id": "25447", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "DA", "label_entity_id": "3", "label_seq_id": "1301",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "230.267", "Cartn_y": "162.588",
            "Cartn_z": "201.139", "occupancy": "1.0", "B_iso_or_equiv": "70.59",
            "pdbx_formal_charge": "0", "auth_seq_id": "1301", "auth_comp_id": "NAG",
            "auth_asym_id": "B", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25460], {
            "group_PDB": "HETATM", "id": "25461", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "EA", "label_entity_id": "3", "label_seq_id": "1302",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "214.299", "Cartn_y": "141.381",
            "Cartn_z": "231.095", "occupancy": "1.0", "B_iso_or_equiv": "79.22",
            "pdbx_formal_charge": "0", "auth_seq_id": "1302", "auth_comp_id": "NAG",
            "auth_asym_id": "B", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25474], {
            "group_PDB": "HETATM", "id": "25475", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "FA", "label_entity_id": "3", "label_seq_id": "1303",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "198.587", "Cartn_y": "143.875",
            "Cartn_z": "195.594", "occupancy": "1.0", "B_iso_or_equiv": "74.72",
            "pdbx_formal_charge": "0", "auth_seq_id": "1303", "auth_comp_id": "NAG",
            "auth_asym_id": "B", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25488], {
            "group_PDB": "HETATM", "id": "25489", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "GA", "label_entity_id": "3", "label_seq_id": "1304",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "227.744", "Cartn_y": "206.039",
            "Cartn_z": "220.892", "occupancy": "1.0", "B_iso_or_equiv": "69.75",
            "pdbx_formal_charge": "0", "auth_seq_id": "1304", "auth_comp_id": "NAG",
            "auth_asym_id": "B", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25502], {
            "group_PDB": "HETATM", "id": "25503", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "HA", "label_entity_id": "3", "label_seq_id": "1305",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "209.034", "Cartn_y": "197.05",
            "Cartn_z": "243.993", "occupancy": "1.0", "B_iso_or_equiv": "67.23",
            "pdbx_formal_charge": "0", "auth_seq_id": "1305", "auth_comp_id": "NAG",
            "auth_asym_id": "B", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25516], {
            "group_PDB": "HETATM", "id": "25517", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "IA", "label_entity_id": "3", "label_seq_id": "1306",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "208.711", "Cartn_y": "156.753",
            "Cartn_z": "177.03", "occupancy": "1.0", "B_iso_or_equiv": "73.18",
            "pdbx_formal_charge": "0", "auth_seq_id": "1306", "auth_comp_id": "NAG",
            "auth_asym_id": "B", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25530], {
            "group_PDB": "HETATM", "id": "25531", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "JA", "label_entity_id": "3", "label_seq_id": "1307",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "228.437", "Cartn_y": "190.86",
            "Cartn_z": "182.552", "occupancy": "1.0", "B_iso_or_equiv": "70.05",
            "pdbx_formal_charge": "0", "auth_seq_id": "1307", "auth_comp_id": "NAG",
            "auth_asym_id": "B", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25544], {
            "group_PDB": "HETATM", "id": "25545", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "KA", "label_entity_id": "3", "label_seq_id": "1308",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "232.273", "Cartn_y": "179.738",
            "Cartn_z": "166.555", "occupancy": "1.0", "B_iso_or_equiv": "63.51",
            "pdbx_formal_charge": "0", "auth_seq_id": "1308", "auth_comp_id": "NAG",
            "auth_asym_id": "B", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25558], {
            "group_PDB": "HETATM", "id": "25559", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "LA", "label_entity_id": "3", "label_seq_id": "1309",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "216.571", "Cartn_y": "196.044",
            "Cartn_z": "129.428", "occupancy": "1.0", "B_iso_or_equiv": "42.57",
            "pdbx_formal_charge": "0", "auth_seq_id": "1309", "auth_comp_id": "NAG",
            "auth_asym_id": "B", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25572], {
            "group_PDB": "HETATM", "id": "25573", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "MA", "label_entity_id": "3", "label_seq_id": "1310",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "217.329", "Cartn_y": "181.072",
            "Cartn_z": "134.101", "occupancy": "1.0", "B_iso_or_equiv": "32.92",
            "pdbx_formal_charge": "0", "auth_seq_id": "1310", "auth_comp_id": "NAG",
            "auth_asym_id": "B", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25586], {
            "group_PDB": "HETATM", "id": "25587", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "NA", "label_entity_id": "3", "label_seq_id": "1301",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "186.159", "Cartn_y": "236.14",
            "Cartn_z": "202.719", "occupancy": "1.0", "B_iso_or_equiv": "72.3",
            "pdbx_formal_charge": "0", "auth_seq_id": "1301", "auth_comp_id": "NAG",
            "auth_asym_id": "C", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25600], {
            "group_PDB": "HETATM", "id": "25601", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "OA", "label_entity_id": "3", "label_seq_id": "1302",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "211.695", "Cartn_y": "234.217",
            "Cartn_z": "231.695", "occupancy": "1.0", "B_iso_or_equiv": "77.68",
            "pdbx_formal_charge": "0", "auth_seq_id": "1302", "auth_comp_id": "NAG",
            "auth_asym_id": "C", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25614], {
            "group_PDB": "HETATM", "id": "25615", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "PA", "label_entity_id": "3", "label_seq_id": "1303",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "218.505", "Cartn_y": "218.599",
            "Cartn_z": "195.565", "occupancy": "1.0", "B_iso_or_equiv": "74.3",
            "pdbx_formal_charge": "0", "auth_seq_id": "1303", "auth_comp_id": "NAG",
            "auth_asym_id": "C", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25628], {
            "group_PDB": "HETATM", "id": "25629", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "QA", "label_entity_id": "3", "label_seq_id": "1304",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "150.033", "Cartn_y": "212.547",
            "Cartn_z": "220.786", "occupancy": "1.0", "B_iso_or_equiv": "67.57",
            "pdbx_formal_charge": "0", "auth_seq_id": "1304", "auth_comp_id": "NAG",
            "auth_asym_id": "C", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25642], {
            "group_PDB": "HETATM", "id": "25643", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "RA", "label_entity_id": "3", "label_seq_id": "1305",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "167.22", "Cartn_y": "201.11",
            "Cartn_z": "244.004", "occupancy": "1.0", "B_iso_or_equiv": "65.33",
            "pdbx_formal_charge": "0", "auth_seq_id": "1305", "auth_comp_id": "NAG",
            "auth_asym_id": "C", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25656], {
            "group_PDB": "HETATM", "id": "25657", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "SA", "label_entity_id": "3", "label_seq_id": "1306",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "202.543", "Cartn_y": "221.089",
            "Cartn_z": "177.163", "occupancy": "1.0", "B_iso_or_equiv": "74.29",
            "pdbx_formal_charge": "0", "auth_seq_id": "1306", "auth_comp_id": "NAG",
            "auth_asym_id": "C", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25670], {
            "group_PDB": "HETATM", "id": "25671", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "TA", "label_entity_id": "3", "label_seq_id": "1307",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "162.958", "Cartn_y": "220.963",
            "Cartn_z": "182.552", "occupancy": "1.0", "B_iso_or_equiv": "69.79",
            "pdbx_formal_charge": "0", "auth_seq_id": "1307", "auth_comp_id": "NAG",
            "auth_asym_id": "C", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25684], {
            "group_PDB": "HETATM", "id": "25685", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "UA", "label_entity_id": "3", "label_seq_id": "1308",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "170.895", "Cartn_y": "230.002",
            "Cartn_z": "166.592", "occupancy": "1.0", "B_iso_or_equiv": "56.9",
            "pdbx_formal_charge": "0", "auth_seq_id": "1308", "auth_comp_id": "NAG",
            "auth_asym_id": "C", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25698], {
            "group_PDB": "HETATM", "id": "25699", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "VA", "label_entity_id": "3", "label_seq_id": "1309",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "164.329", "Cartn_y": "208.107",
            "Cartn_z": "129.441", "occupancy": "1.0", "B_iso_or_equiv": "40.7",
            "pdbx_formal_charge": "0", "auth_seq_id": "1309", "auth_comp_id": "NAG",
            "auth_asym_id": "C", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25712], {
            "group_PDB": "HETATM", "id": "25713", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "WA", "label_entity_id": "3", "label_seq_id": "1310",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "177.014", "Cartn_y": "216.238",
            "Cartn_z": "134.114", "occupancy": "1.0", "B_iso_or_equiv": "32.94",
            "pdbx_formal_charge": "0", "auth_seq_id": "1310", "auth_comp_id": "NAG",
            "auth_asym_id": "C", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25726], {
            "group_PDB": "HETATM", "id": "25727", "type_symbol": "C",
            "label_atom_id": "C1", "label_alt_id": ".", "label_comp_id": "NAG",
            "label_asym_id": "XA", "label_entity_id": "3", "label_seq_id": "1311",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "194.17", "Cartn_y": "223.354",
            "Cartn_z": "243.545", "occupancy": "1.0", "B_iso_or_equiv": "77.23",
            "pdbx_formal_charge": "0", "auth_seq_id": "1311", "auth_comp_id": "NAG",
            "auth_asym_id": "C", "auth_atom_id": "C1", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25740], {
            "group_PDB": "HETATM", "id": "25741", "type_symbol": "O",
            "label_atom_id": "O", "label_alt_id": ".", "label_comp_id": "HOH",
            "label_asym_id": "YA", "label_entity_id": "4", "label_seq_id": "1401",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "182.603", "Cartn_y": "190.854",
            "Cartn_z": "142.721", "occupancy": "1.0", "B_iso_or_equiv": "15.17",
            "pdbx_formal_charge": "0", "auth_seq_id": "1401", "auth_comp_id": "HOH",
            "auth_asym_id": "A", "auth_atom_id": "O", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25823], {
            "group_PDB": "HETATM", "id": "25824", "type_symbol": "O",
            "label_atom_id": "O", "label_alt_id": ".", "label_comp_id": "HOH",
            "label_asym_id": "ZA", "label_entity_id": "4", "label_seq_id": "1401",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "210.8", "Cartn_y": "207.888",
            "Cartn_z": "210.711", "occupancy": "1.0", "B_iso_or_equiv": "40.9",
            "pdbx_formal_charge": "0", "auth_seq_id": "1401", "auth_comp_id": "HOH",
            "auth_asym_id": "B", "auth_atom_id": "O", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][25905], {
            "group_PDB": "HETATM", "id": "25906", "type_symbol": "O",
            "label_atom_id": "O", "label_alt_id": ".", "label_comp_id": "HOH",
            "label_asym_id": "AB", "label_entity_id": "4", "label_seq_id": "1401",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "157.066", "Cartn_y": "197.2",
            "Cartn_z": "210.725", "occupancy": "1.0", "B_iso_or_equiv": "37.89",
            "pdbx_formal_charge": "0", "auth_seq_id": "1401", "auth_comp_id": "HOH",
            "auth_asym_id": "C", "auth_atom_id": "O", "pdbx_PDB_model_num": "1",
        })
        self.assertEqual(d["atom_site"][-1], {
            "group_PDB": "HETATM", "id": "25990", "type_symbol": "O",
            "label_atom_id": "O", "label_alt_id": ".", "label_comp_id": "HOH",
            "label_asym_id": "AB", "label_entity_id": "4", "label_seq_id": "1485",
            "pdbx_PDB_ins_code": "?", "Cartn_x": "195.435", "Cartn_y": "177.65",
            "Cartn_z": "181.957", "occupancy": "1.0", "B_iso_or_equiv": "21.55",
            "pdbx_formal_charge": "0", "auth_seq_id": "1485", "auth_comp_id": "HOH",
            "auth_asym_id": "C", "auth_atom_id": "O", "pdbx_PDB_model_num": "1",
        })



class DictToFileTests(TestCase):

    def setUp(self):
        if os.path.exists("tests/integration/files/output/"):
            shutil.rmtree("tests/integration/files/output/")
        os.mkdir("tests/integration/files/output")


    def tearDown(self):
        if os.path.exists("tests/integration/files/output/"):
            shutil.rmtree("tests/integration/files/output/")
    

    def save(self, code):
        original = atomium.open(f"tests/integration/files/{code}.mmtf", dictionary=True)
        atomium.save_dictionary(original, f"tests/integration/files/output/{code}.mmtf")
        saved = atomium.open(f"tests/integration/files/output/{code}.mmtf", dictionary=True)
        self.assertEqual(original, saved)
    

    def test_1lol(self):
        self.save("1lol")
    

    def test_1xda(self):
        self.save("1xda")
    

    def test_5xme(self):
        self.save("5xme")


    def test_1cbn(self):
        self.save("1cbn")
    

    def test_1m4x(self):
        self.save("1m4x")
    

    def test_6xlu(self):
        self.save("6xlu")



'''class ParsingTests(TestCase):

    def test_1lol(self):
        pdb = atomium.open("tests/integration/files/1lol.mmtf")
        self.assertEqual(pdb.name, "1LOL")
        self.assertEqual(pdb.source["entry"][0]["id"], "1LOL")
        self.assertEqual(pdb.entry__id, "1LOL")
        self.assertEqual(str(pdb.model), "<Model (2 polymers, 4 non-polymers)>")
        self.assertEqual(len(pdb.model.molecules()), 8)
        self.assertEqual(len(pdb.model.polymers()), 2)
        self.assertEqual(len(pdb.model.branched_polymers()), 0)
        self.assertEqual(len(pdb.model.non_polymers()), 4)
        self.assertEqual(len(pdb.model.waters()), 2)
        self.assertEqual(len(pdb.model.residues()), 418)
        self.assertEqual(len(pdb.model.atoms()), 3431)

        mol_a = pdb.model.polymer("A")
        self.assertEqual(str(mol_a), "<Polymer A (204 residues)>")
        self.assertEqual(mol_a.id, "A")
        self.assertEqual(mol_a.auth_id, "A")
        self.assertEqual(len(mol_a), 204)
        self.assertEqual(len(mol_a.residues()), 204)
        for res in mol_a: self.assertIn(res, mol_a)
        self.assertIs(mol_a[2], mol_a.residues()[2])
        self.assertIs(mol_a.model, pdb.model)
        self.assertEqual(len(mol_a.atoms()), 1557)
        self.assertTrue(list(mol_a.atoms())[0] in mol_a)
        self.assertEqual(len(mol_a.helices), 11)
        self.assertEqual(mol_a.helices[0], mol_a.residues()[1:4])
        self.assertEqual(len(mol_a.helices[-1]), 10)
        self.assertEqual(len(mol_a.strands), 8)
        self.assertEqual(mol_a.strands[0], mol_a.residues()[4:9])
        self.assertEqual(mol_a.entity_id, "1")
        self.assertEqual(mol_a.entity_name, "orotidine 5'-monophosphate decarboxylase")
        self.assertEqual(len(mol_a.sequence), 229)
        self.assertTrue(mol_a.sequence.startswith("LRSRRVDVMDVMNRLILAMDL"))
        self.assertTrue(mol_a.sequence.endswith("LADNPAAAAAGIIESIKDLLIPE"))
        self.assertIn(mol_a, pdb.model)

        mol_c = pdb.model.non_polymer("C")
        self.assertEqual(str(mol_c), "<NonPolymer C (BU2)>")
        self.assertEqual(mol_c.id, "C")
        self.assertEqual(mol_c.auth_id, "A")
        self.assertEqual(mol_c.name, "BU2")
        self.assertIs(mol_c.model, pdb.model)
        self.assertEqual(len(mol_c.atoms()), 6)
        for atom in mol_c.atoms(): self.assertIn(atom, mol_c)
        self.assertEqual(mol_c.entity_id, "2")
        self.assertEqual(mol_c.entity_name, "1,3-BUTANEDIOL")
        self.assertIn(mol_c, pdb.model)

        mol_g = pdb.model.water("G")
        self.assertEqual(str(mol_g), "<Water G (HOH)>")
        self.assertEqual(mol_g.id, "G")
        self.assertEqual(mol_g.auth_id, "A")
        self.assertEqual(mol_g.name, "HOH")
        self.assertIs(mol_g.model, pdb.model)
        self.assertEqual(len(mol_g.atoms()), 96)
        for atom in mol_g.atoms(): self.assertIn(atom, mol_g)
        self.assertEqual(mol_g.entity_id, "4")
        self.assertEqual(mol_g.entity_name, "water")
        self.assertIn(mol_g, pdb.model)
    

    def test_1lol_mmtf_compressed(self):
        d = atomium.open("tests/integration/files/1lol.mmtf.gz", dictionary=True)
        self.assertEqual(d["entry"], [{"id": "1LOL"}])
    

    def test_5xme(self):
        # Parse
        pdb = atomium.open("tests/integration/files/5xme.mmtf")
        self.assertIsNone(pdb.resolution)

        # Models are correct
        self.assertEqual(len(pdb.models), 10)
        self.assertIs(pdb.model, pdb.models[0])
        x_values = [
            33.969, 34.064, 37.369, 36.023, 35.245,
            35.835, 37.525, 35.062, 36.244, 37.677
        ]
        all_atoms = set()
        for x, model in zip(x_values, pdb.models):
            self.assertEqual(len(model.atoms()), 1827)
            all_atoms.update(model.atoms())
            res = model.polymer()[0]
            atom = res.atom(name="N")
            self.assertEqual(atom.location[0], x)
        self.assertEqual(len(all_atoms), 18270)
    

    def test_1cbn(self):
        # Multiple occupancy is handled
        pdb = atomium.open("tests/integration/files/1cbn.mmtf")
        polymer = pdb.model.polymer()
        residue1, residue2, residue3 = polymer[:3]
        self.assertEqual(len(residue1.atoms()), 16)
        self.assertEqual(len(residue2.atoms()), 14)
        self.assertEqual(len(residue3.atoms()), 10)
        for residue in polymer[:3]:
            for name in ["N", "C", "CA", "CB"]:
                self.assertEqual(len(residue.atoms(name=name)), 1)
    

    def test_6xlu(self):
        # Carbs parsed
        pdb = atomium.open("tests/integration/files/6xlu.mmtf")

        # Model
        self.assertEqual(len(pdb.model.polymers()), 3)
        self.assertEqual(len(pdb.model.branched_polymers()), 15)
        self.assertEqual(len(pdb.model.residues()), 3208)
        self.assertEqual(len(pdb.model.non_polymers()), 32)
        self.assertEqual(len(pdb.model.atoms()), 25948)

        # Carb
        carb = pdb.model.branched_polymer("D")
        self.assertIn(carb, pdb.model)
        self.assertEqual(str(carb), "<BranchedPolymer D (2 residues)>")
        self.assertEqual(carb.id, "D")
        self.assertEqual(carb.auth_id, "D")
        self.assertEqual(len(carb.residues()), 2)
        for residue in carb.residues(): self.assertIn(residue, carb)
        self.assertIs(carb.model, pdb.model)
        self.assertEqual(len(carb.atoms()), 28)
        self.assertTrue(list(carb.atoms())[0] in carb)
        self.assertEqual(carb.entity_id, "2")
        self.assertEqual(carb.entity_name, "2-acetamido-2-deoxy-beta-D-glucopyranose-(1-4)-2-acetamido-2-deoxy-beta-D-glucopyranose")
    

    def test_1msh_data_dict_model(self):
        pdb = atomium.open("tests/integration/files/1msh.mmtf")
        self.assertEqual(len(pdb.models), 30)
        for model in pdb.models[:-1]:
            self.assertEqual(len(model.polymers()), 2)
        self.assertEqual(len(pdb.models[-1].polymers()), 1)



class NetworkTests(TestCase):

    def test_fetching_3nir_mmtf_rcsb(self):
        pdb = atomium.fetch("3nir.mmtf")
        self.assertEqual(pdb.title, "Crystal structure of small protein crambin at 0.48 A resolution")



class AssemblyTests(TestCase):

    def test_1xda(self):
        pdb = atomium.open("tests/integration/files/1xda.mmtf")

        model = pdb.generate_assembly(1)
        self.assertEqual(len(model.polymers()), 2)
        self.assertEqual(set([p.id for p in model.polymers()]), {"A", "B"})
        self.assertEqual(len(model.non_polymers()), 4)

        model = pdb.generate_assembly(2)
        self.assertEqual(len(model.polymers()), 2)
        self.assertEqual(set([p.id for p in model.polymers()]), {"C", "D"})
        self.assertEqual(len(model.non_polymers()), 4)

        model = pdb.generate_assembly(7)
        self.assertEqual(len(model.polymers()), 6)
        self.assertEqual(set([p.id for p in model.polymers()]), {"A", "B"})
        self.assertEqual(len(model.non_polymers()), 12)
        zn = model.atom(element="Zn")
        liganding_residues = zn.nearby_residues(3, is_metal=False, element__ne="CL")
        self.assertEqual(len(liganding_residues), 3)
        self.assertEqual(set([r.id for r in liganding_residues]), {"B.10"})
        self.assertEqual(set([r.name for r in liganding_residues]), {"HIS"})
        res1, res2, res3 = liganding_residues
        self.assertGreater(res1.atom(name="N").distance_to(res2.atom(name="N")), 10)
    

    def test_complex_operations(self):
        pdb = atomium.open("tests/integration/files/1m4x.mmtf")
        model = pdb.generate_assembly(6)
        self.assertEqual(len(model.polymers()), 198)
        self.assertEqual([round(c, 1) for c in model.center_of_mass], [433.1, 433.1, 433.1])

'''