import os
import shutil
from unittest import TestCase
import atomium

class FileToDictTests(TestCase):

    def test_1lol(self):
        d = atomium.open("tests/integration/files/1lol.mmtf", dictionary=True)
        self.assertEqual(len(d.keys()), 17)

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
        self.assertEqual(len(d.keys()), 17)

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
        for r in original["struct_conf"]: print(r)
        atomium.save_dictionary(original, f"tests/integration/files/output/{code}.mmtf")
        saved = atomium.open(f"tests/integration/files/output/{code}.mmtf", dictionary=True)
        self.assertEqual(original, saved)
    

    def test_1lol(self):
        self.save("1lol")
    

    def test_1xda(self):
        self.save("1xda")



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