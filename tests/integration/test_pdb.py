from unittest import TestCase
import atomium

class DictParsingTests(TestCase):

    def test_1LOL_pdb(self):
        d = atomium.open("tests/integration/files/1lol.pdb", dictionary=True)
        
        self.assertEqual(d["entry"], [{"id": "1LOL"}])
        self.assertEqual(d["struct_keywords"], [
            {"entry_id": "1LOL", "pdbx_keywords": "LYASE", "text": "TIM BARREL, LYASE"}
        ])
        self.assertEqual(d["pdbx_database_status"], [{
            "status_code": "REL", "entry_id": "1LOL",
            "recvd_initial_deposition_date": "2002-05-06"
        }])

        self.assertEqual(d["struct"], [{
            "entry_id": "1LOL", "pdbx_CASP_flag": "?", "pdbx_descriptor": "?",
            "pdbx_model_details": "?", "pdbx_model_type_details": "?",
            "title": "CRYSTAL STRUCTURE OF OROTIDINE MONOPHOSPHATE DECARBOXYLASE COMPLEX WITH XMP"
        }])

        self.assertEqual(d["exptl"], [{"entry_id": "1LOL", "method": "X-RAY DIFFRACTION"}])

        self.assertEqual(d["audit_author"], [
            {"name": "Wu, N", "pdbx_ordinal": "1"},
            {"name": "Pai, E.F", "pdbx_ordinal": "2"},
        ])

        self.assertEqual(d["pdbx_audit_revision_history"], [{
            "ordinal": "1", "data_content_type": "Structure model",
            "major_revision": "1", "minor_revision": "0", "revision_date": "2002-08-07"
        }, {
            "ordinal": "2", "data_content_type": "Structure model", 
            "major_revision": "2","minor_revision": "0", "revision_date": "2002-08-14"
        }, {
            "ordinal": "3", "data_content_type": "Structure model",
            "major_revision": "3", "minor_revision": "0", "revision_date": "2003-04-01"
        }, {
            "ordinal": "4", "data_content_type": "Structure model",
            "major_revision": "4", "minor_revision": "0", "revision_date": "2009-02-24"
        }])

        self.assertEqual(d["citation"], [{
            "id": "primary",
            "title": "CRYSTAL STRUCTURES OF INHIBITOR COMPLEXES REVEAL AN ALTERNATE BINDING MODE IN OROTIDINE-5'-MONOPHOSPHATE DECARBOXYLASE.",
            "journal_abbrev": "J.Biol.Chem.", "journal_volume": "277",
            "page_first": "28080", "page_last": "?", "year": "2002", "journal_id_ASTM": "?",
            "country": "?", "journal_id_ISSN": "0021-9258", "journal_id_CSD": "?",
            "book_publisher": "?", "pdbx_database_id_PubMed": "12011084",
            "pdbx_database_id_DOI": "10.1074/jbc.m202362200"
        }])

        self.assertEqual(d["citation_author"], [
            {"citation_id": "primary", "name": "Wu, N", "pdbx_ordinal": "1"},
            {"citation_id": "primary", "name": "Pai, E.F", "pdbx_ordinal": "2"}
        ])
        self.assertEqual(d["reflns"], [{
            "B_iso_Wilson_estimate": "?", "entry_id": "1LOL", "data_reduction_details": "?",
            "data_reduction_method": "?", "d_resolution_high": "1.90", "d_resolution_low": "27.07",
            "details": "?", "limit_h_max": "?", "limit_h_min": "?", "limit_k_max": "?",
            "limit_k_min": "?", "limit_l_max": "?", "limit_l_min": "?", "number_all": "?",
            "number_obs": "?", "observed_criterion": "?", "observed_criterion_F_max": "?",
            "observed_criterion_F_min": "?", "observed_criterion_I_max": "?",
            "observed_criterion_I_min": "?", "observed_criterion_sigma_F": "?",
            "observed_criterion_sigma_I": "?", "percent_possible_obs": "?", "R_free_details": "?",
            "Rmerge_F_all": "?", "Rmerge_F_obs": "?", "Friedel_coverage": "?", "number_gt": "?",
            "threshold_expression": "?", "pdbx_redundancy": "?", "pdbx_Rmerge_I_obs": "?",
            "pdbx_Rmerge_I_all": "?", "pdbx_Rsym_value": "?", "pdbx_netI_over_av_sigmaI": "?",
            "pdbx_netI_over_sigmaI": "?", "pdbx_res_netI_over_av_sigmaI_2": "?",
            "pdbx_res_netI_over_sigmaI_2": "?", "pdbx_chi_squared": "?", "pdbx_scaling_rejects": "?",
            "pdbx_d_res_high_opt": "?", "pdbx_d_res_low_opt": "?", "pdbx_d_res_opt_method": "?",
            "phase_calculation_details": "?", "pdbx_Rrim_I_all": "?", "pdbx_Rpim_I_all": "?",
            "pdbx_d_opt": "?", "pdbx_number_measured_all": "?", "pdbx_diffrn_id": "?",
            "pdbx_ordinal": "?", "pdbx_CC_half": "?", "pdbx_R_split": "?",
            "ls_R_factor_R_work": "0.193", "ls_R_factor_R_free": "0.229",
            "ls_percent_reflns_R_free": "4.900", "ls_number_reflns_R_free": "1583"
        }])

        self.assertEqual(d["diffrn"], [
            {"id": "1", "crystal_id": "1", "ambient_temp": "100", "ambient_temp_details": "?"}
        ])
        self.assertEqual(d["diffrn_detector"], [{
            "details": "?", "detector": "CCD", "diffrn_id": "1",
            "pdbx_collection_date": "2001-09-12", "type": "ADSC QUANTUM 4"
        }])

        self.assertEqual(len(d["pdbx_unobs_or_zero_occ_residues"]), 38)
        self.assertEqual(d["pdbx_unobs_or_zero_occ_residues"][0], {
            "id": "1", "PDB_model_num": "1", "polymer_flag": "Y", "occupancy_flag": "1",
            "auth_asym_id": "A", "auth_comp_id": "SER", "auth_seq_id": "3", "PDB_ins_code": "?",
            "label_asym_id": "A", "label_comp_id": "SER", "label_seq_id": "3"
        })
        self.assertEqual(d["pdbx_unobs_or_zero_occ_residues"][-1], {
            "id": "38", "PDB_model_num": "1", "polymer_flag": "Y", "occupancy_flag": "1",
            "auth_asym_id": "B", "auth_comp_id": "GLY", "auth_seq_id": "1186", "PDB_ins_code": "?",
            "label_asym_id": "B", "label_comp_id": "GLY", "label_seq_id": "1186"
        })

        self.assertEqual(d["struct_site"], [{
            "id": "AC1", "pdbx_evidence_code": "Software", "pdbx_auth_asym_id": "?",
            "pdbx_auth_comp_id": "?", "pdbx_auth_seq_id": "?", "pdbx_auth_ins_code": "?",
            "pdbx_num_residues": "?", "details": "BINDING SITE FOR RESIDUE BU2 A 5001"
        }, {
            "id": "AC2", "pdbx_evidence_code": "Software", "pdbx_auth_asym_id": "?",
            "pdbx_auth_comp_id": "?", "pdbx_auth_seq_id": "?", "pdbx_auth_ins_code": "?",
            "pdbx_num_residues": "?", "details": "BINDING SITE FOR RESIDUE BU2 B 5002"
        }, {
            "id": "AC3", "pdbx_evidence_code": "Software", "pdbx_auth_asym_id": "?",
            "pdbx_auth_comp_id": "?", "pdbx_auth_seq_id": "?", "pdbx_auth_ins_code": "?",
            "pdbx_num_residues": "?", "details": "BINDING SITE FOR RESIDUE XMP A 2001"
        }, {
            "id": "AC4", "pdbx_evidence_code": "Software", "pdbx_auth_asym_id": "?",
            "pdbx_auth_comp_id": "?", "pdbx_auth_seq_id": "?", "pdbx_auth_ins_code": "?",
            "pdbx_num_residues": "?", "details": "BINDING SITE FOR RESIDUE XMP B 2002"
        }])

        self.assertEqual(d["cell"], [{
            "entry_id": "1LOL", "length_a": "57.570", "length_b": "55.482", 
            "length_c": "66.129", "angle_alpha": "90.00", "angle_beta": "94.28",
            "angle_gamma": "90.00", "Z_pdb": "4", "pdbx_unique_axis": "?"
        }])
        self.assertEqual(d["symmetry"], [{
            "entry_id": "1LOL", "space_group_name_H-M": "P 1 21 1", 
            "pdbx_full_space_group_name_H-M": "?",
            "cell_setting": "?", "Int_Tables_number": "?"
        }])
        self.assertEqual(d["database_PDB_matrix"], [{
            "entry_id": "1LOL", "origx[1][1]": "1.000000", "origx[2][1]": "0.000000",
            "origx[3][1]": "0.000000", "origx[1][2]": "0.000000", "origx[2][2]": "1.000000",
            "origx[3][2]": "0.000000", "origx[1][3]": "0.000000", "origx[2][3]": "0.000000",
            "origx[3][3]": "1.000000", "origx_vector[1]": "0.00000",
            "origx_vector[2]": "0.00000", "origx_vector[3]": "0.00000"
        }])
        self.assertEqual(d["atom_sites"], [{
            "entry_id": "1LOL", "fract_transf_matrix[1][1]": "0.017370",
            "fract_transf_matrix[2][1]": "0.000000", "fract_transf_matrix[3][1]": "0.000000",
            "fract_transf_matrix[1][2]": "0.000000", "fract_transf_matrix[2][2]": "0.018024",
            "fract_transf_matrix[3][2]": "0.000000", "fract_transf_matrix[1][3]": "0.001301",
            "fract_transf_matrix[2][3]": "0.000000", "fract_transf_matrix[3][3]": "0.015164",
            "fract_transf_vector[1]": "0.00000", "fract_transf_vector[2]": "0.00000",
            "fract_transf_vector[3]": "0.00000", "fract_transf_matrix[P][1]": "?",
            "fract_transf_matrix[P][2]": "?", "fract_transf_matrix[P][3]": "?",
            "fract_transf_vector[P]": "?"
        }])

        self.assertEqual(d["entity"], [{
            "id": "1", "type": "polymer", "src_method": "man",
            "pdbx_description": "OROTIDINE 5'-MONOPHOSPHATE DECARBOXYLASE",
            "formula_weight": "?", "pdbx_number_of_molecules": "2",
            "pdbx_ec": "4.1.1.23", "pdbx_mutation": "?", "pdbx_fragment": "?",
            "details": "?"
        }, {
            "id": "2", "type": "non-polymer", "src_method": "syn",
            "pdbx_description": "1,3-BUTANEDIOL",
            "formula_weight": "?", "pdbx_number_of_molecules": "2",
            "pdbx_ec": "?", "pdbx_mutation": "?", "pdbx_fragment": "?",
            "details": "?"
        }, {
            "id": "3", "type": "non-polymer", "src_method": "syn",
            "pdbx_description": "XANTHOSINE-5'-MONOPHOSPHATE",
            "formula_weight": "?", "pdbx_number_of_molecules": "2",
            "pdbx_ec": "?", "pdbx_mutation": "?", "pdbx_fragment": "?",
            "details": "?"
        }, {
            "id": "4", "type": "water", "src_method": "nat",
            "pdbx_description": "water",
            "formula_weight": "?", "pdbx_number_of_molecules": "180",
            "pdbx_ec": "?", "pdbx_mutation": "?", "pdbx_fragment": "?",
            "details": "?"
        }])

        self.assertEqual(d["entity_name_com"], [{
            "entity_id": "1", "name": "OMP DECARBOXYLASE, OMPDCASE, OMPDECASE"
        }])

        self.assertEqual(d["entity_poly"], [{
            "entity_id": "1", "type": "polypeptide(L)", "nstd_linkage": "no", "nstd_monomer": "no",
            "pdbx_seq_one_letter_code": "LRSRRVDVMDVMNRLILAMDLMNRDDALRVTGEVREYIDTVKIGYPLVLSEGMDIIAEFRKRFGCRIIADFKVADIPETNEKICRATFKAGADAIIVHGFPGADSVRACLNVAEEMGREVFLLTEMSHPGAEMFIQGAADEIARMGVDLGVKNYVGPSTRPERLSRLREIIGQDSFLISPGVGAQGGDPGETLRFADAIIVGRSIYLADNPAAAAAGIIESIKDLLIPE",
            "pdbx_seq_one_letter_code_can": "LRSRRVDVMDVMNRLILAMDLMNRDDALRVTGEVREYIDTVKIGYPLVLSEGMDIIAEFRKRFGCRIIADFKVADIPETNEKICRATFKAGADAIIVHGFPGADSVRACLNVAEEMGREVFLLTEMSHPGAEMFIQGAADEIARMGVDLGVKNYVGPSTRPERLSRLREIIGQDSFLISPGVGAQGGDPGETLRFADAIIVGRSIYLADNPAAAAAGIIESIKDLLIPE",
            "pdbx_strand_id": "A,B", "pdbx_target_identifier": "?"
        }])

        self.assertEqual(d["entity_poly_seq"][0], {
            "entity_id": "1", "num": "1", "mon_id": "LEU", "hetero": "n"
        })
        self.assertEqual(d["entity_poly_seq"][1], {
            "entity_id": "1", "num": "2", "mon_id": "ARG", "hetero": "n"
        })
        self.assertEqual(d["entity_poly_seq"][-1], {
            "entity_id": "1", "num": "229", "mon_id": "GLU", "hetero": "n"
        })

        self.assertEqual(d["struct_ref"][0], {
            "id": "1", "db_name": "UNP", "db_code": "PYRF_METTH", "entity_id": "1",
            "pdbx_seq_one_letter_code": "?", "pdbx_align_begin": "1",
            "pdbx_db_accession": "O26232", "pdbx_db_isoform": "?"
        })
        
        self.assertEqual(d["struct_ref_seq"], [{
            "align_id": "1", "ref_id": "1", "pdbx_PDB_id_code": "1LOL", "pdbx_strand_id": "A",
            "seq_align_beg": "?", "pdbx_seq_align_beg_ins_code": "?", "seq_align_end": "?",
            "pdbx_seq_align_end_ins_code": "?", "pdbx_db_accession": "O26232", "db_align_beg": "1",
            "pdbx_db_align_beg_ins_code": "?", "db_align_end": "228", "pdbx_db_align_end_ins_code": "?",
            "pdbx_auth_seq_align_beg": "1", "pdbx_auth_seq_align_end": "229"
        }, {
            "align_id": "2", "ref_id": "1", "pdbx_PDB_id_code": "1LOL", "pdbx_strand_id": "B",
            "seq_align_beg": "?", "pdbx_seq_align_beg_ins_code": "?", "seq_align_end": "?",
            "pdbx_seq_align_end_ins_code": "?", "pdbx_db_accession": "O26232", "db_align_beg": "1",
            "pdbx_db_align_beg_ins_code": "?", "db_align_end": "228", "pdbx_db_align_end_ins_code": "?",
            "pdbx_auth_seq_align_beg": "1001", "pdbx_auth_seq_align_end": "1229"
        }])

        self.assertEqual(len(d["struct_ref_seq_dif"]), 8)
        self.assertEqual(d["struct_ref_seq_dif"][0], {
            "align_id": "1", "pdbx_pdb_id_code": "1LOL", "mon_id": "LEU",
            "pdbx_pdb_strand_id": "A", "seq_num": "?", "pdbx_pdb_ins_code": "?",
            "pdbx_seq_db_name": "UNP", "pdbx_seq_db_accession_code": "O26232",
            "db_mon_id": "MET", "pdbx_seq_db_seq_num": "1",
            "details": "SEE REMARK 999", "pdbx_auth_seq_num": "1", "pdbx_ordinal": "1"
        })
        self.assertEqual(d["struct_ref_seq_dif"][-1], {
            "align_id": "2", "pdbx_pdb_id_code": "1LOL", "mon_id": "GLU",
            "pdbx_pdb_strand_id": "B", "seq_num": "?", "pdbx_pdb_ins_code": "?",
            "pdbx_seq_db_name": "UNP", "pdbx_seq_db_accession_code": "O26232",
            "db_mon_id": "?", "pdbx_seq_db_seq_num": "?", "details": "INSERTION",
            "pdbx_auth_seq_num": "1229", "pdbx_ordinal": "8"
        })

        self.assertEqual(d["pdbx_entity_nonpoly"], [
            {"entity_id": "2", "name": "1,3-BUTANEDIOL", "comp_id": "BU2"},
            {"entity_id": "3", "name": "XANTHOSINE-5'-MONOPHOSPHATE", "comp_id": "XMP"},
            {"entity_id": "4", "name": "water", "comp_id": "HOH"},
        ])

        self.assertEqual(d["chem_comp"][0], {
            "id": "ALA", "type": "L-peptide linking", "mon_nstd_flag": "y", "name": "ALANINE",
            "pdbx_synonyms": "?", "formula": "C3 H7 N O2", "formula_weight": "89.093"
        })
        self.assertEqual(d["chem_comp"][1], {
            "id": "ARG", "type": "L-peptide linking", "mon_nstd_flag": "y", "name": "ARGININE",
            "pdbx_synonyms": "?", "formula": "C6 H15 N4 O2 1", "formula_weight": "175.209"
        })
        self.assertEqual(d["chem_comp"][4], {
            "id": "BU2", "type": "non-polymer", "mon_nstd_flag": ".", "name": "1,3-BUTANEDIOL",
            "pdbx_synonyms": "?", "formula": "C4 H10 O2", "formula_weight": "90.121"
        })
        self.assertEqual(d["chem_comp"][10], {
            "id": "HOH", "type": "non-polymer", "mon_nstd_flag": ".", "name": "WATER",
            "pdbx_synonyms": "?", "formula": "H2 O", "formula_weight": "18.015"
        })
        self.assertEqual(d["chem_comp"][21], {
            "id": "XMP", "type": "non-polymer", "mon_nstd_flag": ".", "name": "XANTHOSINE-5'-MONOPHOSPHATE",
            "pdbx_synonyms": "5-MONOPHOSPHATE-9-BETA-D-RIBOFURANOSYL XANTHINE",
            "formula": "C10 H14 N4 O9 P 1+", "formula_weight": "365.213"
        })

        self.assertEqual(d["atom_type"], [
            {"symbol": "C"}, {"symbol": "N"}, {"symbol": "O"}, {"symbol": "P"}, {"symbol": "S"}
        ])

        self.assertEqual(d["atom_site"][0], {
            "group_pdb": "ATOM", "id": "1", "type_symbol": "N", "label_atom_id": "N",
            "label_alt_id": ".", "label_comp_id": "VAL", "label_asym_id": "A",
            "label_entity_id": "1", "label_seq_id": "11", "pdbx_PDB_ins_code": "?",
            "Cartn_x": "3.696", "Cartn_y": "33.898", "Cartn_z": "63.219",
            "occupancy": "1.00", "B_iso_or_equiv": "21.50", "pdbx_formal_charge": "?",
            "auth_seq_id": "11", "auth_comp_id": "VAL", "auth_asym_id": "A",
            "auth_atom_id": "N", "pdbx_PDB_model_num": "1"
        })
        self.assertEqual(d["atom_site"][1557], {
            "group_pdb": "ATOM", "id": "1558", "type_symbol": "N", "label_atom_id": "N",
            "label_alt_id": ".", "label_comp_id": "VAL", "label_asym_id": "B",
            "label_entity_id": "1", "label_seq_id": "11", "pdbx_PDB_ins_code": "?",
            "Cartn_x": "-26.384", "Cartn_y": "61.433", "Cartn_z": "36.898",
            "occupancy": "1.00", "B_iso_or_equiv": "39.30", "pdbx_formal_charge": "?",
            "auth_seq_id": "1011", "auth_comp_id": "VAL", "auth_asym_id": "B",
            "auth_atom_id": "N", "pdbx_PDB_model_num": "1"
        })
        self.assertEqual(d["atom_site"][3191], {
            "group_pdb": "HETATM", "id": "3192", "type_symbol": "C", "label_atom_id": "C1",
            "label_alt_id": ".", "label_comp_id": "BU2", "label_asym_id": "C",
            "label_entity_id": "2", "label_seq_id": ".", "pdbx_PDB_ins_code": "?",
            "Cartn_x": "2.646", "Cartn_y": "45.112", "Cartn_z": "48.995",
            "occupancy": "1.00", "B_iso_or_equiv": "43.24", "pdbx_formal_charge": "?",
            "auth_seq_id": "5001", "auth_comp_id": "BU2", "auth_asym_id": "A",
            "auth_atom_id": "C1", "pdbx_PDB_model_num": "1"
        })
        self.assertEqual(d["atom_site"][3197], {
            "group_pdb": "HETATM", "id": "3198", "type_symbol": "P", "label_atom_id": "P",
            "label_alt_id": ".", "label_comp_id": "XMP", "label_asym_id": "D",
            "label_entity_id": "3", "label_seq_id": ".", "pdbx_PDB_ins_code": "?",
            "Cartn_x": "3.293", "Cartn_y": "36.948", "Cartn_z": "44.605",
            "occupancy": "1.00", "B_iso_or_equiv": "20.47", "pdbx_formal_charge": "?",
            "auth_seq_id": "2001", "auth_comp_id": "XMP", "auth_asym_id": "A",
            "auth_atom_id": "P", "pdbx_PDB_model_num": "1"
        })
        self.assertEqual(d["atom_site"][3221], {
            "group_pdb": "HETATM", "id": "3222", "type_symbol": "C", "label_atom_id": "C1",
            "label_alt_id": ".", "label_comp_id": "BU2", "label_asym_id": "E",
            "label_entity_id": "2", "label_seq_id": ".", "pdbx_PDB_ins_code": "?",
            "Cartn_x": "-14.563", "Cartn_y": "61.208", "Cartn_z": "49.005",
            "occupancy": "1.00", "B_iso_or_equiv": "45.50", "pdbx_formal_charge": "?",
            "auth_seq_id": "5002", "auth_comp_id": "BU2", "auth_asym_id": "B",
            "auth_atom_id": "C1", "pdbx_PDB_model_num": "1"
        })
        self.assertEqual(d["atom_site"][3227], {
            "group_pdb": "HETATM", "id": "3228", "type_symbol": "P", "label_atom_id": "P",
            "label_alt_id": ".", "label_comp_id": "XMP", "label_asym_id": "F",
            "label_entity_id": "3", "label_seq_id": ".", "pdbx_PDB_ins_code": "?",
            "Cartn_x": "-23.846", "Cartn_y": "63.007", "Cartn_z": "53.020",
            "occupancy": "1.00", "B_iso_or_equiv": "33.88", "pdbx_formal_charge": "?",
            "auth_seq_id": "2002", "auth_comp_id": "XMP", "auth_asym_id": "B",
            "auth_atom_id": "P", "pdbx_PDB_model_num": "1"
        })
        self.assertEqual(d["atom_site"][3251], {
            "group_pdb": "HETATM", "id": "3252", "type_symbol": "O", "label_atom_id": "O",
            "label_alt_id": ".", "label_comp_id": "HOH", "label_asym_id": "G",
            "label_entity_id": "4", "label_seq_id": ".", "pdbx_PDB_ins_code": "?",
            "Cartn_x": "-7.624", "Cartn_y": "43.710", "Cartn_z": "47.691",
            "occupancy": "1.00", "B_iso_or_equiv": "15.72", "pdbx_formal_charge": "?",
            "auth_seq_id": "3005", "auth_comp_id": "HOH", "auth_asym_id": "A",
            "auth_atom_id": "O", "pdbx_PDB_model_num": "1"
        })
        self.assertEqual(d["atom_site"][3252], {
            "group_pdb": "HETATM", "id": "3253", "type_symbol": "O", "label_atom_id": "O",
            "label_alt_id": ".", "label_comp_id": "HOH", "label_asym_id": "G",
            "label_entity_id": "4", "label_seq_id": ".", "pdbx_PDB_ins_code": "?",
            "Cartn_x": "-8.138", "Cartn_y": "37.776", "Cartn_z": "46.819",
            "occupancy": "1.00", "B_iso_or_equiv": "31.82", "pdbx_formal_charge": "?",
            "auth_seq_id": "3007", "auth_comp_id": "HOH", "auth_asym_id": "A",
            "auth_atom_id": "O", "pdbx_PDB_model_num": "1"
        })
        self.assertEqual(d["atom_site"][3347], {
            "group_pdb": "HETATM", "id": "3348", "type_symbol": "O", "label_atom_id": "O",
            "label_alt_id": ".", "label_comp_id": "HOH", "label_asym_id": "H",
            "label_entity_id": "4", "label_seq_id": ".", "pdbx_PDB_ins_code": "?",
            "Cartn_x": "-9.220", "Cartn_y": "50.800", "Cartn_z": "49.092",
            "occupancy": "1.00", "B_iso_or_equiv": "17.69", "pdbx_formal_charge": "?",
            "auth_seq_id": "3001", "auth_comp_id": "HOH", "auth_asym_id": "B",
            "auth_atom_id": "O", "pdbx_PDB_model_num": "1"
        })
        self.assertEqual(d["atom_site"][3348], {
            "group_pdb": "HETATM", "id": "3349", "type_symbol": "O", "label_atom_id": "O",
            "label_alt_id": ".", "label_comp_id": "HOH", "label_asym_id": "H",
            "label_entity_id": "4", "label_seq_id": ".", "pdbx_PDB_ins_code": "?",
            "Cartn_x": "-8.391", "Cartn_y": "48.561", "Cartn_z": "47.191",
            "occupancy": "1.00", "B_iso_or_equiv": "23.47", "pdbx_formal_charge": "?",
            "auth_seq_id": "3002", "auth_comp_id": "HOH", "auth_asym_id": "B",
            "auth_atom_id": "O", "pdbx_PDB_model_num": "1"
        })
    

    def test_1CK8_pdb(self):
        # Tests OBSLTE, CAVEAT, AUTHOR
        d = atomium.open("tests/integration/files/1ck8.pdb", dictionary=True)
        self.assertEqual(d["pdbx_database_PDB_obs_spr"], [{
            "id": "OBSLTE", "date": "2006-07-25", "pdb_id": "1T0K",
            "replace_pdb_id": "1CK8", "details": "?"
        }])
        self.assertEqual(d["database_PDB_caveat"], [{
            "text": "THERE ARE CHIRALITY ERRORS IN C-ALPHA CENTERS"
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
    

    def test_4y60_pdb(self):
        # Test ANISOU
        d = atomium.open("tests/integration/files/4y60.pdb", dictionary=True)
        
        self.assertEqual(d["entry"], [{"id": "4Y60"}])

        self.assertEqual(d["atom_site"][0], {
            "group_pdb": "ATOM", "id": "1", "type_symbol": "N", "label_atom_id": "N",
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
            "group_pdb": "ATOM", "id": "1", "type_symbol": "N", "label_atom_id": "N",
            "label_alt_id": ".", "label_comp_id": "ALA", "label_asym_id": "A",
            "label_entity_id": "1", "label_seq_id": "14", "pdbx_PDB_ins_code": "?",
            "Cartn_x": "33.969", "Cartn_y": "-8.430", "Cartn_z": "-0.271",
            "occupancy": "1.00", "B_iso_or_equiv": "0.00", "pdbx_formal_charge": "?",
            "auth_seq_id": "199", "auth_comp_id": "ALA", "auth_asym_id": "A",
            "auth_atom_id": "N", "pdbx_PDB_model_num": "1"
        })

        self.assertEqual(d["atom_site"][1827], {
            "group_pdb": "ATOM", "id": "1828", "type_symbol": "N", "label_atom_id": "N",
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
        self.assertEqual(str(pdb.model), "<Model (2 polymers, 4 non-polymers)>")