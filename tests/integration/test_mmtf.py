from datetime import date
from unittest import TestCase
import atomium

class DictParsingTests(TestCase):

    def test_1lol_mmtf(self):
        d = atomium.open("tests/integration/files/1lol.mmtf", dictionary=True)

        self.assertEqual(d["entry"], [{"id": "1LOL"}])
        self.assertEqual(d["struct"], [{
            "entry_id": "1LOL", "title": "Crystal structure of orotidine monophosphate decarboxylase complex with XMP"
        }])
        self.assertEqual(d["pdbx_database_status"], [{
            "entry_id": "1LOL", "recvd_initial_deposition_date": "2002-05-06"
        }])
        self.assertEqual(d["exptl"], [{
            "entry_id": "1LOL", "method": "X-RAY DIFFRACTION"
        }])
        self.assertEqual(d["refine"][0]["ls_d_res_high"], "1.9")
        self.assertEqual(d["refine"][0]["ls_R_factor_R_work"], "0.193")
        self.assertEqual(d["refine"][0]["ls_R_factor_R_free"], "0.229")
        self.assertEqual(d["symmetry"], [{
            "entry_id": "1LOL", "space_group_name_H-M": "P 1 21 1"
        }])
        self.assertEqual(d["cell"][0]["length_a"], "57.57")
        self.assertEqual(d["cell"][0]["length_alpha"], "90.0")

        self.assertEqual(d["pdbx_struct_assembly"], [{"id": "1"}])
        self.assertEqual(d["pdbx_struct_assembly_gen"], [{
            "assembly_id": "1", "oper_expression": "1", "asym_id_list": "A,B,C,D,E,F,G,H"
        }])
        self.assertEqual(d["pdbx_struct_oper_list"][0]["pdbx_struct_oper_list.matrix[1][1]"], "1.0000000000")
        self.assertEqual(d["pdbx_struct_oper_list"][0]["pdbx_struct_oper_list.vector[3]"], "0.0000000000")

        self.assertEqual(len(d["entity"]), 4)       
        self.assertEqual(d["entity"][0]["pdbx_description"], "orotidine 5'-monophosphate decarboxylase")
        self.assertEqual(d["entity"][1]["pdbx_description"], "1,3-BUTANEDIOL")
        self.assertEqual(d["entity"][2]["pdbx_description"], "XANTHOSINE-5'-MONOPHOSPHATE")
        self.assertEqual(d["entity"][3]["pdbx_description"], "water")
        self.assertEqual(d["entity"][0]["type"], "polymer")
        self.assertEqual(d["entity"][1]["type"], "non-polymer")
        self.assertEqual(d["entity"][2]["type"], "non-polymer")
        self.assertEqual(d["entity"][3]["type"], "water")
        self.assertEqual(d["entity"][0]["pdbx_number_of_molecules"], "2")
        self.assertEqual(d["entity"][1]["pdbx_number_of_molecules"], "2")
        self.assertEqual(d["entity"][2]["pdbx_number_of_molecules"], "2")
        self.assertEqual(
            d["entity_poly"][0]["pdbx_seq_one_letter_code"],
            "LRSRRVDVMDVMNRLILAMDLMNRDDALRVTGEVREYIDTVKIGYPLVLSEGMDIIAEFRKRFG" +
            "CRIIADFKVADIPETNEKICRATFKAGADAIIVHGFPGADSVRACLNVAEEMGREVFLLTEMSH" +
            "PGAEMFIQGAADEIARMGVDLGVKNYVGPSTRPERLSRLREIIGQDSFLISPGVGAQGGDPGET" +
            "LRFADAIIVGRSIYLADNPAAAAAGIIESIKDLLIPE"
        )

        self.assertEqual(len(d["atom_site"]), 3431)
        self.assertEqual(d["atom_site"][0], {
            "group_PDB": "ATOM", "id": "1", "type_symbol": "N", "label_atom_id": "N",
            "label_alt_id": ".", "label_comp_id": "VAL", "label_asym_id": "A",
            "label_entity_id": "1", "label_seq_id": "11", "pdbx_PDB_ins_code": "?",
            "Cartn_x": "3.696", "Cartn_y": "33.898", "Cartn_z": "63.219", "occupancy": "1.0",
            "B_iso_or_equiv": "21.5", "pdbx_formal_charge": "0", "auth_seq_id": "11",
            "auth_comp_id": "VAL", "auth_asym_id": "A", "auth_atom_id": "N",
            "pdbx_PDB_model_num": "1"
        })
        self.assertEqual(d["atom_site"][-1], {
            "group_PDB": "HETATM", "id": "3431", "type_symbol": "O", "label_atom_id": "O",
            "label_alt_id": ".", "label_comp_id": "HOH", "label_asym_id": "H",
            "label_entity_id": "4", "label_seq_id": "3180", "pdbx_PDB_ins_code": "?",
            "Cartn_x": "-39.239", "Cartn_y": "51.357", "Cartn_z": "40.064", "occupancy": "1.0",
            "B_iso_or_equiv": "37.93", "pdbx_formal_charge": "0", "auth_seq_id": "3180",
            "auth_comp_id": "HOH", "auth_asym_id": "B", "auth_atom_id": "O", "pdbx_PDB_model_num": "1"
        })


class ParsingTests(TestCase):

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