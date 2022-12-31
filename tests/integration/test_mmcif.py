import os
import shutil
from unittest import TestCase
import atomium

class FileToDictTests(TestCase):

    def test_1lol(self):
        d = atomium.open("tests/integration/files/1lol.cif", dictionary=True)
        self.assertEqual(len(d.keys()), 66)
        self.assertEqual(d["entry"], [{"id": "1LOL"}])
        self.assertEqual(d["database_2"], [
            {"database_id": "PDB", "database_code": "1LOL"},
            {"database_id": "RCSB", "database_code": "RCSB016137"},
            {"database_id": "WWPDB", "database_code": "D_1000016137"},
        ])
        self.assertEqual(
            d["pdbx_database_related"][0]["details"],
            "orotidine monophosphate decarboxylase complexed with UMP"
        )
        self.assertEqual(d["entity"][2]["pdbx_description"], "XANTHOSINE-5'-MONOPHOSPHATE")
        self.assertEqual(
            d["entity_poly"][0]["pdbx_seq_one_letter_code"],
            "LRSRRVDVMDVMNRLILAMDLMNRDDALRVTGEVREYIDTVKIGYPLVLSEGMDIIAEFRKRFG" +
            "CRIIADFKVADIPETNEKICRATFKAGADAIIVHGFPGADSVRACLNVAEEMGREVFLLTEMSH" +
            "PGAEMFIQGAADEIARMGVDLGVKNYVGPSTRPERLSRLREIIGQDSFLISPGVGAQGGDPGET" +
            "LRFADAIIVGRSIYLADNPAAAAAGIIESIKDLLIPE"
        )
        self.assertEqual(
            d["struct_biol"][0]["details"],
            "The biologically functional unit is a dimer composed of the two monomers in the asymmetric unit."
        )
        self.assertEqual(
            d["pdbx_database_remark"][0]["text"], "\n".join([
                "sequence",
                "Author states that although residues 1 and 1001 are MET",
                "and residues 101 and 1101 are Arg according to the",
                "SwissProt entry, residues 1 and 1001 were LEU and residues",
                "101 and 1101 were Pro in the original construct cloned",
                "of MT genomic dna.",
            ])
        )
    

    def test_1lol_compressed(self):
        d = atomium.open("tests/integration/files/1lol.cif.gz", dictionary=True)
        self.assertEqual(len(d.keys()), 66)
        self.assertEqual(d["entry"], [{"id": "1LOL"}])
    

    def test_2bfb(self):
        d = atomium.open("tests/integration/files/2bfb.cif", dictionary=True)
        self.assertEqual(d["pdbx_database_related"][0], {
            "db_name": "PDB", "db_id": "1DTW", "content_type": "unspecified",
            "details": "HUMAN BRANCHED-CHAIN ALPHA-KETO ACID DEHYDROGENASE"
        })
        self.assertEqual(d["pdbx_database_related"][1], {
            "db_name": "PDB", "db_id": "1OLS", "content_type": "unspecified",
            "details": "ROLES OF HIS291-ALPHA AND HIS146-BETA' IN THE REDUCTIVE"
            " ACYLATION REACTION CATALYZED BY HUMAN BRANCHED-CHAIN ALPHA-KETOACID DEHYDROGENASE"
        })
    

    def test_6uaj(self):
        d = atomium.open("tests/integration/files/6uaj.cif", dictionary=True)
        self.assertEqual(d["pdbx_struct_assembly_gen"], [{
            "assembly_id": "1", "oper_expression": "1",
            "asym_id_list": "A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z"
            ",AA,BA,CA,DA,EA,FA,GA,HA,IA,JA,KA,LA,MA,NA,OA,PA,QA,RA,SA,TA,UA,VA"
        }])
    

    def test_4jtd(self):
        d = atomium.open("tests/integration/files/4jtd.cif", dictionary=True)
        self.assertEqual(d["refine"][0]["details"], (
            "THERE IS A TOXIN MOLECULE BOUND TO THE CHANNEL TETRAMER GENERATED BY\n"
            "FOUR COPIES OF A TOGETHER WITH FOUR COPIES OF B. HOWEVER IT WAS NOT BUILT\n"
            "BECAUSE IT WAS NOT SUFFICIENTLY WELL ORDERED.\n"
            "\n"
            "RESIDUES 133-144 IN CHAIN B WAS BUILT AS A POLYGLYCINE CHAIN BECAUSE\n"
            "OF LACK OF ADEQUATE ELECTRON DENSITY FOR THE SIDE CHAINS."
        ))
    

    def test_2kg4(self):
        # Handle spaces
        d = atomium.open("tests/integration/files/2kg4.cif", dictionary=True)
        self.assertEqual(d["pdbx_nmr_sample_details"], [{
            "contents": "0.45mM [U-99% 13C; U-99% 15N] gadd45a-1, 20mM sodium phosphate-2, 100mM potassium chloride-3, 2mM DTT-4, 2mM EDTA-5, 0.03% sodium azide-6, 93% H2O/7% D2O",
            "solution_id": "1", "solvent_system": "93% H2O/7% D2O"
        }, {
            "contents": "0.45mM [U-99% 15N] gadd45a-7, 20mM sodium phosphate-8, 100mM potassium chloride-9, 2mM DTT-10, 2mM EDTA-11, 0.03% sodium azide-12, 93% H2O/7% D2O",
            "solution_id": "2", "solvent_system": "93% H2O/7% D2O"
        }, {
            "contents": "0.195mM gadd45a-13, 20mM sodium phosphate-14, 100mM potassium chloride-15, 2mM DTT-16, 2mM EDTA-17, 0.03% sodium azide-18, 100% D2O",
            "solution_id": "3", "solvent_system": "100% D2O"
        }, {
            "contents": "0.45mM [U-15N]-Leu gadd45a-19, 20mM sodium phosphate-20, 100mM potassium chloride-21, 2mM DTT-22, 2mM EDTA-23, 0.03% sodium azide-24, 93% H2O/7% D2O",
            "solution_id": "4", "solvent_system": "93% H2O/7% D2O"
        }, {
            "contents": "0.45mM [U-15N]-Val gadd45a-25, 20mM sodium phosphate-26, 100mM potassium chloride-27, 2mM DTT-28, 2mM EDTA-29, 0.03% sodium azide-30, 93% H2O/7% D2O",
            "solution_id": "5", "solvent_system": "93% H2O/7% D2O"
        }, {
            "contents": "0.45mM [U-15N]-Ile gadd45a-31, 20mM sodium phosphate-32, 100mM potassium chloride-33, 2mM DTT-34, 2mM EDTA-35, 0.03% sodium azide-36, 93% H2O/7% D2O",
            "solution_id": "6", "solvent_system": "93% H2O/7% D2O"
        }, {
            "contents": "0.45mM [U-15N]-Phe/Tyr gadd45a-37, 20mM sodium phosphate-38, 100mM potassium chloride-39, 2mM DTT-40, 2mM EDTA-41, 0.03% sodium azide-42, 93% H2O/7% D2O",
            "solution_id": "7", "solvent_system": "93% H2O/7% D2O"
        }, {
            "contents": "0.45mM [U-99% 15N] gadd45a-43, 20mM sodium phosphate-44, 100mM potassium chloride-45, 2mM DTT-46, 2mM EDTA-47, 0.03% sodium azide-48, 10 mg/mL Pf1 phage-49, 93% H2O/7% D2O",
            "solution_id": "8", "solvent_system": "93% H2O/7% D2O"
        }])



class DictToFileTests(TestCase):

    def setUp(self):
        if os.path.exists("tests/integration/files/output/"):
            shutil.rmtree("tests/integration/files/output/")
        os.mkdir("tests/integration/files/output")


    def tearDown(self):
        if os.path.exists("tests/integration/files/output/"):
            shutil.rmtree("tests/integration/files/output/")
    

    def save(self, code):
        original = atomium.open(f"tests/integration/files/{code}.cif", dictionary=True)
        atomium.save_dictionary(original, f"tests/integration/files/output/{code}.cif")
        saved = atomium.open(f"tests/integration/files/output/{code}.cif", dictionary=True)
        self.assertEqual(original, saved)
    

    def test_1lol(self):
        self.save("1lol")
    

    def test_1xda(self):
        self.save("1xda")
    

    def test_1m4x(self):
        self.save("1m4x")
    

    def test_12ca(self):
        self.save("12ca")
    

    def test_1cbn(self):
        self.save("1cbn")
    

    def test_2bfb(self):
        self.save("2bfb")
    

    def test_2igd(self):
        self.save("2igd")
    

    def test_3jbp(self):
        self.save("3jbp")
    

    def test_6xlu(self):
        self.save("6xlu")


    def test_6uaj(self):
        self.save("6uaj")
    
    
    def test_4jtd(self):
        self.save("4jtd")
    

    def test_3cvx(self):
        self.save("3cvx")
    

    def test_2fur(self):
        self.save("2fur")
    

    def test_2fur(self):
        self.save("2kg4")