import os
import shutil
from unittest import TestCase
import atomium

class FileToDictTests(TestCase):

    def test_1lol(self):
        d = atomium.open("tests/integration/files/1lol.bcif", dictionary=True)
        self.assertEqual(len(d.keys()), 66)
        self.assertEqual(d["entry"], [{"id": "1LOL"}])
        self.assertEqual(d["database_2"], [
            {"database_id": "PDB", "database_code": "1LOL"},
            {"database_id": "RCSB", "database_code": "RCSB016137"},
            {"database_id": "WWPDB", "database_code": "D_1000016137"},
        ])
        self.assertEqual(d["pdbx_database_status"][0]["SG_entry"], "?")
        self.assertEqual(
            d["pdbx_database_related"][0]["details"],
            "orotidine monophosphate decarboxylase complexed with UMP"
        )

        self.assertEqual(d["entity"][0], {
            "id": "1", "type": "polymer", "src_method": "man",
            "pdbx_description": "orotidine 5'-monophosphate decarboxylase",
            "formula_weight": "24994.83", "pdbx_number_of_molecules": "2",
            "pdbx_ec": "4.1.1.23", "pdbx_mutation": "?", "pdbx_fragment": "?",
            "details": "?"
        })
        self.assertEqual(d["entity"][1], {
            "id": "2", "type": "non-polymer", "src_method": "syn",
            "pdbx_description": "1,3-BUTANEDIOL", "formula_weight": "90.121",
            "pdbx_number_of_molecules": "2", "pdbx_ec": "?", "pdbx_mutation": "?",
            "pdbx_fragment": "?", "details": "?"
        })

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
    

    def test_5o5u(self):
        # Floating point errors
        d = atomium.open("tests/integration/files/5o5u.bcif", dictionary=True)
        self.assertEqual(len(d.keys()), 76)
        self.assertEqual(d["entry"], [{"id": "5O5U"}])
        self.assertEqual(d["citation"], [{
            "abstract": "?", "abstract_id_CAS": "?", "book_id_ISBN": "?", "book_publisher": "?",
            "book_publisher_city": "?", "book_title": "?", "coordinate_linkage": "?",
            "country": "?", "database_id_Medline": "?", "details": "?", "id": "primary",
            "journal_abbrev": "To Be Published", "journal_id_ASTM": "?", "journal_id_CSD": "353",
            "journal_id_ISSN": "?", "journal_full": "?", "journal_issue": "?", "journal_volume": "?",
            "language": "?", "page_first": "?", "page_last": "?",
            "title": "X-ray structure of human glutamate carboxypeptidase II (GCPII) in complex with a urea based inhibitor PSMA 1027",
            "year": "?", "database_id_CSD": "?", "pdbx_database_id_DOI": "?",
            "pdbx_database_id_PubMed": "?", "unpublished_flag": "?"
        }])
        self.assertEqual(d["atom_site_anisotrop"][6232], {
            "id": "6233", "type_symbol": "O", "pdbx_label_atom_id": "OAM",
            "pdbx_label_alt_id": ".", "pdbx_label_comp_id": "9LZ", "pdbx_label_asym_id": "M",
            "pdbx_label_seq_id": "0", "pdbx_PDB_ins_code": "?", "U[1][1]": "1.7614",
            "U[2][2]": "1.7716", "U[3][3]": "1.852", "U[1][2]": "-0.0093", "U[1][3]": "0.0111",
            "U[2][3]": "-0.0342", "pdbx_auth_seq_id": "824", "pdbx_auth_comp_id": "9LZ",
            "pdbx_auth_asym_id": "A", "pdbx_auth_atom_id": "OAM"
        })
    

    def test_5qu6(self):
        # Unicode characters
        d = atomium.open("tests/integration/files/5qu6.bcif", dictionary=True)
        self.assertEqual(len(d.keys()), 65)
        self.assertEqual(d["entry"], [{"id": "5QU6"}])
        self.assertEqual(d["citation"], [{
            "id": "primary", "journal_abbrev": "J.Biol.Chem.",
            "title": "Small molecule AX-024 reduces T cell proliferation independently of CD3ε/Nck1 interaction, which is governed by a domain swap in the Nck1-SH3.1 domain",
            "year": "2020", "journal_volume": "295", "page_first": "7849",
            "page_last": "7864", "journal_id_ASTM": "JBCHA3", "country": "US",
            "journal_id_ISSN": "1083-351X", "journal_id_CSD": "71", "book_publisher": "?",
            "pdbx_database_id_PubMed": "32317279", "pdbx_database_id_DOI": "10.1074/jbc.RA120.012788"
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
        original = atomium.open(f"tests/integration/files/{code}.bcif", dictionary=True)
        atomium.save_dictionary(original, f"tests/integration/files/output/{code}.bcif")
        saved = atomium.open(f"tests/integration/files/output/{code}.bcif", dictionary=True)
        self.assertEqual(original.keys(), saved.keys())
        for key in original:
            self.assertEqual(original[key], saved[key])
        self.assertEqual(original, saved)
    

    def test_1lol(self):
        self.save("1lol")
    

    def test_5xme(self):
        self.save("5xme")
    

    def test_5o5u(self):
        # Floating point errors
        self.save("5o5u")
    

    def test_5qu6(self):
        # Unicode characters
        self.save("5qu6")