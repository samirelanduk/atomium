from datetime import date
from unittest import TestCase
import atomium

class DictParsingTests(TestCase):

    def test_1lol_mmtf(self):
        d = atomium.open("tests/integration/files/1lol.mmtf", dictionary=True)
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



class ParsingTests(TestCase):

    def test_1lol(self):
        # Open file
        pdb = atomium.open("tests/integration/files/1lol.mmtf")
        self.assertEqual(pdb.name, "1LOL")
        self.assertEqual(len(pdb.source.keys()), 66)
        self.assertEqual(pdb.source["entry"][0]["id"], "1LOL")
        self.assertEqual(pdb.entry__id, "1LOL")
        missing_residues = [{"id": id, "name": name} for id, name in zip([
            "A.1", "A.2", "A.3", "A.4", "A.5", "A.6", "A.7", "A.8", "A.9", "A.10",
            "A.182", "A.183", "A.184", "A.185", "A.186", "A.187", "A.188", "A.189",
            "A.223", "A.224", "A.225", "A.226", "A.227", "A.228", "A.229", "B.1",
            "B.2", "B.3", "B.4", "B.5", "B.6", "B.7", "B.8",
            "B.9", "B.10", "B.182", "B.183", "B.184", "B.185", "B.186"
        ], [
            "LEU", "ARG", "SER", "ARG", "ARG", "VAL", "ASP", "VAL", "MET", "ASP",
            "VAL", "GLY", "ALA", "GLN", "GLY", "GLY", "ASP", "PRO", "LYS", "ASP",
            "LEU", "LEU", "ILE", "PRO", "GLU", "LEU", "ARG", "SER", "ARG", "ARG",
            "VAL", "ASP", "VAL", "MET", "ASP", "VAL", "GLY", "ALA", "GLN", "GLY"
        ])]
        self.assertEqual(pdb.missing_residues, missing_residues)
        self.assertEqual(str(pdb.model), "<Model (2 polymers, 4 non-polymers)>")
    

    def test_1lol_mmtf_compressed(self):
        d = atomium.open("tests/integration/files/1lol.mmtf.gz", dictionary=True)
        self.assertEqual(len(d.keys()), 66)
        self.assertEqual(d["entry"], [{"id": "1LOL"}])



class NetworkTests(TestCase):

    def test_fetching_3nir_mmtf_rcsb(self):
        pdb = atomium.fetch("3nir.mmtf")
        self.assertEqual(pdb.title, "Crystal structure of small protein crambin at 0.48 A resolution")