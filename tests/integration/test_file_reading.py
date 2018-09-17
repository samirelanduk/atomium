from datetime import date
import atomium
from unittest import TestCase

class MmcifReadingTests(TestCase):

    def test_1lol_file_dict(self):
        d = atomium.open("tests/integration/files/1lol.cif", file_dict=True)
        self.assertEqual(d["entry"], [{"id": "1LOL"}])
        self.assertEqual(d["audit_author"], [
         {"name": "Wu, N.", "pdbx_ordinal": "1"},
         {"name": "Pai, E.F.", "pdbx_ordinal": "2"}
        ])
        self.assertEqual(d["struct"][0]["title"][:7], "Crystal")
        entity = d["entity"]
        self.assertEqual(len(entity), 4)
        self.assertEqual(
         entity[0]["pdbx_description"], "orotidine 5'-monophosphate decarboxylase"
        )
        self.assertTrue(
         d["entity_poly"][0]["pdbx_seq_one_letter_code"].startswith("LRSRRVDVM")
        )
        self.assertTrue(
         d["entity_poly"][0]["pdbx_seq_one_letter_code"].endswith("IESIKDLLIPE")
        )
        self.assertEqual(entity[1]["type"], "non-polymer")
        self.assertTrue(d["citation"][0]["title"].startswith("Crystal"))
        self.assertTrue(d["citation"][0]["title"].endswith("decarboxylase."))


    def test_1lol_data_dict(self):
        d = atomium.open("tests/integration/files/1lol.cif", data_dict=True)
        self.assertEqual(set(d.keys()), {
         "description", "experiment", "quality"
        })
        self.assertEqual(d["description"], {
         "code": "1LOL",
         "title": "Crystal structure of orotidine monophosphate decarboxylase complex with XMP",
         "deposition_date": date(2002, 5, 6),
         "classification": "LYASE",
         "keywords": ["TIM barrel", "LYASE"],
         "authors": ["Wu, N.", "Pai, E.F."]
        })
        self.assertEqual(d["experiment"], {
         "technique": "X-RAY DIFFRACTION",
         "source_organism": "Methanothermobacter thermautotrophicus str. Delta H",
         "expression_system": "Escherichia coli"
        })
        self.assertEqual(d["quality"], {
         "resolution": 1.9, "rvalue": 0.193, "rfree": 0.229
        })


    def test_5xme_data_dict(self):
        d = atomium.open("tests/integration/files/5xme.cif", data_dict=True)
        self.assertEqual(d["experiment"]["technique"], "SOLUTION NMR")
        self.assertEqual(d["quality"], {
         "resolution": None, "rvalue": None, "rfree": None
        })



class MmtfReadingTests(TestCase):

    def test_1lol_file_dict(self):
        d = atomium.open("tests/integration/files/1lol.mmtf", file_dict=True)
        self.assertEqual(d["mmtfVersion"], "1.0.0")
        self.assertEqual(len(d["unitCell"]), 6)
        self.assertAlmostEqual(d["unitCell"][0], 57.57, delta=0.00005)
        self.assertAlmostEqual(d["resolution"], 1.9, delta=0.00005)
        self.assertEqual(d["numAtoms"], 3431)
        self.assertEqual(len(d["secStructList"]), 602)
        self.assertEqual(d["secStructList"][:5], (7, 4, 4, 4, 3))
        self.assertEqual(len(d["bondAtomList"]), 828)
        self.assertEqual(d["bondAtomList"][:6], (7, 2, 15, 9, 23, 17))
        self.assertEqual(d["chainIdList"], list("ABCDEFGH"))
        self.assertEqual(d["insCodeList"], [""] * 602)
        self.assertEqual(d["sequenceIndexList"][:6], [10, 11, 12, 13, 14, 15])
        self.assertEqual(d["occupancyList"], [1.0] * 3431)
        self.assertEqual(d["xCoordList"][:3], [3.696, 3.198, 3.914])
        self.assertEqual(d["groupList"][0]["groupName"], "ASN")
        self.assertEqual(d["groupList"][0]["atomNameList"][:3], ["N", "CA", "C"])


    def test_1igt_file_dict(self):
        d = atomium.open("tests/integration/files/1igt.mmtf", file_dict=True)
        self.assertEqual(d["mmtfVersion"], "1.0.0")
        self.assertEqual(d["insCodeList"][266], "A")



class PdbReadingTests(TestCase):

    def test_1lol_file_dict(self):
        d = atomium.open("tests/integration/files/1lol.pdb", file_dict=True)
        self.assertEqual(d["HEADER"], [
         "HEADER    LYASE                                   06-MAY-02   1LOL"
        ])
        self.assertEqual(d["HETNAM"], [
         "HETNAM     BU2 1,3-BUTANEDIOL", "HETNAM     XMP XANTHOSINE-5'-MONOPHOSPHATE"
        ])
        self.assertEqual(d["CONECT"][0], "CONECT 3194 3195 3196")
        self.assertEqual(len(d["REMARK"].keys()), 16)
        self.assertEqual(d["REMARK"]["2"][1], "REMARK   2 RESOLUTION.    1.90 ANGSTROMS.")
        self.assertEqual(len(d["MODEL"]), 1)
        self.assertEqual(len(d["MODEL"][0]), 3433)
        self.assertEqual(
         d["MODEL"][0][0],
         "ATOM      1  N   VAL A  11       3.696  33.898  63.219  1.00 21.50           N"
        )


    def test_5xme_file_dict(self):
        d = atomium.open("tests/integration/files/5xme.pdb", file_dict=True)
        self.assertEqual(d["HEADER"], [
         "HEADER    APOPTOSIS                               15-MAY-17   5XME"
        ])
        self.assertEqual(len(d["MODEL"]), 10)
        self.assertEqual(len(d["MODEL"][0]), 1828)
        self.assertEqual(
         d["MODEL"][1][4],
         "ATOM      5  CB  ALA A 199      36.093  -8.556  -1.452  1.00  0.00           C"
        )
