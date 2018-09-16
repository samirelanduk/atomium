from unittest import TestCase
import atomium

class MmcifReadingTests(TestCase):

    def test_1lol_file_dict(self):
        d = atomium.open("tests/integration/files/1lol.cif", file_dict=True)
        self.assertEqual(d["entry"], [{"id": "1LOL"}])
        self.assertEqual(d["audit_author"], [
         {"name": "Wu, N.", "pdbx_ordinal": "1"},
         {"name": "Pai, E.F.", "pdbx_ordinal": "2"}
        ])
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
