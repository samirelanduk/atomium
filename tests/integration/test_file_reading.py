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
         "description", "experiment", "quality", "geometry"
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
        self.assertEqual(d["geometry"], {"assemblies": [{
         "id": 1,
         "software": "PISA",
         "delta_energy": -31.0,
         "buried_surface_area": 5230,
         "surface_area": 16550,
         "transformations": [{
          "chains": ["A", "B", "C", "D", "E", "F", "G", "H"],
          "matrix": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
          "vector": [0.0, 0.0, 0.0]
         }]
        }]})


    def test_5xme_data_dict(self):
        d = atomium.open("tests/integration/files/5xme.cif", data_dict=True)
        self.assertEqual(d["experiment"]["technique"], "SOLUTION NMR")
        self.assertEqual(d["quality"], {
         "resolution": None, "rvalue": None, "rfree": None
        })


    def test_1xda_data_dict(self):
        d = atomium.open("tests/integration/files/1xda.cif", data_dict=True)
        self.assertEqual(len(d["geometry"]["assemblies"]), 12)
        self.assertEqual(d["geometry"]["assemblies"][0], {
         "id": 1,
         "software": "PISA",
         "delta_energy": -7.0,
         "buried_surface_area": 1720.0,
         "surface_area": 3980.0,
         "transformations": [{
          "chains": ["A", "B", "I", "J", "K", "L", "Y", "Z"],
          "matrix": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
          "vector": [0.0, 0.0, 0.0]
         }]
        })
        self.assertEqual(d["geometry"]["assemblies"][4], {
         "id": 5,
         "software": "PISA",
         "delta_energy": -332.0,
         "buried_surface_area": 21680.0,
         "surface_area": 12240.0,
         "transformations": [{
          "chains": ["E", "F", "G", "H", "Q", "R", "S", "T", "U", "V", "W", "X", "CA", "DA", "EA", "FA"],
          "matrix": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
          "vector": [0.0, 0.0, 0.0]
         }, {
          "chains": ["E", "F", "G", "H", "Q", "R", "S", "T", "U", "V", "W", "X", "CA", "DA", "EA", "FA"],
          "matrix": [[-0.5, -0.8660254038, 0.0], [0.8660254038, -0.5, 0.0], [0.0, 0.0, 1.0]],
          "vector": [0.0, 0.0, 0.0]
         }, {
          "chains": ["E", "F", "G", "H", "Q", "R", "S", "T", "U", "V", "W", "X", "CA", "DA", "EA", "FA"],
          "matrix": [[-0.5, 0.8660254038, 0.0], [-0.8660254038, -0.5, 0.0], [0.0, 0.0, 1.0]],
          "vector": [0.0, 0.0, 0.0]
         }]
        })


    def test_1m4x_data_dict(self):
        d = atomium.open("tests/integration/files/1m4x.cif", data_dict=True)
        self.assertEqual(len(d["geometry"]["assemblies"]), 7)
        self.assertEqual(len(d["geometry"]["assemblies"][0]["transformations"]), 1680)
        self.assertEqual(len(d["geometry"]["assemblies"][2]["transformations"]), 140)
        self.assertEqual(len(d["geometry"]["assemblies"][3]["transformations"]), 168)
        self.assertEqual(len(d["geometry"]["assemblies"][4]["transformations"]), 30)
        self.assertEqual(len(d["geometry"]["assemblies"][5]["transformations"]), 66)
        self.assertEqual(
         d["geometry"]["assemblies"][0]["transformations"][29]["chains"],
         ["A", "B", "C"]
        )
        self.assertAlmostEqual(
         d["geometry"]["assemblies"][0]["transformations"][29]["vector"][0],
         -18.95, delta=0.005
        )
        self.assertAlmostEqual(
         d["geometry"]["assemblies"][0]["transformations"][29]["matrix"][0][0],
         0.812, delta=0.005
        )
        self.assertAlmostEqual(
         d["geometry"]["assemblies"][0]["transformations"][29]["matrix"][-1][-1],
         0.286, delta=0.005
        )



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


    def test_1lol_data_dict(self):
        d = atomium.open("tests/integration/files/1lol.mmtf", data_dict=True)
        self.assertEqual(set(d.keys()), {
         "description", "experiment", "quality", "geometry"
        })
        self.assertEqual(d["description"], {
         "code": "1LOL",
         "title": "Crystal structure of orotidine monophosphate decarboxylase complex with XMP",
         "deposition_date": date(2002, 5, 6),
         "classification": None,
         "keywords": [],
         "authors": []
        })
        self.assertEqual(d["experiment"], {
         "technique": "X-RAY DIFFRACTION",
         "source_organism": None,
         "expression_system": None
        })
        self.assertEqual(d["quality"], {
         "resolution": 1.9, "rvalue": 0.193, "rfree": 0.229
        })
        self.assertEqual(d["geometry"], {"assemblies": [{
         "id": 1, "software": None, "delta_energy": None,
         "buried_surface_area": None, "surface_area": None,
         "transformations": [{
          "chains": ["A", "B", "C", "D", "E", "F", "G", "H"],
          "matrix": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
          "vector": [0.0, 0.0, 0.0]
         }]
        }]})


    def test_1m4x_data_dict(self):
        d = atomium.open("tests/integration/files/1m4x.mmtf", data_dict=True)
        self.assertEqual(len(d["geometry"]["assemblies"]), 7)
        self.assertEqual(len(d["geometry"]["assemblies"][0]["transformations"]), 1680)
        self.assertEqual(len(d["geometry"]["assemblies"][2]["transformations"]), 140)
        self.assertEqual(len(d["geometry"]["assemblies"][3]["transformations"]), 168)
        self.assertEqual(len(d["geometry"]["assemblies"][4]["transformations"]), 30)
        self.assertEqual(len(d["geometry"]["assemblies"][5]["transformations"]), 66)
        self.assertEqual(d["geometry"]["assemblies"][0]["transformations"][29], {
         "chains": ["A", "B", "C"],
         "matrix": [
          [0.81187091297769, 0.33948028980191997, 0.47499304403271003],
          [0.28918995207974996, 0.47292247575946, -0.83229391203448],
          [-0.5071820391020601, 0.81307881404246, 0.28577895606718995]
         ], "vector": [-18.950923036069995, -34.41379177212, 47.03986873605]
        })



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


    def test_1lol_data_dict(self):
        d = atomium.open("tests/integration/files/1lol.pdb", data_dict=True)
        self.assertEqual(set(d.keys()), {
         "description", "experiment", "quality", "geometry"
        })
        self.assertEqual(d["description"], {
         "code": "1LOL",
         "title": "CRYSTAL STRUCTURE OF OROTIDINE MONOPHOSPHATE DECARBOXYLASE COMPLEX WITH XMP",
         "deposition_date": date(2002, 5, 6),
         "classification": "LYASE",
         "keywords": ["TIM BARREL", "LYASE"],
         "authors": ["N.WU", "E.F.PAI"]
        })
        self.assertEqual(d["experiment"], {
         "technique": "X-RAY DIFFRACTION",
         "source_organism": "METHANOTHERMOBACTER THERMAUTOTROPHICUS STR. DELTA H",
         "expression_system": "ESCHERICHIA COLI"
        })
        self.assertEqual(d["quality"], {
         "resolution": 1.9, "rvalue": 0.193, "rfree": 0.229
        })
        self.assertEqual(d["geometry"], {"assemblies": [{
         "id": 1,
         "software": "PISA",
         "delta_energy": -31.0,
         "buried_surface_area": 5230,
         "surface_area": 16550,
         "transformations": [{
          "chains": ["A", "B"],
          "matrix": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
          "vector": [0.0, 0.0, 0.0]
         }]
        }]})


    def test_1xda_data_dict(self):
        d = atomium.open("tests/integration/files/1xda.pdb", data_dict=True)
        self.assertEqual(len(d["geometry"]["assemblies"]), 12)
        self.assertEqual(d["geometry"]["assemblies"][0], {
         "id": 1,
         "software": "PISA",
         "delta_energy": -7.0,
         "buried_surface_area": 1720.0,
         "surface_area": 3980.0,
         "transformations": [{
          "chains": ["A", "B"],
          "matrix": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
          "vector": [0.0, 0.0, 0.0]
         }]
        })
        self.assertEqual(d["geometry"]["assemblies"][4], {
         "id": 5,
         "software": "PISA",
         "delta_energy": -332.0,
         "buried_surface_area": 21680.0,
         "surface_area": 12240.0,
         "transformations": [{
          "chains": ["E", "F", "G", "H"],
          "matrix": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
          "vector": [0.0, 0.0, 0.0]
         }, {
          "chains": ["E", "F", "G", "H"],
          "matrix": [[-0.5, -0.866025, 0.0], [0.866025, -0.5, 0.0], [0.0, 0.0, 1.0]],
          "vector": [0.0, 0.0, 0.0]
         }, {
          "chains": ["E", "F", "G", "H"],
          "matrix": [[-0.5, 0.866025, 0.0], [-0.866025, -0.5, 0.0], [0.0, 0.0, 1.0]],
          "vector": [0.0, 0.0, 0.0]
         }]
        })


    def test_1m4x_data_dict(self):
        self.maxDiff=None
        d = atomium.open("tests/integration/files/1m4x.pdb", data_dict=True)
        self.assertEqual(len(d["geometry"]["assemblies"]), 1)
        self.assertEqual(len(d["geometry"]["assemblies"][0]["transformations"]), 1680)
        self.assertEqual(d["geometry"]["assemblies"][0]["transformations"][29], {
         "chains": ["A", "B", "C"],
         "matrix": [
          [0.811871,  0.339480,  0.474993],
          [0.289190,  0.472922, -0.832294],
          [-0.507182,  0.813079,  0.285779]
         ], "vector": [-18.95092, -34.41379, 47.03987]
        })
