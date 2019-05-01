from datetime import date
import atomium
from unittest import TestCase

class MmcifFileDictReadingTests(TestCase):

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



class MmtfFileDictReadingTests(TestCase):

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
        self.assertEqual(d["bFactorList"][:3], [21.5, 19.76, 19.29])
        self.assertEqual(d["groupList"][0]["groupName"], "ASN")
        self.assertEqual(d["groupList"][0]["atomNameList"][:3], ["N", "CA", "C"])


    def test_1igt_file_dict(self):
        d = atomium.open("tests/integration/files/1igt.mmtf", file_dict=True)
        self.assertEqual(d["mmtfVersion"], "1.0.0")
        self.assertEqual(d["insCodeList"][266], "A")



class PdbFileDictReadingTests(TestCase):

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



class DataDictReadingTests(TestCase):

    def test_1lol_data_dict(self):
        for e in ["cif", "mmtf", "pdb"]:
            d = atomium.open("tests/integration/files/1lol." + e, data_dict=True)
            self.assertEqual(set(d.keys()), {
             "description", "experiment", "quality", "geometry", "models"
            })
            if e == "pdb":
                self.assertEqual(d["description"], {
                 "code": "1LOL",
                 "title": "CRYSTAL STRUCTURE OF OROTIDINE MONOPHOSPHATE DECARBOXYLASE COMPLEX WITH XMP",
                 "deposition_date": date(2002, 5, 6),
                 "classification": "LYASE",
                 "keywords": ["TIM BARREL", "LYASE"],
                 "authors": ["N.WU", "E.F.PAI"]
                })
            else:
                self.assertEqual(d["description"], {
                 "code": "1LOL",
                 "title": "Crystal structure of orotidine monophosphate decarboxylase complex with XMP",
                 "deposition_date": date(2002, 5, 6),
                 "classification": None if e == "mmtf" else "LYASE",
                 "keywords": [] if e == "mmtf" else ["TIM barrel", "LYASE"],
                 "authors": [] if e == "mmtf" else ["Wu, N.", "Pai, E.F."]
                })
            missing_residues = [{"id": id, "name": name} for id, name in zip([
             "A.1", "A.2", "A.3", "A.4", "A.5", "A.6", "A.7", "A.8", "A.9", "A.10",
             "A.182", "A.183", "A.184", "A.185", "A.186", "A.187", "A.188", "A.189",
             "A.223", "A.224", "A.225", "A.226", "A.227", "A.228", "A.229", "B.1001",
             "B.1002", "B.1003", "B.1004", "B.1005", "B.1006", "B.1007", "B.1008",
             "B.1009", "B.1010", "B.1182", "B.1183", "B.1184", "B.1185", "B.1186"
            ], [
             "LEU", "ARG", "SER", "ARG", "ARG", "VAL", "ASP", "VAL", "MET", "ASP",
             "VAL", "GLY", "ALA", "GLN", "GLY", "GLY", "ASP", "PRO", "LYS", "ASP",
             "LEU", "LEU", "ILE", "PRO", "GLU", "LEU", "ARG", "SER", "ARG", "ARG",
             "VAL", "ASP", "VAL", "MET", "ASP", "VAL", "GLY", "ALA", "GLN", "GLY"
            ])]
            if e == "pdb":
                self.assertEqual(d["experiment"], {
                 "technique": "X-RAY DIFFRACTION",
                 "source_organism": "METHANOTHERMOBACTER THERMAUTOTROPHICUS STR. DELTA H",
                 "expression_system": "ESCHERICHIA COLI",
                 "missing_residues": missing_residues
                })
            else:
                self.assertEqual(d["experiment"], {
                 "technique": "X-RAY DIFFRACTION",
                 "source_organism": None if e == "mmtf" else "Methanothermobacter thermautotrophicus str. Delta H",
                 "expression_system": None if e == "mmtf" else "Escherichia coli",
                 "missing_residues": missing_residues if e == "cif" else []
                })
            self.assertEqual(d["quality"], {
             "resolution": 1.9, "rvalue": 0.193, "rfree": 0.229
            })
            self.assertEqual(d["geometry"], {"assemblies": [{
             "id": 1,
             "software": None if e == "mmtf" else "PISA",
             "delta_energy": None if e == "mmtf" else -31.0,
             "buried_surface_area": None if e == "mmtf" else 5230,
             "surface_area": None if e == "mmtf" else 16550,
             "transformations": [{
              "chains": ["A", "B"] if e == "pdb" else ["A", "B", "C", "D", "E", "F", "G", "H"],
              "matrix": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
              "vector": [0.0, 0.0, 0.0]
             }]
            }], "crystallography": {
             "space_group": "P 1 21 1", "unit_cell": [
              57.57, 55.482, 66.129, 90, 94.28, 90
             ]
            }})
            self.assertEqual(len(d["models"]), 1)
            self.assertEqual(len(d["models"][0]["polymer"]), 2)
            self.assertEqual(len(d["models"][0]["polymer"]["A"]["sequence"]), 229)
            self.assertEqual(len(d["models"][0]["polymer"]["A"]["residues"]), 204)
            self.assertEqual(d["models"][0]["polymer"]["A"]["sequence"][:2], "LR")
            self.assertEqual(d["models"][0]["polymer"]["A"]["residues"]["A.11"]['number'], 1)
            self.assertEqual(d["models"][0]["polymer"]["A"]["residues"]["A.13"]['number'], 3)
            self.assertEqual(d["models"][0]["polymer"]["B"]["residues"]["B.1011"]['number'], 1)
            self.assertEqual(d["models"][0]["polymer"]["B"]["residues"]["B.1013"]['number'], 3)
            self.assertEqual(len(d["models"][0]["polymer"]["A"]["residues"]["A.11"]["atoms"]), 7)
            self.assertEqual(d["models"][0]["polymer"]["A"]["residues"]["A.11"]["atoms"][1], {
             "element": "N", "name": "N", "x": 3.696, "y": 33.898, "z": 63.219,
             "bvalue": 21.5, "charge": 0.0, "occupancy": 1.0, "alt_loc": None,
             "anisotropy": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            })
            self.assertEqual(d["models"][0]["polymer"]["A"]["residues"]["A.11"]["name"], "VAL")
            self.assertEqual(d["models"][0]["polymer"]["B"]["residues"]["B.1011"]["atoms"][1558 + (e == "pdb")], {
             "element": "N", "name": "N", "x": -26.384, "y": 61.433, "z": 36.898,
             "bvalue": 39.3, "charge": 0.0, "occupancy": 1.0, "alt_loc": None,
             "anisotropy": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            })
            self.assertEqual(len(d["models"][0]["non-polymer"]), 4)
            self.assertEqual(d["models"][0]["non-polymer"]["A.5001"]["name"], "BU2")
            self.assertEqual(d["models"][0]["non-polymer"]["A.5001"]["internal_id"], "A" if e == "pdb" else "C")
            self.assertEqual(len(d["models"][0]["non-polymer"]["A.5001"]["atoms"]), 6)
            self.assertEqual(d["models"][0]["non-polymer"]["A.5001"]["polymer"], "A")
            self.assertEqual(d["models"][0]["non-polymer"]["B.2002"]["name"], "XMP")
            self.assertEqual(d["models"][0]["non-polymer"]["B.2002"]["internal_id"], "B" if e == "pdb" else "F")
            self.assertEqual(len(d["models"][0]["non-polymer"]["B.2002"]["atoms"]), 24)
            self.assertEqual(d["models"][0]["non-polymer"]["B.2002"]["polymer"], "B")
            self.assertEqual(len(d["models"][0]["water"]), 180)
            self.assertEqual(d["models"][0]["water"]["A.3005"]["name"], "HOH")
            self.assertEqual(d["models"][0]["water"]["A.3005"]["internal_id"], "A" if e == "pdb" else "G")
            self.assertEqual(d["models"][0]["water"]["A.3005"]["polymer"], "A")
            self.assertEqual(d["models"][0]["water"]["B.3180"]["name"], "HOH")
            self.assertEqual(d["models"][0]["water"]["B.3180"]["internal_id"], "B" if e == "pdb" else "H")
            self.assertEqual(d["models"][0]["water"]["B.3180"]["polymer"], "B")


    def test_1xda_data_dict(self):
        for e in ["cif", "mmtf", "pdb"]:
            d = atomium.open("tests/integration/files/1xda." + e, data_dict=True)
            self.assertEqual(len(d["geometry"]["assemblies"]), 12)
            self.assertEqual(d["geometry"]["assemblies"][0], {
             "id": 1,
             "software": None if e == "mmtf" else "PISA",
             "delta_energy": None if e == "mmtf" else -7.0,
             "buried_surface_area": None if e == "mmtf" else 1720.0,
             "surface_area": None if e == "mmtf" else 3980.0,
             "transformations": [{
              "chains": ["A", "B"] if e == "pdb" else ["A", "B", "I", "J", "K", "L", "Y", "Z"],
              "matrix": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
              "vector": [0.0, 0.0, 0.0]
             }]
            })
            if e == "pdb":
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
            else:
                self.assertEqual(d["geometry"]["assemblies"][4], {
                 "id": 5,
                 "software": None if e == "mmtf" else "PISA",
                 "delta_energy": None if e == "mmtf" else -332.0,
                 "buried_surface_area": None if e == "mmtf" else 21680.0,
                 "surface_area": None if e == "mmtf" else 12240.0,
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
        for e in ["cif", "mmtf", "pdb"]:
            d = atomium.open("tests/integration/files/1m4x." + e, data_dict=True)
            self.assertEqual(len(d["geometry"]["assemblies"]), 1 if e == "pdb" else 7)
            self.assertEqual(len(d["geometry"]["assemblies"][0]["transformations"]), 1680)
            if e != "pdb":
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


    def test_4opj_data_dict(self):
        for e in ["cif", "mmtf", "pdb"]:
            d = atomium.open("tests/integration/files/4opj." + e, data_dict=True)
            self.assertEqual(len(d["geometry"]["assemblies"]), 2)
            self.assertEqual(d["geometry"]["assemblies"][0]["transformations"][0]["vector"], [42.387, 0, 0])


    def test_5xme_data_dict(self):
        for e in ["cif", "mmtf", "pdb"]:
            d = atomium.open("tests/integration/files/5xme." + e, data_dict=True)
            self.assertEqual(d["experiment"]["technique"], "SOLUTION NMR")
            self.assertEqual(d["quality"], {
             "resolution": None, "rvalue": None, "rfree": None
            })
            self.assertEqual(d["geometry"]["crystallography"], {
             "space_group": "P 1", "unit_cell": [1, 1, 1, 90, 90, 90]
            } if e == "pdb" else {})
            self.assertEqual(len(d["models"]), 10)
            for model in d["models"][1:]:
                self.assertEqual(len(model["polymer"]), len(d["models"][0]["polymer"]))
                self.assertEqual(len(model["non-polymer"]), len(d["models"][0]["non-polymer"]))
                self.assertEqual(len(model["water"]), len(d["models"][0]["water"]))
            self.assertEqual(d["models"][0]["polymer"]["A"]["residues"]["A.199"]["atoms"][1]["x"], 33.969)
            self.assertEqual(d["models"][1]["polymer"]["A"]["residues"]["A.199"]["atoms"][1 if e == "pdb" else 1828]["x"], 34.064)


    def test_1msh_data_dict(self):
        for e in ["cif", "mmtf", "pdb"]:
            d = atomium.open("tests/integration/files/1msh." + e, data_dict=True)
        self.assertEqual(len(d["models"]), 30)
        for m in d["models"][:-1]:
            self.assertEqual(len(m["polymer"]), 2)
        self.assertEqual(len(d["models"][-1]["polymer"]), 1)


    def test_4y60_data_dict(self):
        for e in ["cif", "mmtf", "pdb"]:
            d = atomium.open("tests/integration/files/4y60.cif", data_dict=True)
            self.assertEqual(d["models"][0]["polymer"]["A"]["sequence"][:2], "CA")
            if e != "mmtf":
                self.assertEqual(d["models"][0]["polymer"]["C"]["residues"]["C.0"]["atoms"][1], {
                 "element": "N", "name": "N", "x": 43.447, "y": -56.622, "z": -20.561,
                 "bvalue": 56.53, "charge": 0.0, "occupancy": 1.0, "alt_loc": None,
                 "anisotropy": [0.9838, 0.7489, 0.4152, -0.1159, -0.0115, -0.2655],
                })


    def test_1cbn_data_dict(self):
        for e in ["cif", "mmtf", "pdb"]:
            d = atomium.open("tests/integration/files/1cbn.cif", data_dict=True)
            self.assertEqual(d["models"][0]["polymer"]["A"]["residues"]["A.1"]["atoms"][1], {
             "element": "N", "name": "N", "x": 16.864, "y": 14.059, "z": 3.442,
             "bvalue": 6.22, "charge": 0.0, "occupancy": 0.8, "alt_loc": "A",
             "anisotropy": [0, 0, 0, 0, 0, 0],
            })



class FileReadingTests(TestCase):

    def test_1lol(self):
        for e in ["cif", "mmtf", "pdb"]:
            f = atomium.open("tests/integration/files/1lol." + e)
            self.assertEqual(f.filetype, e)
            self.assertEqual(f.code, "1LOL")
            if e == "pdb":
                self.assertEqual(
                 f.title, "CRYSTAL STRUCTURE OF OROTIDINE MONOPHOSPHATE DECARBOXYLASE COMPLEX WITH XMP"
                )
            else:
                self.assertEqual(
                 f.title, "Crystal structure of orotidine monophosphate decarboxylase complex with XMP"
                )
            self.assertEqual(f.deposition_date, date(2002, 5, 6))
            self.assertEqual(f.classification, None if e == "mmtf" else "LYASE")
            self.assertEqual(f.keywords, [] if e == "mmtf" else ["TIM BARREL", "LYASE"] if e == "pdb" else ["TIM barrel", "LYASE"])
            self.assertEqual(f.authors, [] if e == "mmtf" else ["N.WU", "E.F.PAI"] if e == "pdb" else ["Wu, N.", "Pai, E.F."])
            self.assertEqual(f.technique, "X-RAY DIFFRACTION")
            missing_residues = [{"id": id, "name": name} for id, name in zip([
             "A.1", "A.2", "A.3", "A.4", "A.5", "A.6", "A.7", "A.8", "A.9", "A.10",
             "A.182", "A.183", "A.184", "A.185", "A.186", "A.187", "A.188", "A.189",
             "A.223", "A.224", "A.225", "A.226", "A.227", "A.228", "A.229", "B.1001",
             "B.1002", "B.1003", "B.1004", "B.1005", "B.1006", "B.1007", "B.1008",
             "B.1009", "B.1010", "B.1182", "B.1183", "B.1184", "B.1185", "B.1186"
            ], [
             "LEU", "ARG", "SER", "ARG", "ARG", "VAL", "ASP", "VAL", "MET", "ASP",
             "VAL", "GLY", "ALA", "GLN", "GLY", "GLY", "ASP", "PRO", "LYS", "ASP",
             "LEU", "LEU", "ILE", "PRO", "GLU", "LEU", "ARG", "SER", "ARG", "ARG",
             "VAL", "ASP", "VAL", "MET", "ASP", "VAL", "GLY", "ALA", "GLN", "GLY"
            ])]
            if e == "pdb":
                self.assertEqual(f.source_organism, "METHANOTHERMOBACTER THERMAUTOTROPHICUS STR. DELTA H")
                self.assertEqual(f.expression_system, "ESCHERICHIA COLI")
                self.assertEqual(f.missing_residues, missing_residues)
            else:
                self.assertEqual(
                 f.source_organism,
                 None if e == "mmtf" else "Methanothermobacter thermautotrophicus str. Delta H"
                )
                self.assertEqual(
                 f.expression_system, None if e == "mmtf" else "Escherichia coli"
                )
                self.assertEqual(f.missing_residues, [] if e == "mmtf" else missing_residues)
            self.assertEqual(f.resolution, 1.9)
            self.assertEqual(f.rvalue, 0.193)
            self.assertEqual(f.rfree, 0.229)
            self.assertEqual(f.assemblies, [{
             "id": 1,
             "software": None if e == "mmtf" else "PISA",
             "delta_energy": None if e == "mmtf" else -31.0,
             "buried_surface_area": None if e == "mmtf" else 5230,
             "surface_area": None if e == "mmtf" else 16550,
             "transformations": [{
              "chains": ["A", "B"] if e == "pdb" else ["A", "B", "C", "D", "E", "F", "G", "H"],
              "matrix": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
              "vector": [0.0, 0.0, 0.0]
             }]
            }])

            self.assertEqual(len(f.models), 1)
            model = f.model
            self.assertEqual(len(model.chains()), 2)
            self.assertIsInstance(model.chains(), set)
            self.assertEqual(len(model.ligands()), 4)
            self.assertIsInstance(model.ligands(), set)
            self.assertEqual(len(model.waters()), 180)
            self.assertIsInstance(model.waters(), set)
            self.assertEqual(len(model.molecules()), 186)
            self.assertIsInstance(model.molecules(), set)
            self.assertEqual(len(model.residues()), 418)
            self.assertIsInstance(model.residues(), set)
            self.assertEqual(len(model.atoms()), 3431)
            self.assertIsInstance(model.atoms(), set)
            self.assertEqual(len(model.chains(length__gt=200)), 2)
            self.assertEqual(len(model.chains(length__gt=210)), 1)
            self.assertEqual(len(model.ligands(name="XMP")), 2)
            self.assertEqual(len(model.residues(name="VAL")), 28)
            self.assertEqual(len(model.residues(name="CYS")), 6)
            self.assertEqual(len(model.residues(name__regex="CYS|VAL")), 34)
            self.assertAlmostEqual(
             model.mass, 46018.5, delta=0.005
            )

            chaina = model.chain("A")
            chainb = model.chain(id="B")
            self.assertIs(chaina.model, model)
            self.assertIs(chainb.model, model)
            self.assertEqual(chaina.id, "A")
            self.assertEqual(chaina.length, 204)
            self.assertEqual(chainb.length, 214)
            self.assertTrue(chaina.sequence.startswith("LRSRRVDVMDVMNRLILAMDL"))
            self.assertTrue(chaina.sequence.endswith("LADNPAAAAAGIIESIKDLLIPE"))
            self.assertTrue(chainb.sequence.startswith("LRSRRVDVMDVMNRLILAMDL"))
            self.assertTrue(chainb.sequence.endswith("LADNPAAAAAGIIESIKDLLIPE"))
            for res in chaina: self.assertIn(res, chaina)
            self.assertEqual(len(chaina.residues()), 204)
            self.assertIsInstance(chaina.residues(), tuple)
            self.assertEqual(len(chaina.ligands()), 2)
            self.assertIsInstance(chaina.ligands(), set)
            self.assertEqual(len(chaina.atoms()), 1557)
            self.assertIsInstance(chaina.atoms(), set)
            self.assertEqual(len(chainb.atoms()), 1634)
            self.assertIsInstance(chainb.atoms(), set)
            res = chaina.residue("A.15")
            self.assertIs(res.chain, chaina)
            self.assertIs(res.model, model)
            self.assertEqual(res.name, "LEU")
            self.assertEqual(res.code, "L")
            self.assertEqual(res.full_name, "leucine")
            self.assertEqual(len(res.atoms()), 8)
            self.assertIsInstance(chaina.atoms(), set)
            self.assertEqual(len(res.atoms(element="C")), 6)
            self.assertEqual(len(res.atoms(element__regex="C|O")), 7)
            self.assertEqual(len(res.atoms(name__regex="^CD")), 2)
            self.assertIs(chaina[0], chaina.residue("A.11"))
            self.assertIs(res.next, chaina[5])
            self.assertIn(chaina.residue(name="GLN"), [chaina.residue("A.136"), chaina.residue("A.173")])

            lig = model.ligand(name="XMP")
            self.assertIs(lig.model, model)
            self.assertEqual(len(lig.atoms()), 24)
            self.assertEqual(lig.formula, {"C": 10, "O": 9, "N": 4, "P": 1})
            lig = model.ligand("A.5001")
            self.assertIs(lig.model, model)
            self.assertIs(lig.chain, chaina)
            self.assertEqual(len(lig.atoms()), 6)
            self.assertEqual(lig.mass, 80.0416)
            pairs = list(lig.pairwise_atoms())
            self.assertEqual(len(pairs), 15)
            for pair in pairs:
                pair = list(pair)
                self.assertTrue(0 < pair[0].distance_to(pair[1]), 5)
            hoh = model.water("A.3005")
            self.assertEqual(hoh.name, "HOH")
            self.assertIs(lig.model, model)
            self.assertIs(lig.chain, chaina)
            lig1, lig2 = model.ligands(name="XMP")
            self.assertAlmostEqual(lig1.rmsd_with(lig2), 0.133, delta=0.001)
            self.assertAlmostEqual(lig2.rmsd_with(lig1), 0.133, delta=0.001)

            atom = model.atom(934)
            self.assertEqual(atom.anisotropy, [0, 0, 0, 0, 0, 0])
            self.assertEqual(atom.element, "C")
            self.assertEqual(atom.name, "CA")
            self.assertEqual(atom.location, (4.534, 53.864, 43.326))
            self.assertEqual(atom.bvalue, 17.14)
            self.assertEqual(atom.charge, 0)
            self.assertAlmostEqual(atom.mass, 12, delta=0.1)
            self.assertIs(atom.chain, chaina)
            self.assertIs(atom.model, model)
            self.assertIs(atom.structure, model.residue("A.131"))

            self.assertEqual(model.molecule("A"), chaina)
            self.assertEqual(model.molecule("A.5001"), lig)
            self.assertEqual(model.molecule("A.3005"), hoh)
            self.assertEqual(len(model.molecules(mass__gt=18)), 6)
            self.assertEqual(len(model.molecules(mass__gt=90)), 4)
            self.assertEqual(len(model.molecules(mass__gt=1000)), 2)
            self.assertEqual(len(model.molecules(mass__gt=90, mass__lt=1000)), 2)

            atom = model.atom(1587 if e == "pdb" else 1586)
            four_angstrom = atom.nearby_atoms(cutoff=4)
            self.assertEqual(len(four_angstrom), 10)
            self.assertEqual(
             sorted([atom.id for atom in four_angstrom]),
             [n - (e != "pdb") for n in [1576, 1582, 1583, 1584, 1586, 1588, 1589, 1590, 1591, 2957]]
            )
            self.assertEqual(len(atom.nearby_atoms(cutoff=4, element="O")), 1)
            four_angstrom = model.atoms_in_sphere(atom.location, 4)
            self.assertEqual(len(four_angstrom), 11)
            self.assertEqual(
             sorted([atom.id for atom in four_angstrom]),
             [n - (e != "pdb") for n in [1576, 1582, 1583, 1584, 1586, 1587, 1588, 1589, 1590, 1591, 2957]]
            )
            self.assertEqual(len(model.atoms_in_sphere(atom.location, 4, element="O")), 1)

            atom = model.atom(905)
            self.assertEqual(len(atom.nearby_structures(5)), 8)
            self.assertEqual(len(atom.nearby_structures(5, ligands=False)), 7)
            self.assertEqual(len(atom.nearby_structures(5, waters=True)), 9)
            self.assertEqual(len(atom.nearby_structures(5, residues=False)), 1)
            self.assertEqual(len(atom.nearby_structures(5, element="O")), 3)

            model.dehydrate()
            self.assertEqual(model.waters(), set())


    def test_5xme(self):
        for e in ["cif", "mmtf", "pdb"]:
            f = atomium.open("tests/integration/files/5xme." + e)
            self.assertEqual(f.resolution, None)
            models = f.models
            self.assertEqual(len(models), 10)
            self.assertIs(f.model, f.models[0])
            x_values = [
             33.969, 34.064, 37.369, 36.023, 35.245,
             35.835, 37.525, 35.062, 36.244, 37.677
            ]
            all_atoms = set()
            for x, model in zip(x_values, models):
                self.assertEqual(len(model.atoms()), 1827)
                all_atoms.update(model.atoms())
                atom = model.chain()[0].atom(name="N")
                self.assertEqual(atom.x, x)
            self.assertEqual(len(all_atoms), 18270)


    def test_1cbn(self):
        for e in ["cif", "mmtf", "pdb"]:
            f = atomium.open("tests/integration/files/1cbn." + e)
            chain = f.model.chain()
            residue1, residue2, residue3 = chain[:3]
            self.assertEqual(len(residue1.atoms()), 16)
            self.assertEqual(len(residue2.atoms()), 14)
            self.assertEqual(len(residue3.atoms()), 10)
            for residue in chain[:3]:
                for name in ["N", "C", "CA", "CB"]:
                    self.assertEqual(len(residue.atoms(name=name)), 1)


    def test_1xda(self):
        for e in ["cif"]:#, "mmtf", "pdb"]:
            f = atomium.open("tests/integration/files/1xda." + e)
            self.assertEqual(len(f.model.atoms()), 1842)
            self.assertEqual(len(f.model.atoms(is_metal=True)), 4)
            self.assertEqual(len(f.model.atoms(is_metal=False)), 1838)

            model = f.model
            self.assertEqual(len(model.atoms()), 1842)
            self.assertEqual(len(model.chains()), 8)
            self.assertEqual(len(model.ligands()), 16)

            model = f.generate_assembly(1)
            self.assertEqual(len(model.chains()), 2)
            self.assertEqual(set([c.id for c in model.chains()]), {"A", "B"})
            self.assertEqual(len(model.ligands()), 4)

            model = f.generate_assembly(2)
            self.assertEqual(len(model.chains()), 2)
            self.assertEqual(set([c.id for c in model.chains()]), {"C", "D"})
            self.assertEqual(len(model.ligands()), 4)

            model = f.generate_assembly(3)
            self.assertEqual(len(model.chains()), 2)
            self.assertEqual(set([c.id for c in model.chains()]), {"E", "F"})
            self.assertEqual(len(model.ligands()), 4)

            model = f.generate_assembly(4)
            self.assertEqual(len(model.chains()), 2)
            self.assertEqual(set([c.id for c in model.chains()]), {"G", "H"})
            self.assertEqual(len(model.ligands()), 4)

            model = f.generate_assembly(7)
            self.assertEqual(len(model.chains()), 6)
            self.assertEqual(set([c.id for c in model.chains()]), {"A", "B"})
            self.assertEqual(len(model.ligands()), 12)
            zn = model.atom(element="ZN")
            liganding_residues = zn.nearby_structures(3, is_metal=False, element__ne="CL")
            self.assertEqual(len(liganding_residues), 3)
            self.assertEqual(set([r.id for r in liganding_residues]), {"B.10"})
            self.assertEqual(set([r.name for r in liganding_residues]), {"HIS"})
            res1, res2, res3 = liganding_residues

            self.assertGreater(res1.atom(name="N").distance_to(res2.atom(name="N")), 10)
