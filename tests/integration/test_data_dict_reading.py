from datetime import date
import re
import atomium
from unittest import TestCase

class DataDictTest(TestCase):

    def open(self, name):
        return {e: atomium.open(
         f"tests/integration/files/{name}.{e}", data_dict=True
        ) for e in ("cif", "mmtf", "pdb")}
    

    def check(self, dicts, subdict, values):
        for ext, d in dicts.items():
            d = d[subdict]
            for key, value in values.items():
                # Is this key specific to one filetype that isn't this one?
                if re.search("_(cif|mmtf|pdb)$", key) and key.split("_")[-1] != ext:
                    continue
                # Is there a more specific version of this key?
                if any(k == f"{key}_{ext}" for k in values.keys()):
                    continue
                try:
                    self.assertEqual(d[key.replace(f"_{ext}", "")], value)
                except AssertionError as e:
                    raise AssertionError(f"({ext.upper()}) {str(e)}")



class DescriptionDictTests(DataDictTest):

    def test_1lol_data_dict_description(self):
        data_dicts = self.open("1lol")
        self.check(data_dicts, "description", {
         "code": "1LOL",
         "title": "Crystal structure of orotidine monophosphate decarboxylase complex with XMP",
         "title_pdb": "CRYSTAL STRUCTURE OF OROTIDINE MONOPHOSPHATE DECARBOXYLASE COMPLEX WITH XMP",
         "deposition_date": date(2002, 5, 6),
         "classification": "LYASE", "classification_mmtf": None,
         "keywords_cif": ["TIM barrel", "LYASE"],
         "keywords_mmtf": [],
         "keywords_pdb": ["TIM BARREL", "LYASE"],
         "authors_cif": ["Wu, N.", "Pai, E.F."],
         "authors_mmtf": [],
         "authors_pdb": ["N.WU", "E.F.PAI"]
        })



class ExperimentDictTests(DataDictTest):

    def test_1lol_data_dict_experiment(self):
        data_dicts = self.open("1lol")
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
        self.check(data_dicts, "experiment", {
         "technique": "X-RAY DIFFRACTION",
         "source_organism_cif": "Methanothermobacter thermautotrophicus str. Delta H",
         "source_organism_mmtf": None,
         "source_organism_pdb": "METHANOTHERMOBACTER THERMAUTOTROPHICUS STR. DELTA H",
         "expression_system_cif": "Escherichia coli",
         "expression_system_mmtf": None,
         "expression_system_pdb": "ESCHERICHIA COLI",
         "missing_residues": missing_residues, "missing_residues_mmtf": []
        })
        

    def test_1cbn_data_dict_experiment(self):
        data_dicts = self.open("1cbn")
        self.check(data_dicts, "experiment", {
         "technique": "X-RAY DIFFRACTION",
         "source_organism_cif": "Crambe hispanica subsp. abyssinica",
         "source_organism_mmtf": None,
         "source_organism_pdb": "CRAMBE HISPANICA SUBSP. ABYSSINICA",
         "expression_system": None,
        })



class QualityDictTests(DataDictTest):

    def test_1lol_data_dict_quality(self):
        data_dicts = self.open("1lol")
        self.check(data_dicts, "quality", {
         "resolution": 1.9, "rvalue": 0.193, "rfree": 0.229
        })
    

    def test_5xme_data_dict_quality(self):
        data_dicts = self.open("5xme")
        self.check(data_dicts, "quality", {
         "resolution": None, "rvalue": None, "rfree": None
        })



class GeometryDictTests(DataDictTest):

    def test_1lol_data_dict_geometry(self):
        data_dicts = self.open("1lol")
        self.check(data_dicts, "geometry", {
         "crystallography": {"space_group": "P 1 21 1", "unit_cell": [
          57.57, 55.482, 66.129, 90, 94.28, 90
         ]},
         "assemblies_cif": [{
          "id": 1, "software": "PISA", "delta_energy": -31.0,
          "buried_surface_area": 5230, "surface_area": 16550,
          "transformations": [{
           "chains": ["A", "B", "C", "D", "E", "F", "G", "H"],
           "matrix": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
           "vector": [0.0, 0.0, 0.0]
          }]
         }],
         "assemblies_mmtf": [{
          "id": 1, "software": None, "delta_energy": None,
          "buried_surface_area": None, "surface_area": None,
          "transformations": [{
           "chains": ["A", "B", "C", "D", "E", "F", "G", "H"],
           "matrix": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
           "vector": [0.0, 0.0, 0.0]
          }]
         }],
         "assemblies_pdb": [{
          "id": 1, "software": "PISA", "delta_energy": -31.0,
          "buried_surface_area": 5230, "surface_area": 16550,
          "transformations": [{
           "chains": ["A", "B"],
           "matrix": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
           "vector": [0.0, 0.0, 0.0]
          }]
         }]
        })


    def test_5xme_data_dict_geometry(self):
        data_dicts = self.open("5xme")
        self.check(data_dicts, "geometry", {
         "crystallography": {},
         "crystallography_pdb": {
          "space_group": "P 1", "unit_cell": [1, 1, 1, 90, 90, 90]
         }
        })


    def test_1xda_data_dict_geometry(self):
        # Multiple assemblies with different chains
        data_dicts = self.open("1xda")
        for d in data_dicts.values():
            self.assertEqual(len(d["geometry"]["assemblies"]), 12)
        self.assertEqual(data_dicts["cif"]["geometry"]["assemblies"][0], {
         "id": 1, "software": "PISA", "delta_energy": -7,
         "buried_surface_area": 1720, "surface_area": 3980,
         "transformations": [{
          "chains": ["A", "B", "I", "J", "K", "L", "Y", "Z"],
          "matrix": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
          "vector": [0.0, 0.0, 0.0]
         }]
        })
        self.assertEqual(data_dicts["mmtf"]["geometry"]["assemblies"][0], {
         "id": 1, "software": None, "delta_energy": None,
         "buried_surface_area": None, "surface_area": None,
         "transformations": [{
          "chains": ["A", "B", "I", "J", "K", "L", "Y", "Z"],
          "matrix": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
          "vector": [0.0, 0.0, 0.0]
         }]
        })
        self.assertEqual(data_dicts["pdb"]["geometry"]["assemblies"][0], {
         "id": 1, "software": "PISA", "delta_energy": -7,
         "buried_surface_area": 1720, "surface_area": 3980,
         "transformations": [{
          "chains": ["A", "B"],
          "matrix": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
          "vector": [0.0, 0.0, 0.0]
         }]
        })
        self.assertEqual(data_dicts["cif"]["geometry"]["assemblies"][4], {
         "id": 5, "software": "PISA", "delta_energy": -332.0,
         "buried_surface_area": 21680.0, "surface_area": 12240.0,
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
        self.assertEqual(data_dicts["mmtf"]["geometry"]["assemblies"][4], {
         "id": 5, "software": None, "delta_energy": None,
         "buried_surface_area": None, "surface_area": None,
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
        self.assertEqual(data_dicts["pdb"]["geometry"]["assemblies"][4], {
         "id": 5, "software": "PISA", "delta_energy": -332.0,
         "buried_surface_area": 21680.0, "surface_area": 12240.0,
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


    def test_1m4x_data_dict_geometry(self):
        # Assemblies with lots of transformations to create virus
        data_dicts = self.open("1m4x")
        self.assertEqual(len(data_dicts["cif"]["geometry"]["assemblies"]), 7)
        self.assertEqual(len(data_dicts["mmtf"]["geometry"]["assemblies"]), 7)
        self.assertEqual(len(data_dicts["pdb"]["geometry"]["assemblies"]), 1)
        self.assertEqual(len(data_dicts["cif"]["geometry"]["assemblies"][2]["transformations"]), 140)
        self.assertEqual(len(data_dicts["cif"]["geometry"]["assemblies"][3]["transformations"]), 168)
        self.assertEqual(len(data_dicts["cif"]["geometry"]["assemblies"][4]["transformations"]), 30)
        self.assertEqual(len(data_dicts["cif"]["geometry"]["assemblies"][5]["transformations"]), 66)
        self.assertEqual(len(data_dicts["mmtf"]["geometry"]["assemblies"][2]["transformations"]), 140)
        self.assertEqual(len(data_dicts["mmtf"]["geometry"]["assemblies"][3]["transformations"]), 168)
        self.assertEqual(len(data_dicts["mmtf"]["geometry"]["assemblies"][4]["transformations"]), 30)
        self.assertEqual(len(data_dicts["mmtf"]["geometry"]["assemblies"][5]["transformations"]), 66)
        for d in data_dicts.values():
            self.assertEqual(len(d["geometry"]["assemblies"][0]["transformations"]), 1680)
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
        

    def test_4opj_data_dict_geometry(self):
        # Weird assemblies
        data_dicts = self.open("4opj")
        for d in data_dicts.values():
            self.assertEqual(len(d["geometry"]["assemblies"]), 2)
            self.assertEqual(d["geometry"]["assemblies"][0]["transformations"][0]["vector"], [42.387, 0, 0])



class ModelDictTests(DataDictTest):

    def test_1lol_data_dict_model(self):
        data_dicts = self.open("1lol")
        for e, d in data_dicts.items():
            self.assertEqual(len(d["models"]), 1)
            self.assertEqual(len(d["models"][0]["polymer"]), 2)
            self.assertEqual(len(d["models"][0]["polymer"]["A"]["sequence"]), 229)
            self.assertEqual(len(d["models"][0]["polymer"]["A"]["residues"]), 204)
            self.assertEqual(d["models"][0]["polymer"]["A"]["sequence"][:2], "LR")
            self.assertEqual(d["models"][0]["polymer"]["A"]["residues"]["A.11"]['number'], 1)
            self.assertEqual(d["models"][0]["polymer"]["A"]["residues"]["A.13"]['number'], 3)
            self.assertEqual(d["models"][0]["polymer"]["A"]["residues"]["A.13"]['full_name'], None)
            self.assertEqual(d["models"][0]["polymer"]["B"]["residues"]["B.1011"]['number'], 1)
            self.assertEqual(d["models"][0]["polymer"]["B"]["residues"]["B.1013"]['number'], 3)
            self.assertEqual(len(d["models"][0]["polymer"]["A"]["residues"]["A.11"]["atoms"]), 7)
            self.assertEqual(d["models"][0]["polymer"]["A"]["residues"]["A.11"]["atoms"][1], {
             "element": "N", "name": "N", "x": 3.696, "y": 33.898, "z": 63.219,
             "bvalue": 21.5, "charge": 0.0, "occupancy": 1.0, "alt_loc": None,
             "anisotropy": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0], "is_hetatm": False
            })
            self.assertEqual(d["models"][0]["polymer"]["A"]["residues"]["A.11"]["name"], "VAL")
            self.assertEqual(d["models"][0]["polymer"]["B"]["residues"]["B.1011"]["atoms"][1558 + (e == "pdb")], {
             "element": "N", "name": "N", "x": -26.384, "y": 61.433, "z": 36.898,
             "bvalue": 39.3, "charge": 0.0, "occupancy": 1.0, "alt_loc": None,
             "anisotropy": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0], "is_hetatm": False
            })
            if e == "mmtf":
                self.assertEqual(d["models"][0]["polymer"]["A"]["helices"][0], [f"A.{n}" for n in range(12, 15)])
                self.assertEqual(d["models"][0]["polymer"]["A"]["strands"][0], [f"A.{n}" for n in range(15, 20)])
                self.assertEqual(d["models"][0]["polymer"]["B"]["helices"][0], [f"B.{n}" for n in range(1012, 1015)])
                self.assertEqual(d["models"][0]["polymer"]["B"]["strands"][0], [f"B.{n}" for n in range(1015, 1019)])
            else:
                self.assertEqual(d["models"][0]["polymer"]["A"]["helices"][0], [f"A.{n}" for n in range(11, 14)])
                self.assertEqual(d["models"][0]["polymer"]["A"]["strands"][0], [f"A.{n}" for n in range(15, 20)])
                self.assertEqual(d["models"][0]["polymer"]["B"]["helices"][0], [f"B.{n}" for n in range(1011, 1014)])
                self.assertEqual(d["models"][0]["polymer"]["B"]["strands"][0], [f"B.{n}" for n in range(1015, 1019)])
            self.assertEqual(len(d["models"][0]["non-polymer"]), 4)
            self.assertEqual(d["models"][0]["non-polymer"]["A.5001"]["name"], "BU2")
            self.assertEqual(d["models"][0]["non-polymer"]["A.5001"]["internal_id"], "A" if e == "pdb" else "C")
            self.assertEqual(len(d["models"][0]["non-polymer"]["A.5001"]["atoms"]), 6)
            self.assertEqual(d["models"][0]["non-polymer"]["A.5001"]["polymer"], "A")
            self.assertEqual(d["models"][0]["non-polymer"]["B.2002"]["name"], "XMP")
            self.assertEqual(d["models"][0]["non-polymer"]["B.2002"]["full_name"], "XANTHOSINE-5'-MONOPHOSPHATE")
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


    def test_5xme_data_dict_model(self):
        data_dicts = self.open("5xme")
        for e, d in data_dicts.items():
            self.assertEqual(len(d["models"]), 10)
            for model in d["models"][1:]:
                self.assertEqual(len(model["polymer"]), len(d["models"][0]["polymer"]))
                self.assertEqual(len(model["non-polymer"]), len(d["models"][0]["non-polymer"]))
                self.assertEqual(len(model["water"]), len(d["models"][0]["water"]))
            self.assertEqual(d["models"][0]["polymer"]["A"]["residues"]["A.199"]["atoms"][1]["x"], 33.969)
            self.assertEqual(d["models"][1]["polymer"]["A"]["residues"]["A.199"]["atoms"][1 if e == "pdb" else 1828]["x"], 34.064)


    def test_1msh_data_dict_model(self):
        data_dicts = self.open("1msh")
        for d in data_dicts.values():
            self.assertEqual(len(d["models"]), 30)
            for m in d["models"][:-1]:
                self.assertEqual(len(m["polymer"]), 2)
            self.assertEqual(len(d["models"][-1]["polymer"]), 1)
    

    def test_4y60_data_dict_model(self):
        data_dicts = self.open("4y60")
        for d in data_dicts.values():
            self.assertEqual(d["models"][0]["polymer"]["A"]["sequence"][:2], "CA")
        self.assertEqual(data_dicts["cif"]["models"][0]["polymer"]["C"]["residues"]["C.0"]["atoms"][1], {
         "element": "N", "name": "N", "x": 43.447, "y": -56.622, "z": -20.561,
         "bvalue": 56.53, "charge": 0.0, "occupancy": 1.0, "alt_loc": None,
         "anisotropy": [0.9838, 0.7489, 0.4152, -0.1159, -0.0115, -0.2655], "is_hetatm": False
        })
        self.assertEqual(data_dicts["mmtf"]["models"][0]["polymer"]["C"]["residues"]["C.0"]["atoms"][1], {
         "element": "N", "name": "N", "x": 43.447, "y": -56.622, "z": -20.561,
         "bvalue": 56.53, "charge": 0.0, "occupancy": 1.0, "alt_loc": None,
         "anisotropy": [0, 0, 0, 0, 0, 0], "is_hetatm": False
        })
        self.assertEqual(data_dicts["pdb"]["models"][0]["polymer"]["C"]["residues"]["C.0"]["atoms"][1], {
         "element": "N", "name": "N", "x": 43.447, "y": -56.622, "z": -20.561,
         "bvalue": 56.53, "charge": 0.0, "occupancy": 1.0, "alt_loc": None,
         "anisotropy": [0.9838, 0.7489, 0.4152, -0.1159, -0.0115, -0.2655], "is_hetatm": False
        })
    

    def test_1cbn_data_dict_model(self):
        data_dicts = self.open("1cbn")
        for e, d in data_dicts.items():
            if e == "mmtf":
                self.assertEqual(d["models"][0]["polymer"]["A"]["helices"], [
                 [f"A.{n}" for n in range(7, 18)],
                 [f"A.{n}" for n in range(23, 31)],
                ])
                self.assertEqual(d["models"][0]["polymer"]["A"]["strands"], [
                 ['A.2', 'A.3'], ['A.33', 'A.34']
                ])
            else:
                self.assertEqual(d["models"][0]["polymer"]["A"]["helices"], [
                 [f"A.{n}" for n in range(7, 20)],
                 [f"A.{n}" for n in range(23, 31)],
                ])
                self.assertEqual(d["models"][0]["polymer"]["A"]["strands"], [
                 [f"A.{n}" for n in range(32, 36)],
                 [f"A.{n}" for n in range(1, 5)],
                 ["A.46"],
                 [f"A.{n}" for n in range(39, 42)],
                ])
            self.assertEqual(d["models"][0]["polymer"]["A"]["residues"]["A.1"]["atoms"][1], {
             "element": "N", "name": "N", "x": 16.864, "y": 14.059, "z": 3.442,
             "bvalue": 6.22, "charge": 0.0, "occupancy": 0.8, "alt_loc": "A",
             "anisotropy": [0, 0, 0, 0, 0, 0], "is_hetatm": False
            })
    

    def test_4opj_data_dict_model(self):
        data_dicts = self.open("4opj")
        for e ,d in data_dicts.items():
            if e == "cif":
                self.assertEqual(
                 d["models"][0]["polymer"]["B"]["residues"]["B.6"]["full_name"],
                 "(2R,3aS,4aR,5aR,5bS)-2-(6-amino-9H-purin-9-yl)-3a-hydroxyhexahydrocyclopropa[4,5]cyclopenta[1,2-b]furan-5a(4H)-yl dihydrogen phosphate"
                )
            elif e =="mmtf":
                self.assertEqual(
                 d["models"][0]["polymer"]["B"]["residues"]["B.6"]["full_name"],
                 None
                )
            else:
                self.assertEqual(
                 d["models"][0]["polymer"]["B"]["residues"]["B.6"]["full_name"],
                 "(2R,3AS,4AR,5AR,5BS)-2-(6-AMINO-9H-PURIN-9-YL)-3A-HYDROXYHEXAHYDROCYCLOPROPA[4,5]CYCLOPENTA[1,2-B]FURAN-5A(4H)-YL DIHYDROGEN PHOSPHATE"
                )
    

    def test_4pgp_data_dict_model(self):
        data_dicts = self.open("4gpg")
        for d in data_dicts.values():
            self.assertEqual(len(d["models"]), 1)
            self.assertEqual(len(d["models"][0]["polymer"]), 1)
            self.assertEqual(len(d["models"][0]["non-polymer"]), 0)
            self.assertEqual(len(d["models"][0]["water"]), 140)
            for water in d["models"][0]["water"].values():
                self.assertEqual(water["name"], "DOD")