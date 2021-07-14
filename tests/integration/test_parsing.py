from datetime import date
from unittest import TestCase
import atomium

class DictParsingTests(TestCase):

    def test_1lol_mmcif(self):
        d = atomium.open("tests/integration/files/1lol.cif", dictionary=True)
        self.assertEqual(len(d.keys()), 66)



class ParsingTests(TestCase):

    def test_1lol_mmcif(self):
        # Open file
        pdb = atomium.open("tests/integration/files/1lol.cif")

        # File information
        self.assertEqual(pdb.name, "1LOL")
        self.assertEqual(len(pdb.source.keys()), 66)
        self.assertEqual(pdb.source["entry"][0]["id"], "1LOL")
        self.assertEqual(pdb.entry__id, "1LOL")
        self.assertEqual(pdb.title, "Crystal structure of orotidine monophosphate decarboxylase complex with XMP")
        self.assertEqual(pdb.deposition_date, date(2002, 5, 6))
        self.assertEqual(pdb.classification, "LYASE")
        self.assertEqual(pdb.keywords, ["TIM barrel", "LYASE"])
        self.assertEqual(pdb.authors, ["Wu, N.", "Pai, E.F."])
        self.assertEqual(pdb.technique, "X-RAY DIFFRACTION")
        self.assertEqual(pdb.source_organism, "Methanothermobacter thermautotrophicus str. Delta H")
        self.assertEqual(pdb.expression_system, "Escherichia coli")
        self.assertEqual(pdb.resolution, 1.9)
        self.assertEqual(pdb.rvalue, 0.193)
        self.assertEqual(pdb.rfree, 0.229)

        # More complex file information
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
        self.assertEqual(pdb.missing_residues, missing_residues)
        self.assertEqual(pdb.assemblies, [{
            "id": 1, "software": "PISA", "delta_energy": -31.0,
            "buried_surface_area":5230, "surface_area": 16550,
            "transformations": [{
                "chains": ["A", "B", "C", "D", "E", "F", "G", "H"],
                "matrix": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                "vector": [0.0, 0.0, 0.0]
            }]
        }])
