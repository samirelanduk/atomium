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


        