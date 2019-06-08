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
        