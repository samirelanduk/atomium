from unittest import TestCase
import unittest.mock
from molecupy.pdb import Pdb
from molecupy.pdbdatafile import PdbDataFile

class PdbCreationTests(TestCase):

    def test_can_create_pdb(self):
        data_file = unittest.mock.Mock(spec=PdbDataFile)
        pdb = Pdb(data_file)
        self.assertIs(pdb.data_file(), data_file)


    def test_can_get_data_attributes(self):
        data_file = unittest.mock.Mock(spec=PdbDataFile)
        data_file.classification.return_value = "val1",
        data_file.deposition_date.return_value = "val2",
        data_file.pdb_code.return_value = "val3",
        data_file.is_obsolete.return_value = "val4",
        data_file.obsolete_date.return_value = "val5",
        data_file.replacement_code.return_value = "val6",
        data_file.title.return_value = "val7",
        data_file.split_codes.return_value = "val8",
        data_file.caveat.return_value = "val9",
        data_file.keywords.return_value = "val10",
        data_file.experimental_techniques.return_value = "val11",
        data_file.model_count.return_value = "val12",
        data_file.model_annotations.return_value = "val13",
        data_file.authors.return_value = "val14",
        data_file.revisions.return_value = "val15",
        data_file.supercedes.return_value = "val16",
        data_file.supercede_date.return_value = "val17",
        data_file.journal.return_value = "val18"
        pdb = Pdb(data_file)
        self.assertIs(
         pdb.classification(),
         data_file.classification()
        )
        self.assertIs(
         pdb.deposition_date(),
         data_file.deposition_date()
        )
        self.assertIs(
         pdb.pdb_code(),
         data_file.pdb_code(),
        )
        self.assertIs(
         pdb.is_obsolete(),
         data_file.is_obsolete()
        )
        self.assertIs(
         pdb.obsolete_date(),
         data_file.obsolete_date()
        )
        self.assertIs(
         pdb.replacement_code(),
         data_file.replacement_code()
        )
        self.assertIs(
         pdb.title(),
         data_file.title()
        )
        self.assertIs(
         pdb.split_codes(),
         data_file.split_codes()
        )
        self.assertIs(
         pdb.caveat(),
         data_file.caveat()
        )
        self.assertIs(
         pdb.keywords(),
         data_file.keywords()
        )
        self.assertIs(
         pdb.experimental_techniques(),
         data_file.experimental_techniques()
        )
        self.assertIs(
         pdb.model_count(),
         data_file.model_count()
        )
        self.assertIs(
         pdb.model_annotations(),
         data_file.model_annotations()
        )
        self.assertIs(
         pdb.revisions(),
         data_file.revisions()
        )
        self.assertIs(
         pdb.supercedes(),
         data_file.supercedes()
        )
        self.assertIs(
         pdb.supercede_date(),
         data_file.supercede_date()
        )
        self.assertIs(
         pdb.journal(),
         data_file.journal()
        )
