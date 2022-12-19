from pathlib import Path
from unittest import TestCase
from unittest.mock import patch, Mock
from atomium.io import *

class OpenTests(TestCase):

    def setUp(self):
        self.patch1 = patch("gzip.open")
        self.patch2 = patch("builtins.open")
        self.patch3 = patch("atomium.io.parse_filestring")
        self.mock_gzopen = self.patch1.start()
        self.mock_open = self.patch2.start()
        self.mock_parse = self.patch3.start()
        self.mock_gzread = self.mock_gzopen.return_value.__enter__.return_value.read
        self.mock_read = self.mock_open.return_value.__enter__.return_value.read
    

    def tearDown(self):
        self.patch1.stop()
        self.patch2.stop()
        self.patch3.stop()


    def test_str_dict_compressed_binary(self):
        self.mock_gzread.side_effect = [UnicodeDecodeError, b"xxx"]
        f = open("/path/file.bcif.gz", dictionary=True)
        self.assertEqual(f, self.mock_parse.return_value)
        self.mock_gzopen.assert_any_call("/path/file.bcif.gz")
        self.mock_gzopen.assert_any_call("/path/file.bcif.gz", "rb")
        self.mock_parse.assert_called_with(b"xxx", "/path/file.bcif", dictionary=True)
    

    def test_str_dict_compressed_plain(self):
        self.mock_gzread.return_value = b"xxx"
        f = open("/path/file.pdb.gz", dictionary=True)
        self.assertEqual(f, self.mock_parse.return_value)
        self.mock_gzopen.assert_called_with("/path/file.pdb.gz")
        self.mock_parse.assert_called_with("xxx", "/path/file.pdb", dictionary=True)
    

    def test_str_dict_uncompressed_binary(self):
        self.mock_read.side_effect = [UnicodeDecodeError, b"xxx"]
        f = open("/path/file.bcif", dictionary=True)
        self.assertEqual(f, self.mock_parse.return_value)
        self.mock_open.assert_any_call("/path/file.bcif")
        self.mock_open.assert_any_call("/path/file.bcif", "rb")
        self.mock_parse.assert_called_with(b"xxx", "/path/file.bcif", dictionary=True)
    

    def test_str_dict_uncompressed_plain(self):
        self.mock_read.return_value = "xxx"
        f = open("/path/file.pdb", dictionary=True)
        self.assertEqual(f, self.mock_parse.return_value)
        self.mock_open.assert_called_with("/path/file.pdb")
        self.mock_parse.assert_called_with("xxx", "/path/file.pdb", dictionary=True)
    

    def test_str_full_compressed_binary(self):
        self.mock_gzread.side_effect = [UnicodeDecodeError, b"xxx"]
        f = open("/path/file.bcif.gz")
        self.assertEqual(f, self.mock_parse.return_value)
        self.mock_gzopen.assert_any_call("/path/file.bcif.gz")
        self.mock_gzopen.assert_any_call("/path/file.bcif.gz", "rb")
        self.mock_parse.assert_called_with(b"xxx", "/path/file.bcif", dictionary=False)
    

    def test_str_full_compressed_plain(self):
        self.mock_gzread.return_value = b"xxx"
        f = open("/path/file.pdb.gz")
        self.assertEqual(f, self.mock_parse.return_value)
        self.mock_gzopen.assert_called_with("/path/file.pdb.gz")
        self.mock_parse.assert_called_with("xxx", "/path/file.pdb", dictionary=False)
    

    def test_str_full_uncompressed_binary(self):
        self.mock_read.side_effect = [UnicodeDecodeError, b"xxx"]
        f = open("/path/file.bcif")
        self.assertEqual(f, self.mock_parse.return_value)
        self.mock_open.assert_any_call("/path/file.bcif")
        self.mock_open.assert_any_call("/path/file.bcif", "rb")
        self.mock_parse.assert_called_with(b"xxx", "/path/file.bcif", dictionary=False)
    

    def test_str_full_uncompressed_plain(self):
        self.mock_read.return_value = "xxx"
        f = open("/path/file.pdb")
        self.assertEqual(f, self.mock_parse.return_value)
        self.mock_open.assert_called_with("/path/file.pdb")
        self.mock_parse.assert_called_with("xxx", "/path/file.pdb", dictionary=False)
    

    def test_path_dict_compressed_binary(self):
        self.mock_gzread.side_effect = [UnicodeDecodeError, b"xxx"]
        f = open(Path("/path/file.bcif.gz"), dictionary=True)
        self.assertEqual(f, self.mock_parse.return_value)
        self.mock_gzopen.assert_any_call("/path/file.bcif.gz")
        self.mock_gzopen.assert_any_call("/path/file.bcif.gz", "rb")
        self.mock_parse.assert_called_with(b"xxx", "/path/file.bcif", dictionary=True)
    

    def test_path_dict_compressed_plain(self):
        self.mock_gzread.return_value = b"xxx"
        f = open(Path("/path/file.pdb.gz"), dictionary=True)
        self.assertEqual(f, self.mock_parse.return_value)
        self.mock_gzopen.assert_called_with("/path/file.pdb.gz")
        self.mock_parse.assert_called_with("xxx", "/path/file.pdb", dictionary=True)
    

    def test_path_dict_uncompressed_binary(self):
        self.mock_read.side_effect = [UnicodeDecodeError, b"xxx"]
        f = open(Path("/path/file.bcif"), dictionary=True)
        self.assertEqual(f, self.mock_parse.return_value)
        self.mock_open.assert_any_call("/path/file.bcif")
        self.mock_open.assert_any_call("/path/file.bcif", "rb")
        self.mock_parse.assert_called_with(b"xxx", "/path/file.bcif", dictionary=True)
    

    def test_path_dict_uncompressed_plain(self):
        self.mock_read.return_value = "xxx"
        f = open(Path("/path/file.pdb"), dictionary=True)
        self.assertEqual(f, self.mock_parse.return_value)
        self.mock_open.assert_called_with("/path/file.pdb")
        self.mock_parse.assert_called_with("xxx", "/path/file.pdb", dictionary=True)
    

    def test_path_full_compressed_binary(self):
        self.mock_gzread.side_effect = [UnicodeDecodeError, b"xxx"]
        f = open(Path("/path/file.bcif.gz"))
        self.assertEqual(f, self.mock_parse.return_value)
        self.mock_gzopen.assert_any_call("/path/file.bcif.gz")
        self.mock_gzopen.assert_any_call("/path/file.bcif.gz", "rb")
        self.mock_parse.assert_called_with(b"xxx", "/path/file.bcif", dictionary=False)
    

    def test_path_full_compressed_plain(self):
        self.mock_gzread.return_value = b"xxx"
        f = open(Path("/path/file.pdb.gz"))
        self.assertEqual(f, self.mock_parse.return_value)
        self.mock_gzopen.assert_called_with("/path/file.pdb.gz")
        self.mock_parse.assert_called_with("xxx", "/path/file.pdb", dictionary=False)
    

    def test_path_full_uncompressed_binary(self):
        self.mock_read.side_effect = [UnicodeDecodeError, b"xxx"]
        f = open(Path("/path/file.bcif"))
        self.assertEqual(f, self.mock_parse.return_value)
        self.mock_open.assert_any_call("/path/file.bcif")
        self.mock_open.assert_any_call("/path/file.bcif", "rb")
        self.mock_parse.assert_called_with(b"xxx", "/path/file.bcif", dictionary=False)
    

    def test_path_full_uncompressed_plain(self):
        self.mock_read.return_value = "xxx"
        f = open(Path("/path/file.pdb"))
        self.assertEqual(f, self.mock_parse.return_value)
        self.mock_open.assert_called_with("/path/file.pdb")
        self.mock_parse.assert_called_with("xxx", "/path/file.pdb", dictionary=False)



class FetchTests(TestCase):

    def setUp(self):
        self.patch1 = patch("requests.get")
        self.patch2 = patch("atomium.io.parse_filestring")
        self.mock_get = self.patch1.start()
        self.mock_parse = self.patch2.start()
        self.response = Mock(status_code=200)
        self.mock_get.return_value = self.response
    

    def tearDown(self):
        self.patch1.stop()
        self.patch2.stop()


    def test_url_bcif_dictionary(self):
        f = fetch("https://files.com/file.bcif", dictionary=True)
        self.assertEqual(f, self.mock_parse.return_value)
        self.mock_get.assert_called_with("https://files.com/file.bcif", stream=True)
        self.mock_parse.assert_called_with(self.response.content, "https://files.com/file.bcif", dictionary=True)
    

    def test_url_mmtf_dictionary(self):
        f = fetch("https://files.com/file.mmtf", dictionary=True)
        self.assertEqual(f, self.mock_parse.return_value)
        self.mock_get.assert_called_with("https://files.com/file.mmtf", stream=True)
        self.mock_parse.assert_called_with(self.response.content, "https://files.com/file.mmtf", dictionary=True)
    

    def test_url_pdb_dictionary(self):
        f = fetch("https://files.com/file.pdb", dictionary=True)
        self.assertEqual(f, self.mock_parse.return_value)
        self.mock_get.assert_called_with("https://files.com/file.pdb", stream=True)
        self.mock_parse.assert_called_with(self.response.text, "https://files.com/file.pdb", dictionary=True)
    

    def test_url_bcif_full(self):
        f = fetch("https://files.com/file.bcif")
        self.assertEqual(f, self.mock_parse.return_value)
        self.mock_get.assert_called_with("https://files.com/file.bcif", stream=True)
        self.mock_parse.assert_called_with(self.response.content, "https://files.com/file.bcif", dictionary=False)
    

    def test_url_mmtf_full(self):
        f = fetch("https://files.com/file.mmtf")
        self.assertEqual(f, self.mock_parse.return_value)
        self.mock_get.assert_called_with("https://files.com/file.mmtf", stream=True)
        self.mock_parse.assert_called_with(self.response.content, "https://files.com/file.mmtf", dictionary=False)
    

    def test_url_pdb_full(self):
        f = fetch("https://files.com/file.pdb")
        self.assertEqual(f, self.mock_parse.return_value)
        self.mock_get.assert_called_with("https://files.com/file.pdb", stream=True)
        self.mock_parse.assert_called_with(self.response.text, "https://files.com/file.pdb", dictionary=False)
    

    def test_code_bcif_dictionary(self):
        f = fetch("1XXX.bcif", dictionary=True)
        self.assertEqual(f, self.mock_parse.return_value)
        self.mock_get.assert_called_with("https://models.rcsb.org/1xxx.bcif", stream=True)
        self.mock_parse.assert_called_with(self.response.content, "1XXX.bcif", dictionary=True)
    

    def test_code_bcif_full(self):
        f = fetch("1XXX.bcif")
        self.assertEqual(f, self.mock_parse.return_value)
        self.mock_get.assert_called_with("https://models.rcsb.org/1xxx.bcif", stream=True)
        self.mock_parse.assert_called_with(self.response.content, "1XXX.bcif", dictionary=False)
    

    def test_code_mmtf_dictionary(self):
        f = fetch("1XXX.mmtf", dictionary=True)
        self.assertEqual(f, self.mock_parse.return_value)
        self.mock_get.assert_called_with("https://mmtf.rcsb.org/v1.0/full/1xxx", stream=True)
        self.mock_parse.assert_called_with(self.response.content, "1XXX.mmtf", dictionary=True)
    

    def test_code_mmtf_full(self):
        f = fetch("1XXX.mmtf")
        self.assertEqual(f, self.mock_parse.return_value)
        self.mock_get.assert_called_with("https://mmtf.rcsb.org/v1.0/full/1xxx", stream=True)
        self.mock_parse.assert_called_with(self.response.content, "1XXX.mmtf", dictionary=False)
    

    def test_code_pdb_dictionary(self):
        f = fetch("1XXX.pdb", dictionary=True)
        self.assertEqual(f, self.mock_parse.return_value)
        self.mock_get.assert_called_with("https://files.rcsb.org/view/1xxx.pdb", stream=True)
        self.mock_parse.assert_called_with(self.response.text, "1XXX.pdb", dictionary=True)
    

    def test_code_pdb_full(self):
        f = fetch("1XXX.pdb")
        self.assertEqual(f, self.mock_parse.return_value)
        self.mock_get.assert_called_with("https://files.rcsb.org/view/1xxx.pdb", stream=True)
        self.mock_parse.assert_called_with(self.response.text, "1XXX.pdb", dictionary=False)
    

    def test_code_cif_dictionary(self):
        f = fetch("1XXX.cif", dictionary=True)
        self.assertEqual(f, self.mock_parse.return_value)
        self.mock_get.assert_called_with("https://files.rcsb.org/view/1xxx.cif", stream=True)
        self.mock_parse.assert_called_with(self.response.text, "1XXX.cif", dictionary=True)
    

    def test_code_cif_full(self):
        f = fetch("1XXX.cif")
        self.assertEqual(f, self.mock_parse.return_value)
        self.mock_get.assert_called_with("https://files.rcsb.org/view/1xxx.cif", stream=True)
        self.mock_parse.assert_called_with(self.response.text, "1XXX.cif", dictionary=False)
    

    def test_code_no_ext_cif_dictionary(self):
        f = fetch("1XXX", dictionary=True)
        self.assertEqual(f, self.mock_parse.return_value)
        self.mock_get.assert_called_with("https://files.rcsb.org/view/1xxx.cif", stream=True)
        self.mock_parse.assert_called_with(self.response.text, "1XXX.cif", dictionary=True)
    

    def test_code_no_ext_cif_full(self):
        f = fetch("1XXX")
        self.assertEqual(f, self.mock_parse.return_value)
        self.mock_get.assert_called_with("https://files.rcsb.org/view/1xxx.cif", stream=True)
        self.mock_parse.assert_called_with(self.response.text, "1XXX.cif", dictionary=False)
    

    def test_failed_request(self):
        self.response.status_code = 404
        with self.assertRaises(ValueError) as err:
            fetch("1XXX")
        self.assertEqual(
            str(err.exception),
            "Could not find anything at https://files.rcsb.org/view/1xxx.cif"
        )
    


class ParseFilestringTests(TestCase):

    @patch("atomium.io.determine_filetype")
    @patch("atomium.io.mmcif_string_to_mmcif_dict")
    def test_cif_to_dict(self, mock_convert, mock_type):
        mock_type.return_value = "mmcif"
        d = parse_filestring("filestring", "file.cif", dictionary=True)
        self.assertEqual(d, mock_convert.return_value)
        mock_type.assert_called_with("filestring", "file.cif")
        mock_convert.assert_called_with("filestring")
    

    @patch("atomium.io.determine_filetype")
    @patch("atomium.io.mmcif_string_to_mmcif_dict")
    @patch("atomium.io.File")
    def test_cif_to_file(self, mock_file, mock_convert, mock_type):
        mock_type.return_value = "mmcif"
        d = parse_filestring("filestring", "file.cif")
        self.assertEqual(d, mock_file.return_value)
        mock_type.assert_called_with("filestring", "file.cif")
        mock_convert.assert_called_with("filestring")
        mock_file.assert_called_with(mock_convert.return_value)


    @patch("atomium.io.determine_filetype")
    @patch("atomium.io.pdb_string_to_mmcif_dict")
    def test_pdb_to_dict(self, mock_convert, mock_type):
        mock_type.return_value = "pdb"
        d = parse_filestring("filestring", "file.pdb", dictionary=True)
        self.assertEqual(d, mock_convert.return_value)
        mock_type.assert_called_with("filestring", "file.pdb")
        mock_convert.assert_called_with("filestring")
    

    @patch("atomium.io.determine_filetype")
    @patch("atomium.io.pdb_string_to_mmcif_dict")
    @patch("atomium.io.File")
    def test_pdb_to_file(self, mock_file, mock_convert, mock_type):
        mock_type.return_value = "pdb"
        d = parse_filestring("filestring", "file.pdb")
        self.assertEqual(d, mock_file.return_value)
        mock_type.assert_called_with("filestring", "file.pdb")
        mock_convert.assert_called_with("filestring")
        mock_file.assert_called_with(mock_convert.return_value)
    

    @patch("atomium.io.determine_filetype")
    @patch("atomium.io.bcif_string_to_mmcif_dict")
    def test_bcif_to_dict(self, mock_convert, mock_type):
        mock_type.return_value = "bcif"
        d = parse_filestring("filestring", "file.bcif", dictionary=True)
        self.assertEqual(d, mock_convert.return_value)
        mock_type.assert_called_with("filestring", "file.bcif")
        mock_convert.assert_called_with("filestring")
    

    @patch("atomium.io.determine_filetype")
    @patch("atomium.io.bcif_string_to_mmcif_dict")
    @patch("atomium.io.File")
    def test_bcif_to_file(self, mock_file, mock_convert, mock_type):
        mock_type.return_value = "bcif"
        d = parse_filestring("filestring", "file.bcif")
        self.assertEqual(d, mock_file.return_value)
        mock_type.assert_called_with("filestring", "file.bcif")
        mock_convert.assert_called_with("filestring")
        mock_file.assert_called_with(mock_convert.return_value)
    

    @patch("atomium.io.determine_filetype")
    @patch("atomium.io.mmtf_string_to_mmcif_dict")
    def test_mmtf_to_dict(self, mock_convert, mock_type):
        mock_type.return_value = "mmtf"
        d = parse_filestring("filestring", "file.mmtf", dictionary=True)
        self.assertEqual(d, mock_convert.return_value)
        mock_type.assert_called_with("filestring", "file.mmtf")
        mock_convert.assert_called_with("filestring")
    

    @patch("atomium.io.determine_filetype")
    @patch("atomium.io.mmtf_string_to_mmcif_dict")
    @patch("atomium.io.File")
    def test_mmtf_to_file(self, mock_file, mock_convert, mock_type):
        mock_type.return_value = "mmtf"
        d = parse_filestring("filestring", "file.mmtf")
        self.assertEqual(d, mock_file.return_value)
        mock_type.assert_called_with("filestring", "file.mmtf")
        mock_convert.assert_called_with("filestring")
        mock_file.assert_called_with(mock_convert.return_value)



class DetermineFiletypeTests(TestCase):

    def test_can_get_cif(self):
        self.assertEqual(determine_filetype("data_", "file.cif"), "mmcif")
    

    def test_can_get_pdb(self):
        self.assertEqual(determine_filetype("ATOM", "file.pdb"), "pdb")
    

    def test_can_get_bcif(self):
        self.assertEqual(determine_filetype(b"data_", "file.bcif"), "bcif")
    

    def test_can_get_mmtf(self):
        self.assertEqual(determine_filetype(b"data_", "file.mmtf"), "mmtf")



class SaveDictionaryTests(TestCase):

    @patch("atomium.io.save_mmcif_dict_to_mmcif")
    def test_str_path_cif(self, mock_save):
        save_dictionary({"mmcif": 1}, "/path.cif")
        mock_save.assert_called_with({"mmcif": 1}, "/path.cif")
    

    @patch("atomium.io.save_mmcif_dict_to_bcif")
    def test_str_path_bcif(self, mock_save):
        save_dictionary({"mmcif": 1}, "/path.bcif")
        mock_save.assert_called_with({"mmcif": 1}, "/path.bcif")
    

    @patch("atomium.io.save_mmcif_dict_to_pdb")
    def test_str_path_pdb(self, mock_save):
        save_dictionary({"mmcif": 1}, "/path.pdb")
        mock_save.assert_called_with({"mmcif": 1}, "/path.pdb")
    

    @patch("atomium.io.save_mmcif_dict_to_mmcif")
    def test_str_path_cif(self, mock_save):
        save_dictionary({"mmcif": 1}, Path("/path.cif"))
        mock_save.assert_called_with({"mmcif": 1}, Path("/path.cif"))
    

    @patch("atomium.io.save_mmcif_dict_to_bcif")
    def test_str_path_bcif(self, mock_save):
        save_dictionary({"mmcif": 1}, Path("/path.bcif"))
        mock_save.assert_called_with({"mmcif": 1}, Path("/path.bcif"))
    

    @patch("atomium.io.save_mmcif_dict_to_pdb")
    def test_str_path_pdb(self, mock_save):
        save_dictionary({"mmcif": 1}, Path("/path.pdb"))
        mock_save.assert_called_with({"mmcif": 1}, Path("/path.pdb"))