from unittest import TestCase
from unittest.mock import Mock, patch, PropertyMock, MagicMock
from atomium.files.utilities import *

class DetermineFileTypeTests(TestCase):

    def test_can_recognise_endings(self):
        self.assertEqual(determine_file_type("/a/d.e/c.pdb", ""), "pdb")
        self.assertEqual(determine_file_type("/a/d.e/c.cif", ""), "cif")
        self.assertEqual(determine_file_type("/a/d.e/c.xyz", ""), "xyz")


    def test_can_recognise_pdb_content(self):
        self.assertEqual(determine_file_type(
         "/a/d.e/c", "HEADER\nATOM\nHETATM"
        ), "pdb")


    def test_can_recognise_mmcif_content(self):
        self.assertEqual(determine_file_type(
         "/a/d.e/c", "loop_\n_atom\nATOM\nHETATM"
        ), "cif")


    def test_can_recognise_xyz_content(self):
        self.assertEqual(determine_file_type(
         "/a/d.e/c", "header\n1 2 3\n3 4 5"
        ), "xyz")



class OpeningTests(TestCase):

    def setUp(self):
        self.patch1 = patch("builtins.open")
        self.patch2 = patch("atomium.files.utilities.parse_string")
        self.mock_open = self.patch1.start()
        self.mock_parse = self.patch2.start()
        open_return = MagicMock()
        mock_file = Mock()
        open_return.__enter__.return_value = mock_file
        mock_file.read.return_value = "returnstring"
        self.mock_open.return_value = open_return


    def tearDown(self):
        self.patch1.stop()
        self.patch2.stop()


    def test_can_open_file(self):
        f = open("path/to/file", 1, a=2)
        self.mock_open.assert_called_with("path/to/file")
        self.mock_parse.assert_called_with("returnstring", "path/to/file", 1, a=2)
        self.assertIs(f, self.mock_parse.return_value)



class FetchingTests(TestCase):

    def setUp(self):
        self.patch1 = patch("atomium.files.utilities.get")
        self.mock_get = self.patch1.start()
        self.response = Mock(status_code=200, text="FILE")
        self.mock_get.return_value = self.response
        self.patch2 = patch("atomium.files.utilities.parse_string")
        self.mock_parse = self.patch2.start()


    def tearDown(self):
        self.patch1.stop()
        self.patch2.stop()


    def test_can_fetch_pdb_code(self):
        f = fetch("1XXX", 1, a=2)
        self.mock_get.assert_called_with("https://files.rcsb.org/view/1xxx.pdb")
        self.mock_parse.assert_called_with("FILE", "1XXX.pdb", 1, a=2)
        self.assertIs(f, self.mock_parse.return_value)


    def test_can_fetch_mmcif(self):
        f = fetch("1XXX.cif", 1, a=2)
        self.mock_get.assert_called_with("https://files.rcsb.org/view/1xxx.cif")
        self.mock_parse.assert_called_with("FILE", "1XXX.cif", 1, a=2)
        self.assertIs(f, self.mock_parse.return_value)


    def test_can_fetch_url(self):
        f = fetch("https://url/", 1, a=2)
        self.mock_get.assert_called_with("https://url/")
        self.mock_parse.assert_called_with("FILE", "https://url/", 1, a=2)
        self.assertIs(f, self.mock_parse.return_value)


    def test_fetching_throws_value_error_if_404(self):
        self.response.status_code = 404
        with self.assertRaises(ValueError):
            fetch("1XXX")
        self.assertFalse(self.mock_parse.called)



class FetchingOverSshTests(TestCase):

    def setUp(self):
        self.patch1 = patch("paramiko.SSHClient")
        self.mock_ssh = self.patch1.start()
        self.mock_client = Mock()
        self.mock_ssh.return_value = self.mock_client
        self.mock_client.exec_command.return_value = (Mock(), Mock(), Mock())
        self.mock_client.exec_command.return_value[1].read.return_value = b"STRING"
        self.patch2 = patch("paramiko.AutoAddPolicy")
        self.mock_policy = self.patch2.start()
        self.mock_policy.return_value = "POLICY"
        self.patch3 = patch("atomium.files.utilities.parse_string")
        self.mock_parse = self.patch3.start()


    def tearDown(self):
        self.patch1.stop()
        self.patch2.stop()
        self.patch3.stop()


    def test_can_get_filestring_over_ssh_with_keys(self):
        f = fetch_over_ssh("HOST", "USER", "/path/", 1, a=2)
        self.mock_client.set_missing_host_key_policy.assert_called_with("POLICY")
        self.mock_client.load_system_host_keys.assert_called_with()
        self.mock_client.connect.assert_called_with(hostname="HOST", username="USER")
        self.mock_client.exec_command.assert_called_with("less /path/")
        self.mock_client.close.assert_called_with()
        self.mock_parse.assert_called_with("STRING", "/path/", 1, a=2)
        self.assertIs(f, self.mock_parse.return_value)


    def test_can_get_filestring_over_ssh_with_password(self):
        f = fetch_over_ssh("HOST", "USER", "/path/", 1, password="xxx", a=2)
        self.mock_client.set_missing_host_key_policy.assert_called_with("POLICY")
        self.assertFalse(self.mock_client.load_system_host_keys.called)
        self.mock_client.connect.assert_called_with(
         hostname="HOST", username="USER", password="xxx"
        )
        self.mock_client.exec_command.assert_called_with("less /path/")
        self.mock_client.close.assert_called_with()
        self.mock_parse.assert_called_with("STRING", "/path/", 1, a=2)
        self.assertIs(f, self.mock_parse.return_value)


    def test_connection_is_always_closed(self):
        self.mock_client.set_missing_host_key_policy.side_effect = Exception
        try:
            fetch_over_ssh("HOST", "USER", "/path/")
        except: pass
        self.mock_client.close.assert_called_with()



class StringParsingTests(TestCase):

    def setUp(self):
        self.patch1 = patch("atomium.files.utilities.determine_file_type")
        self.patch2 = patch("atomium.files.utilities.data_dict_to_file")
        self.mock_type = self.patch1.start()
        self.mock_cont = self.patch2.start()


    @patch("atomium.files.utilities.pdb_string_to_pdb_dict")
    @patch("atomium.files.utilities.pdb_dict_to_data_dict")
    def test_can_parse_pdb(self, mock_data, mock_pdb):
        self.mock_type.return_value = "pdb"
        mock_pdb.return_value = {"PDB": "DICT"}
        mock_data.return_value = {"DATA": "DICT"}
        pdb = parse_string("returnstring", "path/to/file")
        self.mock_type.assert_called_with("path/to/file", "returnstring")
        mock_pdb.assert_called_with("returnstring")
        mock_data.assert_called_with({"PDB": "DICT"})
        self.mock_cont.assert_called_with({"DATA": "DICT"})
        self.assertIs(pdb, self.mock_cont.return_value)
        self.assertEqual(pdb._filetype, "pdb")


    @patch("atomium.files.utilities.pdb_string_to_pdb_dict")
    @patch("atomium.files.utilities.pdb_dict_to_data_dict")
    def test_can_parse_pdb_data(self, mock_data, mock_pdb):
        self.mock_type.return_value = "pdb"
        mock_pdb.return_value = {"PDB": "DICT"}
        mock_data.return_value = {"DATA": "DICT"}
        pdb = parse_string("returnstring", "path/to/file", data_dict=True)
        self.mock_type.assert_called_with("path/to/file", "returnstring")
        mock_pdb.assert_called_with("returnstring")
        mock_data.assert_called_with({"PDB": "DICT"})
        self.assertEqual(pdb, {"DATA": "DICT"})


    @patch("atomium.files.utilities.pdb_string_to_pdb_dict")
    def test_can_parse_pdb_file_dict(self, mock_pdb):
        self.mock_type.return_value = "pdb"
        mock_pdb.return_value = {"PDB": "DICT"}
        pdb = parse_string("returnstring", "path/to/file", file_dict=True)
        self.mock_type.assert_called_with("path/to/file", "returnstring")
        mock_pdb.assert_called_with("returnstring")
        self.assertEqual(pdb, {"PDB": "DICT"})


    @patch("atomium.files.utilities.mmcif_string_to_mmcif_dict")
    @patch("atomium.files.utilities.mmcif_dict_to_data_dict")
    def test_can_parse_mmcif(self, mock_data, mock_mmcif):
        self.mock_type.return_value = "cif"
        mock_mmcif.return_value = {"MMCIF": "DICT"}
        mock_data.return_value = {"DATA": "DICT"}
        mmcif = parse_string("returnstring", "path/to/file")
        self.mock_type.assert_called_with("path/to/file", "returnstring")
        mock_mmcif.assert_called_with("returnstring")
        mock_data.assert_called_with({"MMCIF": "DICT"})
        self.mock_cont.assert_called_with({"DATA": "DICT"})
        self.assertIs(mmcif, self.mock_cont.return_value)
        self.assertEqual(mmcif._filetype, "cif")


    @patch("atomium.files.utilities.mmcif_string_to_mmcif_dict")
    @patch("atomium.files.utilities.mmcif_dict_to_data_dict")
    def test_can_parse_mmcif_data(self, mock_data, mock_mmcif):
        self.mock_type.return_value = "cif"
        mock_mmcif.return_value = {"MMCIF": "DICT"}
        mock_data.return_value = {"DATA": "DICT"}
        mmcif = parse_string("returnstring", "path/to/file", data_dict=True)
        self.mock_type.assert_called_with("path/to/file", "returnstring")
        mock_mmcif.assert_called_with("returnstring")
        mock_data.assert_called_with({"MMCIF": "DICT"})
        self.assertEqual(mmcif, {"DATA": "DICT"})


    @patch("atomium.files.utilities.mmcif_string_to_mmcif_dict")
    def test_can_parse_mmcif_file_dict(self, mock_mmcif):
        self.mock_type.return_value = "cif"
        mock_mmcif.return_value = {"MMCIF": "DICT"}
        mmcif = parse_string("returnstring", "path/to/file", file_dict=True)
        self.mock_type.assert_called_with("path/to/file", "returnstring")
        mock_mmcif.assert_called_with("returnstring")
        self.assertEqual(mmcif, {"MMCIF": "DICT"})


    @patch("atomium.files.utilities.xyz_string_to_xyz_dict")
    @patch("atomium.files.utilities.xyz_dict_to_data_dict")
    def test_can_parse_xyz(self, mock_data, mock_xyz):
        self.mock_type.return_value = "xyz"
        mock_xyz.return_value = {"XYZ": "DICT"}
        mock_data.return_value = {"DATA": "DICT"}
        xyz = parse_string("returnstring", "path/to/file")
        self.mock_type.assert_called_with("path/to/file", "returnstring")
        mock_xyz.assert_called_with("returnstring")
        mock_data.assert_called_with({"XYZ": "DICT"})
        self.mock_cont.assert_called_with({"DATA": "DICT"})
        self.assertIs(xyz, self.mock_cont.return_value)
        self.assertEqual(xyz._filetype, "xyz")


    @patch("atomium.files.utilities.xyz_string_to_xyz_dict")
    @patch("atomium.files.utilities.xyz_dict_to_data_dict")
    def test_can_parse_xyz_data(self, mock_data, mock_xyz):
        self.mock_type.return_value = "xyz"
        mock_xyz.return_value = {"XYZ": "DICT"}
        mock_data.return_value = {"DATA": "DICT"}
        xyz = parse_string("returnstring", "path/to/file", data_dict=True)
        self.mock_type.assert_called_with("path/to/file", "returnstring")
        mock_xyz.assert_called_with("returnstring")
        mock_data.assert_called_with({"XYZ": "DICT"})
        self.assertEqual(xyz, {"DATA": "DICT"})


@patch("atomium.files.utilities.xyz_string_to_xyz_dict")
def test_can_parse_xyz_dfile_dict(self, mock_xyz):
    self.mock_type.return_value = "xyz"
    mock_xyz.return_value = {"XYZ": "DICT"}
    xyz = parse_string("returnstring", "path/to/file", file_dict=True)
    self.mock_type.assert_called_with("path/to/file", "returnstring")
    mock_xyz.assert_called_with("returnstring")
    self.assertEqual(xyz, {"XYZ": "DICT"})
