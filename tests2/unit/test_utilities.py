from unittest import TestCase
from unittest.mock import Mock, patch, PropertyMock, MagicMock
from atomium.utilities import *

class OpeningTests(TestCase):

    def setUp(self):
        self.patch1 = patch("builtins.open")
        self.patch2 = patch("atomium.utilities.parse_string")
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


    def test_can_open_string(self):
        self.assertEqual(open("path/to/file", 1, a=2), self.mock_parse.return_value)
        self.mock_open.assert_called_with("path/to/file")
        self.mock_parse.assert_called_with("returnstring", "path/to/file", 1, a=2)


    def test_can_open_bytestring(self):
        self.mock_open.side_effect = [Exception, self.mock_open.return_value]
        self.assertEqual(open("path/to/file", 1, a=2), self.mock_parse.return_value)
        self.mock_open.assert_called_with("path/to/file", "rb")
        self.mock_parse.assert_called_with("returnstring", "path/to/file", 1, a=2)



class FetchingTests(TestCase):

    def setUp(self):
        self.patch1 = patch("atomium.utilities.get")
        self.mock_get = self.patch1.start()
        self.mock_get.return_value = Mock(status_code=200, text="ABC")
        self.patch2 = patch("atomium.utilities.parse_string")
        self.mock_parse = self.patch2.start()


    def tearDown(self):
        self.patch1.stop()
        self.patch2.stop()


    def test_can_fetch_cif(self):
        f = fetch("1ABC", 1, b=2)
        self.mock_get.assert_called_with("https://files.rcsb.org/view/1abc.cif", stream=True)
        self.mock_parse.assert_called_with("ABC", "1ABC.cif", 1, b=2)
        self.assertEqual(f, self.mock_parse.return_value)


    def test_can_fetch_pdb(self):
        f = fetch("1ABC.pdb", 1, b=2)
        self.mock_get.assert_called_with("https://files.rcsb.org/view/1abc.pdb", stream=True)
        self.mock_parse.assert_called_with("ABC", "1ABC.pdb", 1, b=2)
        self.assertEqual(f, self.mock_parse.return_value)


    def test_can_fetch_mmtf(self):
        self.mock_get.return_value.content = b"ABC"
        f = fetch("1ABC.mmtf", 1, b=2)
        self.mock_get.assert_called_with("https://mmtf.rcsb.org/v1.0/full/1abc", stream=True)
        self.mock_parse.assert_called_with(b"ABC", "1ABC.mmtf", 1, b=2)
        self.assertEqual(f, self.mock_parse.return_value)


    def test_can_fetch_by_url(self):
        f = fetch("https://website.com/1ABC", 1, b=2)
        self.mock_get.assert_called_with("https://website.com/1ABC", stream=True)
        self.mock_parse.assert_called_with("ABC", "https://website.com/1ABC", 1, b=2)
        self.assertEqual(f, self.mock_parse.return_value)


    def test_can_handle_no_results(self):
        self.mock_get.return_value.status_code = 400
        with self.assertRaises(ValueError):
            f = fetch("1ABC", 1, b=2)


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
        self.patch3 = patch("atomium.utilities.parse_string")
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

    @patch("atomium.utilities.get_parse_functions")
    def test_can_get_file_dict(self, mock_get):
        mock_get.return_value = [MagicMock(), MagicMock()]
        f = parse_string("ABCD", "file.xyz", file_dict=True)
        mock_get.assert_called_with("ABCD", "file.xyz")
        mock_get.return_value[0].assert_called_with("ABCD")
        self.assertEqual(f, mock_get.return_value[0].return_value)


    @patch("atomium.utilities.get_parse_functions")
    def test_can_get_data_dict(self, mock_get):
        mock_get.return_value = [MagicMock(), MagicMock()]
        f = parse_string("ABCD", "file.xyz", data_dict=True)
        mock_get.assert_called_with("ABCD", "file.xyz")
        mock_get.return_value[0].assert_called_with("ABCD")
        mock_get.return_value[1].assert_called_with(mock_get.return_value[0].return_value)
        self.assertEqual(f, mock_get.return_value[1].return_value)


    @patch("atomium.utilities.get_parse_functions")
    @patch("atomium.utilities.data_dict_to_file")
    def test_can_get_file(self, mock_data, mock_get):
        mock_get.return_value = [MagicMock(), MagicMock()]
        mock_get.return_value[1].__name__ = "mmcif_Z"
        f = parse_string("ABCD", "file.cif")
        mock_get.assert_called_with("ABCD", "file.cif")
        mock_get.return_value[0].assert_called_with("ABCD")
        mock_get.return_value[1].assert_called_with(mock_get.return_value[0].return_value)
        mock_data.assert_called_with(mock_get.return_value[1].return_value, "cif")
        self.assertEqual(f, mock_data.return_value)



class ParseFunctionGettingTests(TestCase):

    def test_can_get_cif_functions(self):
        f1, f2 = get_parse_functions("ABC", "x.cif")
        self.assertIs(f1, mmcif_string_to_mmcif_dict)
        self.assertIs(f2, mmcif_dict_to_data_dict)


    def test_can_get_mmtf_functions(self):
        f1, f2 = get_parse_functions("ABC", "x.mmtf")
        self.assertIs(f1, mmtf_bytes_to_mmtf_dict)
        self.assertIs(f2, mmtf_dict_to_data_dict)


    def test_can_get_pdb_functions(self):
        f1, f2 = get_parse_functions("ABC", "x.pdb")
        self.assertIs(f1, pdb_string_to_pdb_dict)
        self.assertIs(f2, pdb_dict_to_data_dict)


    def test_bytes_mean_mmtf(self):
        f1, f2 = get_parse_functions(b"ABC", "x.xxx")
        self.assertIs(f1, mmtf_bytes_to_mmtf_dict)
        self.assertIs(f2, mmtf_dict_to_data_dict)


    def test_can_identify_cif(self):
        f1, f2 = get_parse_functions("ABC_atom_sites", "x.xxx")
        self.assertIs(f1, mmcif_string_to_mmcif_dict)
        self.assertIs(f2, mmcif_dict_to_data_dict)


    def test_can_identify_pdb(self):
        f1, f2 = get_parse_functions("ABC", "x.xxx")
        self.assertIs(f1, pdb_string_to_pdb_dict)
        self.assertIs(f2, pdb_dict_to_data_dict)



class Saving(TestCase):

    @patch("builtins.open")
    def test_saves_string_to_file(self, mock_open):
        open_return = MagicMock()
        mock_file = Mock()
        mock_write = MagicMock()
        mock_file.write = mock_write
        open_return.__enter__.return_value = mock_file
        mock_open.return_value = open_return
        save("filestring", "filename")
        mock_open.assert_called_once_with("filename", "w")
        mock_write.assert_called_once_with("filestring")


    @patch("builtins.open")
    def test_saves_bytestring_to_file(self, mock_open):
        open_return = MagicMock()
        mock_file = Mock()
        mock_write = MagicMock()
        mock_file.write = mock_write
        open_return.__enter__.return_value = mock_file
        mock_open.return_value = open_return
        mock_write.side_effect = [Exception, None]
        save(b"filestring", "filename")
        mock_open.assert_called_with("filename", "wb")
        mock_write.assert_called_with(b"filestring")
