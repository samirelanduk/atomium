from collections import deque
import struct
from datetime import date
from unittest import TestCase
from unittest.mock import Mock, patch, PropertyMock, MagicMock
from atomium.mmtf import *

class MmtfBytesToMmtfDictTests(TestCase):

    @patch("msgpack.unpackb")
    @patch("atomium.mmtf.decode_dict")
    def test_can_turn_mmtf_bytes_to_mmtf_dict(self, mock_dec, mock_unpack):
        d = mmtf_bytes_to_mmtf_dict(b"1234")
        mock_unpack.assert_called_with(b"1234")
        mock_dec.assert_called_with(mock_unpack.return_value)
        self.assertIs(d, mock_dec.return_value)



class BytesDictDecodingTests(TestCase):

    def test_can_handle_normal_dictionary_with_byte_keys(self):
        d = decode_dict({
         b"key": "value", b"key2": ["a", "b"], b"key3": 19.5
        })
        self.assertEqual(d, {
         "key": "value", "key2": ["a", "b"], "key3": 19.5
        })


    def test_can_handle_dictionary_with_byte_lists(self):
        d = decode_dict({
         b"key": "value", b"key2": [b"a", b"b"], b"key3": 19.5
        })
        self.assertEqual(d, {
         "key": "value", "key2": ["a", "b"], "key3": 19.5
        })


    @patch("atomium.mmtf.parse_binary_field")
    def test_can_handle_dictionary_with_binary_fields(self, mock_parse):
        mock_parse.side_effect = [1, 1]
        d = decode_dict({
         b"key": b"\x00\x08", b"key2": [b"a", b"b"], b"key3": "\x00123"
        })
        mock_parse.assert_any_call(b"\x00\x08")
        mock_parse.assert_any_call(b"\x00123")
        self.assertEqual(d, {
         "key": 1, "key2": ["a", "b"], "key3": 1
        })


    def test_can_handle_dictionary_with_recursive_fields(self):
        d = decode_dict({
         b"key": "value", b"key2": [{b"keya": [b"a", b"b"], b"key3": 19.5}]
        })
        self.assertEqual(d, {
         "key": "value", "key2": [{"keya": ["a", "b"], "key3": 19.5}]
        })



class BinaryFieldParsingTests(TestCase):

    def test_can_raise_exception_on_invalid_codec(self):
        with self.assertRaises(ValueError):
            parse_binary_field(struct.pack(">iii", 20, 0, 0))


    def test_can_parse_with_codec_1(self):
        b = parse_binary_field(
         struct.pack(">iii", 1, 3, 0) + struct.pack("fff", 4.5, 5.5, 6.5)
        )
        self.assertEqual(b, (4.5, 5.5, 6.5))


    def test_can_parse_with_codec_2(self):
        b = parse_binary_field(
         struct.pack(">iii", 2, 3, 0) + b"\x04\x05\x06"
        )
        self.assertEqual(b, (4, 5, 6))


    def test_can_parse_with_codec_3(self):
        b = parse_binary_field(
         struct.pack(">iii", 3, 3, 0) + struct.pack(">hhh", 4, 5, 6)
        )
        self.assertEqual(b, (4, 5, 6))


    def test_can_parse_with_codec_4(self):
        b = parse_binary_field(
         struct.pack(">iii", 4, 3, 0) + struct.pack(">iii", 4, 5, 6)
        )
        self.assertEqual(b, (4, 5, 6))


    def test_can_parse_with_codec_5(self):
        b = parse_binary_field(
         struct.pack(">iii", 5, 3, 0) + bytes([65, 0, 0, 0, 66, 0, 0, 0, 67, 0, 0, 0])
        )
        self.assertEqual(b, ["A", "B", "C"])
        b = parse_binary_field(
         struct.pack(">iii", 5, 2, 0) + bytes([65, 0, 0, 0, 68, 65, 0, 0])
        )
        self.assertEqual(b, ["A", "DA"])


    @patch("atomium.mmtf.run_length_decode")
    def test_can_parse_with_codec_6(self, mock_rldec):
        mock_rldec.return_value = [100, 0, 0, 105, 0, 0]
        b = parse_binary_field(
         struct.pack(">iii", 6, 3, 0) + struct.pack(">iii", 0, 5, 65)
        )
        mock_rldec.assert_called_with((0, 5, 65))
        self.assertEqual(b, ["d", "", "", "i", "", ""])


    @patch("atomium.mmtf.run_length_decode")
    def test_can_parse_with_codec_7(self, mock_rldec):
        mock_rldec.return_value = [100, 0, 0, 105, 0, 0]
        b = parse_binary_field(
         struct.pack(">iii", 7, 3, 0) + struct.pack(">iii", 0, 5, 65)
        )
        mock_rldec.assert_called_with((0, 5, 65))
        self.assertEqual(b, [100, 0, 0, 105, 0, 0])


    @patch("atomium.mmtf.run_length_decode")
    @patch("atomium.mmtf.delta_decode")
    def test_can_parse_with_codec_8(self, mock_deldec, mock_rldec):
        mock_rldec.return_value = [100, 0, 0, 105, 0, 0]
        mock_deldec.return_value = [1, 2, 3]
        b = parse_binary_field(
         struct.pack(">iii", 8, 3, 0) + struct.pack(">iii", 0, 5, 65)
        )
        mock_rldec.assert_called_with((0, 5, 65))
        mock_deldec.assert_called_with([100, 0, 0, 105, 0, 0])
        self.assertEqual(b, [1, 2, 3])


    @patch("atomium.mmtf.run_length_decode")
    def test_can_parse_with_codec_9(self, mock_rldec):
        mock_rldec.return_value = [1000, 0, 0, 1050, 0, 0]
        b = parse_binary_field(
         struct.pack(">iii", 9, 3, 0) + struct.pack(">iii", 0, 5, 65)
        )
        mock_rldec.assert_called_with((0, 5, 65))
        self.assertEqual(b, [10, 0, 0, 10.5, 0, 0])


    @patch("atomium.mmtf.recursive_decode")
    @patch("atomium.mmtf.delta_decode")
    def test_can_parse_with_codec_8(self, mock_deldec, mock_recdec):
        mock_recdec.return_value = [100, 0, 0, 105, 0, 0]
        mock_deldec.return_value = [1000, 2000, 3000]
        b = parse_binary_field(
         struct.pack(">iii", 10, 3, 0) + struct.pack(">hhh", 0, 5, 65)
        )
        mock_recdec.assert_called_with((0, 5, 65))
        mock_deldec.assert_called_with([100, 0, 0, 105, 0, 0])
        self.assertEqual(b, [1, 2, 3])



class RunLengthDecodingTests(TestCase):

    def test_can_run_length_decode(self):
        self.assertEqual(
         run_length_decode([1, 10, 2, 1, 1, 4]),
         [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1]
        )



class DeltaDecodingTests(TestCase):

    def test_can_delta_decode(self):
        self.assertEqual(
         delta_decode([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1]),
         [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16]
        )



class RecursiveDecodingTests(TestCase):

    def test_can_recursive_decode(self):
        self.assertEqual(
         recursive_decode([32767, 32767, 32767, 6899, 0, 2, -1, 100, -3, 5]),
         [105200, 0, 2, -1, 100, -3, 5]
        )


    def test_can_recursive_decode_other_powers(self):
        self.assertEqual(recursive_decode([
         127, 41, 34, 1, 0, -50, -128, 0, 7, 127, 0, 127, 127, 14
        ], bits=8), [168, 34, 1, 0, -50, -128, 7, 127, 268])



class MmtfDictToDataDictTests(TestCase):

    @patch("atomium.mmtf.mmtf_to_data_transfer")
    def test_can_convert_mmtf_dict_to_data_dict(self, mock_trans):
        m = {"A": "B", "bioAssemblyList": [{
         "name": "1", "transformList": [{
          "chainIndexList": [1, 2], "matrix": "abcdefghijklmnop"
         }, {
          "chainIndexList": [0, 2], "matrix": "ABCDEFGHIJKLMNOP"
         }]
        }], "chainIdList": ["A", "B", "C"]}
        d = mmtf_dict_to_data_dict(m)
        mock_trans.assert_any_call(m, d, "description", "code", "structureId")
        mock_trans.assert_any_call(m, d, "description", "title", "title")
        mock_trans.assert_any_call(m, d, "description", "deposition_date", "depositionDate", date=True)
        mock_trans.assert_any_call(m, d, "experiment", "technique", "experimentalMethods", first=True)
        mock_trans.assert_any_call(m, d, "quality", "resolution", "resolution", trim=3)
        mock_trans.assert_any_call(m, d, "quality", "rvalue", "rWork", trim=3)
        mock_trans.assert_any_call(m, d, "quality", "rfree", "rFree", trim=3)
        self.assertEqual(d, {
         "description": {
          "code": None, "title": None, "deposition_date": None,
          "classification": None, "keywords": [], "authors": []
         }, "experiment": {
          "technique": None, "source_organism": None, "expression_system": None
         }, "quality": {"resolution": None, "rvalue": None, "rfree": None},
         "geometry": {"assemblies": [{
          "id": 1, "software": None, "delta_energy": None,
          "buried_surface_area": None, "surface_area": None, "transformations": [{
           "chains": ["B", "C"], "matrix": ["abc", "efg", "ijk"], "vector": "dhl"
          }, {
           "chains": ["A", "C"], "matrix": ["ABC", "EFG", "IJK"], "vector": "DHL"
          }]
         }]}
        })



class MmtfDictTransferTests(TestCase):

    def setUp(self):
        self.d = {"A": {1: None, 2: None, 3: None}, "B": {4: None, 5: None}}
        self.m = {"M": 10.127, "C": 100, "D": "2018-09-17", "L": [8, 9]}


    def test_can_do_nothing(self):
        mmtf_to_data_transfer(self.m, self.d, "A", 1, "X", 10)
        self.assertEqual(self.d["A"][1], None)
        mmtf_to_data_transfer(self.m, self.d, "A", 1, "M", 100)
        self.assertEqual(self.d["A"][1], None)


    def test_can_transfer_from_mmtf_to_data_dict(self):
        mmtf_to_data_transfer(self.m, self.d, "A", 1, "M")
        self.assertEqual(self.d["A"][1], 10.127)
        mmtf_to_data_transfer(self.m, self.d, "B", 5, "C")
        self.assertEqual(self.d["B"][5], 100)


    def test_can_transfer_from_mmtf_to_data_dict_date(self):
        mmtf_to_data_transfer(self.m, self.d, "B", 5, "D", date=True)
        self.assertEqual(self.d["B"][5], date(2018, 9, 17))


    def test_can_transfer_from_mmtf_to_data_dict_first(self):
        mmtf_to_data_transfer(self.m, self.d, "B", 5, "L", first=True)
        self.assertEqual(self.d["B"][5], 8)


    def test_can_transfer_from_mmtf_to_data_dict_round(self):
        mmtf_to_data_transfer(self.m, self.d, "B", 5, "M", trim=1)
        self.assertEqual(self.d["B"][5], 10.1)
        mmtf_to_data_transfer(self.m, self.d, "B", 5, "M", trim=2)
        self.assertEqual(self.d["B"][5], 10.13)
