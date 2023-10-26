import msgpack
from unittest import TestCase
from unittest.mock import patch
from atomium.bcif import *

class BcifStringToMmcifDictTests(TestCase):

    @patch("atomium.bcif.category_to_table")
    def test_can_get_dict(self, mock_cat):
        filestring = msgpack.packb({"dataBlocks": [{b"categories": [
            {"name": "_cat1"}, {"name": "_cat2"}, {"name": "_cat3"}
        ]}]})
        mock_cat.side_effect = [1, 2, 3]
        mmcif = bcif_string_to_mmcif_dict(filestring)
        self.assertEqual(mmcif, {"cat1": 1, "cat2": 2, "cat3": 3})
        mock_cat.assert_any_call({b"name": b"_cat1"})
        mock_cat.assert_any_call({b"name": b"_cat2"})
        mock_cat.assert_any_call({b"name": b"_cat3"})



class CategoryToTableTests(TestCase):

    @patch("atomium.bcif.parse_column_data")
    def test_can_convert_category_to_table(self, mock_col):
        mock_col.side_effect = [["0", "1", "2"], [""], ["10", "20"], ["30", "40"], ["50", "60"]]
        category = {b"columns": [
            {b"mask": "M1", b"data": 1, b"name": b"prop1"},
            {b"mask": None, b"data": 2, b"name": b"prop2"},
            {b"mask": "M2", b"data": 3, b"name": b"prop3"},
        ], b"rowCount": 2}
        table = category_to_table(category)
        mock_col.assert_any_call("M1")
        mock_col.assert_any_call("M2")
        mock_col.assert_any_call(1, ["", ".", "?"])
        mock_col.assert_any_call(2, None)
        mock_col.assert_any_call(3, None)
        self.assertEqual(table, [
            {"prop1": "10", "prop2": "30", "prop3": "50"},
            {"prop1": "20", "prop2": "40", "prop3": "60"},
        ])



class ColumnParsingTests(TestCase):

    @patch("atomium.bcif.decode")
    def test_can_get_simple_data(self, mock_decode):
        mock_decode.side_effect = ["stage1", "stage2", [1, 2, 3]]
        column = parse_column_data({b"data": "initial", b"encoding": ["enc3", "enc2", "enc1"]})
        mock_decode.assert_any_call("initial", "enc1")
        mock_decode.assert_any_call("stage1", "enc2")
        mock_decode.assert_any_call("stage2", "enc3")
        self.assertEqual(column, ["1", "2", "3"])
    

    @patch("atomium.bcif.decode")
    def test_can_get_data_with_line_breaks(self, mock_decode):
        mock_decode.side_effect = ["stage1", "stage2", [1, "abc\ndef\nghi", "line 1\line 2\nline3"]]
        column = parse_column_data({b"data": "initial", b"encoding": ["enc3", "enc2", "enc1"]})
        mock_decode.assert_any_call("initial", "enc1")
        mock_decode.assert_any_call("stage1", "enc2")
        mock_decode.assert_any_call("stage2", "enc3")
        self.assertEqual(column, ["1", "abcdefghi", "line 1\line 2\nline3"])
    

    @patch("atomium.bcif.decode")
    def test_can_get_data_with_mask(self, mock_decode):
        mock_decode.side_effect = ["stage1", "stage2", [1, "abc\ndef\nghi", ""]]
        column = parse_column_data({
            b"data": "initial", b"encoding": ["enc3", "enc2", "enc1"]
        }, masks=["", "", "?"])
        mock_decode.assert_any_call("initial", "enc1")
        mock_decode.assert_any_call("stage1", "enc2")
        mock_decode.assert_any_call("stage2", "enc3")
        self.assertEqual(column, ["1", "abcdefghi", "?"])



class DataDecodingTests(TestCase):

    @patch("atomium.bcif.decode_byte_array")
    def test_can_decode_byte_array(self, mock_decode):
        self.assertEqual(
            decode(b"abc", {b"kind": b"ByteArray"}),
            mock_decode.return_value
        )
    

    @patch("atomium.bcif.decode_fixed_point")
    def test_can_decode_fixed_point(self, mock_decode):
        self.assertEqual(
            decode(b"abc", {b"kind": b"FixedPoint"}),
            mock_decode.return_value
        )
    

    @patch("atomium.bcif.decode_interval_quantization")
    def test_can_decode_interval_quantization(self, mock_decode):
        self.assertEqual(
            decode(b"abc", {b"kind": b"IntervalQuantization"}),
            mock_decode.return_value
        )
    

    @patch("atomium.bcif.decode_run_length")
    def test_can_decode_run_length(self, mock_decode):
        self.assertEqual(
            decode(b"abc", {b"kind": b"RunLength"}),
            mock_decode.return_value
        )
    

    @patch("atomium.bcif.decode_delta")
    def test_can_decode_delta(self, mock_decode):
        self.assertEqual(
            decode(b"abc", {b"kind": b"Delta"}),
            mock_decode.return_value
        )
    

    @patch("atomium.bcif.decode_integer_packing")
    def test_can_decode_integer_packing(self, mock_decode):
        self.assertEqual(
            decode(b"abc", {b"kind": b"IntegerPacking"}),
            mock_decode.return_value
        )
    

    @patch("atomium.bcif.decode_string_array")
    def test_can_decode_string_array(self, mock_decode):
        self.assertEqual(
            decode(b"abc", {b"kind": b"StringArray"}),
            mock_decode.return_value
        )
    

    def test_can_handle_unknown_kind(self):
        self.assertEqual(decode(b"abc", {b"kind": b"XXX"}), b"abc")



class ByteArrayDecodingTests(TestCase):

    def test_int8_decoding(self):
        decoded = decode_byte_array(b"\x7F\x80\x81", {b"type": 1})
        self.assertEqual(decoded, (127, -128, -127))


    def test_int16_decoding(self):
        decoded = decode_byte_array(b"\x80\x01\x81\x82", {b"type": 2})
        self.assertEqual(decoded, (384, -32127))
    

    def test_int32_decoding(self):
        decoded = decode_byte_array(b"\x80\x01\x81\x82", {b"type": 3})
        self.assertEqual(decoded, (-2105474688,))
    

    def test_uint8_decoding(self):
        decoded = decode_byte_array(b"\x7F\x80\x81", {b"type": 4})
        self.assertEqual(decoded, (127, 128, 129))
    

    def test_uint16_decoding(self):
        decoded = decode_byte_array(b"\x80\x01\x81\x82", {b"type": 5})
        self.assertEqual(decoded, (384, 33409))
    

    def test_int32_decoding(self):
        decoded = decode_byte_array(b"\x80\x01\x81\x82", {b"type": 6})
        self.assertEqual(decoded, (2189492608,))
    

    def test_float32_decoding(self):
        decoded = decode_byte_array(b"AbcD", {b"type": 32})
        self.assertEqual(decoded, (909.535217285156,))
    

    def test_ufloat64_decoding(self):
        decoded = decode_byte_array(b"AbcDAbcA", {b"type": 33})
        self.assertEqual(decoded, (10162698.137131812,))



class FixedPointDecodingTests(TestCase):

    def test_fixed_decoding(self):
        decoded = decode_fixed_point([120, 123, 12], {b"factor": 100})
        self.assertEqual(decoded, [1.2, 1.23, 0.12])



class IntervalQuantizationDecodingTests(TestCase):

    def test_can_decode_interval_quantization(self):
        decoded = decode_interval_quantization(
            [0, 0, 1, 2, 2, 1], {b"min": 1, b"max": 2, b"numSteps": 3}
        )
        self.assertEqual(decoded, [1, 1, 1.5, 2, 2, 1.5])



class RunLengthDecodingTests(TestCase):

    def test_can_decode_run_length(self):
        decoded = decode_run_length([1, 3, 2, 1, 3, 2], {b"srcSize": 6})
        self.assertEqual(decoded, [1, 1, 1, 2, 3, 3])



class DeltaDecodingTests(TestCase):

    def test_can_decode_delta(self):
        decoded = decode_delta([0, 3, 2, 1], {b"origin": 1000})
        self.assertEqual(decoded, [1000, 1003, 1005, 1006])



class IntegerPackingDecodingTests(TestCase):

    def test_can_decode_signed_single_byte(self):
        decoded = decode_integer_packing([1, 2, -3, 127, 1, -128, -128, -4], {
            b"isUnsigned": False, b"byteCount": 1
        })
        self.assertEqual(decoded, [1, 2, -3, 128, -260])
    

    def test_can_decode_unsigned_single_byte(self):
        decoded = decode_integer_packing([1, 2, 5, 255, 9, 3], {
            b"isUnsigned": True, b"byteCount": 1
        })
        self.assertEqual(decoded, [1, 2, 5, 264, 3])
    

    def test_can_decode_signed_double_byte(self):
        decoded = decode_integer_packing([1, 2, -3, 127, 1, -128, -4, 32767, 9, -32768, -1], {
            b"isUnsigned": False, b"byteCount": 2
        })
        self.assertEqual(decoded, [1, 2, -3, 127, 1, -128, -4, 32776, -32769])
    

    def test_can_decode_unsigned_double_byte(self):
        decoded = decode_integer_packing([1, 2, 3, 127, 1, 65535, 9], {
            b"isUnsigned": True, b"byteCount": 2
        })
        self.assertEqual(decoded, [1, 2, 3, 127, 1, 65544])



class StringArrayDecodingTests(TestCase):

    def test(self):
        decoded = decode_string_array([0, 1, 0], {
            b"stringData": b"aAB", b"dataEncoding": [],
            b"offsets": [0, 1, 3], b"offsetEncoding": []
        })
        self.assertEqual(decoded, ["a", "AB", "a"])



class MmcifDictToBcifTests(TestCase):

    @patch("atomium.bcif.encode_column")
    def test_can_save_dict(self, mock_enc):
        mmcif = {"entry": [{"id": "1XXX"}], "cat2": [{"a": "b", "c": "d"}, {"a": "1", "c": "2"}]}
        mock_enc.side_effect = ["col1", "col2", "col3"]
        bytestring = mmcif_dict_to_bcif_filestring(mmcif)
        mock_enc.assert_any_call("id", ["1XXX"])
        mock_enc.assert_any_call("a", ["b", "1"])
        mock_enc.assert_any_call("c", ["d", "2"])
        self.assertEqual(bytestring, msgpack.packb({
                b"encoder": b"atomium 2.0.0",
                b"version": b"0.3.0",
                b"dataBlocks": [{
                    b"header": b"1XXX",
                    b"categories": [{
                        b"name": b"_entry",
                        b"rowCount": 1,
                        b"columns": ["col1"],
                    }, {
                        b"name": b"_cat2",
                        b"rowCount": 2,
                        b"columns": ["col2", "col3"],
                    }]
                }]
            })
        )



class ColumnEncodingTests(TestCase):

    @patch("atomium.bcif.encode_delta")
    def test_can_encode_int_column(self, mock_delta):
        mock_delta.return_value = ["data", "encoding"]
        encoding = encode_column("cat1", ["1", "-20", "4444"])
        self.assertEqual(encoding, {
            b"name": b"cat1",
            b"data": {b"encoding": ["encoding"], b"data": "data"},
            b"mask": None
        })
        mock_delta.assert_called_with([1, -20, 4444])
    

    @patch("atomium.bcif.encode_fixed_point")
    @patch("atomium.bcif.encode_delta")
    def test_can_encode_float_column(self, mock_delta, mock_fixed):
        mock_fixed.return_value = ["data1", "encoding1"]
        mock_delta.return_value = ["data2", "encoding2"]
        encoding = encode_column("cat1", ["1.0", "-20.3", "0.999"])
        self.assertEqual(encoding, {
            b"name": b"cat1",
            b"data": {b"encoding": ["encoding1", "encoding2"], b"data": "data2"},
            b"mask": None
        })
        mock_fixed.assert_called_with([1.0, -20.3, 0.999])
        mock_delta.assert_called_with("data1")
    

    @patch("atomium.bcif.encode_string_array")
    def test_can_encode_string_column(self, mock_string):
        mock_string.return_value = ["data", "encoding"]
        encoding = encode_column("cat1", ["1", "-20", "A", "4444"])
        self.assertEqual(encoding, {
            b"name": b"cat1",
            b"data": {b"encoding": ["encoding"], b"data": "data"},
            b"mask": None
        })
        mock_string.assert_called_with(["1", "-20", "A", "4444"])



class ByteArrayEncodingTests(TestCase):

    def test_int8_encoding(self):
        encoded, encoding = encode_byte_array((127, -128, -127))
        self.assertEqual(encoding, {b"type": 1, b"kind": b"ByteArray"})
        self.assertEqual(encoded, b"\x7F\x80\x81")


    def test_int16_encoding(self):
        encoded, encoding = encode_byte_array((384, -32127))
        self.assertEqual(encoding, {b"type": 2, b"kind": b"ByteArray"})
        self.assertEqual(encoded, b"\x80\x01\x81\x82")
    

    def test_int32_encoding(self):
        encoded, encoding = encode_byte_array((-2105474688,))
        self.assertEqual(encoding, {b"type": 3, b"kind": b"ByteArray"})
        self.assertEqual(encoded, b"\x80\x01\x81\x82")
    

    def test_uint8_encoding(self):
        encoded, encoding = encode_byte_array((127, 128, 129))
        self.assertEqual(encoding, {b"type": 4, b"kind": b"ByteArray"})
        self.assertEqual(encoded, b"\x7F\x80\x81")
    

    def test_uint16_encoding(self):
        encoded, encoding = encode_byte_array((384, 33409))
        self.assertEqual(encoding, {b"type": 5, b"kind": b"ByteArray"})
        self.assertEqual(encoded, b"\x80\x01\x81\x82")
    

    def test_int32_encoding(self):
        encoded, encoding = encode_byte_array((2189492608,))
        self.assertEqual(encoding, {b"type": 6, b"kind": b"ByteArray"})
        self.assertEqual(encoded, b"\x80\x01\x81\x82")
    

    def test_float32_encoding(self):
        encoded, encoding = encode_byte_array((909.5352172851562,))
        self.assertEqual(encoding, {b"type": 32, b"kind": b"ByteArray"})
        self.assertEqual(encoded, b"AbcD")
    

    def test_float64_encoding(self):
        encoded, encoding = encode_byte_array((10162698.137131812, 101622423698.137131812))
        self.assertEqual(encoding, {b"type": 33, b"kind": b"ByteArray"})
        self.assertEqual(encoded, b"AbcDAbcA\x1b#\x92 +\xa97B")



class FixedPointEncodingTests(TestCase):

    def test_fixed_encoding(self):
        encoded, encoding = encode_fixed_point([1.2, 1.23, 0.12])
        self.assertEqual(encoding, {b"factor": 100, b"kind": b"FixedPoint", b"srcType": 3})
        self.assertEqual(encoded, [120, 123, 12])



class IntervalQuantizationEncodingTests(TestCase):

    def test_can_encode_interval_quantization(self):
        encoded, encoding = encode_interval_quantization([0.5, 1, 1.5, 2, 3, 1.345], 1, 2, 3)
        self.assertEqual(encoding, {
            b"min": 1, b"max": 2, b"numSteps": 3, b"kind": b"IntervalQuantization", b"srcType": 3
        })
        self.assertEqual(encoded, [0, 0, 1, 2, 2, 1])


    def test_can_encode_interval_quantization_thirds(self):
        encoded, encoding = encode_interval_quantization([3, 2, 4, 8, 5, 4.5], 3, 5, 7)
        self.assertEqual(encoding, {
            b"min": 3, b"max": 5, b"numSteps": 7, b"kind": b"IntervalQuantization", b"srcType": 3
        })
        self.assertEqual(encoded, [0, 0, 3, 6, 6, 4])



class RunLengthEncodingTests(TestCase):

    def test_can_decode_run_length(self):
        encoded, encoding = encode_run_length([1, 1, 1, 2, 3, 3])
        self.assertEqual(encoding, {b"srcSize": 6, b"kind": b"RunLength", b"srcType": 3})
        self.assertEqual(encoded, [1, 3, 2, 1, 3, 2])
    

    def test_can_decode_run_length_lone_final(self):
        encoded, encoding = encode_run_length([1, 1, 1, 2, 3, 3, 4])
        self.assertEqual(encoding, {b"srcSize": 7, b"kind": b"RunLength", b"srcType": 3})
        self.assertEqual(encoded, [1, 3, 2, 1, 3, 2, 4, 1])



class DeltaEncodingTests(TestCase):

    def test_can_decode_delta(self):
        encoded, encoding = encode_delta([1000, 1003, 1005, 1006])
        self.assertEqual(encoding, {b"origin": 1000, b"kind": b"Delta", b"srcType": 3})
        self.assertEqual(encoded, [0, 3, 2, 1])
    

    def test_can_decode_single_value(self):
        encoded, encoding = encode_delta([1000])
        self.assertEqual(encoding, {b"origin": 1000, b"kind": b"Delta", b"srcType": 3})
        self.assertEqual(encoded, [0])



class IntegerPackingEncodingTests(TestCase):

    def test_can_encode_signed_single_byte(self):
        encoded, encoding = encode_integer_packing([1, 2, -3, 128, -260])
        self.assertEqual(encoded, [1, 2, -3, 127, 1, -128, -128, -4])
        self.assertEqual(encoding,{
            b"isUnsigned": False, b"byteCount": 1,
            b"kind": b"IntegerPacking", b"srcType": 3
        })
    

    def test_can_encode_unsigned_single_byte(self):
        encoded, encoding = encode_integer_packing([1, 2, 5, 264, 3])
        self.assertEqual(encoded, [1, 2, 5, 255, 9, 3])
        self.assertEqual(encoding,{
            b"isUnsigned": True, b"byteCount": 1,
            b"kind": b"IntegerPacking", b"srcType": 3
        })
    

    def test_can_encode_signed_double_byte(self):
        encoded, encoding = encode_integer_packing([1, 2, -3, 127, 1, -128, -4, 32776, -32769], byte_count=2)
        self.assertEqual(encoded, [1, 2, -3, 127, 1, -128, -4, 32767, 9, -32768, -1])
        self.assertEqual(encoding,{
            b"isUnsigned": False, b"byteCount": 2,
            b"kind": b"IntegerPacking", b"srcType": 3
        })
    

    def test_can_encode_unsigned_double_byte(self):
        encoded, encoding = encode_integer_packing([1, 2, 3, 127, 1, 65544], byte_count=2)
        self.assertEqual(encoded, [1, 2, 3, 127, 1, 65535, 9])
        self.assertEqual(encoding,{
            b"isUnsigned": True, b"byteCount": 2,
            b"kind": b"IntegerPacking", b"srcType": 3
        })



class StringArrayEncodingTests(TestCase):

    @patch("atomium.bcif.encode_run_length")
    def test_can_encode_string_array(self, mock_run):
        mock_run.return_value = ("encoded_indices", "indices_encoding")
        encoded, encoding = encode_string_array(["a", "AB", "a"])
        self.assertEqual(encoding, {
            b"stringData": b"aAB", b"dataEncoding": ["indices_encoding"], b"kind": b"StringArray",
            b"offsets": [0, 1, 3], b"offsetEncoding": []
        })
        self.assertEqual(encoded, "encoded_indices")


    @patch("atomium.bcif.encode_run_length")
    def test_can_encode_longer_list(self, mock_run):
        strings = ["sam", "joe", "sam", "ABC", "a", "xxxx", "a"]
        mock_run.return_value = ("encoded_indices", "indices_encoding")
        encoded, encoding = encode_string_array(strings)
        self.assertEqual(encoding, {
            b"stringData": b"samjoeABCaxxxx", b"dataEncoding": ["indices_encoding"], b"kind": b"StringArray",
            b"offsets": [0, 3, 6, 9, 10, 14], b"offsetEncoding": []
        })
        mock_run.assert_called_with([0, 1, 0, 2, 3, 4, 3])
        self.assertEqual(encoded, "encoded_indices")
    

    @patch("atomium.bcif.encode_run_length")
    def test_can_encode_unicode_characters(self, mock_run):
        mock_run.return_value = ("encoded_indices", "indices_encoding")
        encoded, encoding = encode_string_array(["a", "A😬B", "a"])
        self.assertEqual(encoding, {
            b"stringData": b"aA\xf0\x9f\x98\xacB", b"dataEncoding": ["indices_encoding"], b"kind": b"StringArray",
            b"offsets": [0, 1, 7], b"offsetEncoding": []
        })
        self.assertEqual(encoded, "encoded_indices")