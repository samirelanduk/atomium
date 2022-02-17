import struct
import msgpack

def bcif_string_to_mmcif_dict(filestring):
    data = msgpack.unpackb(filestring)
    categories = data[b"dataBlocks"][0][b"categories"]
    return {
        category[b"name"][1:].decode(): category_to_table(category)
        for category in categories
    }


def category_to_table(category):
    masks = [parse_column_data(col[b"mask"]) if col[b"mask"] else [""]
        for col in category[b"columns"]]   
    masks = [{"2": "?", "1": "."}.get(mask[0], "") for mask in masks]
    columns = [parse_column_data(col[b"data"], mask)
        for col, mask in zip(category[b"columns"], masks)]        
    return [{
        col[b"name"].decode(): col_data[n]
        for col, col_data in zip(category[b"columns"], columns)
    } for n in range(category[b"rowCount"])]


def parse_column_data(column_data, mask=""):
    data = column_data[b"data"]
    for encoding in column_data[b"encoding"][::-1]:
        data = decode(data, encoding)
    data = list(data)
    for n in range(len(data)):
        value = str(data[n])
        if "\n" in value:
            if " " not in value:
                data[n] = value.replace("\n", "")
            else:
                data[n] = "\n".join([
                    line.strip() for line in value.splitlines()
                ])
        else:
            data[n] = value
        if value == "": data[n] = mask
    return list(data)


def decode(data, encoding):
    if encoding[b"kind"] == b"ByteArray":
        return decode_byte_array(data, encoding)
    if encoding[b"kind"] == b"FixedPoint":
        return decode_fixed_point(data, encoding)
    if encoding[b"kind"] == b"IntervalQuantization":
        return decode_interval_quantization(data, encoding)
    if encoding[b"kind"] == b"RunLength":
        return decode_run_length(data, encoding)
    if encoding[b"kind"] == b"Delta":
        return decode_delta(data, encoding)
    if encoding[b"kind"] == b"IntegerPacking":
        return decode_integer_packing(data, encoding)
    if encoding[b"kind"] == b"StringArray":
        return decode_string_array(data, encoding)
    return data


def decode_byte_array(data, encoding):
    buffer, bytesize = {
        1: ["b", 1], 2: ["h", 2], 3: ["i", 4], 4: ["B", 1],
        5: ["H", 2], 6: ["I", 4], 32: ["f", 4], 33: ["d", 8]
    }[encoding[b"type"]]
    buffer = buffer * int(len(data) / bytesize)
    return struct.unpack(buffer, data)


def decode_fixed_point(data, encoding):
    return [n / encoding[b"factor"] for n in data]


def decode_interval_quantization(data, encoding):
    delta = float(encoding[b"max"] - encoding[b"min"]) / float(encoding[b"numSteps"] - 1.0)
    minVal = encoding[b"min"]
    values = []
    for val in data:
        values.append(minVal + delta * val)
    return values


def decode_run_length(data, encoding):
    values, counts = data[:-1:2], data[1::2]
    decoded = []
    for value, count in zip(values, counts):
        decoded += [value] * count
    return decoded


def decode_delta(data, encoding):
    values = [encoding[b"origin"]]
    for diff in data[1:]:
        values.append(values[-1] + diff)
    return values


def decode_integer_packing(data, encoding):

    if encoding[b"isUnsigned"]:
        upperLimit = 0xFF if encoding[b"byteCount"] == 1 else 0xFFFF
        i = 0
        colDataList = list(data)
        integers = []
        while i < len(colDataList):
            value = 0
            tVal = colDataList[i]
            while tVal == upperLimit:
                value += tVal
                i += 1
                tVal = colDataList[i]
            integers.append(value + tVal)
            i += 1
        assert len(integers) == encoding[b"srcSize"]
        return integers
    else:
        upperLimit = 0x7F if encoding[b"byteCount"] == 1 else 0x7FFF
        lowerLimit = -0x7F - 1 if encoding[b"byteCount"] == 1 else -0x7FFF - 1
        i = 0
        colDataList = list(data)
        integers = []
        while i < len(colDataList):
            value = 0
            tVal = colDataList[i]
            while tVal == upperLimit or tVal == lowerLimit:
                value += tVal
                i += 1
                tVal = colDataList[i]
            integers.append(value + tVal)
            i += 1
        assert len(integers) == encoding[b"srcSize"]
        return integers


def decode_string_array(data, encoding):
    for sub_encoding in encoding[b"dataEncoding"][::-1]:
        data = decode(data, sub_encoding)
    offsets = encoding[b"offsets"]
    for sub_encoding in encoding[b"offsetEncoding"][::-1]:
        offsets = decode(offsets, sub_encoding)
    unique_strings = []
    for start, end in zip(offsets, offsets[1:]):
        unique_strings.append(encoding[b"stringData"][start:end].decode())
    return [unique_strings[n] if n >= 0 else None for n in data]