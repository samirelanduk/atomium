import re
import struct
import msgpack

def bcif_string_to_mmcif_dict(filestring):
    """Converts a BinaryCif bytestring to an mmCIF dictionary.
    
    :param bytes filestring: the raw BinaryCif filestring.
    :rtype: ``dict``"""
    
    data = msgpack.unpackb(filestring)
    categories = data[b"dataBlocks"][0][b"categories"]
    return {
        category[b"name"][1:].decode(): category_to_table(category)
        for category in categories
    }


def category_to_table(category):
    """Converts a category from its compressed bcif representation to a full
    mmCIF list of rows.
    
    :param dict category: the category data to unpack:
    :rtype: ``list``"""

    masks = [parse_column_data(col[b"mask"]) if col[b"mask"] else [""]
        for col in category[b"columns"]]       
    masks = [[
        {"2": "?", "1": "."}.get(m, "") for m in mask_values
    ] if mask_values[0] != "" else None for mask_values in masks]
    columns = [parse_column_data(col[b"data"], mask)
        for col, mask in zip(category[b"columns"], masks)]
    return [{
        col[b"name"].decode(): col_data[n]
        for col, col_data in zip(category[b"columns"], columns)
    } for n in range(category[b"rowCount"])]


def parse_column_data(column_data, masks=None):
    """Decodes column data into a list of strings. You can provide a list of
    mask values to replace any empty values. The column should contain the
    encoded data, and a list of encodings to undo.
    
    :param dict column_data: a description of a column.
    :param list masks: the empty values.
    :rtype: ``list``"""
    
    data = column_data[b"data"]
    for encoding in column_data[b"encoding"][::-1]:
        data = decode(data, encoding)
    data = list(data)
    for n in range(len(data)):
        value = str(data[n])
        if "\n" in value:
            if " " in value:
                data[n] = "\n".join([l.strip() for l in value.splitlines()])
            else:
                data[n] = value.replace("\n", "")
        else:
            data[n] = value
        if value == "" and masks: data[n] = masks[n]
    return list(data)


def decode(data, encoding):
    """Decodes some encoded data using a description of the encoding that was
    used to create it.
    
    :param data: the encoded data - usually a list, can be a bytestring.
    :param dict encoding: a description of the encoding used.
    :rtype: ``list``"""

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
    """Decodes a bytestring back into the original list of integers that encode
    it, using one of a variety of encodings.
    
    :param list data: the integers to unpack.
    :param dict encoding: a description of the encoding used.
    :rtype: ``list``"""

    buffer, bytesize = {
        1: ["b", 1], 2: ["h", 2], 3: ["i", 4], 4: ["B", 1],
        5: ["H", 2], 6: ["I", 4], 32: ["f", 4], 33: ["d", 8]
    }[encoding[b"type"]]
    buffer = buffer * int(len(data) / bytesize)
    return tuple(round(i, 13) for i in struct.unpack(buffer, data))


def decode_fixed_point(data, encoding):
    """Decodes a list of integers back into the list of floats it was before
    encoding.
    
    :param list data: the integers to unpack.
    :param dict encoding: a description of the encoding used.
    :rtype: ``list``"""

    return [n / encoding[b"factor"] for n in data]


def decode_interval_quantization(data, encoding):
    """Decodes a list of interal indexes back into the original values, with a
    loss of some of the original precision.
    
    :param list data: the integers to decode.
    :param dict encoding: a description of the encoding used.
    :rtype: ``list``"""

    delta = float(
        encoding[b"max"] - encoding[b"min"]) / float(encoding[b"numSteps"] - 1
    )
    return [encoding[b"min"] + delta * value for value in data]


def decode_run_length(data, encoding):
    """Decodes a list of value, count pairs back into the original list of
    numbers.
    
    :param list data: the integers to decode.
    :param dict encoding: a description of the encoding used.
    :rtype: ``list``"""

    values, counts = data[:-1:2], data[1::2]
    decoded = []
    for value, count in zip(values, counts):
        decoded += [value] * count
    return decoded


def decode_delta(data, encoding):
    """Decodes a list of nteger differences back into the original list of
    numbers.
    
    :param list data: the integers to decode.
    :param dict encoding: a description of the encoding used.
    :rtype: ``list``"""

    values = [encoding[b"origin"]]
    for diff in data[1:]:
        values.append(round(values[-1] + diff))
    return values


def decode_integer_packing(data, encoding):
    """Decodes integer packed integers into the original list, which can contain
    large numbers as a result.
    
    :param list data: the integers to unpack.
    :param dict encoding: a description of the encoding used.
    :rtype: ``list``"""

    i, column, integers, unsigned = 0, list(data), [], encoding[b"isUnsigned"]
    if unsigned:
        upper = (0x100 ** encoding[b"byteCount"]) - 1
    else:
        upper = ((0x100 ** encoding[b"byteCount"]) // 2) - 1
        lower = -(0x100 ** encoding[b"byteCount"]) // 2
    while i < len(column):
        value = 0
        t_value = column[i]
        while t_value == upper or t_value == (upper if unsigned else lower):
            value += t_value
            i += 1
            t_value = column[i]
        integers.append(value + t_value)
        i += 1
    return integers


def decode_string_array(data, encoding):
    """Decodes a list of indices into a long concatenated string back into a
    list of strings.
    
    :param list data: the integers to unpack.
    :param dict encoding: a description of the encoding used.
    :rtype: ``list``"""

    for sub_encoding in encoding[b"dataEncoding"][::-1]:
        data = decode(data, sub_encoding)
    offsets = encoding[b"offsets"]
    for sub_encoding in encoding[b"offsetEncoding"][::-1]:
        offsets = decode(offsets, sub_encoding)
    unique_strings = []
    for start, end in zip(offsets[:-1], offsets[1:]):
        unique_strings.append(encoding[b"stringData"][start:end].decode())
    return [unique_strings[n] if n >= 0 else None for n in data]


def save_mmcif_dict(mmcif_dict, path):
    """Saves an mmCIF dictionary to a .bcif file.

    :param dict mmcif_dict: the dictionary to save.
    :param Path path: the location to save to."""


    from atomium import __version__
    code = mmcif_dict["entry"][0]["id"]
    bcif = {
        b"encoder": b"atomium " + __version__.encode(),
        b"version": b"0.3.0",
        b"dataBlocks": [{
            b"header": code.encode(),
            b"categories": [{
                b"name": b"_" + name.encode(),
                b"rowCount": len(rows),
                b"columns": [encode_column(
                    key, [row[key] for row in rows]
                ) for key in rows[0].keys()],
            } for name, rows in mmcif_dict.items()]
        }]
    }
    data = msgpack.packb(bcif)
    with open(path, "wb") as f:
        f.write(data)


def encode_column(name, values):
    """Encodes a column of values, with the specific encoding used depending on
    the kinds of data provided. The aim is to find the most efficient method
    possible.
    
    :param str name: the name of the column.
    :param list values: the list of string values to encode.
    :rtype: ``dict``"""

    encodings = []
    if all(v.lstrip("-").isdigit() for v in values):
        data, encoding = encode_delta([int(v) for v in values])
        encodings.append(encoding)
    elif all(re.search(r"^[+-]?[0-9]+\.[0-9]+$", v) for v in values):
        data, encoding = encode_fixed_point([float(v) for v in values])
        encodings.append(encoding)
        data, encoding = encode_delta(data)
        encodings.append(encoding)
    else:
        data, encoding = encode_string_array(values)
        encodings.append(encoding)
    return {
        b"name": name.encode(), b"mask": None,
        b"data": {b"encoding": encodings, b"data": data},
    }


def encode_byte_array(values):
    """Encodes a list of integers as a bytestring, using one of a variety of 
    encodings.
    
    :param list values: a list of integers.
    :rtype: ``tuple``"""

    is_float = isinstance(values[0], float)
    largest, smallest = max(values), min(values)
    if is_float:
        type = 33 if abs(largest) >= 2 ** 24 else 32
    else:
        if any(value < 0 for value in values):
            if smallest >= -128 and largest < 128:
                type = 1
            else:
                type = 2 if smallest >= -32768 and largest < 32768 else 3
        else:
            type = 4 if largest <= 255 else 5 if largest <= 65535 else 6
    encoding = {b"type": type, b"kind": b"ByteArray"}
    buffer = ".bhiBHI.........................fd"[type]
    return struct.pack(buffer * len(values), *values), encoding
    

def encode_fixed_point(values):
    """Encodes a list of floats as a list of integers, using the smallest
    possible multiple.
    
    :param list values: a list of floats.
    :rtype: ``tuple``"""

    dp = max([len(str(v).split(".")[1]) for v in values])
    encoding = {b"factor": 10 ** dp, b"kind": b"FixedPoint", b"srcType": 3}
    return [round(n * encoding[b"factor"]) for n in values], encoding


def encode_interval_quantization(values, min, max, steps):
    """Encodes a list of numbers as a list of integers, where each integer
    refers to the index of some evenly spaced points between two bounds. It is a
    lossy way of encoding high precision values.
    
    :param list values: a list of numbers.
    :rtype: ``tuple``"""

    encoding = {
        b"min": min, b"max": max, b"numSteps": steps,
        b"kind": b"IntervalQuantization", b"srcType": 3
    }
    step = (max - min) / (steps - 1)
    encoded = []
    for value in values:
        if value <= min:
            encoded.append(0)
            continue
        for i in range(steps):
            point = min + (i * step)
            if (point - step / 2) <= value <= (point + step / 2):
                encoded.append(i)
                break
        else:
            encoded.append(i)
    return encoded, encoding


def encode_run_length(values):
    """Encodes a list of numbers as a list of integers, where each pair is a
    value and the number of times it should be repeated.
    
    :param list values: a list of numbers.
    :rtype: ``tuple``"""

    encoding = {b"srcSize": len(values), b"kind": b"RunLength", b"srcType": 3}
    encoded, n = [], 0
    while True:
        value = values[n]
        count = 1
        while n < len(values) - 1 and values[n + 1] == value:
            count += 1
            n += 1
        encoded += [value, count]
        if n == len(values) - 1: break
        n += 1
    return encoded, encoding


def encode_delta(values):
    """Takes a list of numbers and encodes them as a list of differences.
    
    :param list values: a list of integers.
    :rtype: ``tuple``"""

    encoding = {b"kind": b"Delta", b"origin": values[0], b"srcType": 3}
    result = [0]
    for val1, val2 in zip(values[:-1], values[1:]):
        result.append(round(val2 - val1))
    return result, encoding


def encode_integer_packing(values, byte_count=1):
    """Takes a list of integers and encodes them with the integr packing method,
    breaking larger numbers into smaller numbers to sum. The relevant encoding
    object describing it is produced also.
    
    :param list values: a list of integers.
    :param int byte_count: the number of bytes the output numbers should have.
    :rtype: ``tuple``"""

    signed = any(v < 0 for v in values)
    if signed:
        upper = ((0x100 ** byte_count) // 2) - 1
        lower = -(0x100 ** byte_count) // 2
    else:
        upper = (0x100 ** byte_count) - 1
    encoding = {
        b"kind": b"IntegerPacking", b"byteCount": byte_count,
        b"srcType": 3, b"isUnsigned": not signed
    }
    result = []
    for value in values:
        if value <= upper and (not signed or lower <= value):
            result.append(value)
        elif value > 0:
            mult, remain = divmod(value, upper)
            result += ([upper] * mult)
            result.append(remain)
        elif signed and value < 0:
            mult, remain = divmod(value, lower)
            result += ([lower] * mult)
            result.append(remain)
    return result, encoding


def encode_string_array(values):
    """Takes a list of strings and encodes them as indices within a long string.
    
    :param list values: a list of strings.
    :rtype: ``tuple``"""

    strings, offsets, indices = [], [0], []
    for value in values:
        if value in strings:
            indices.append(strings.index(value))
        else:
            strings.append(value)
            offsets.append(offsets[-1] + len(value.encode()))
            indices.append(len(strings) - 1)
    indices, indices_encoding = encode_run_length(indices)
    encoding = {
        b"stringData": "".join(strings).encode(),
        b"dataEncoding": [indices_encoding], b"kind": b"StringArray",
        b"offsets": offsets, b"offsetEncoding": []
    }
    return indices, encoding