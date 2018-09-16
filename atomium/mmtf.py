"""Contains functions for dealing with the .mmtf file format."""

import msgpack
import struct
from collections import deque

def mmtf_bytes_to_mmtf_dict(bytestring):
    """Takes the raw bytestring of a .mmtf file and turns it into a normal,
    fully decoded JSON dictionary.

    :patam bytes bytestring: the .mmtf filestring.
    :rtype: ``dict``"""

    raw = msgpack.unpackb(bytestring)
    return decode_dict(raw)


def decode_dict(d):
    """Takes a dictionary that might have bytestring keys, lists of bytestring
    values, .mmtf binary values, or other weirdness, and returns a dully decoded
    version of it which is a JSON-valid dictionary.

    :param dict d: the dictionary to read.
    :rtype: ``dict``"""

    new = {}
    for key, value in d.items():
        try:
            new_value = value.decode()
        except: new_value = value
        if isinstance(new_value, str) and new_value and new_value[0] == "\x00":
            new_value = new_value.encode()
        if isinstance(new_value, bytes):
            new_value = parse_binary_field(new_value)
        if isinstance(new_value, list) and new_value:
            if isinstance(new_value[0], dict):
                new_value = [decode_dict(x) for x in new_value]
            elif isinstance(new_value[0], bytes):
                new_value = [x.decode() for x in new_value]
        new[key.decode()] = new_value
    return new


def parse_binary_field(b):
    """Some fields in a .mmtf file cannot be unpacked by msgpack and have
    special .mmtf encoding, as specified in its documentation. This function
    takes such a field and decodes it.

    :param bytestring b: the field to parse.
    :returns: the parsed result (type varies)."""

    codec, length, params = struct.unpack(">iii", b[:12])
    len4 = lambda b: int(len(b[12:]) / 4)
    if codec == 1: return struct.unpack("f" * length, b[12:])
    elif codec == 2: return struct.unpack("b" * length, b[12:])
    elif codec == 3: return struct.unpack(">" + "h" * length, b[12:])
    elif codec == 4: return struct.unpack(">" + "i" * length, b[12:])
    elif codec == 5:
        chars = struct.unpack("c" * (length * 4), b[12:])
        return [b"".join([
         c for c in chars[i * 4: (i + 1) * 4] if c != b"\x00"
        ]).decode() for i in range(length)]
    elif codec == 6:
        integers = struct.unpack(">" + ("i" * len4(b)), b[12:])
        return [chr(c) if c != 0 else "" for c in run_length_decode(integers)]
    elif codec == 7:
        integers = struct.unpack(">" + ("i" * len4(b)), b[12:])
        return run_length_decode(integers)
    elif codec == 8:
        integers = struct.unpack(">" + ("i" * len4(b)), b[12:])
        return delta_decode(run_length_decode(integers))
    elif codec == 9:
        integers = struct.unpack(">" + ("i" * len4(b)), b[12:])
        return [n / 100 for n in run_length_decode(integers)]
    elif codec == 10:
        integers = struct.unpack(">" + ("h" * int(len(b[12:]) / 2)), b[12:])
        return [n / 1000 for n in delta_decode(recursive_decode(integers))]
    else: raise ValueError(".mmtf error: {} is invalid codec".format(codec))


def run_length_decode(integers):
    """Expands a list of integers where every second integer is a count of the
    integer before it.

    :param list integers: the integers to decode.
    :rtype: ``list``"""

    x = []
    for index, val in enumerate(integers[::2]):
        x += [val] * integers[1::2][index]
    return x


def delta_decode(integers):
    """Turns a list of integers into a new list of integers where the values in
    the first are treated as deltas to be applied to the previous value.

    :param list integers: the integers to decode.
    :rtype: ``list``"""

    array, last = [], 0
    for i in integers:
        last += i
        array.append(last)
    return array


def recursive_decode(integers, bits=16):
    """Turns a list of integers into a new list of integers where the values in
    the first are merged if it looks like a higher order integer split over two
    integers.

    (Code here adapted from the official python-mmtf package.)

    :param list integers: the integers to decode.
    :rtype: ``list``"""

    new = []
    power = 2 ** (bits - 1)
    cutoff = [power - 1, 0 - power]
    index = 0
    while index < len(integers):
        value = 0
        while integers[index] in cutoff:
            value += integers[index]
            index += 1
            if integers[index] == 0: break
        value += integers[index]
        index += 1
        new.append(value)
    return new
