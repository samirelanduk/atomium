"""Contains functions for dealing with the .mmtf file format."""

import msgpack
import struct
from collections import deque
from datetime import datetime

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
        return [n / params for n in run_length_decode(integers)]
    elif codec == 10:
        integers = struct.unpack(">" + ("h" * int(len(b[12:]) / 2)), b[12:])
        return [n / params for n in delta_decode(recursive_decode(integers))]
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


def mmtf_dict_to_data_dict(mmtf_dict):
    """Converts an .mmtf dictionary into an atomium data dictionary, with the
    same standard layout that the other file formats get converted into.

    :param dict mmtf_dict: the .mmtf dictionary.
    :rtype: ``dict``"""

    data_dict = {
     "description": {
      "code": None, "title": None, "deposition_date": None,
      "classification": None, "keywords": [], "authors": []
     }, "experiment": {
      "technique": None, "source_organism": None, "expression_system": None
     }, "quality": {"resolution": None, "rvalue": None, "rfree": None},
     "geometry": {"assemblies": []}, "models": []
    }
    mmtf_to_data_transfer(mmtf_dict, data_dict,
     "description", "code", "structureId")
    mmtf_to_data_transfer(mmtf_dict, data_dict,
     "description", "title", "title")
    mmtf_to_data_transfer(mmtf_dict, data_dict,
     "description", "deposition_date", "depositionDate", date=True)
    mmtf_to_data_transfer(mmtf_dict, data_dict,
     "experiment", "technique", "experimentalMethods", first=True)
    mmtf_to_data_transfer(mmtf_dict, data_dict,
     "quality", "resolution", "resolution", trim=3)
    mmtf_to_data_transfer(mmtf_dict, data_dict,
     "quality", "rvalue", "rWork", trim=3)
    mmtf_to_data_transfer(mmtf_dict, data_dict,
     "quality", "rfree", "rFree", trim=3)
    data_dict["geometry"]["assemblies"] = [{
     "id": int(a["name"]), "software": None, "delta_energy": None,
     "buried_surface_area": None, "surface_area": None, "transformations": [{
      "chains": [mmtf_dict["chainIdList"][i] for i in t["chainIndexList"]],
      "matrix": [t["matrix"][n * 4: (n * 4) + 3] for n in range(3)],
      "vector": t["matrix"][3:-4:4]} for t in a["transformList"]
     ]
    } for a in mmtf_dict.get("bioAssemblyList", [{"transformList": [{}]}])]
    update_models_list(mmtf_dict, data_dict)
    return data_dict


def update_models_list(mmtf_dict, data_dict):
    """Takes a data dictionary and updates its models list with
    information from a .mmtf dictionary.

    :param dict mmtf_dict: the .mmtf dictionary to read.
    :param dict data_dict: the data dictionary to update."""

    gc, atom_count = 0, 0
    atoms = get_atoms_list(mmtf_dict)
    for model_num in range(mmtf_dict["numModels"]):
        model = {"polymer": {}, "non-polymer": {}, "water": {}}
        for chain_num in range(mmtf_dict["chainsPerModel"][model_num]):
            type_ = get_type(mmtf_dict, chain_num)
            if type_ == "polymer": chain = make_chain(mmtf_dict, chain_num)
            for g in range(gc, gc + mmtf_dict["groupsPerChain"][chain_num]):
                group = mmtf_dict["groupList"][mmtf_dict["groupTypeList"][g]]
                group_id = make_group_id(mmtf_dict, chain_num, gc)
                group_name = group["groupName"]
                mol = make_molecule(type_, group_name, mmtf_dict, chain_num)
                atoms_in_group = len(group["atomNameList"])
                for atom_num in range(atoms_in_group):
                    add_atom(mol, atoms, group, atom_num, atom_count)
                atom_count += atoms_in_group
                d = chain["residues"] if type_ == "polymer" else model[type_]
                if type_ == "polymer": mol["number"] = len(d) + 1
                d[group_id] = mol
                gc += 1
            if type_ == "polymer":
                model["polymer"][mmtf_dict["chainNameList"][chain_num]] = chain
        data_dict["models"].append(model)


def get_atoms_list(mmtf_dict):
    """Creates a list of atom dictionaries from a .mmtf dictionary by zipping
    together some of its fields.

    :param dict mmtf_dict: the .mmtf dictionary to read.
    :rtype: ``list``"""

    return [{
     "x": x, "y": y, "z": z, "alt_loc": a or None, "bvalue": b, "occupancy": o,
     "id": i
    } for x, y, z, a, b, i, o in zip(
     mmtf_dict["xCoordList"], mmtf_dict["yCoordList"], mmtf_dict["zCoordList"],
     mmtf_dict["altLocList"], mmtf_dict["bFactorList"], mmtf_dict["atomIdList"],
     mmtf_dict["occupancyList"]
    )]


def get_type(mmtf_dict, chain_num):
    """Takes an .mmtf dictionary and a chain number, and works out what type of
    molecule it corresponds to.

    :param dict mmtf_dict: the .mmtf dictionary to read.
    :param int chain_num: the chain number.
    :rtype: ``str``"""

    for entity in mmtf_dict["entityList"]:
        if chain_num in entity["chainIndexList"]:
            return entity["type"]


def make_chain(mmtf_dict, chain_num):
    """Creates a chain dictionary

    :param dict mmtf_dict: the .mmtf dictionary to read.
    :param int chain_num: the chain number.
    :rtype: ``dict``"""

    return {
     "internal_id": mmtf_dict["chainIdList"][chain_num],
     "residues": {}, "sequence": get_chain_sequence(mmtf_dict, chain_num)
    }


def get_chain_sequence(mmtf_dict, chain_num):
    """Gets the chain sequence for a given chain number.

    :param dict mmtf_dict: the .mmtf dictionary to read.
    :param int chain_num: the chain number.
    :rtype: ``str``"""

    for entity in mmtf_dict["entityList"]:
        if chain_num in entity["chainIndexList"]:
            return entity["sequence"]
    return ""


def make_group_id(mmtf_dict, chain_num, g_count):
    """Creates an atomium residue/molecule ID for a given group.

    :param dict mmtf_dict: the .mmtf dictionary to read.
    :param int chain_num: the chain number.
    :param int g_count: the group number.
    :rtype: ``str``"""

    return "{}.{}{}".format(
     mmtf_dict["chainNameList"][chain_num], mmtf_dict["groupIdList"][g_count],
     mmtf_dict["insCodeList"][g_count]
    )


def make_molecule(type_, group_name, mmtf_dict, chain_num):
    """Creates a residue or ligand dictionary.

    :param str type_: 'polymer' or 'non-polymer'.
    :param str group_name: the name of the molecule.
    :param dict mmtf_dict: the .mmtf dictionary to read.
    :param int chain_num: the chain number.
    :rtype: ``dict``"""

    return {"name": group_name, "atoms": {}} if type_ == "polymer" else {
     "name": group_name,
     "internal_id": mmtf_dict["chainIdList"][chain_num],
     "polymer": mmtf_dict["chainNameList"][chain_num],
     "atoms": {}
    }


def add_atom(d, atoms, group, atom_num, atom_count):
    """Adds an atom dictionary to a molecule dictionary.

    :param dict d: the molecule dictionary to update.
    :param list atoms: the list of half-created atom dictionaries.
    :param dict group: a .mmtf group dictionary.
    :param int atom_num: the atom number within the group.
    :param int atom_count: the total number of processed atoms."""

    atom = atoms[atom_num + atom_count]
    atom["name"] = group["atomNameList"][atom_num]
    atom["element"] = group["elementList"][atom_num]
    atom["charge"] = group["formalChargeList"][atom_num]
    atom["anisotropy"] = [0, 0, 0, 0, 0, 0]
    d["atoms"][atom["id"]] = atom
    del atom["id"]


def mmtf_to_data_transfer(mmtf_dict, data_dict, d_cat, d_key, m_key,
                           date=False, first=False, trim=False):
    """A function for transfering a bit of data from a .mmtf dictionary to a
    data dictionary, or doing nothing if the data doesn't exist.

    :param dict mmtf_dict: the .mmtf dictionary to read.
    :param dict data_dict: the data dictionary to update.
    :param str d_cat: the top-level key in the data dictionary.
    :param str d_key: the data dictionary field to update.
    :param str m_key: the .mmtf field to read.
    :param bool date: if True, the value will be converted to a date.
    :param bool first: if True, the value's first item will be split used.
    :param int trim: if given, the value will be rounded by this amount."""

    try:
        value = mmtf_dict[m_key]
        if date: value = datetime.strptime(value, "%Y-%m-%d").date()
        if first: value = value[0]
        if trim: value = round(value, trim)
        data_dict[d_cat][d_key] = value
    except: pass
