import msgpack
import struct

def mmtf_string_to_mmcif_dict(bytestring):
    """Takes the raw bytestring of a .mmtf file and turns it into a normal,
    fully decoded JSON dictionary.

    :patam bytes bytestring: the .mmtf filestring.
    :rtype: ``dict``"""

    raw = msgpack.unpackb(bytestring)
    mmtf = decode_dict(raw)
    return mmtf_dict_to_mmcif_dict(mmtf)


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
        new[key.decode() if isinstance(key, bytes) else key] = new_value
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


def mmtf_dict_to_mmcif_dict(mmtf_dict):
    mmcif_dict = {}
    parse_mmtf_header(mmtf_dict, mmcif_dict)
    parse_mmtf_quality(mmtf_dict, mmcif_dict)
    parse_mmtf_crystal(mmtf_dict, mmcif_dict)
    parse_mmtf_assemblies(mmtf_dict, mmcif_dict)
    parse_mmtf_entities(mmtf_dict, mmcif_dict)
    parse_mmtf_models(mmtf_dict, mmcif_dict)
    return mmcif_dict


def parse_mmtf_header(mmtf_dict, mmcif_dict):
    mmcif_dict["entry"] = [{"id": mmtf_dict["structureId"]}]
    mmcif_dict["struct"] = [{
        "entry_id": mmtf_dict["structureId"], "title": mmtf_dict["title"]
    }]
    mmcif_dict["pdbx_database_status"] = [{
        "entry_id": mmtf_dict["structureId"],
        "recvd_initial_deposition_date": mmtf_dict["depositionDate"]
    }]
    mmcif_dict["exptl"] = [{
        "entry_id": mmtf_dict["structureId"],
        "method": ",".join(mmtf_dict["experimentalMethods"])
    }]


def parse_mmtf_quality(mmtf_dict, mmcif_dict):
    if any(k in mmtf_dict for k in ["resolution", "rWork", "rFree"]):
        get = lambda k: str(round(mmtf_dict[k], 3)) if k in mmtf_dict else "?"
        mmcif_dict["refine"] = [{
            "entry_id": mmtf_dict["structureId"],
            "ls_d_res_high": get("resolution"),
            "ls_R_factor_obs": get("rWork"),
            "ls_R_factor_all": get("rWork"),
            "ls_R_factor_R_work": get("rWork"),
            "ls_R_factor_R_free": get("rFree"),
        }]


def parse_mmtf_crystal(mmtf_dict, mmcif_dict):
    mmcif_dict["symmetry"] = [{
        "entry_id": mmtf_dict["structureId"],
        "space_group_name_H-M": mmtf_dict["spaceGroup"]
    }]
    if "unitCell" in mmtf_dict:
        mmcif_dict["cell"] = [{
            "entry_id": mmtf_dict["structureId"],
            "length_a": str(round(mmtf_dict["unitCell"][0], 3)),
            "length_b": str(round(mmtf_dict["unitCell"][1], 3)),
            "length_c": str(round(mmtf_dict["unitCell"][2], 3)),
            "length_alpha": str(round(mmtf_dict["unitCell"][3], 3)),
            "length_beta": str(round(mmtf_dict["unitCell"][4], 3)),
            "length_gamma": str(round(mmtf_dict["unitCell"][5], 3))
        }]


def parse_mmtf_assemblies(mmtf_dict, mmcif_dict):
    if mmtf_dict["bioAssemblyList"]:
        mmcif_dict["pdbx_struct_assembly"] = []
        mmcif_dict["pdbx_struct_assembly_gen"] = []
        mmcif_dict["pdbx_struct_oper_list"] = []
        for assembly in mmtf_dict["bioAssemblyList"]:
            mmcif_dict["pdbx_struct_assembly"].append({"id": assembly["name"]})
            for transformation in assembly["transformList"]:
                op_id = str(len(mmcif_dict["pdbx_struct_oper_list"]) + 1)
                mmcif_dict["pdbx_struct_assembly_gen"].append({
                    "assembly_id": assembly["name"],
                    "oper_expression": op_id,
                    "asym_id_list": ",".join([
                        mmtf_dict["chainIdList"][index]
                        for index in transformation["chainIndexList"]
                    ])
                })
                get = lambda i: f"{round(transformation['matrix'][i], 10):.10f}"
                mmcif_dict["pdbx_struct_oper_list"].append({
                    "id": op_id, "matrix[1][1]": get(0),
                    "matrix[1][2]": get(1), "matrix[1][3]": get(2),
                    "vector[1]": get(3), "matrix[2][1]": get(4),
                    "matrix[2][2]": get(5), "matrix[2][3]": get(6),
                    "vector[2]": get(7), "matrix[3][1]": get(8),
                    "matrix[3][2]": get(9), "matrix[3][3]": get(10),
                    "vector[3]": get(11),
                })


def parse_mmtf_entities(mmtf_dict, mmcif_dict):
    mmcif_dict["entity"] = []
    for entity in mmtf_dict["entityList"]:
        mmcif_dict["entity"].append({
            "id": str(len(mmcif_dict["entity"]) + 1),
            "type": entity["type"],
            "pdbx_description": entity["description"],
            "pdbx_number_of_molecules": str(len(entity["chainIndexList"]))
        })
        if entity["type"] == "polymer":
            if "entity_poly" not in mmcif_dict: mmcif_dict["entity_poly"] = []
            mmcif_dict["entity_poly"].append({
                "entity_id": str(len(mmcif_dict["entity"])),
                "pdbx_seq_one_letter_code": entity["sequence"],
                "pdbx_seq_one_letter_code_can": entity["sequence"],
            })


def parse_mmtf_models(mmtf_dict, mmcif_dict):
    mmcif_dict["atom_site"] = []
    mmcif_dict["struct_conf"] = []
    mmcif_dict["struct_sheet_range"] = []
    atoms = parse_atoms(mmtf_dict)
    for m_num, c_count in enumerate(mmtf_dict["chainsPerModel"], start=1):
        for c_index, g_count in enumerate(mmtf_dict["groupsPerChain"][:c_count]):
            parse_chain(mmtf_dict, mmcif_dict, m_num, c_index, g_count, atoms)
    return mmcif_dict


def parse_atoms(mmtf_dict):
    return list(zip(
        mmtf_dict["atomIdList"], mmtf_dict["xCoordList"],
        mmtf_dict["yCoordList"], mmtf_dict["zCoordList"],
        mmtf_dict["altLocList"], mmtf_dict["occupancyList"],
        mmtf_dict["bFactorList"]
    ))


def parse_chain(mmtf_dict, mmcif_dict, model_num, chain_index, group_count, atoms):
    entity_id = get_chain_entity_id(mmtf_dict, chain_index)
    groups = get_chain_groups(mmtf_dict, group_count)
    for group_type, group_id, insert, _ in groups:
        g = mmtf_dict["groupList"][group_type]
        for atom_name, element, charge in zip(
            g["atomNameList"], g["elementList"], g["formalChargeList"]
        ):
            id, x, y, z, alt, occ, b = atoms.pop(0)
            mmcif_dict["atom_site"].append({
                "group_PDB": "HETATM" if g["singleLetterCode"] == "?" else "ATOM",
                "id": str(id),
                "type_symbol": element,
                "label_atom_id": atom_name,
                "label_alt_id": alt or ".",
                "label_comp_id": g["groupName"],
                "label_asym_id": mmtf_dict["chainIdList"][chain_index],
                "label_entity_id": str(entity_id),
                "label_seq_id": str(group_id),
                "pdbx_PDB_ins_code": insert or "?",
                "Cartn_x": str(x), "Cartn_y": str(y), "Cartn_z": str(z),
                "occupancy": str(occ), "B_iso_or_equiv": str(b),
                "pdbx_formal_charge": str(charge),
                "auth_seq_id": str(group_id),
                "auth_comp_id": g["groupName"],
                "auth_asym_id": mmtf_dict["chainNameList"][chain_index],
                "auth_atom_id": atom_name,
                "pdbx_PDB_model_num": str(model_num)
            })


    helices, strands = get_chain_secondary_structure(mmtf_dict, groups, chain_index)
    add_secondary_structure_to_mmcif(helices, strands, mmcif_dict)
    

def get_chain_entity_id(mmtf_dict, chain_index):
    entity_id = "?"
    for num, entity in enumerate(mmtf_dict["entityList"], start=1):
        if chain_index in entity["chainIndexList"]:
            entity_id = num
            break
    return entity_id


def get_chain_groups(mmtf_dict, group_count):
    group_ids = mmtf_dict["groupIdList"][:group_count]
    group_types = mmtf_dict["groupTypeList"][:group_count]
    insert_codes = mmtf_dict["insCodeList"][:group_count]
    ss_codes = mmtf_dict["secStructList"][:group_count]
    mmtf_dict["groupIdList"] = mmtf_dict["groupIdList"][group_count:]
    mmtf_dict["groupTypeList"] = mmtf_dict["groupTypeList"][group_count:]
    mmtf_dict["insCodeList"] = mmtf_dict["insCodeList"][group_count:]
    mmtf_dict["secStructList"] = mmtf_dict["secStructList"][group_count:]
    return list(zip(group_types, group_ids, insert_codes, ss_codes))


def get_chain_secondary_structure(mmtf_dict, groups, chain_index):
    in_helix, helices, helix = False, [], []
    in_strand, strands, strand = False, [], []
    for group_type, group_id, insert, ss_code in groups:
        group = mmtf_dict["groupList"][group_type]
        ss_type = [
            "helices", None, "helices", "strands",
            "helices", "strands", None, None
        ][ss_code]
        if ss_type == "helices":
            helix.append([
                mmtf_dict["chainIdList"][chain_index],
                mmtf_dict["chainNameList"][chain_index],
                group["groupName"], group_id, insert
            ])
        elif in_helix:
            helices.append(helix)
            helix = []
        in_helix = ss_type == "helices"
        if ss_type == "strands":
            strand.append([
                mmtf_dict["chainIdList"][chain_index],
                mmtf_dict["chainNameList"][chain_index],
                group["groupName"], group_id, insert
            ])
        elif in_strand:
            strands.append(strand)
            strand = []
        in_strand = ss_type == "strands"
    if helix: helices.append(helix)
    if strand: strands.append(strand)
    return helices, strands


def add_secondary_structure_to_mmcif(helices, strands, mmcif_dict):
    for helix in helices:
        mmcif_dict["struct_conf"].append({
            "id": f"HELIX_P{len(mmcif_dict['struct_conf'])}",
            "pdbx_PDB_helix_id": str(len(mmcif_dict["struct_conf"])),
            "beg_label_comp_id": helix[0][2],
            "beg_label_asym_id": helix[0][0],
            "beg_label_seq_id": str(helix[0][3]),
            "pdbx_beg_PDB_ins_code": helix[0][4] or "?",
            "end_label_comp_id": helix[-1][2],
            "end_label_asym_id": helix[-1][0],
            "end_label_seq_id": str(helix[-1][3]),
            "pdbx_end_PDB_ins_code": helix[-1][4] or "?",
            "beg_auth_comp_id": helix[0][2],
            "beg_auth_asym_id": helix[0][1],
            "beg_auth_seq_id": str(helix[0][3]),
            "end_auth_comp_id": helix[-1][2],
            "end_auth_asym_id": helix[-1][1],
            "end_auth_seq_id": str(helix[-1][3]),
            "pdbx_PDB_helix_class": "?",
            "details": "?",
            "pdbx_PDB_helix_length": str(len(helix))
        })
    for strand in strands:
        mmcif_dict["struct_sheet_range"].append({
            "sheet_id": strand[0][1],
            "id": str(len(mmcif_dict["struct_sheet_range"])),
            "beg_label_comp_id": strand[0][2],
            "beg_label_asym_id": strand[0][0],
            "beg_label_seq_id": str(strand[0][3]),
            "pdbx_beg_PDB_ins_code": strand[0][4] or "?",
            "end_label_comp_id": strand[-1][2],
            "end_label_asym_id": strand[-1][0],
            "end_label_seq_id": str(strand[-1][3]),
            "pdbx_end_PDB_ins_code": strand[-1][4] or "?",
            "beg_auth_comp_id": strand[0][2],
            "beg_auth_asym_id": strand[0][1],
            "beg_auth_seq_id": str(strand[0][3]),
            "end_auth_comp_id": strand[-1][2],
            "end_auth_asym_id": strand[-1][1],
            "end_auth_seq_id": str(strand[-1][3]),
        })