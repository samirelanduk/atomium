import msgpack
import struct
from atomium.data import FULL_NAMES, CODES

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
    values, .mmtf binary values, or other weirdness, and returns a fully decoded
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
    parse_mmtf_compounds(mmtf_dict, mmcif_dict)
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
        get = lambda k: str(round(float(mmtf_dict[k]), 3)) if k in mmtf_dict and mmtf_dict[k] != "?" else "?"
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
            "length_a": str(round(float(mmtf_dict["unitCell"][0]), 3)),
            "length_b": str(round(float(mmtf_dict["unitCell"][1]), 3)),
            "length_c": str(round(float(mmtf_dict["unitCell"][2]), 3)),
            "length_alpha": str(round(float(mmtf_dict["unitCell"][3]), 3)),
            "length_beta": str(round(float(mmtf_dict["unitCell"][4]), 3)),
            "length_gamma": str(round(float(mmtf_dict["unitCell"][5]), 3))
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
    mmcif_dict["struct_asym"] = []
    for entity in mmtf_dict["entityList"]:
        mmcif_dict["entity"].append({
            "id": str(len(mmcif_dict["entity"]) + 1),
            "type": entity["type"],
            "pdbx_description": entity["description"],
            "pdbx_number_of_molecules": str(len(entity["chainIndexList"]))
        })
        for id in entity["chainIndexList"]:
            mmcif_dict["struct_asym"].append({
                "id": mmtf_dict["chainIdList"][id], "pdbx_blank_PDB_chainid_flag": "N",
                "pdbx_modified": "N", "entity_id": mmcif_dict["entity"][-1]["id"],
                "details": "?"
            })
        mmcif_dict["struct_asym"].sort(key=lambda row: mmtf_dict["chainIdList"].index(row["id"]))
        if entity["type"] == "polymer":
            if "entity_poly" not in mmcif_dict: mmcif_dict["entity_poly"] = []
            mmcif_dict["entity_poly"].append({
                "entity_id": str(len(mmcif_dict["entity"])),
                "pdbx_seq_one_letter_code": entity["sequence"],
                "pdbx_seq_one_letter_code_can": entity["sequence"],
            })


def parse_mmtf_compounds(mmtf_dict, mmcif_dict):
    names, mmcif_dict["chem_comp"] = [], []
    for group in mmtf_dict["groupList"]:
        if group["groupName"] in names: continue
        mmcif_dict["chem_comp"].append({
            "id": group["groupName"],
            "type": group["chemCompType"],
            "mon_nstd_flag": "y" if "peptide linking" in group["chemCompType"].lower() else ".",
            "name": FULL_NAMES.get(group["groupName"], "?"),
            "pdbx_synonyms": "?",
            "formula": "?",
            "formula_weight": "?",
        })
        names.append(group["groupName"])
    mmcif_dict["chem_comp"].sort(key=lambda c: c["id"])


def parse_mmtf_models(mmtf_dict, mmcif_dict):
    mmcif_dict["atom_type"] = []
    mmcif_dict["atom_site"] = []
    mmcif_dict["struct_conf"] = []
    mmcif_dict["struct_sheet_range"] = []
    atoms = parse_atoms(mmtf_dict)
    for m_num, c_count in enumerate(mmtf_dict["chainsPerModel"], start=1):
        for c_index, g_count in enumerate(mmtf_dict["groupsPerChain"][:c_count]):
            parse_chain(mmtf_dict, mmcif_dict, m_num, c_index, g_count, atoms)
    for name in sorted(set([a["type_symbol"] for a in mmcif_dict["atom_site"]])):
        mmcif_dict["atom_type"].append({"symbol": name})
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


def mmcif_dict_to_mmtf_filestring(mmcif_dict):
    """Converts an mmCIF dictionary to a .mmtf filestring.

    :param dict mmcif_dict: the dictionary to save.
    :rtype: ``bytes``"""

    # Make basic dict
    mmtf = make_base_mmtf_dict(mmcif_dict)

    # Get 'chain' info
    chain_ids, chains_per_model, chain_names = get_chain_info(mmcif_dict)
    mmtf["chainIdList"] = chain_ids
    mmtf["chainNameList"] = chain_names
    mmtf["chainsPerModel"] = chains_per_model

    # Add assembly info
    mmtf["bioAssemblyList"] = make_bio_assembly_list(mmcif_dict, chain_ids)

    # Add entity info
    mmtf["entityList"] = make_entities(mmcif_dict, chain_ids)

    # Need to:
    # (1) Get the 7 atom lists
    # (2) Get the 3 group lists
    # (3) Get the group type list
    # (4) Get groupsPerChain
    atom_ids, atom_x, atom_y, atom_z, atom_alts, atom_occupancy, atom_b = (
        [], [], [], [], [], [], []
    )
    group_ids, group_types, group_ins = [], [], []
    group_atoms = []
    groups = []
    groups_per_chain = []
    group_count = 0
    for i, atom in enumerate(mmcif_dict["atom_site"]):
        # Get atom values
        atom_ids.append("" if atom["id"] in ".?" else int(atom["id"]))
        atom_x.append("" if atom["Cartn_x"] in ".?" else float(atom["Cartn_x"]))
        atom_y.append("" if atom["Cartn_y"] in ".?" else float(atom["Cartn_y"]))
        atom_z.append("" if atom["Cartn_z"] in ".?" else float(atom["Cartn_z"]))
        atom_alts.append("" if atom["label_alt_id"] in ".?" else atom["label_alt_id"])
        atom_occupancy.append("" if atom["occupancy"] in ".?" else float(atom["occupancy"]))
        atom_b.append("" if atom["B_iso_or_equiv"] in ".?" else float(atom["B_iso_or_equiv"]))

        # Add atom to group_atoms
        group_atoms.append(atom)

        # Is this the last atom in the group?
        this_sig = atom_sig(atom)
        last_atom = i == len(mmcif_dict["atom_site"]) - 1
        next_sig = None if last_atom else atom_sig(mmcif_dict["atom_site"][i + 1])
        is_last = last_atom or this_sig != next_sig

        if is_last:
            # Make new group
            group = make_mmtf_group(group_atoms, mmcif_dict)
            for group_index, existing in enumerate(groups):
                if group["atomNameList"] == existing["atomNameList"] and \
                    group["elementList"] == existing["elementList"] and \
                        group["formalChargeList"] == existing["formalChargeList"]:
                    break
            else:
                groups.append(group)
                group_index = len(groups) - 1
            
            # Get group values
            group_ids.append(this_sig[1])
            group_ins.append(this_sig[2])
            group_types.append(group_index)

            # Last of chain?
            group_count += 1
            if last_atom or this_sig[0] != next_sig[0]:
                groups_per_chain.append(group_count)
                group_count = 0



            # Reset group
            group_atoms = []
            

    
    # Update atom lists
    mmtf["atomIdList"] = atom_ids
    mmtf["xCoordList"] = atom_x
    mmtf["yCoordList"] = atom_y
    mmtf["zCoordList"] = atom_z
    mmtf["altLocList"] = atom_alts
    mmtf["occupancyList"] = atom_occupancy
    mmtf["bFactorList"] = atom_b
    mmtf["numAtoms"] = len(atom_ids)

    # Update group lists
    mmtf["groupIdList"] = group_ids
    mmtf["insCodeList"] = group_ins
    mmtf["groupTypeList"] = group_types
    mmtf["secStructList"] = [-1 for _ in group_ids]
    mmtf["groupList"] = groups
    mmtf["groupsPerChain"] = groups_per_chain
    mmtf["numGroups"] = len(group_ids)
    mmtf["numChains"] = len(chain_ids)
    mmtf["numModels"] = len(chains_per_model)



    # Save
    return msgpack.packb(mmtf)


def make_base_mmtf_dict(mmcif):
    from atomium import __version__
    code = mmcif["entry"][0]["id"]
    get = lambda t, c: mmcif[t][0][c] if t in mmcif and c in mmcif[t][0] else ""
    return {
        "mmtfVersion": b"1.0.0",
        "mmtfProducer": b"atomium " + __version__.encode(),
        "structureId": code,
        "title": mmcif["struct"][0]["title"],
        "depositionDate": get("pdbx_database_status", "recvd_initial_deposition_date"),
        "experimentalMethods":  [get("exptl", "method")],
        "spaceGroup":  get("symmetry", "space_group_name_H-M"),
        "resolution": get("refine", "ls_d_res_high"),
        "rFree": get("refine", "ls_R_factor_R_free"), 
        "rWork": get("refine", "ls_R_factor_R_work"),
        "unitCell": [
            get("cell", "length_a"),
            get("cell", "length_b"),
            get("cell", "length_c"),
            get("cell", "length_alpha"),
            get("cell", "length_beta"),
            get("cell", "length_gamma"),
        ],
        "numBonds": 0,
        "bondAtomList": [],
        "bondOrderList": [],
        "sequenceIndexList": [],
        "ncsOperatorList": [],
    }


def make_bio_assembly_list(mmcif, chain_ids):
    assemblies = []
    if "pdbx_struct_assembly" not in mmcif: return []
    for assembly in mmcif["pdbx_struct_assembly"]:
        gens = [gen for gen in mmcif["pdbx_struct_assembly_gen"]
            if gen["assembly_id"] == assembly["id"]]
        assembly = {"name": assembly["id"], "transformList": []}
        for gen in gens:
            ops = [op for op in mmcif["pdbx_struct_oper_list"]
                if op["id"] in gen["oper_expression"].split(",")]
            for op in ops:
                assembly["transformList"].append({
                    "chainIndexList": [chain_ids.index(n)
                        for n in gen["asym_id_list"].split(",")],
                    "matrix": [
                        float(op["matrix[1][1]"]), float(op["matrix[1][2]"]),
                        float(op["matrix[1][3]"]), float(op["vector[1]"]),
                        float(op["matrix[2][1]"]), float(op["matrix[2][2]"]),
                        float(op["matrix[2][3]"]), float(op["vector[2]"]),
                        float(op["matrix[3][1]"]), float(op["matrix[3][2]"]),
                        float(op["matrix[3][3]"]), float(op["vector[3]"]),
                    ]
                })
        assemblies.append(assembly)
    return assemblies


def make_entities(mmcif, chain_ids):
    entities = []
    for e in mmcif["entity"]:
        index_list = [chain_ids.index(asym["id"])
            for asym in mmcif["struct_asym"] if asym["entity_id"] == e["id"]]
        poly_seqs = [
            pol for pol in mmcif["entity_poly"] if pol["entity_id"] == e["id"]
        ]
        sequence = poly_seqs[0]["pdbx_seq_one_letter_code"] if poly_seqs else ""
        entity = {
            "description": e["pdbx_description"], "type": e["type"],
            "chainIndexList": index_list, "sequence": sequence
        }
        entities.append(entity)
    return entities


def get_chain_info(mmcif):
    chain_ids = [asym["id"] for asym in mmcif["struct_asym"]]
    model_ids = sorted(set(a["pdbx_PDB_model_num"] for a in mmcif["atom_site"]))
    chains_per_model = [len(set(
        a["label_asym_id"] for a in mmcif["atom_site"]
            if a["pdbx_PDB_model_num"] == model
    )) for model in model_ids]
    chain_names = []
    for asym_id in chain_ids:
        for atom in mmcif["atom_site"]:
            if atom["label_asym_id"] == asym_id:
                chain_names.append(atom["auth_asym_id"])
                break
    return chain_ids, chains_per_model, chain_names


def atom_sig(atom):
    return (
        atom["label_asym_id"], atom["auth_seq_id"],
        atom["pdbx_PDB_ins_code"].replace("?", ""), atom["auth_comp_id"]
    )


def make_mmtf_group(atoms, mmcif):
    lookup = {c["id"]: c["type"] for c in mmcif["chem_comp"]}
    names, elements, charges = [], [], []
    for atom in atoms:
        names.append(atom["label_atom_id"])
        elements.append(atom["type_symbol"])
        charges.append(atom["pdbx_formal_charge"])
    group = {
        "groupName": atoms[0]["label_comp_id"],
        "atomNameList": names,
        "elementList": elements,
        "bondOrderList": [],
        "bondAtomList": [],
        "formalChargeList": charges,
        "singleLetterCode": CODES.get(atoms[0]["label_comp_id"], "?"),
        "chemCompType": lookup[atoms[0]["label_comp_id"]]
    }
    return group