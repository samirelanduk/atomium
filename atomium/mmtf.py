import msgpack
import struct

def mmtf_string_to_mmcif_dict(bytestring):
    """Takes the raw bytestring of a .mmtf file and turns it into a normal,
    fully decoded JSON dictionary.

    :patam bytes bytestring: the .mmtf filestring.
    :rtype: ``dict``"""

    raw = msgpack.unpackb(bytestring)
    mmtf = decode_dict(raw)
    mmcif_dict = {

    }
    #print(mmtf["entityList"])

    mmcif_dict["entry"] = [{"id": mmtf["structureId"]}]

    mmcif_dict["struct"] = [{"entry_id": mmtf["structureId"], "title": mmtf["title"]}]

    mmcif_dict["pdbx_database_status"] = [{
        "entry_id": mmtf["structureId"], "recvd_initial_deposition_date": mmtf["depositionDate"]
    }]

    mmcif_dict["exptl"] = [{
        "entry_id": mmtf["structureId"], "method": ",".join(mmtf["experimentalMethods"])
    }]

    if "resolution" in mmtf or "rWork" in mmtf or "rFree" in mmtf:
        mmcif_dict["refine"] = [{
            "entry_id": mmtf["structureId"],
            "ls_d_res_high": str(round(mmtf["resolution"], 3)) if "resolution" in mmtf else "?",
            "ls_R_factor_obs":str(round( mmtf["rWork"], 3)) if "rWork" in mmtf else "?",
            "ls_R_factor_all": str(round(mmtf["rWork"], 3)) if "rWork" in mmtf else "?",
            "ls_R_factor_R_work": str(round(mmtf["rWork"], 3)) if "rWork" in mmtf else "?",
            "ls_R_factor_R_free": str(round(mmtf["rFree"], 3)) if "rFree" in mmtf else "?",
        }]

    mmcif_dict["symmetry"] = [{
        "entry_id": mmtf["structureId"],
        "space_group_name_H-M": mmtf["spaceGroup"]
    }]

    if "unitCell" in mmtf:
        mmcif_dict["cell"] = [{
            "entry_id": mmtf["structureId"],
            "length_a": str(round(mmtf["unitCell"][0], 3)),
            "length_b": str(round(mmtf["unitCell"][1], 3)),
            "length_c": str(round(mmtf["unitCell"][2], 3)),
            "length_alpha": str(round(mmtf["unitCell"][3], 3)),
            "length_beta": str(round(mmtf["unitCell"][4], 3)),
            "length_gamma": str(round(mmtf["unitCell"][5], 3))
        }]

    if mmtf["bioAssemblyList"]:
        mmcif_dict["pdbx_struct_assembly"] = []
        mmcif_dict["pdbx_struct_assembly_gen"] = []
        mmcif_dict["pdbx_struct_oper_list"] = []
        for assembly in mmtf["bioAssemblyList"]:
            mmcif_dict["pdbx_struct_assembly"].append({"id": assembly["name"]})
            for transformation in assembly["transformList"]:
                op_id = str(len(mmcif_dict["pdbx_struct_oper_list"]) + 1)
                mmcif_dict["pdbx_struct_assembly_gen"].append({
                    "assembly_id": assembly["name"],
                    "oper_expression": op_id,
                    "asym_id_list": ",".join([mmtf["chainIdList"][index] for index in transformation["chainIndexList"]])
                })
                mmcif_dict["pdbx_struct_oper_list"].append({
                    "id": op_id,
                    "matrix[1][1]": f"{round(transformation['matrix'][0], 10):.10f}",
                    "matrix[1][2]": f"{round(transformation['matrix'][1], 10):.10f}",
                    "matrix[1][3]": f"{round(transformation['matrix'][2], 10):.10f}",
                    "vector[1]": f"{round(transformation['matrix'][3], 10):.10f}",
                    "matrix[2][1]": f"{round(transformation['matrix'][4], 10):.10f}",
                    "matrix[2][2]": f"{round(transformation['matrix'][5], 10):.10f}",
                    "matrix[2][3]": f"{round(transformation['matrix'][6], 10):.10f}",
                    "vector[2]": f"{round(transformation['matrix'][7], 10):.10f}",
                    "matrix[3][1]": f"{round(transformation['matrix'][8], 10):.10f}",
                    "matrix[3][2]": f"{round(transformation['matrix'][9], 10):.10f}",
                    "matrix[3][3]": f"{round(transformation['matrix'][10], 10):.10f}",
                    "vector[3]": f"{round(transformation['matrix'][11], 10):.10f}",
                })
    
    mmcif_dict["entity"] = []
    for entity in mmtf["entityList"]:
        mmcif_dict["entity"].append({
            "id": str(len(mmcif_dict["entity"]) + 1),
            "type": entity["type"],
            "pdbx_description": entity["description"],
            "pdbx_number_of_molecules": str(len(entity["chainIndexList"]))
        })
        if entity["type"] == "polymer":
            if "entity_poly" not in mmcif_dict:
                mmcif_dict["entity_poly"] = []
            mmcif_dict["entity_poly"].append({
                "entity_id": str(len(mmcif_dict["entity"])),
                "pdbx_seq_one_letter_code": entity["sequence"],
                "pdbx_seq_one_letter_code_can": entity["sequence"],
            })
        
    mmcif_dict["atom_site"] = []
    mmcif_dict["struct_conf"] = []
    mmcif_dict["struct_sheet_range"] = []
    atoms = list(zip(
        mmtf["atomIdList"], mmtf["xCoordList"], mmtf["yCoordList"],
        mmtf["zCoordList"], mmtf["altLocList"],
        mmtf["occupancyList"], mmtf["bFactorList"]
    ))
    
    for model_num, chain_count in enumerate(mmtf["chainsPerModel"], start=1):
        for chain_index, group_count in enumerate(mmtf["groupsPerChain"][:chain_count]):

            in_helix = False
            helices = []
            helix = []

            in_strand = False
            strands = []
            strand = []

            entity_id = "?"
            for num, entity in enumerate(mmtf["entityList"], start=1):
                if chain_index in entity["chainIndexList"]:
                    entity_id = num
                    break
            group_ids = mmtf["groupIdList"][:group_count]
            group_types = mmtf["groupTypeList"][:group_count]
            insert_codes = mmtf["insCodeList"][:group_count]
            ss_codes = mmtf["secStructList"][:group_count]
            mmtf["groupIdList"] = mmtf["groupIdList"][group_count:]
            mmtf["groupTypeList"] = mmtf["groupTypeList"][group_count:]
            mmtf["insCodeList"] = mmtf["insCodeList"][group_count:]
            mmtf["secStructList"] = mmtf["secStructList"][group_count:]
            for group_type, group_id, insert, ss_code in zip(group_types, group_ids, insert_codes, ss_codes):
                group = mmtf["groupList"][group_type]
                ss_type = [
                    "helices", None, "helices", "strands", "helices", "strands", None, None
                ][ss_code]

                if ss_type == "helices":
                    helix.append([
                        mmtf["chainIdList"][chain_index], mmtf["chainNameList"][chain_index],
                        group["groupName"], group_id, insert
                    ])
                elif in_helix:
                    helices.append(helix)
                    helix = []
                in_helix = ss_type == "helices"

                if ss_type == "strands":
                    strand.append([
                        mmtf["chainIdList"][chain_index], mmtf["chainNameList"][chain_index],
                        group["groupName"], group_id, insert
                    ])
                elif in_strand:
                    strands.append(strand)
                    strand = []
                in_strand = ss_type == "strands"

                for atom_name, element, charge in zip(group["atomNameList"], group["elementList"], group["formalChargeList"]):
                    id, x, y, z, alt, occ, b = atoms.pop(0)
                    mmcif_dict["atom_site"].append({
                        "group_PDB": "HETATM" if group["singleLetterCode"] == "?" else "ATOM",
                        "id": str(id),
                        "type_symbol": element,
                        "label_atom_id": atom_name,
                        "label_alt_id": alt or ".",
                        "label_comp_id": group["groupName"],
                        "label_asym_id": mmtf["chainIdList"][chain_index],
                        "label_entity_id": str(entity_id),
                        "label_seq_id": str(group_id),
                        "pdbx_PDB_ins_code": insert or "?",
                        "Cartn_x": str(x),
                        "Cartn_y": str(y),
                        "Cartn_z": str(z),
                        "occupancy": str(occ),
                        "B_iso_or_equiv": str(b),
                        "pdbx_formal_charge": str(charge),
                        "auth_seq_id": str(group_id),
                        "auth_comp_id": group["groupName"],
                        "auth_asym_id": mmtf["chainNameList"][chain_index],
                        "auth_atom_id": atom_name,
                        "pdbx_PDB_model_num": str(model_num)
                    })
    
            if helix: helices.append(helix)
            if strand: strands.append(strand)
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

    return mmcif_dict



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

"""
dict_keys([
    'mmtfVersion', 'mmtfProducer',
    'numBonds', 'numAtoms', 'numGroups', 'numChains', 'numModels', 'structureId',
    'title', 'chainsPerModel', 'groupsPerChain', 'chainNameList', 'chainIdList',
    'spaceGroup', 'unitCell', 'bioAssemblyList', 'bondAtomList', 'bondOrderList', 
    'groupList', 'xCoordList', 'yCoordList', 'zCoordList', 'bFactorList', 'secStructList',
    'occupancyList', 'altLocList', 'insCodeList', 'groupTypeList', 'groupIdList',
    'atomIdList', 'sequenceIndexList', 'experimentalMethods', 'resolution',
    'rFree', 'rWork', 'entityList', 'depositionDate', 'releaseDate', 'ncsOperatorList'
])"""
'''
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
      "technique": None, "source_organism": None, "expression_system": None,
      "missing_residues": []
     }, "quality": {"resolution": None, "rvalue": None, "rfree": None},
     "geometry": {"assemblies": [], "crystallography": {}}, "models": []
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
    mmtf_to_data_transfer(mmtf_dict, data_dict["geometry"],
     "crystallography", "space_group", "spaceGroup")
    mmtf_to_data_transfer(mmtf_dict, data_dict["geometry"],
     "crystallography", "unit_cell", "unitCell", trim=3)
    if data_dict["geometry"]["crystallography"].get("space_group") == "NA":
        data_dict["geometry"]["crystallography"] = {}
    data_dict["geometry"]["assemblies"] = [{
     "id": int(a["name"]), "software": None, "delta_energy": None,
     "buried_surface_area": None, "surface_area": None, "transformations": [{
      "chains": [mmtf_dict["chainIdList"][i] for i in t["chainIndexList"]],
      'matrix': [t['matrix'][n * 4: (n * 4) + 3] for n in range(3)],
      "vector": t['matrix'][3:-4:4]} for t in a.get("transformList", [])
     ]
    } for a in mmtf_dict.get("bioAssemblyList", [])]
    update_models_list(mmtf_dict, data_dict)
    return data_dict


def update_models_list(mmtf_dict, data_dict):
    """Takes a data dictionary and updates its models list with
    information from a .mmtf dictionary.

    :param dict mmtf_dict: the .mmtf dictionary to read.
    :param dict data_dict: the data dictionary to update."""

    atoms = get_atoms_list(mmtf_dict)
    group_definitions = get_group_definitions_list(mmtf_dict)
    groups = get_groups_list(mmtf_dict, group_definitions)
    chains = get_chains_list(mmtf_dict, groups)
    for model_num in range(mmtf_dict["numModels"]):
        model = {"polymer": {}, "non-polymer": {}, "water": {}, "branched": {}}
        for chain_num in range(mmtf_dict["chainsPerModel"][model_num]):
            chain = chains[chain_num]
            add_chain_to_model(chain, model, atoms)
        data_dict["models"].append(model)


def get_atoms_list(mmtf_dict):
    """Creates a list of atom dictionaries from a .mmtf dictionary by zipping
    together some of its fields.

    :param dict mmtf_dict: the .mmtf dictionary to read.
    :rtype: ``list``"""

    return [{
     "x": x, "y": y, "z": z, "alt_loc": a or None, "bvalue": b, "occupancy": o,
     "id": i, "is_hetatm": False
    } for x, y, z, a, b, i, o in zip(
     mmtf_dict["xCoordList"], mmtf_dict["yCoordList"], mmtf_dict["zCoordList"],
     mmtf_dict["altLocList"], mmtf_dict["bFactorList"], mmtf_dict["atomIdList"],
     mmtf_dict["occupancyList"]
    )]


def get_group_definitions_list(mmtf_dict):
    """Gets a list of group definitions from the .mmtf dict and packs its atom
    attributes into atoms dicts.

    :param dict mmtf_dict: the .mmtf dictionary to read.
    :rtype: ``list``"""


    group_definitions = []
    for group in mmtf_dict["groupList"]:
        atoms = [{
         "name": name, "element": element.upper(), "charge": charge
        } for name, element, charge in zip(
         group["atomNameList"], group["elementList"], group["formalChargeList"],
        )]
        group_definitions.append({
         "name": group["groupName"], "atoms": atoms
        })
    return group_definitions


def get_groups_list(mmtf_dict, group_definitions):
    """Creates a list of group dictionaries from a .mmtf dictionary by zipping
    together some of its fields.

    :param dict mmtf_dict: the .mmtf dictionary to read.
    :rtype: ``list``"""

    sec_struct = [
     "helices", None, "helices", "strands", "helices", "strands", None, None
    ]
    return [{
     "number": id, "insert": insert, "secondary_structure": sec_struct[ss],
     **group_definitions[type_]
    } for id, insert, ss, type_, in zip(
     mmtf_dict["groupIdList"], mmtf_dict["insCodeList"],
     mmtf_dict.get("secStructList", [-1] * len(mmtf_dict["groupIdList"])),
     mmtf_dict["groupTypeList"]
    )]


def get_chains_list(mmtf_dict, groups):
    """Creates a list of chain dictionaries from a .mmtf dictionary by zipping
    together some of its fields.

    :param dict mmtf_dict: the .mmtf dictionary to read.
    :rtype: ``list``"""

    chains = []
    for i_id, id, group_num in zip(mmtf_dict["chainIdList"],
     mmtf_dict["chainNameList"], mmtf_dict["groupsPerChain"]):
        chain = {"id": id, "internal_id": i_id, "groups": groups[:group_num]}
        del groups[:group_num]
        for entity in mmtf_dict["entityList"]:
            if len(chains) in entity["chainIndexList"]:
                chain["type"] = entity["type"]
                chain["sequence"] = entity.get("sequence", "")
                chain["full_name"] = entity.get("description", None)
                break
        chains.append(chain)
    return chains


def add_chain_to_model(chain, model, atoms):
    """Adds a 'chain' to a model - a chain in the .mmtf dict, which can also be
    a non-polymer.

    :param dict chain: the 'chain' to add.
    :param dict model: the model to add it to.
    :param list atoms: the atoms list to work through."""

    if chain["type"] == "polymer" or chain["type"] == "branched":
        polymer = {
         "internal_id": chain["internal_id"], "sequence": chain["sequence"],
         "helices": [], "strands": [], "residues": {}
        }
        for i, group in enumerate(chain["groups"], start=1):
            add_het_to_dict(group, chain, atoms, polymer["residues"], number=i)
        add_ss_to_chain(polymer)
        model["polymer"][chain["id"]] = polymer
    else:
        for group in chain["groups"]:
            add_het_to_dict(group, chain, atoms, model[chain["type"]])
         

def add_het_to_dict(group, chain, atoms, d, number=None):
    """Adds a ligand or water or residue to the appropriate dict. An ID and name
    will be generated for it.

    :param dict group: the group template the het should be based on.
    :param dict chain: the chain (in the real sense) the het is associated with.
    :param list atoms: the atoms list to work through.
    :param dict d: the dictionary to add to.
    :param int number: if given, the residue number to use."""

    het_id = f"{chain['id']}.{group['number']}{group['insert']}"
    het_atoms = atoms[:len(group["atoms"])]
    del atoms[:len(het_atoms)]
    het_atoms = {a["id"]: {
     "anisotropy": [0] * 6, **a, **g_a
    } for a, g_a in zip(het_atoms, group["atoms"])}
    for a in het_atoms.values(): del a["id"]
    het = {
     "name": group["name"], "atoms": het_atoms, "full_name": None,
     "secondary_structure": group["secondary_structure"]
    }
    if number is None:
        het["internal_id"] = chain["internal_id"]
        het["polymer"] = chain["id"]
        het["full_name"] = chain["full_name"]
    else:
        het["number"] = number
    d[het_id] = het


def add_ss_to_chain(chain):
    """Updates polymer dictionary with secondary structure information, from
    information temporarily stored in its residue dicts.

    :param dict chain: the chain to update."""
    
    in_ss = {"helices": False, "strands": False}
    for res_id, res in chain["residues"].items():
        ss = res["secondary_structure"]
        if ss:
            if not in_ss[ss]:
                chain[ss].append([])
            in_ss[ss] = True
            chain[ss][-1].append(res_id)
        else:
            if in_ss["helices"]: in_ss["helices"] = False
            if in_ss["strands"]: in_ss["strands"] = False
        del res["secondary_structure"]'''