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

    mmcif_dict["refine"] = [{
        "entry_id": mmtf["structureId"],
        "ls_d_res_high": mmtf["resolution"],
        "ls_R_factor_obs": mmtf["rWork"],
        "ls_R_factor_all": mmtf["rWork"],
        "ls_R_factor_R_work": mmtf["rWork"],
        "ls_R_factor_R_free": mmtf["rFree"],
    }]

    mmcif_dict["symmetry"] = [{
        "entry_id": mmtf["structureId"],
        "space_group_name_H-M": mmtf["spaceGroup"]
    }]

    mmcif_dict["cell"] = [{
        "entry_id": mmtf["structureId"],
        "length_a": mmtf["unitCell"][0],
        "length_b": mmtf["unitCell"][1],
        "length_c": mmtf["unitCell"][2],
        "length_alpha": mmtf["unitCell"][3],
        "length_beta": mmtf["unitCell"][4],
        "length_gamme": mmtf["unitCell"][5],
    }]

    if mmtf["bioAssemblyList"]:
        mmcif_dict["pdbx_struct_assembly"] = []
        mmcif_dict["pdbx_struct_assembly_gen"] = []
        mmcif_dict["pdbx_struct_assembly_prop"] = []
        for assembly in mmtf["bioAssemblyList"]:
            mmcif_dict["pdbx_struct_assembly"].append({"id": assembly["name"]})
            for transformation in assembly["transformList"]:
                op_id = len(mmcif_dict["pdbx_struct_assembly_prop"]) + 1
                mmcif_dict["pdbx_struct_assembly_gen"].append({
                    "assembly_id": assembly["name"],
                    "oper_expression": op_id,
                    "asym_id_list": ",".join([mmtf["chainIdList"][index] for index in transformation["chainIndexList"]])
                })
                mmcif_dict["pdbx_struct_assembly_gen"].append({
                    "id": op_id,
                    "pdbx_struct_oper_list.matrix[1][1]": transformation["matrix"][0],
                    "pdbx_struct_oper_list.matrix[1][2]": transformation["matrix"][1],
                    "pdbx_struct_oper_list.matrix[1][3]": transformation["matrix"][2],
                    "pdbx_struct_oper_list.matrix[1]": transformation["matrix"][3],
                    "pdbx_struct_oper_list.matrix[2][1]": transformation["matrix"][4],
                    "pdbx_struct_oper_list.matrix[2][2]": transformation["matrix"][5],
                    "pdbx_struct_oper_list.matrix[2][3]": transformation["matrix"][6],
                    "pdbx_struct_oper_list.matrix[2]": transformation["matrix"][7],
                    "pdbx_struct_oper_list.matrix[3][1]": transformation["matrix"][8],
                    "pdbx_struct_oper_list.matrix[3][2]": transformation["matrix"][9],
                    "pdbx_struct_oper_list.matrix[3][3]": transformation["matrix"][10],
                    "pdbx_struct_oper_list.matrix[3]": transformation["matrix"][11],
                })
    
    mmcif_dict["entity"] = []
    for entity in mmtf["entityList"]:
        mmcif_dict["entity"].append({
            "id": len(mmcif_dict["entity"]) + 1,
            "type": entity["type"],
            "pdbx_description": entity["description"],
            "pdbx_number_of_molecules": len(entity["chainIndexList"])
        })
        if entity["type"] == "polymer":
            if "entity_poly" not in mmcif_dict:
                mmcif_dict["entity_poly"] = []
            mmcif_dict["entity_poly"].append({
                "entity_id": len(mmcif_dict["entity"]),
                "pdbx_seq_one_letter_code": entity["sequence"],
                "pdbx_seq_one_letter_code_can": entity["sequence"],
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
      "matrix": [t["matrix"][n * 4: (n * 4) + 3] for n in range(3)],
      "vector": t["matrix"][3:-4:4]} for t in a.get("transformList", [])
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