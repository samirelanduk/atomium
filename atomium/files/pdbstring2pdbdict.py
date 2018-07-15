"""This module handles the conversion of .pdb filestrings to PDB data
dictionaries."""

from datetime import datetime
import re

def pdb_string_to_pdb_dict(filestring):
    """Converts the string of a .pdb file to a parsed data ``dict``.

    :param str filestring: The filestring to parse.
    :rtype: ``dict``"""

    from .utilities import string_to_lines
    lines = string_to_lines(filestring, width=80)
    pdb_dict = {}
    extract_annotation(pdb_dict, lines)
    extract_structure(pdb_dict, lines)
    return pdb_dict


def extract_annotation(pdb_dict, lines):
    """Takes a ``dict`` and adds header information to it by parsing file lines.

    :param dict pdb_dict: the ``dict`` to update.
    :param list lines: the file lines to read from."""

    extract_header(pdb_dict, lines)
    extract_title(pdb_dict, lines)
    extract_resolution(pdb_dict, lines)
    extract_rfactor(pdb_dict, lines)
    extract_source(pdb_dict, lines)
    extract_technique(pdb_dict, lines)
    extract_keywords(pdb_dict, lines)
    extract_sequence(pdb_dict, lines)
    extract_biomolecules(pdb_dict, lines)


def extract_header(pdb_dict, lines):
    """Takes a ``dict`` and adds header information to it by parsing the HEADER
    line.

    :param dict pdb_dict: the ``dict`` to update.
    :param list lines: the file lines to read from."""

    headline = get_line("HEADER", lines)
    if headline:
        pdb_dict["deposition_date"] = datetime.strptime(
         headline[50:59], "%d-%b-%y"
        ).date() if headline[50:59].strip() else None
        pdb_dict["code"] = headline[62:66] if headline[62:66].strip() else None
        pdb_dict["classification"] = (
         headline[10:50].strip() if headline[10:50].strip() else None
        )
    else:
        pdb_dict["deposition_date"] = None
        pdb_dict["code"] = None
        pdb_dict["classification"] = None


def extract_title(pdb_dict, lines):
    """Takes a ``dict`` and adds title information to it by parsing the TITLE
    line.

    :param dict pdb_dict: the ``dict`` to update.
    :param list lines: the file lines to read from."""

    title_lines = get_lines("TITLE", lines)
    pdb_dict["title"] = merge_lines(title_lines, 10) if title_lines else None


def extract_resolution(pdb_dict, lines):
    """Takes a ``dict`` and adds resolution information to it by parsing
    REMARK 2.

    :param dict pdb_dict: the ``dict`` to update.
    :param list lines: the file lines to read from."""

    remark_lines = get_lines("REMARK", lines)
    for remark in remark_lines:
        if int(remark[7:10]) == 2 and remark[10:].strip():
            try:
                pdb_dict["resolution"] = float(remark[10:].strip().split()[1])
            except ValueError: pdb_dict["resolution"] = None
            break
    else:
        pdb_dict["resolution"] = None


def extract_rfactor(pdb_dict, lines):
    """Takes a ``dict`` and adds rfactor information to it by parsing
    REMARK 3.

    :param dict pdb_dict: the ``dict`` to update.
    :param list lines: the file lines to read from."""

    remark_lines = get_lines("REMARK", lines)
    patterns = {
     "rfactor": r"R VALUE[ ]{2,}\(WORKING SET\) : (.+)",
     "rfree": r"FREE R VALUE[ ]{2,}: (.+)",
     "rcount": r"FREE R VALUE TEST SET COUNT[ ]{2,}: (.+)"
    }
    for attribute, pattern in patterns.items():
        for remark in remark_lines:
            if int(remark[7:10]) == 3 and remark[10:].strip():
                matches = re.findall(pattern, remark[10:].strip())
                if matches:
                    try:
                        pdb_dict[attribute] = float(matches[0].strip())
                        break
                    except: pass
    for attr in patterns:
        if attr not in pdb_dict:
            pdb_dict[attr] = None


def extract_source(pdb_dict, lines):
    """Takes a ``dict`` and adds source information to it by parsing
    SOURCE.

    :param dict pdb_dict: the ``dict`` to update.
    :param list lines: the file lines to read from."""

    pdb_dict["organism"], pdb_dict["expression_system"] = None, None
    source_lines = get_lines("SOURCE", lines)
    if source_lines:
        data = merge_lines(source_lines, 10)
        pattern = r"ORGANISM_SCIENTIFIC\: (.+?);"
        matches = re.findall(pattern, data)
        if matches:
            pdb_dict["organism"] = matches[0]
        pattern = r"EXPRESSION_SYSTEM\: (.+?);"
        matches = re.findall(pattern, data)
        if matches:
            pdb_dict["expression_system"] = matches[0]


def extract_technique(pdb_dict, lines):
    """Takes a ``dict`` and adds technique information to it by parsing file
    lines.

    :param dict pdb_dict: the ``dict`` to update.
    :param list lines: the file lines to read from."""

    line = get_line("EXPDTA", lines)
    if line:
        pdb_dict["technique"] = line[6:].strip()
        if pdb_dict["technique"]: return
    pdb_dict["technique"] = None


def extract_keywords(pdb_dict, lines):
    """Takes a ``dict`` and adds keyword information to it by parsing file
    lines.

    :param dict pdb_dict: the ``dict`` to update.
    :param list lines: the file lines to read from."""

    lines = get_lines("KEYWDS", lines)
    pdb_dict["keywords"] = []
    if lines:
        text = merge_lines(lines, 10)
        pdb_dict["keywords"] = [word.strip() for word in text.split(",")]


def extract_sequence(pdb_dict, lines):
    """Takes a ``dict`` and adds sequence information to it by parsing file
    lines.

    :param dict pdb_dict: the ``dict`` to update.
    :param list lines: the file lines to read from."""

    lines = get_lines("SEQRES", lines)
    pdb_dict["sequences"] = {}
    if lines:
        text = merge_lines(lines, 10)
        blocks, code = text.split(), None
        blocks = [b for b in blocks if not b.isdigit()]
        for block in blocks:
            if len(block) == 1 and block not in pdb_dict["sequences"]:
                pdb_dict["sequences"][block], code = [], block
            elif code and not len(block) == 1:
                pdb_dict["sequences"][code].append(block)


def extract_biomolecules(pdb_dict, lines):
    """Takes a ``dict`` and adds sequence information to it by parsing file
    lines - specifically REMARK 350.

    :param dict pdb_dict: the ``dict`` to update.
    :param list lines: the file lines to read from."""

    lines = [l for l in get_lines("REMARK", lines) if int(l[7:10]) == 350]
    biomolecules = break_lines_on_regex(lines, r"REMARK 350 BIOMOLECULE: (\d+)")
    pdb_dict["biomolecules"] = []
    for biomolecule_lines in biomolecules:
        biomolecule = {
         "transformations": [], "software": None, "buried_surface_area": None,
         "surface_area": None, "delta_energy": None
        }
        patterns = [
         [r"(.+)SOFTWARE USED: (.+)", "software", lambda x: x],
         [r"(.+)BIOMOLECULE: (.+)", "id", int],
         [r"(.+)BURIED SURFACE AREA: (.+) [A-Z]", "buried_surface_area", float],
         [r"(.+)AREA OF THE COMPLEX: (.+) [A-Z]", "surface_area", float],
         [r"(.+)FREE ENERGY: (.+) [A-Z]", "delta_energy", float],
        ]
        for line in biomolecule_lines:
            for pattern in patterns:
                matches = re.findall(pattern[0], line)
                if matches:
                    biomolecule[pattern[1]] = pattern[2](matches[0][1].strip())
        operations = break_lines_on_regex(
         biomolecule_lines, r"(.+)APPLY THE FOLLOWING TO CHAINS"
        )
        for operation_lines in operations:
            chain_lines = [l.strip().split(": ")[1] for l in operation_lines
             if "CHAINS:" in l]
            chain_lines = [l[:-1] if l.endswith(",") else l for l in chain_lines]
            chains = ", ".join(chain_lines)
            chains = chains.split(", ")
            for matrix in break_lines_on_regex(operation_lines, r"(.+)BIOMT1"):
                matrix = [l for l in matrix if "BIOMT" in l]
                vector = [float(line.strip().split()[-1]) for line in matrix]
                matrix = [[float(x) for x in line.split()[-4:-1]] for line in matrix]
                biomolecule["transformations"].append({
                 "chains": chains, "vector": vector, "matrix": matrix
                })
        pdb_dict["biomolecules"].append(biomolecule)



def extract_structure(pdb_dict, lines):
    """Takes a ``dict`` and adds structure information to it by parsing file
    lines.

    :param dict pdb_dict: the ``dict`` to update.
    :param list lines: the file lines to read from."""

    model_lines = get_lines("MODEL", lines, number=True)
    atom_lines = get_lines("ATOM", lines, number=bool(model_lines))
    hetatm_lines = get_lines("HETATM", lines, number=bool(model_lines))
    anisou_lines = get_lines("ANISOU", lines, number=bool(model_lines))
    conect_lines = get_lines("CONECT", lines)
    pdb_dict["models"] = []
    if model_lines:
        for index, (line, line_number) in enumerate(model_lines):
            next_ = model_lines[index + 1][1] if (
             index < len(model_lines) - 1
            ) else len(lines)
            m_atoms = [line for line, number in atom_lines
             if line_number < number < next_]
            m_het = [line for line, number in hetatm_lines
             if line_number < number < next_]
            m_anisou = [line for line, number in anisou_lines
             if line_number < number < next_]
            pdb_dict["models"].append(lines_to_model(m_atoms, m_het, m_anisou))
    else:
        pdb_dict["models"].append(
         lines_to_model(atom_lines, hetatm_lines, anisou_lines)
        )
    extract_connections(pdb_dict, conect_lines)


def lines_to_model(atom_lines, hetatm_lines, anisou_lines):
    """Creates a model ``dict`` from ATOM, ANISOU and HETATM lines.

    :param list atom_lines: the ATOM lines.
    :param list hetatm_lines: the HETATM lines.
    :param list anisou_lines: the ANISOU lines.
    :rtype: ``dict``"""

    atoms = [atom_line_to_atom_dict(line) for line in atom_lines]
    heteroatoms = [atom_line_to_atom_dict(line) for line in hetatm_lines]
    assign_anisou(anisou_lines, atoms, heteroatoms)
    chains = atoms_to_chains(atoms, heteroatoms)
    return chains


def atom_line_to_atom_dict(line):
    """Takes an ATOM or HETATM line and converts it to an atom ``dict``.

    :param str line: the atom record to parse.
    :rtype: ``dict``"""

    charge = line[78:80].strip() if line[78:80].strip() else 0
    try:
        charge = float(charge)
    except:
        charge = float(charge[::-1])
    d = {
     "atom_id": int(line[6:11].strip()) if line[6:11].strip() else None,
     "atom_name": line[12:16].strip() if line[12:16].strip() else None,
     "alt_loc": line[16] if line[16].strip() else None,
     "residue_name": line[17:20].strip() if line[17:20].strip() else None,
     "chain_id": line[21] if line[21].strip() else "",
     "residue_id": int(line[22:26].strip()) if line[22:26].strip() else None,
     "insert_code": line[26] if line[26].strip() else "",
     "x": float(line[30:38].strip()) if line[30:38].strip() else None,
     "y": float(line[38:46].strip()) if line[38:46].strip() else None,
     "z": float(line[46:54].strip()) if line[46:54].strip() else None,
     "occupancy": float(line[54:60].strip()) if line[54:60].strip() else 1,
     "temp_factor": float(line[60:66].strip()) if line[60:66].strip() else None,
     "element": line[76:78].strip() if line[76:78].strip() else None,
     "charge": charge
    }
    d["full_id"] = d["chain_id"] + (
     line[22:26].strip() if line[22:26].strip() else ""
    ) + d["insert_code"]
    return d


def atoms_to_residues(atoms):
    """Takes a list of atoms, clusters them into lists belonging to the
    different residues, and makes residue ``dict`` objects out of them. The
    order of residues will be the order they appear in the atoms.

    :param list atoms: The atom ``dict`` objects to collate.
    :rtype: ``list``"""

    residue_ids = [atom["full_id"] for atom in atoms]
    residue_ids = sorted(set(residue_ids), key=lambda r: residue_ids.index(r))
    residue_dict = {res_id: [] for res_id in residue_ids}
    residues = []
    for atom in atoms:
        residue_dict[atom["full_id"]].append(atom)
    for res_id in residue_ids:
        r_atoms = residue_dict[res_id]
        residues.append({
         "id": res_id, "name": r_atoms[0]["residue_name"], "atoms": r_atoms
        })
    return residues


def atoms_to_chains(atoms, heteroatoms):
    """Takes a list of atoms, clusters them into lists belonging to the
    different chains, and makes residue ``dict`` objects out of them.

    :param list atoms: The atom ``dict`` objects to collate.
    :rtype: ``list``"""

    chain_ids = sorted(set([atom["chain_id"] for atom in atoms + heteroatoms]))
    chains = []
    for id_ in chain_ids:
        c_atoms = [atom for atom in atoms if atom["chain_id"] == id_]
        c_heteroatoms = [atom for atom in heteroatoms if atom["chain_id"] == id_]
        chains.append({
         "chain_id": id_,
         "residues": atoms_to_residues(c_atoms),
         "ligands": atoms_to_residues(c_heteroatoms)
        })
    return chains


def assign_anisou(anisou_lines, atoms, heteroatoms):
    """Takes a list of ANISOU records, and atom and heteroatom lists, and
    assigns anisotropy information the atoms and heteroatoms.

    :param list anisou_lines: The ANISOU records.
    :param list atoms: The atom ``dict`` objects to update.
    :param list heteroatoms: The heteroatom ``dict`` objects to update."""

    anisotropy = {int(line[6:11].strip()): [
     int(line[n * 7 + 28:n * 7 + 35]) / 10000 for n in range(6)
    ] for line in anisou_lines}
    for atom in atoms + heteroatoms:
        atom["anisotropy"] = anisotropy.get(atom["atom_id"], [0, 0, 0, 0, 0, 0])


def extract_connections(pdb_dict, lines):
    """Takes a ``dict`` and adds connection information to it by parsing conect
    lines.

    :param dict pdb_dict: the ``dict`` to update.
    :param list lines: the conect lines to read from."""

    atom_ids = sorted(list(set([int(r[6:11].strip()) for r in lines])))
    pdb_dict["connections"] = [{"atom": id_, "bond_to": []} for id_ in atom_ids]
    for connection in pdb_dict["connections"]:
        relevant_lines = [line[11:] for line in lines
          if int(line[6:11].strip()) == connection["atom"]]
        line_numbers = [[
         line[n * 5:n * 5 + 5].strip() for n in range(4)
        ] for line in relevant_lines]
        for line in line_numbers:
            for number in line:
                if number:
                    connection["bond_to"].append(int(number))


def get_line(name, lines, number=False):
    """Gets a single line by record name from a list of lines. If there are no
    matches, ``None`` will be returned.

    :param str name: The record name to search for.
    :param list lines: The list of lines to search through.
    :rtype: ``str``"""

    for index, line in enumerate(lines):
        if line[:6].strip() == name:
            return [line, index] if number else line


def get_lines(name, lines, number=False):
    """Searches a list of lines and returns those whose record name matches the
    name given.

    :param str name: The record name to search for.
    :param list lines: The list of lines to search through.
    :rtype: ``list``"""

    return [[line, index] if number else line for
     index, line in enumerate(lines) if line[:6].strip() == name]


def merge_lines(lines, start, join=" "):
    """Gets a single continuous string from a sequence of lines.

    :param list lines: The lines to merge.
    :param int start: The start point in each record.
    :param str join: The string to join on.
    :rtype: ``str``"""

    string = join.join([line[start:].strip() for line in lines])
    return string


def break_lines_on_regex(lines, regex):
    """Takes a list of strings, and divides them into smaller lists based on
    which lines match a regex pattern.

    :param list lines: The list of lines.
    :param str regex: The pattern to use.
    :rtype: ``list``"""

    groups, group = [], []
    for line in lines:
        if re.match(regex, line):
            groups.append(group)
            group = []
        group.append(line)
    groups.append(group)
    return filter(lambda g: bool(g) and re.match(regex, g[0]), groups)
