import re
import math
import calendar
import itertools
from collections import Counter
from atomium.file import get_operations
from atomium.sequences import get_alignment_indices
from atomium.data import WATER_NAMES, CODES, FULL_NAMES, FORMULAE
from atomium.data import RESIDUE_MASSES, PERIODIC_TABLE
from atomium.categories import REFLNS, REFINE

def pdb_string_to_mmcif_dict(filestring):
    """Takes the contents of a .pdb file and returns a fully parsed dictionary
    representation of it. First the metadata is parsed to create the initial
    mmCIF dictionary, Then a description of the entities involved is obtained,
    which are then incorporated into the mmCIF. Structural annotations are then
    added, and then IDs are updated globally.
    
    :param str filestring: the contents of the .pdb file.
    :rtype: ``dict``"""

    mmcif = parse_metadata(filestring)
    polymer_entities, non_polymer_entities = get_entities(filestring)
    build_structure_categories(
        filestring, polymer_entities, non_polymer_entities, mmcif
    )
    parse_annotation(filestring, mmcif)
    update_auth_and_label(mmcif)
    return mmcif


def parse_metadata(filestring):
    """Creates an initial mmCIF dictionary with information about metadata from
    thr PDB filestring - that is, descriptions of the experiment etc., rather
    than details of the structure contained.
    
    :param str filestring: the contents of the .pdb file.
    :rtype: ``dict``"""

    mmcif = {}
    parse_header(filestring, mmcif)
    parse_obslte(filestring, mmcif)
    parse_title(filestring, mmcif)
    parse_split(filestring, mmcif)
    parse_caveat(filestring, mmcif)
    parse_keywds(filestring, mmcif)
    parse_expdta(filestring, mmcif)
    parse_mdltyp(filestring, mmcif)
    parse_author(filestring, mmcif)
    parse_revdat(filestring, mmcif)
    parse_sprsde(filestring, mmcif)
    parse_jrnl(filestring, mmcif)
    parse_remarks(filestring, mmcif)
    parse_cryst1(filestring, mmcif)
    parse_origx(filestring, mmcif)
    parse_scale(filestring, mmcif)
    parse_mtrix(filestring, mmcif)
    return mmcif


def get_entities(filestring):
    """Gets a description of the polymer and non-polymer entities in a PDB file,
    based on descriptive records and the actual contents of the ATOM/HETATM
    records. Each entity contains all the descriptive information needed to
    incorporate them into the mmCIF dictionary, and a list of the actual polymer
    molecules that represent the entity.

    :param str filestring: the contents of the .pdb file.
    :rtype: ``tuple``"""

    polymer_entities = parse_compnd_and_source(filestring)
    polymers = parse_dbref(filestring)
    parse_seqadv(filestring, polymers)
    parse_seqres(filestring, polymers)
    parse_modres(filestring, polymers)
    non_polymer_entities, non_polymers = parse_het(filestring)
    parse_hetnam(filestring, non_polymer_entities)
    parse_hetsyn(filestring, non_polymer_entities)
    parse_formul(filestring, non_polymer_entities)
    update_entities_from_atoms(
        filestring, polymer_entities, polymers,
        non_polymer_entities, non_polymers
    )
    finalize_entities(polymer_entities, non_polymer_entities)
    non_polymer_entities = ordered_by_atom(non_polymer_entities, filestring)
    finalize_polymers(polymers)
    add_molecules_to_entities(
        polymer_entities, polymers, non_polymer_entities, non_polymers
    )
    return polymer_entities, non_polymer_entities


def build_structure_categories(filestring, polymer_entities,
                               non_polymer_entities, mmcif):
    """Builds the structure categories in an mmCIF dictionary, based on the PDB
    file contents and the pre-calculated entity information.
    
    :param str filestring: the contents of the .pdb file.
    :param dict polymer_entities: the polymer entity information.
    :param dict non_polymer_entities: the non-polymer entity information.
    :param dict mmcif: the dictionary to update."""
    
    build_entity_category(polymer_entities, non_polymer_entities, mmcif)
    build_entity_name_com(polymer_entities, mmcif)
    build_entity_poly(polymer_entities, mmcif)
    build_entity_poly_seq(polymer_entities, mmcif)
    build_struct_ref(polymer_entities, mmcif)
    build_struct_ref_seq(polymer_entities, mmcif)
    build_struct_ref_seq_dif(polymer_entities, mmcif)
    build_pdbx_struct_mod_residue(polymer_entities, mmcif)
    build_pdbx_entity_nonpoly(non_polymer_entities, mmcif)
    build_chem_comp(non_polymer_entities, mmcif)
    build_atom_type(filestring, mmcif)
    build_atom_site(filestring, polymer_entities, non_polymer_entities, mmcif)
    build_atom_site_anisotrop(filestring, mmcif)
    update_atom_ids(mmcif)
    build_struct_asym(mmcif)


def parse_header(filestring, mmcif):
    """Parses the HEADER record of a PDB and creates the three mmCIF categories
    that contain its information.
    
    :param str filestring: the contents of the .pdb file.
    :param dict mmcif: the dictionary to update."""
    
    header = re.search(f"HEADER.+", filestring, re.M)
    header = header.group(0) if header else None
    code = (header[62:66].strip() or "XXXX") if header else "XXXX"
    date = header[50:59].strip() if header else ""
    date = pdb_date_to_mmcif_date(date) or "?"
    keyword = (header[10:50].strip() or "?") if header else "?"
    mmcif["entry"] = [{"id": code}]
    mmcif["pdbx_database_status"] = [{
        "status_code": "REL", "entry_id": code,
        "recvd_initial_deposition_date": date
    }]
    mmcif["struct_keywords"] = [{
        "entry_id": code, "pdbx_keywords": keyword, "text": "?"
    }]


def parse_obslte(filestring, mmcif):
    """Parses the OBSLTE records and makes the ``pdbx_database_PDB_obs_spr``
    table as a result if needed. Only the first OBSTLE record is used.
    
    :param str filestring: the contents of the .pdb file.
    :param dict mmcif: the dictionary to update."""

    obslte_lines = re.findall(r"^OBSLTE.+", filestring, re.M)
    if not obslte_lines: return
    codes = [code for l in obslte_lines for code in l[31:].strip().split()]
    mmcif["pdbx_database_PDB_obs_spr"] = [{
        "id": "OBSLTE", "details": "?",
        "date": pdb_date_to_mmcif_date(obslte_lines[0][11:20].strip()),
        "pdb_id": " ".join(codes),
        "replace_pdb_id": obslte_lines[0][21:25].strip()
    }]


def parse_title(filestring, mmcif):   
    """Parses the TITLE records and makes the ``struct`` table as a result if
    needed.
    
    :param str filestring: the contents of the .pdb file.
    :param dict mmcif: the dictionary to update."""

    mmcif["struct"] = [{"entry_id": mmcif["entry"][0]["id"], **{
        k: "?" for k in [
            "title", "pdbx_descriptor", "pdbx_model_details", "pdbx_CASP_flag",
            "pdbx_model_type_details"
        ]
    }}]
    title_lines = re.findall(r"^TITLE.+", filestring, re.M)
    title = " ".join([l[10:80].strip() for l in title_lines]).strip()
    mmcif["struct"][0]["title"] = title or "?"


def parse_split(filestring, mmcif):
    """Parses the SPLIT records and creates the ``pdbx_database_related`` table
    as a result if needed. The contents here may be enhanced by the contents of
    REMARK 900.
    
    :param str filestring: the contents of the .pdb file.
    :param dict mmcif: the dictionary to update."""

    split_lines = re.findall(r"^SPLIT.+", filestring, re.M)
    if split_lines:
        codes = [code for l in split_lines for code in l[11:].strip().split()]
        mmcif["pdbx_database_related"] = [{
            "db_name": "PDB", "db_id": code, "content_type": "split",
            "details": f"Split {n}"
        } for n, code in enumerate(codes, start=1)]


def parse_caveat(filestring, mmcif):
    """Parses the CAVEAT records and creates the ``database_PDB_caveat`` table
    as a result if needed.

    :param str filestring: the contents of the .pdb file.
    :param dict mmcif: the dictionary to update."""

    caveat_lines = re.findall(r"^CAVEAT.+", filestring, re.M)
    if not caveat_lines: return
    caveat = " ".join([l[19:80] for l in caveat_lines]).strip()
    mmcif["database_PDB_caveat"] = [{"id": "1", "text": " ".join(caveat.split())}]


def parse_keywds(filestring, mmcif):
    """Parses the KEYWDS records and updates the ``struct_keywords`` table,
    which it expects to exist.

    :param str filestring: the contents of the .pdb file.
    :param dict mmcif: the dictionary to update."""

    lines = re.findall(r"^KEYWDS.+", filestring, re.M)
    if not lines: return
    keywords = " ".join([l[10:].strip() for l in lines]).strip()
    mmcif["struct_keywords"][0]["text"] = keywords


def parse_expdta(filestring, mmcif):
    """Parses the EXPDTA records and creates the ``exptl`` table from it.

    :param str filestring: the contents of the .pdb file.
    :param dict mmcif: the dictionary to update."""

    exptl = re.findall(f"^EXPDTA.+", filestring, re.M)
    if not exptl: return
    mmcif["exptl"] = [{
        "entry_id": mmcif["entry"][0]["id"], "method": method
    } for method in " ".join([l[10:79].strip() for l in exptl]).split("; ")]


def parse_mdltyp(filestring, mmcif):
    """Parses the MDLTYP records. The ``pdbx_coordinate_model`` table will be
    created if there are any comments about specific chains, and anything else
    will be added to the ``struct`` table, which is expected to exist.
    
    :param str filestring: the contents of the .pdb file.
    :param dict mmcif: the dictionary to update."""

    mdltyp = re.findall(f"^MDLTYP.+", filestring, re.M)
    if not mdltyp: return
    mmcif["pdbx_coordinate_model"] = []
    text = " ".join([l[10:80].strip() for l in mdltyp]).strip()
    features = text.split(";")
    sections = []
    for feature in features:
        match = re.match(r"(.+?), CHAIN (.+)", feature)
        if match:
            for chain in match[2].split(","):
                mmcif["pdbx_coordinate_model"].append({
                    "asym_id": chain.strip(), "type": match[1].strip()
                })
        else:
            sections.append(feature.strip())
    if sections:
        mmcif["struct"][0]["pdbx_model_type_details"] = " ; ".join(sections)
    if not mmcif["pdbx_coordinate_model"]: del mmcif["pdbx_coordinate_model"]


def parse_author(filestring, mmcif):
    """Parses the AUTHOR records and creates the ``audit_author`` table from it.
    Names will be lightly reformatted.

    :param str filestring: the contents of the .pdb file.
    :param dict mmcif: the dictionary to update."""

    lines = re.findall(r"^AUTHOR.+", filestring, re.M)
    if not lines: return
    author_lines = [line[10:].strip() for line in lines]
    authors = pdb_names_to_mmcif_names(author_lines)
    mmcif["audit_author"] = [{
        "name": author, "pdbx_ordinal": str(num)
    }  for num, author in enumerate(authors, start=1)]


def parse_revdat(filestring, mmcif):
    """Parses the REVDAT records and creates the basic
    ``pdbx_audit_revision_history`` table. Each is treated as a major release as
    there is no way of distinguishing them.
    
    :param str filestring: the contents of the .pdb file.
    :param dict mmcif: the dictionary to update."""

    lines = re.findall(r"^REVDAT.+", filestring, re.M)
    if not lines: return
    mod_nums = sorted(set([int(line[7:10].strip()) for line in lines]))
    mmcif["pdbx_audit_revision_history"] = [{
        "ordinal": str(num), "data_content_type": "Structure model",
        "major_revision": str(num), "minor_revision": "0",
        "revision_date": pdb_date_to_mmcif_date([
            l for l in lines if int(l[7:10].strip()) == num
        ][0][13:22])
    } for num in mod_nums]


def parse_sprsde(filestring, mmcif):
    """Parses the SPRSDE records and makes or updates the 
    ``pdbx_database_PDB_obs_spr`` table as a result if needed. Only the first
    SPRSDE record is used.
    
    :param str filestring: the contents of the .pdb file.
    :param dict mmcif: the dictionary to update."""

    sprsde_lines = re.findall(r"^SPRSDE.+", filestring, re.M)
    if not sprsde_lines: return
    codes = [code for l in sprsde_lines for code in l[31:].strip().split()]
    sprsde = [{
        "id": "SPRSDE", "details": "?",
        "date": pdb_date_to_mmcif_date(sprsde_lines[0][11:20].strip()),
        "pdb_id": sprsde_lines[0][21:25].strip(),
        "replace_pdb_id": " ".join(codes),
    }]
    mmcif["pdbx_database_PDB_obs_spr"] = [
        *sprsde, *mmcif.get("pdbx_database_PDB_obs_spr", [])
    ]


def parse_jrnl(filestring, mmcif):
    """Parses the JRNL records, and creates the ``citation_author``,
    ``citation_editor`` and ``citation`` tables if needed.
    
    :param str filestring: the contents of the .pdb file.
    :param dict mmcif: the dictionary to update."""

    journal_lines = re.findall(f"^JRNL.+", filestring, re.M)
    if not journal_lines: return
    for record, table in [["AUTH", "author"], ["EDIT", "editor"]]:
        lines = [l[19:].strip() for l in journal_lines if l[12:16] == record]
        names = pdb_names_to_mmcif_names(lines)
        if names:
            mmcif[f"citation_{table}"] = [{
                "citation_id": "primary", "name": name, "pdbx_ordinal": str(num)
            }  for num, name in enumerate(names, start=1)]
    mmcif["citation"] = [{"id": "primary", **{k: "?" for k in [
        "title", "journal_abbrev", "journal_volume", "page_first", "page_last",
        "year", "journal_id_ASTM", "country", "journal_id_ISSN",
        "journal_id_CSD", "book_publisher", "pdbx_database_id_PubMed",
        "pdbx_database_id_DOI"
    ]}}]
    parse_journal_title(journal_lines, mmcif)
    parse_journal_references(journal_lines, mmcif)
    parse_journal_ids(journal_lines, mmcif)


def parse_journal_title(journal_lines, mmcif):
    """Parses JRNL TITL records and updates the ``citation`` table.
    
    :param list journal_lines: lines beginning with JRNL.
    :param dict mmcif: the dictionary to update."""

    title_lines = [l for l in journal_lines if l[12:16] == "TITL"]
    if not title_lines: return
    title = " ".join([l[19:80] for l in title_lines])
    mmcif["citation"][0]["title"] = " ".join(title.split())


def parse_journal_references(journal_lines, mmcif):
    """Parses the JRNL records which pertain to references, updating the
    ``citation`` table.
    
    :param list journal_lines: lines beginning with JRNL.
    :param dict mmcif: the dictionary to update."""

    lines = [l for l in journal_lines if l[12:16].strip() == "REF"]
    if lines:
        mmcif["citation"][0]["journal_abbrev"] = " ".join(
            [l[19:47].strip() for l in lines]
        ).strip().title() or "?"
        mmcif["citation"][0]["journal_volume"] = lines[0][51:55].strip() or "?"
        mmcif["citation"][0]["page_first"] = lines[0][56:61].strip() or "?"
        mmcif["citation"][0]["year"] = lines[0][62:66].strip() or "?"
    lines = [l for l in journal_lines if l[12:16] == "REFN"]
    if lines:
        mmcif["citation"][0]["journal_id_ISSN"] = lines[0][40:65].strip() or "?"
    lines = [l for l in journal_lines if l[12:16] == "PUBL"]
    if lines:
        mmcif["citation"][0]["book_publisher"] = " ".join(
            [l[19:].strip() for l in lines]
        ).strip() or "?"


def parse_journal_ids(journal_lines, mmcif):
    """Parses the JRNL records which pertain to paper IDs, updating the
    ``citation`` table.
    
    :param list journal_lines: lines beginning with JRNL.
    :param dict mmcif: the dictionary to update."""

    pmid_lines = [l for l in journal_lines if l[12:16] == "PMID"]
    if pmid_lines:
        mmcif["citation"][0]["pdbx_database_id_PubMed"] = "".join(
            [l[19:79] for l in pmid_lines]
        ).strip() or "?"
    doi_lines = [l for l in journal_lines if l[12:16].strip() == "DOI"]
    if doi_lines:
        mmcif["citation"][0]["pdbx_database_id_DOI"] = "".join(
            [l[19:79] for l in doi_lines]
        ).strip().lower() or "?"


def parse_remarks(filestring, mmcif):
    """Parses the various REMARK records.
    
    :param str filestring: the contents of the .pdb file.
    :param dict mmcif: the dictionary to update."""

    parse_remark_2(filestring, mmcif)
    parse_remark_3(filestring, mmcif)
    parse_remark_350(filestring, mmcif)
    parse_remark_465(filestring, mmcif)
    parse_remark_800(filestring, mmcif)


def parse_remark_2(filestring, mmcif):
    """Parses REMARK 2 records. This creates the full ``reflns`` category if
    needed.

    :param str filestring: the contents of the .pdb file.
    :param dict mmcif: the dictionary to update."""

    rec = re.search(r"^REMARK   2 RESOLUTION.    (\d+\.\d+)", filestring, re.M)
    if not rec: return
    mmcif["reflns"] = [{
        **{key: "?" for key in REFLNS}, "entry_id": mmcif["entry"][0]["id"]
    }]
    mmcif["reflns"][0]["d_resolution_high"] = rec[1]


def parse_remark_3(filestring, mmcif):
    """Parses REMARK 3 records. The ``reflns`` and ``refine`` table are updated
    from values here.

    :param str filestring: the contents of the .pdb file.
    :param dict mmcif: the dictionary to update."""

    records = re.findall(r"^REMARK   3 .+", filestring, re.M)
    if not records: return
    update_reflns_from_remark_3(records, mmcif)
    update_refine_from_remark_3(records, mmcif)
    
    
def update_reflns_from_remark_3(lines, mmcif):
    """Creates or updates the ``reflns`` table from REMARK 3 contents.
    
    :param list lines: REMARK 3 lines.
    :param dict mmcif: the dictionary to update."""

    string = "\n".join(lines)
    for regex, key in [
        (r"FROM WILSON PLOT.+?\(A\*\*2\).+?\:(.+)", "B_iso_Wilson_estimate")
    ]:
        match = re.search(regex, string)
        if match:
            value = match[1].strip().replace("NULL", "?")
            if not value or value == "?": continue
            if "reflns" not in mmcif:
                mmcif["reflns"] = [{
                    **{key: "?" for key in REFLNS},
                    "entry_id": mmcif["entry"][0]["id"]
                }]
            mmcif["reflns"][0][key] = value


def update_refine_from_remark_3(lines, mmcif):
    """Creates the ``reflns`` table from REMARK 3 contents.
    
    :param list lines: REMARK 3 lines.
    :param dict mmcif: the dictionary to update."""

    string = "\n".join(lines)
    mmcif["refine"] = [{
        **{key: "?" for key in REFINE}, "entry_id": mmcif["entry"][0]["id"]
    }]
    keep = False
    for regex, key in [
        (r"NUMBER OF REFLECTIONS.+?\:(.+)", "ls_number_reflns_obs"),
        (r"DATA CUTOFF HIGH.+?\(ABS\(F\)\).+?\:(.+)", "pdbx_data_cutoff_high_absF"),
        (r"DATA CUTOFF HIGH.+?\(ABS\(F\)\).+?\:(.+)", "pdbx_data_cutoff_high_rms_absF"),
        (r"DATA CUTOFF LOW.+?\(ABS\(F\)\).+?\:(.+)", "pdbx_data_cutoff_low_absF"),
        (r"RESOLUTION RANGE LOW.+?\(ANGSTROMS\).+?\:(.+)", "ls_d_res_low"),
        (r"RESOLUTION RANGE HIGH.+?\(ANGSTROMS\).+?\:(.+)", "ls_d_res_high"),
        (r"COMPLETENESS \(WORKING\+TEST\).+?\(%\).+?\:(.+)", "ls_percent_reflns_obs"),
        (r"R VALUE[ ]+?\(WORKING SET\)[ ]+?\:(.+)", "ls_R_factor_obs"),
        (r"R VALUE[ ]+?\(WORKING SET\)[ ]+?\:(.+)", "ls_R_factor_all"),
        (r"R VALUE[ ]+?\(WORKING SET\)[ ]+?\:(.+)", "ls_R_factor_R_work"),
        (r"FREE R VALUE[ ]+?\:(.+)", "ls_R_factor_R_free"),
        (r"ESTIMATED ERROR OF FREE R VALUE[ ]+?\:(.+)", "ls_R_factor_R_free_error"),
        (r"FREE R VALUE TEST SET SIZE[ ]+?\(%\)[ ]+?\:(.+)", "ls_percent_reflns_R_free"),
        (r"FREE R VALUE TEST SET COUNT[ ]+?\:(.+)", "ls_number_reflns_R_free"),
        (r"MEAN B VALUE[ ]+?\(OVERALL, A\*\*2\)[ ]+?\:(.+)", "B_iso_mean"),
        (r"B11[ ]+?\(A\*\*2\)[ ]+?\:(.+)", "aniso_B[1][1]"),
        (r"B22[ ]+?\(A\*\*2\)[ ]+?\:(.+)", "aniso_B[2][2]"),
        (r"B33[ ]+?\(A\*\*2\)[ ]+?\:(.+)", "aniso_B[3][3]"),
        (r"B12[ ]+?\(A\*\*2\)[ ]+?\:(.+)", "aniso_B[1][2]"),
        (r"B13[ ]+?\(A\*\*2\)[ ]+?\:(.+)", "aniso_B[1][3]"),
        (r"B23[ ]+?\(A\*\*2\)[ ]+?\:(.+)", "aniso_B[2][3]"),
        (r"METHOD USED[ ]+?\:(.+)", "solvent_model_details"),
        (r"KSOL[ ]+?\:(.+)", "solvent_model_param_ksol"),
        (r"BSOL[ ]+?\:(.+)", "solvent_model_param_bsol"),
        (r"CROSS-VALIDATION METHOD[ ]+?\:(.+)", "pdbx_ls_cross_valid_method"),
        (r"ISOTROPIC THERMAL MODEL[ ]+?\:(.+)", "pdbx_isotropic_thermal_model"),
        (r"REFINEMENT TARGET[ ]+?\:(.+)", "pdbx_stereochemistry_target_values"),
        (r"FREE R VALUE TEST SET SELECTION[ ]+?\:(.+)", "pdbx_R_Free_selection_details"),
    ]:
        match = re.search(regex, string)
        value = match[1].strip().replace("NULL", "?") if match else "?"
        if not match or value == "?": continue
        keep = True
        mmcif["refine"][0][key] = value
    if not keep: del mmcif["refine"]


def parse_remark_350(filestring, mmcif):
    """Parses REMARK 350 records to create the four biological assembly
    categories.
    
    :param str filestring: the contents of the .pdb file.
    :param dict mmcif: the dictionary to update."""

    records = re.findall(r"^REMARK 350 .+", filestring, re.M)
    if not records: return
    groups = [list(g) for _, g in itertools.groupby(
        records, lambda x: "MOLECULE:" in x
    )][1:]
    assemblies = [list(itertools.chain(*a))
        for a in zip(groups[::2], groups[1::2])]
    names = [
        "pdbx_struct_assembly", "pdbx_struct_assembly_gen",
        "pdbx_struct_assembly_prop", "pdbx_struct_oper_list"
    ]
    for name in names: mmcif[name] = []
    for assembly in assemblies: parse_assembly(assembly, mmcif)
    for name in names:
        if not mmcif[name]: del mmcif[name]
    

def parse_assembly(lines, mmcif):
    """Takes the REMARK lines of a single biological assembly and adds the
    information to an mmcif dictionary.
    
    :param list lines: the REMARK lines of one assembly.
    :rtype: ``dict``"""

    assembly = create_assembly_dict(lines)
    add_assembly_info_to_mmcif(assembly, mmcif)
    for gen in assembly["gens"]:
        add_assembly_gen_to_mmcif(gen, mmcif)


def create_assembly_dict(lines):
    """Creates a dictionary representation of a biological assembly. As well as
    metadata, it contains an ID and a list of 'gens', each of which is a list of
    transformations and a list of chain IDs to apply to.
    
    :param list lines: the REMARK lines of one assembly.
    :rtype: ``dict``"""

    patterns = [
        [r"(.+)SOFTWARE USED: (.+)", "software"],
        [r"(.+)SURFACE AREA: (.+) [A-Z]", "ABSA (A^2)"],
        [r"(.+)AREA OF THE COMPLEX: (.+) [A-Z]", "SSA (A^2)"],
        [r"(.+)FREE ENERGY: (.+) [A-Z]", "MORE"]
    ]
    assembly, gen = {"gens": []}, None
    for line in lines:
        value = line.split(":")[-1].strip()
        for p in patterns:
            matches = re.findall(p[0], line)
            if matches: assembly[p[1]] = matches[0][1].strip()
        if "APPLY THE FOLLOWING" in line:
            if gen: assembly["gens"].append(gen)
            gen = {
                "chains": [c.strip() for c in value.split(",") if c.strip()],
                "transformations": []
            }
        if "AND CHAINS: " in line:
            gen["chains"] += [c.strip() for c in value.split(",") if c.strip()]
        if "BIOMT" in line:
            if "BIOMT1" in line: transformation = {"matrix": [], "vector": []}
            values = [float(x) for x in line.split()[4:]]
            transformation["matrix"].append(values[:3])
            transformation["vector"].append(values[-1])
            if "BIOMT3" in line: gen["transformations"].append(transformation)
    if gen: assembly["gens"].append(gen)
    return assembly


def add_assembly_info_to_mmcif(assembly, mmcif):
    """Takes an assembly ``dict`` and adds rows for ``pdbx_struct_assembly`` and
    ``pdbx_struct_assembly_prop`` to an mmcif.
    
    :param dict assembly: a dictionary describing an assembly.
    :param dict mmcif: the dictionary to update."""

    assembly_id = str(len(mmcif["pdbx_struct_assembly"]) + 1)
    mmcif["pdbx_struct_assembly"].append({
        "id": assembly_id, "details": "?",
        "method_details": assembly.get("software", "?") or "?",
        "oligomeric_details": "?", "oligomeric_count": "?"
    })
    for prop in ["ABSA (A^2)", "MORE", "SSA (A^2)"]:
        if assembly.get(prop):
            mmcif["pdbx_struct_assembly_prop"].append({
                "biol_id": assembly_id, "type": prop,
                "value": assembly[prop], "details": "?"
            })


def add_assembly_gen_to_mmcif(gen, mmcif):
    """Adds a biological assembly 'gen' to an mmcif dictionary. Any
    transformations not in ``pdbx_struct_oper_list`` will be added, and a
    record for the gen will be added to ``pdbx_struct_assembly_gen``.
    
    :param dict gen: a dictionary describing an assembly.
    :param dict mmcif: the dictionary to update."""

    operation_ids = []
    assembly_id = mmcif["pdbx_struct_assembly"][-1]["id"] 
    for tf in gen["transformations"]:
        op = {
            "type": "?", "type": "?", "name": "?",
            "symmetry_operation": "x,y,z",
        }
        for i in range(3):
            for j in range(3):
                op[f"matrix[{i + 1}][{j + 1}]"] = str(tf["matrix"][i][j])
            op[f"vector[{i + 1}]"] = str(tf["vector"][i])
        for o in mmcif["pdbx_struct_oper_list"]:
            if {k: v for k, v in o.items() if k != "id"} == op:
                operation_ids.append(o["id"])
                break
        else:
            op_id = str(len(mmcif["pdbx_struct_oper_list"]) + 1)
            mmcif["pdbx_struct_oper_list"].append({"id": op_id, **op})
            operation_ids.append(op_id)
    mmcif["pdbx_struct_assembly_gen"].append({
        "assembly_id": assembly_id, "oper_expression": ",".join(operation_ids),
        "asym_id_list": ",".join(gen["chains"])
    })
    

def parse_remark_465(filestring, mmcif):
    """Parses REMARK 465 to produce the ``pdbx_unobs_or_zero_occ_residues``
    category.
    
    :param str filestring: the contents of the .pdb file.
    :param dict mmcif: the dictionary to update."""

    records = re.findall(r"^REMARK 465 .+[A-Z]+.+", filestring, re.M)
    for i, record in enumerate(records):
        if "SSSEQI" in record:
            i += 1
            break
    else:
        i = 5
    if not records[i:]: return
    mmcif["pdbx_unobs_or_zero_occ_residues"] = [{
        "id": str(i), "PDB_model_num": "1", "polymer_flag": "Y",
        "occupancy_flag": "1", "auth_asym_id": rec[19].strip() or "?",
        "auth_comp_id": rec[15:18].strip() or "?",
        "auth_seq_id": rec[21:26].strip() or "?",
        "PDB_ins_code": (rec[26].strip() if len(rec) > 26 else "") or "?",
        "label_asym_id": rec[19].strip() or "?",
        "label_comp_id": rec[15:18].strip() or "?",
        "label_seq_id": rec[21:26].strip() or "?"
    } for i, rec in enumerate(records[i:], start=1)]


def parse_remark_800(filestring, mmcif):
    """Parses REMARK 800 records to produce the ``struct_site`` category.
    
    :param str filestring: the contents of the .pdb file.
    :param dict mmcif: the dictionary to update."""

    records = re.findall(r"^REMARK 800 .+", filestring, re.M)
    records = [r for r in records if r[15:].strip()]
    if not records: return
    site_obj = {
        "id": "?", "pdbx_evidence_code": "?", "pdbx_auth_asym_id": "?",
        "pdbx_auth_comp_id": "?", "pdbx_auth_seq_id": "?",
        "pdbx_auth_ins_code": "?", "pdbx_num_residues": "?",  "details": "?"
    }
    mmcif["struct_site"], site = [], {**site_obj}
    for record in records:
        if "SITE_IDENTIFIER" in record:
            if site != site_obj:
                mmcif["struct_site"].append(site)
                site = {**site_obj}
            site["id"] = record.split(":")[1].strip()
        if "EVIDENCE_CODE" in record:
            site["pdbx_evidence_code"] = record.split(":")[1].strip().title()
        if "SITE_DESCRIPTION" in record:
            site["details"] = record.split(":")[1].strip()
    if site != site_obj: mmcif["struct_site"].append(site)


def parse_cryst1(filestring, mmcif):
    """Parses the CRYST1 record to create the ``symmetry`` and `cell``
    categories.
    
    :param str filestring: the contents of the .pdb file.
    :param dict mmcif: the dictionary to update."""

    cryst1 = re.search(r"^CRYST1.+", filestring, re.M)
    if not cryst1: return
    mmcif["cell"] = [{**{k: "?" for k in [
        "entry_id", "length_a", "length_b", "length_c", "angle_alpha",
        "angle_beta", "angle_gamma", "Z_pdb", "pdbx_unique_axis"
    ]}, "entry_id": mmcif["entry"][0]["id"], }]
    mmcif["symmetry"] = [{
        **{k: "?" for k in [
            "entry_id", "space_group_name_H-M",
            "pdbx_full_space_group_name_H-M", "cell_setting",
            "Int_Tables_number"
        ]},
        "entry_id": mmcif["entry"][0]["id"],
    }]
    rec = cryst1.group(0)
    mmcif["cell"][0]["length_a"] = rec[6:15].strip() or "?"
    mmcif["cell"][0]["length_b"] = rec[15:24].strip() or "?"
    mmcif["cell"][0]["length_c"] = rec[24:33].strip() or "?"
    mmcif["cell"][0]["angle_alpha"] = rec[33:40].strip() or "?"
    mmcif["cell"][0]["angle_beta"] = rec[40:47].strip() or "?"
    mmcif["cell"][0]["angle_gamma"] = rec[47:54].strip() or "?"
    mmcif["cell"][0]["Z_pdb"] = rec[66:70].strip() or "?"
    mmcif["symmetry"][0]["space_group_name_H-M"] = rec[55:66].strip() or "?"


def parse_origx(filestring, mmcif):
    """Parses the ORIGXx records tom produce the ``database_PDB_matrix``
    category.
    
    :param str filestring: the contents of the .pdb file.
    :param dict mmcif: the dictionary to update."""

    lines = re.findall(r"ORIGX.+", filestring, re.M)
    if len(lines) == 0: return
    key = "database_PDB_matrix"
    mmcif[key] = [{"entry_id": mmcif["entry"][0]["id"], **{
        k: "?" for k in [
        "origx[1][1]", "origx[1][2]", "origx[1][3]", "origx[2][1]", 
        "origx[2][2]", "origx[2][3]", "origx[3][1]", "origx[3][2]", 
        "origx[3][3]", "origx_vector[1]", "origx_vector[2]", "origx_vector[3]", 
    ]}}]
    for rec in lines:
        n = rec[5]
        mmcif[key][0][f"origx[{n}][1]"] = rec[10:20].strip() or "?"
        mmcif[key][0][f"origx[{n}][2]"] = rec[20:30].strip() or "?"
        mmcif[key][0][f"origx[{n}][3]"] = rec[30:40].strip() or "?"
        mmcif[key][0][f"origx_vector[{n}]"] = rec[45:55].strip() or "?"


def parse_scale(filestring, mmcif):
    """Parses the SCALEn records to produce the ``atom_sites`` category.
    
    :param str filestring: the contents of the .pdb file.
    :param dict mmcif: the dictionary to update."""

    lines = re.findall(r"^SCALE.+", filestring, re.M)
    if len(lines) == 0: return
    key = "atom_sites"
    mmcif[key] = [{"entry_id": mmcif["entry"][0]["id"], **{
        k: "?" for k in [
        "fract_transf_matrix[1][1]", "fract_transf_matrix[1][2]",
        "fract_transf_matrix[1][3]", "fract_transf_matrix[2][1]", 
        "fract_transf_matrix[2][2]", "fract_transf_matrix[2][3]",
        "fract_transf_matrix[3][1]", "fract_transf_matrix[3][2]", 
        "fract_transf_matrix[3][3]", "fract_transf_vector[1]",
        "fract_transf_vector[2]", "fract_transf_vector[3]", 
    ]}}]
    for r in lines:
        n = r[5]
        if not n.strip(): continue
        mmcif[key][0][f"fract_transf_matrix[{n}][1]"] = r[10:20].strip() or "?"
        mmcif[key][0][f"fract_transf_matrix[{n}][2]"] = r[20:30].strip() or "?"
        mmcif[key][0][f"fract_transf_matrix[{n}][3]"] = r[30:40].strip() or "?"
        mmcif[key][0][f"fract_transf_vector[{n}]"] = r[45:55].strip() or "?"


def parse_mtrix(filestring, mmcif):
    """Parses the MTRIXn records to produce the ``struct_ncs_oper`` category.
    
    :param str filestring: the contents of the .pdb file.
    :param dict mmcif: the dictionary to update."""

    lines = re.findall(r"^MTRIX\d.+", filestring, re.M)
    if len(lines) == 0: return
    matrices = [lines[n * 3:(n + 1) * 3] for n in range(len(lines) // 3)]
    mmcif["struct_ncs_oper"] = [{
        "id": str(n),
        "code": "given" if matrix[0][59] == "1" else "generate",
        "details": "?",
        "matrix[1][1]": matrix[0][10:20].strip() or "?",
        "matrix[1][2]": matrix[0][20:30].strip() or "?",
        "matrix[1][3]": matrix[0][30:40].strip() or "?",
        "matrix[2][1]": matrix[1][10:20].strip() or "?",
        "matrix[2][2]": matrix[1][20:30].strip() or "?",
        "matrix[2][3]": matrix[1][30:40].strip() or "?",
        "matrix[3][1]": matrix[2][10:20].strip() or "?",
        "matrix[3][2]": matrix[2][20:30].strip() or "?",
        "matrix[3][3]": matrix[2][30:40].strip() or "?",
        "vector[1]": matrix[0][45:55].strip() or "?",
        "vector[2]": matrix[1][45:55].strip() or "?",
        "vector[3]": matrix[2][45:55].strip() or "?",
    } for n, matrix in enumerate(matrices, start=1)]


def parse_compnd_and_source(filestring):
    """Gets a dictionary of polymer entities as defined by COMPND and SOURCE
    records. Each entity object should have a CHAINS key, along with other
    information about the entity.
    
    :param str filestring: the contents of the .pdb file.
    :rtype: ``dict``"""

    polymer_entities = {}
    for record_name in ["COMPND", "SOURCE"]:
        records = re.findall(rf"^{record_name}.+", filestring, re.M)
        molecules, molecule = [], ""
        for record in records:
            if "MOL_ID" in record:
                if molecule: molecules.append(molecule.rstrip())
                molecule = ""
            molecule += record[10:].strip() + " "
        if molecule: molecules.append(molecule.rstrip())
        molecules = [parse_entity_string(mol) for mol in molecules]
        for molecule in molecules:
            molecule["id"] = molecule["MOL_ID"]
            del molecule["MOL_ID"]
            if molecule["id"] not in polymer_entities:
                polymer_entities[molecule["id"]] = {}
            polymer_entities[molecule["id"]].update(molecule)
    polymer_entities = {k: v for k, v in polymer_entities.items() if "CHAIN" in v}
    return polymer_entities


def parse_entity_string(string):
    """Takes a molecule string from a COMPND or SOURCE record and parses it into
    a dict, converting values where appropriate.

    :param str string: the string to convert.
    :rtype: ``dict``"""

    molecule = {"molecules": {}}
    entries = [entry.strip() for entry in string.split(";")]
    for entry in entries:
        if not entry.strip(): continue
        key =  entry.split(":")[0].strip()
        value = ":".join([s.strip() for s in entry.split(":")[1:]])
        if value == "YES": value = True
        if value == "NO": value = False
        if key == "CHAIN": value = tuple([
            c.strip() for c in value.split(",")
        ])
        molecule[key] = value
    return molecule


def parse_dbref(filestring):
    """Parses the DBREF records to get external DB info for each polymer
    molecule. A dictionary of polymers is started, with each being given a
    dbrefs list here.

    DBREF1 and DBREF2 records are also parsed.

    :param str filestring: the contents of the .pdb file.
    :rtype: ``dict``"""

    lines = re.findall(r"^DBREF.+", filestring, re.M)
    polymers = {}
    for i, line in enumerate(lines):
        chain_id = line[12]
        if chain_id not in polymers: polymers[chain_id] = {"dbrefs": []}
        if line[5] == "2": continue
        if line[5] == "1":
            line2 = lines[i + 1]
            polymers[chain_id]["dbrefs"].append({
                "start": line[14:18].strip(), "start_insert": line[18:19].strip(),
                "end": line[20:24].strip(), "end_insert": line[24:25].strip(),
                "database": line[26:32].strip(), "accession": line2[18:40].strip(),
                "id": line[47:67].strip(), "db_start": line2[45:55].strip(),
                "db_start_insert": line2[55:56].strip(), "order": i,
                "db_end": line2[57:67].strip(), "db_end_insert": line2[67:68].strip(),
            })
        else:
            polymers[chain_id]["dbrefs"].append({
                "start": line[14:18].strip(), "start_insert": line[18:19].strip(),
                "end": line[20:24].strip(), "end_insert": line[24:25].strip(),
                "database": line[26:32].strip(), "accession": line[33:41].strip(),
                "id": line[42:54].strip(), "db_start": line[55:60].strip(),
                "db_start_insert": line[60:61].strip(), "order": i,
                "db_end": line[62:67].strip(), "db_end_insert": line[67:68].strip(),
            })
    return polymers


def parse_seqadv(filestring, polymers):
    """Parses the SEQADV records and updates a polymers dictionary with any
    sequence modifications for each polymer, in a 'differences' list.
    
    :param str filestring: the contents of the .pdb file.
    :param dict polymers: the polymers dictionary to update."""

    lines = re.findall(r"^SEQADV.+", filestring, re.M)
    for line in lines:
        chain_id = line[16]
        if chain_id not in polymers: polymers[chain_id] = {}
        if "differences" not in polymers[chain_id]:
            polymers[chain_id]["differences"] = []
        polymers[chain_id]["differences"].append({
            "name": line[12:15].strip(), "number": line[18:22].strip(),
            "insert": line[22].strip(), "comment": line[49:70].strip(),
            "database": line[24:28].strip(), "accession": line[29:38].strip(),
            "db_name": line[39:42].strip(), "db_number": line[43:48].strip(),
        })


def parse_seqres(filestring, polymers):
    """Parses the SEQRES records to get the list of residue names for each
    polymer molecule in a polymers dictionary, adding them to a 'residues'
    list.

    :param str filestring: the contents of the .pdb file.
    :param dict polymers: the polymers dictionary to update."""

    lines = re.findall(r"^SEQRES.+", filestring, re.M)
    for line in lines:
        chain_id = line[11]
        if chain_id not in polymers: polymers[chain_id] = {}
        if "residues" not in polymers[chain_id]:
            polymers[chain_id]["residues"] = []
        polymers[chain_id]["residues"] += line[19:70].strip().split()


def parse_modres(filestring, polymers):
    """Parses the MODRES records to get the list of residue modifications for
    each polymer molecule in a polymers dictionary, adding them to a 'modified'
    list.

    :param str filestring: the contents of the .pdb file.
    :param dict polymers: the polymers dictionary to update."""

    lines = re.findall(r"^MODRES.+", filestring, re.M)
    for line in lines:
        chain_id = line[16]
        if chain_id not in polymers: polymers[chain_id] = {}
        if "modified" not in polymers[chain_id]:
            polymers[chain_id]["modified"] = []
        polymers[chain_id]["modified"].append({
            "name": line[12:15].strip(), "number": line[18:22].strip(),
            "insert": line[22].strip(), "standard_name": line[24:27].strip(),
            "comment": line[29:70].strip(),
        })


def parse_het(filestring):
    """Parses the HET records to get the non-polymers in the structure, and the
    entities they belong to. Some of these will only be modified residues, but
    these will be rectified later.

    :param str filestring: the contents of the .pdb file.
    :rtype: ``tuple``"""

    non_polymer_entities, non_polymers = {}, {}
    lines = re.findall(r"^HET .+", filestring, re.M)
    for line in lines:
        name = line[7:10].strip()
        non_polymer_entities[name] = {}
        sig = (line[12], name, line[13:17].strip(), line[17].strip())
        non_polymers[sig] = {}
    return non_polymer_entities, non_polymers


def parse_hetnam(filestring, non_polymer_entities):
    """Parses the HETNAM records to get the names of non-polymers in the
    structure.

    :param str filestring: the contents of the .pdb file.
    :param dict non_polymer_entities: the non-polymers dictionary to update."""

    lines = re.findall(r"^HETNAM.+", filestring, re.M)
    names = [l[11:14].strip() for l in lines]
    for name in names:
        if name not in non_polymer_entities: non_polymer_entities[name] = {}
        non_polymer_entities[name]["name"] = " ".join([
            l[15:70].strip() for l in lines if l[11:14].strip() == name
        ])


def parse_hetsyn(filestring, non_polymer_entities):
    """Parses the HETSYN records to get the synonyms of non-polymers in the
    structure.

    :param str filestring: the contents of the .pdb file.
    :param dict non_polymer_entities: the non-polymers dictionary to update."""

    lines = re.findall(r"^HETSYN.+", filestring, re.M)
    names = [l[11:14].strip() for l in lines]
    for name in names:
        if name not in non_polymer_entities: non_polymer_entities[name] = {}
        non_polymer_entities[name]["synonyms"] = " ".join([
            l[15:70].strip() for l in lines if l[11:14].strip() == name
        ]).split(";")


def parse_formul(filestring, non_polymer_entities):
    """Parses the FORMUL records to get the formulae of non-polymers in the
    structure. This will also determine which non-polymers are water.

    :param str filestring: the contents of the .pdb file.
    :param dict non_polymer_entities: the non-polymers dictionary to update."""

    lines = re.findall(r"^FORMUL.+", filestring, re.M)
    names = [l[12:15].strip() for l in lines]
    for name in names:
        if name not in non_polymer_entities: non_polymer_entities[name] = {}
        match = [l for l in lines if l[12:15].strip() == name]
        non_polymer_entities[name]["formula"] = "".join(
            [l[19:70].strip() for l in match]
        )
        non_polymer_entities[name]["is_water"] = match[0][18] == "*"


def update_entities_from_atoms(filestring, polymer_entities, polymers,
                                non_polymer_entities, non_polymers):
    """Takes a filestring and looks for any entities or molecules that weren't
    in the annotation and adds them. It will also get the 'observed residues'
    for each polymer.
    
    :param dict polymer_entities: the polymer entities dictionary to update.
    :param dict polymers: the polymers dictionary to update.
    :param dict non_polymer_entities: the non-polymer entities to update.
    :param dict non_polymers: the non-polymers dictionary to update."""

    last_sig, sigs, chain_id = None, [], ""
    for rec in re.findall(r"^ATOM.+|^HETATM.+|^TER", filestring, re.M):
        if rec == "TER":
            if chain_id not in polymers: polymers[chain_id] = {}
            polymers[chain_id]["observed_residues"] = [r[1] for r in sigs]
            for entity in polymer_entities.values():
                if chain_id in entity["CHAIN"]: break
            else:
                new_id = int(list(polymer_entities.keys() or [0])[-1]) + 1
                polymer_entities[str(new_id)] = {"CHAIN": (chain_id,)}
            for sig in sigs:
                if sig in non_polymers: del non_polymers[sig]
            sig, sigs = None, []
        else:
            sig = (rec[21], rec[17:20].strip(), rec[22:26].strip(), rec[26].strip())
            if sig != last_sig: sigs.append(sig)
            last_sig = sig
            chain_id = rec[21]
    for sig in sigs:
        if sig[1] not in non_polymer_entities:
            non_polymer_entities[sig[1]] = {}
        if sig not in non_polymers: non_polymers[sig] = {}


def finalize_entities(polymer_entities, non_polymer_entities):
    """Finalizes the determined entities by giving them all the fields they need
    if they haven't been set, and standardising the IDs.
    
    :param dict polymer_entities: the polymer entities to update.
    :param dict non_polymer_entities: the non-polymer entities to update."""

    entity_id = 1
    holding = []
    non_pol_fields = {
        "name": str, "formula": str, "synonyms": list, "molecules": dict
    }
    for entity in polymer_entities.values():
        entity["id"] = str(entity_id)
        if "molecules" not in entity: entity["molecules"] = {}
        entity_id += 1
        holding.append(entity)
    for name, entity in non_polymer_entities.items():
        entity["id"] = str(entity_id)
        entity_id += 1
        for key, func in non_pol_fields.items():
            if key not in entity: entity[key] = func()
            if "is_water" not in entity: entity["is_water"] = name in WATER_NAMES
    for key in list(polymer_entities.keys()):
        del polymer_entities[key]
    for entity in holding:
        polymer_entities[entity["id"]] = entity


def ordered_by_atom(non_polymer_entities, filestring):
    """Orders the non-polymer entities by the order in which they appear in the
    file. This is done by looking at the first ATOM/HETATM record for each.
    
    :param dict non_polymer_entities: the non-polymer entities to update.
    :param str filestring: the contents of the .pdb file."""

    atoms = re.findall(r"^ATOM.+|^HETATM.+", filestring, re.M)
    names = [a[17:20].strip() for a in atoms]
    return {k: non_polymer_entities[k] for k in sorted(
        non_polymer_entities.keys(), key=lambda x: names.index(x)
    )}


def finalize_polymers(polymers):
    """Finalizes the determined polymers by giving them all the fields they need
    if they haven't been set, and giving them all an alignment.
    
    :param dict polymer_entities: the polymer entities to update.
    :param dict non_polymer_entities: the non-polymer entities to update."""

    polymer_fields = {
        "dbrefs": list, "differences": list, "modified": list,
        "residues": list, "observed_residues": list
    }
    for polymer in polymers.values():
        for key, func in polymer_fields.items():
            if key not in polymer: polymer[key] = func()
        if polymer["residues"] == []:
            polymer["residues"] = polymer["observed_residues"]
        polymer["alignment"] = get_alignment_indices(
            polymer["residues"], polymer["observed_residues"]
        )
    

def add_molecules_to_entities(polymer_entities, polymers, non_polymer_entities,
                                non_polymers):
    """Assings molecules to the entity that represents them.
    
    :param dict polymer_entities: the polymer entities dictionary to update.
    :param dict polymers: the polymers dictionary to update.
    :param dict non_polymer_entities: the non-polymer entities to update.
    :param dict non_polymers: the non-polymers dictionary to update."""

    for chain_id, polymer in polymers.items():
        for entity in polymer_entities.values():
            if chain_id in entity["CHAIN"]:
                entity["molecules"][chain_id] = polymer
                break
    for sig, non_polymer in non_polymers.items():
        for name, entity in non_polymer_entities.items():
            if name == sig[1]:
                entity["molecules"][sig] = non_polymer
                break


def build_entity_category(polymer_entities, non_polymer_entities, mmcif):
    """Creates the entity category in an mmCIF dictionary.
    
    :param dict polymer_entities: the polymer entities dictionary.
    :param dict non_polymer_entities: the non-polymer entities.
    :param dict mmcif: the dictionary to update."""

    mmcif["entity"] = []
    entity_template = {
        "id": "?", "type": "?", "src_method": "?", "pdbx_description": "?",
        "formula_weight": "?", "pdbx_number_of_molecules": "1", "pdbx_ec": "?",
        "pdbx_mutation": "?", "pdbx_fragment": "?", "details": "?",
    }
    for entity in polymer_entities.values():
        entity_id = str(len(mmcif["entity"]) + 1)
        mmcif["entity"].append({**entity_template})
        mmcif["entity"][-1]["id"] = entity["id"] = entity_id
        update_polymer_entity(mmcif["entity"][-1], entity)
    for entity in non_polymer_entities.values():
        if not entity["molecules"]: continue
        entity_id = str(len(mmcif["entity"]) + 1)
        mmcif["entity"].append({**entity_template})
        mmcif["entity"][-1]["id"] = entity["id"] = entity_id
        update_non_polymer_entity(mmcif["entity"][-1], entity)
    if not mmcif["entity"]: del mmcif["entity"]


def update_polymer_entity(entity_row, entity_info):
    """Updates a polymer entity with information from the PDB.
    
    :param dict entity_row: the mmCIF entity row.
    :param dict entity_info: the processed PDB info for this entity."""

    entity_row["type"] = "polymer"
    entity_row["pdbx_description"] = entity_info.get("MOLECULE", "?")
    entity_row["pdbx_number_of_molecules"] = str(len(entity_info["CHAIN"]))
    entity_row["pdbx_ec"] = entity_info.get("EC", "?")
    entity_row["pdbx_mutation"] = "?"
    entity_row["pdbx_fragment"] = entity_info.get("FRAGMENT", "?")
    entity_row["details"] = entity_info.get("OTHER_DETAILS", "?")
    entity_row["src_method"] = "syn" if entity_info.get("SYNTHETIC") \
        else "man" if entity_info.get("ENGINEERED") else "nat"


def update_non_polymer_entity(entity_row, entity_info):
    """Updates a non-polymer entity with information from the PDB.
    
    :param dict entity_row: the mmCIF entity row.
    :param dict entity_info: the processed PDB info for this entity."""

    entity_row["type"] = "water" if entity_info["is_water"] else "non-polymer"
    entity_row["pdbx_description"] = entity_info.get("name") or "?"
    if entity_row["pdbx_description"] == "?" and entity_row["type"] == "water":
        entity_row["pdbx_description"] = "water"
    entity_row["pdbx_number_of_molecules"] = str(len(entity_info["molecules"]))
    entity_row["src_method"] = "nat" if entity_info["is_water"] else "syn"


def build_entity_name_com(polymer_entities, mmcif):
    """Creates the entity_name_com category in an mmCIF dictionary.
    
    :param dict polymer_entities: the polymer entities dictionary.
    :param dict non_polymer_entities: the non-polymer entities.
    :param dict mmcif: the dictionary to update."""

    if any(e.get("SYNONYM") for e in polymer_entities.values()):
        mmcif["entity_name_com"] = [{
            "entity_id": entity["id"], "name": entity.get("SYNONYM", "?")
        } for entity in polymer_entities.values() if entity.get("SYNONYM")]


def build_entity_poly(polymer_entities, mmcif):
    """Creates the entity_poly category in an mmCIF dictionary.
    
    :param dict polymer_entities: the polymer entities dictionary.
    :param dict mmcif: the dictionary to update."""

    mmcif["entity_poly"] = []
    for entity in polymer_entities.values():
        if not entity["molecules"]: continue
        polymer = list(entity["molecules"].values())[0]
        nucleocount = len([r for r in polymer["residues"] if r in "ACGTU"])
        is_nucleotide = nucleocount / len(polymer["residues"]) > 0.75
        sequence, can_sequence = get_sequence_strings(polymer)
        mmcif["entity_poly"].append({
            "entity_id": entity["id"],
            "type": "polyribonucleotide" if is_nucleotide else "polypeptide(L)",
            "nstd_linkage": "no",
            "nstd_monomer": "yes" if polymer["modified"] else "no",
            "pdbx_seq_one_letter_code": sequence,
            "pdbx_seq_one_letter_code_can": can_sequence,
            "pdbx_strand_id": ",".join(entity["molecules"].keys()),
            "pdbx_target_identifier": "?"
        })
    if not mmcif["entity_poly"]: del mmcif["entity_poly"]


def get_sequence_strings(polymer):
    """Generates the single-code polymer strings, in both non-canonical and
    canonical forms.
    
    :param dict polymer: a parsed PDB polymer.
    :type: ``tuple``"""

    sequence, can_sequence = [], []
    mod_lookup = {m["name"]: m["standard_name"] for m in polymer["modified"]}
    for res in polymer["residues"]:
        if res in mod_lookup:
            sequence.append(f"({res})")
            can_sequence.append(CODES.get(mod_lookup[res], "X") )
        else:
            sequence.append(CODES.get(res, f"({res})"))
            can_sequence.append(CODES.get(res, "X"))
    return "".join(sequence), "".join(can_sequence)


def build_entity_poly_seq(polymer_entities, mmcif):
    """Creates the entity_poly_seq category in an mmCIF dictionary.
    
    :param dict polymer_entities: the polymer entities dictionary.
    :param dict mmcif: the dictionary to update."""

    mmcif["entity_poly_seq"] = []
    for entity in polymer_entities.values():
        if not entity["molecules"]: continue
        polymer = list(entity["molecules"].values())[0]
        for i, residue in enumerate(polymer["residues"], start=1):
            mmcif["entity_poly_seq"].append({
                "entity_id": entity["id"], "num": str(i),
                "mon_id": residue, "hetero": "n"
            })
    if not mmcif["entity_poly_seq"]: del mmcif["entity_poly_seq"]


def build_struct_ref(polymer_entities, mmcif):
    """Creates the struct_ref category in an mmCIF dictionary.
    
    :param dict polymer_entities: the polymer entities dictionary.
    :param dict mmcif: the dictionary to update."""

    mmcif["struct_ref"] = []
    for entity in polymer_entities.values():
        dbrefs = []
        for molecule in entity["molecules"].values():
            for dbref in molecule["dbrefs"]:
                if dbref["id"] not in [d["id"] for d in dbrefs]:
                    dbrefs.append(dbref)
        for dbref in dbrefs:
            mmcif["struct_ref"].append({
                "id": str(len(mmcif["struct_ref"]) + 1),
                "db_name": dbref["database"] or "?",
                "db_code": dbref["id"] or "?",
                "pdbx_db_accession": dbref["accession"] or "?",
                "pdbx_db_isoform": "?", "entity_id": entity["id"],
                "pdbx_seq_one_letter_code": "?", "pdbx_align_begin": "?",
            })
    if not mmcif["struct_ref"]: del mmcif["struct_ref"]


def build_struct_ref_seq(polymer_entities, mmcif):
    """Creates the struct_ref_seq category in an mmCIF dictionary.
    
    :param dict polymer_entities: the polymer entities dictionary.
    :param dict mmcif: the dictionary to update."""

    mmcif["struct_ref_seq"] = []
    code = mmcif.get("entry", [{"id": ""}])[0]["id"] or "?"
    dbrefs = [{"chain_id": chain_id, **dbref} for e in polymer_entities.values()
        for chain_id, molecule in e["molecules"].items()
            for dbref in molecule["dbrefs"]]
    dbrefs.sort(key=lambda d: d["order"])
    for dbref in dbrefs:
        refs = [r for r in mmcif.get("struct_ref", [])
            if r["pdbx_db_accession"] == dbref["accession"]]
        if not refs: continue
        ref = refs[0]
        mmcif["struct_ref_seq"].append({
            "align_id": str(len(mmcif["struct_ref_seq"]) + 1),
            "ref_id": ref["id"],
            "pdbx_PDB_id_code": code,
            "pdbx_strand_id": dbref["chain_id"], "seq_align_beg": "?",
            "pdbx_seq_align_beg_ins_code": dbref["start_insert"] or "?",
            "seq_align_end": "?",
            "pdbx_seq_align_end_ins_code": dbref["end_insert"] or "?",
            "pdbx_db_accession": dbref["accession"] or "?",
            "db_align_beg": dbref["db_start"] or "?",
            "pdbx_db_align_beg_ins_code": dbref["db_start_insert"] or "?",
            "db_align_end": dbref["db_end"] or "?",
            "pdbx_db_align_end_ins_code": dbref["db_end_insert"] or "?",
            "pdbx_auth_seq_align_beg": dbref["start"] or "?",
            "pdbx_auth_seq_align_end": dbref["end"] or "?"
        })
    if not mmcif["struct_ref_seq"]: del mmcif["struct_ref_seq"]


def build_struct_ref_seq_dif(polymer_entities, mmcif):
    """Creates the struct_ref_seq_dif category in an mmCIF dictionary.
    
    :param dict polymer_entities: the polymer entities dictionary.
    :param dict mmcif: the dictionary to update."""

    mmcif["struct_ref_seq_dif"] = []
    code = mmcif.get("entry", [{"id": ""}])[0]["id"] or "?"
    for entity in polymer_entities.values():
        for mol_id, molecule in entity["molecules"].items():
            for diff in molecule["differences"]:
                align_id = "?"
                for align in mmcif["struct_ref_seq"]:
                    if align["pdbx_strand_id"] == mol_id:
                        if align["pdbx_db_accession"] == diff["accession"]:
                            align_id = align["align_id"]
                            break
                mmcif["struct_ref_seq_dif"].append({
                    "align_id": align_id, "pdbx_pdb_id_code": code,
                    "mon_id": diff["name"] or "?", "pdbx_pdb_strand_id": mol_id,
                    "seq_num": "?", "pdbx_pdb_ins_code": diff["insert"] or "?",
                    "pdbx_seq_db_name": diff["database"] or "?",
                    "pdbx_seq_db_accession_code": diff["accession"] or "?",
                    "db_mon_id": diff["db_name"] or "?",
                    "pdbx_seq_db_seq_num": diff["db_number"] or "?",
                    "details": diff["comment"] or "?",
                    "pdbx_auth_seq_num": diff["number"] or "?",
                    "pdbx_ordinal": str(len(mmcif["struct_ref_seq_dif"]) + 1),
                })
    if not mmcif["struct_ref_seq_dif"]: del mmcif["struct_ref_seq_dif"]


def build_pdbx_struct_mod_residue(polymer_entities, mmcif):
    """Creates the pdbx_struct_mod_residue category in an mmCIF dictionary.
    
    :param dict polymer_entities: the polymer entities dictionary.
    :param dict mmcif: the dictionary to update."""

    mmcif["pdbx_struct_mod_residue"] = []
    for entity in polymer_entities.values():
        for mol_id, molecule in entity["molecules"].items():
            for mod in molecule["modified"]:
                mmcif["pdbx_struct_mod_residue"].append({
                    "id": str(len(mmcif["pdbx_struct_mod_residue"]) + 1),
                    "label_asym_id": "?", "label_comp_id": mod["name"] or "?",
                    "label_seq_id": "?", "auth_asym_id": mol_id,
                    "auth_comp_id": mod["name"] or "?",
                    "auth_seq_id": mod["number"] or "?",
                    "PDB_ins_code": mod["insert"] or "?",
                    "parent_comp_id": mod["standard_name"] or "?",
                    "details": mod["comment"] or "?",
                })
    if not mmcif["pdbx_struct_mod_residue"]:
        del mmcif["pdbx_struct_mod_residue"]


def build_pdbx_entity_nonpoly(non_polymer_entities, mmcif):
    """Creates the pdbx_entity_nonpoly category in an mmCIF dictionary.
    
    :param dict polymer_entities: the polymer entities dictionary.
    :param dict mmcif: the dictionary to update."""

    mmcif["pdbx_entity_nonpoly"] = []
    for name, entity in non_polymer_entities.items():
        if not entity["molecules"]: continue
        entity_row = [e for e in mmcif["entity"] if e["id"] == entity["id"]][0]
        mmcif["pdbx_entity_nonpoly"].append({
            "entity_id": entity["id"], "name": entity_row["pdbx_description"],
            "comp_id": name,
        })
    if not mmcif["pdbx_entity_nonpoly"]: del mmcif["pdbx_entity_nonpoly"]


def build_chem_comp(non_polymer_entities, mmcif):
    """Creates the chem_comp category in an mmCIF dictionary. All non-polymers
    and residues are added.
    
    :param dict polymer_entities: the polymer entities dictionary.
    :param dict mmcif: the dictionary to update."""

    mmcif["chem_comp"] = []
    formulae_lookup, name_lookup = build_ligand_chem_comp(non_polymer_entities, mmcif)
    build_residue_chem_comp(formulae_lookup, name_lookup,  mmcif)
    mmcif["chem_comp"].sort(key=lambda c: c["id"])
    if not mmcif["chem_comp"]: del mmcif["chem_comp"]


def build_ligand_chem_comp(non_polymer_entities, mmcif):
    """Builds the chem_comp items corresponding to ligands. If the 'ligand' is
    actually a modified residue that has no non-polymer representatives, its
    formula is noted and not added here.

    :param dict polymer_entities: the polymer entities dictionary.
    :param dict mmcif: the dictionary to update."""

    formulae_lookup, name_lookup = {}, {}
    for name, entity in non_polymer_entities.items():
        formula = entity["formula"] or "?"
        if re.match(r"^\d+\(", formula):
            formula = formula[formula.find("(") + 1:-1]
        if not entity["molecules"]:
            formulae_lookup[name] = formula
            name_lookup[name] = entity["name"]
            continue
        entity_row = [e for e in mmcif["entity"] if e["id"] == entity["id"]][0]
        mmcif["chem_comp"].append({
            "id": name, "type": "non-polymer", "mon_nstd_flag": ".",
            "name": entity_row["pdbx_description"].upper(),
            "pdbx_synonyms": ",".join(entity.get("synonyms", [])) or "?",
            "formula": formula,
            "formula_weight": f"{formula_to_weight(formula):.3f}"
        })
    return formulae_lookup, name_lookup


def build_residue_chem_comp(formulae_lookup, name_lookup, mmcif):
    """Builds the chem_comp items corresponding to residues. If the residue is
    non-standard, it will use a lookup generated from the ligands to try to
    assign the formula and weight.

    :param dict formulae_lookup: a mapping of modified residue name to formulae.
    :param dict mmcif: the dictionary to update."""

    residues = sorted({r["mon_id"] for r in mmcif.get("entity_poly_seq", [])})
    for res in residues:
        weight = RESIDUE_MASSES.get(res)
        formula = FORMULAE.get(res, formulae_lookup.get(res, "?"))
        if formula != "?" and not weight: weight = formula_to_weight(formula)
        weight = f"{weight:.3f}" if weight else "?"
        name = FULL_NAMES.get(res, name_lookup.get(res, "?")).upper()
        mmcif["chem_comp"].append({
            "id": res, "type": "L-peptide linking",
            "mon_nstd_flag": "y" if res in FORMULAE else "n",
            "name": name, "pdbx_synonyms": "?",
            "formula": formula, "formula_weight": weight
        })


def formula_to_weight(formula):
    """Gets the weight of a formula, ignoring unrecognised symbols.
    
    :param str formula: the formula to calculate the mass for.
    :rtype: ``float``"""

    weight = 0
    atoms = formula.split()
    for atom in atoms:
        symbol, count, in_count = "", "", False
        for char in atom:
            if char in "+-": continue
            if char.isdigit(): in_count = True
            if in_count:
                count += char
            else:
                symbol += char
        mass = PERIODIC_TABLE.get(symbol, 0)
        weight += (mass * int(count or "1"))
    return weight


def build_atom_type(filestring, mmcif):
    """Creates the atom_type category in an mmCIF dictionary.
    
    :param str filestring: the contents of the .pdb file.
    :param dict mmcif: the dictionary to update."""

    lines = re.findall(r"^ATOM.+|^HETATM.+", filestring, re.M)
    elements = sorted(set(line[76:78].strip() or "?" for line in lines))
    if not elements: return
    mmcif["atom_type"] = [{"symbol": element} for element in elements]


def build_atom_site(filestring, polymer_entities, non_polymer_entities, mmcif):
    """Creates the atom_site category in an mmCIF dictionary.
    
    :param str filestring: the contents of the .pdb file.
    :param dict polymer_entities: the polymer entities dictionary.
    :param dict non_polymer_entities: the non-polymer entities.
    :param dict mmcif: the dictionary to update."""

    lines = re.findall(r"^ATOM.+|^HETATM.+|^ENDMDL", filestring, re.M)
    mmcif["atom_site"], labels = [], {}
    model_index, residue_index, current_sig, current_chain = 0, -1, None, None
    for rec in lines:
        if rec == "ENDMDL":
            model_index += 1
            residue_index = -1
            current_chain = None
        else:
            sig = (rec[21], rec[17:20].strip(), rec[22:26].strip(), rec[26].strip())
            if sig != current_sig: residue_index += 1
            current_sig = sig
            if sig[0] != current_chain: residue_index = 0
            current_chain = sig[0]
            polymer_entity, polymer, non_polymer_entity = get_atom_entity(
                sig, polymer_entities, non_polymer_entities
            )
            label = get_atom_label(sig, polymer, non_polymer_entity, labels)
            entity_id = (polymer_entity or non_polymer_entity)["id"]
            seq_id = str(polymer["alignment"][residue_index] + 1) if polymer else "."
            model_num = str(model_index + 1)
            add_atom(mmcif, rec, label, entity_id, seq_id, model_num)
    if not mmcif["atom_site"]: del mmcif["atom_site"]


def get_atom_entity(sig, polymer_entities, non_polymer_entities):
    """Takes an atom signature and the known polymer and non-polymer entities,
    and works out what the polymer entity, polymer and non-polymer entity is for
    that atom (some of which will be None).
    
    :param tuple sig: the atom signature.
    :param dict polymer_entities: the polymer entities.
    :param dict non_polymer_entities: the non-polymer entities.
    :rtype: ``tuple``"""

    polymer_entity, polymer, non_polymer_entity = None, None, None
    for entity in non_polymer_entities.values():
        if sig in entity["molecules"]:
            non_polymer_entity = entity
            break
    else:
        for entity in polymer_entities.values():
            if sig[0] in entity["molecules"]:
                polymer_entity = entity
                polymer = entity["molecules"][sig[0]]
                break
    return polymer_entity, polymer, non_polymer_entity


def get_atom_label(sig, polymer, non_polymer, labels):
    """Works out what the label_asym_id of an atom should be. To do this we need
    to know the signature of the atom, whether it is a polymer, non-polymer or
    water, and what previous label_asym_ids have been set for atoms.
    
    If the atom is a polymer, look to see if the chain ID has been assigned a
    label_asym_id before, and if so use that. Otherwise generate a new one based
    on the current biggest label_asym_id.

    If the atom is a non-polymer, look to see if the signature has been assigned
    label_asym_id before, and if so use that. Otherwise generate a new one based
    on the current biggest label_asym_id.

    If the atom is a water, look to see if the chain ID has been assigned a
    label_asym_id before in the context of water, and if so use that. Otherwise
    generate a new one based on the current biggest label_asym_id.

    The lookup of atom identifier to label_asym_id is modified when a new one is
    created.
    
    :param tuple sig: the atom signature.
    :param dict polymer: the polymer the atom belongs to.
    :param dict non_polymer: the non-polymer the atom belongs to.
    :rtype: ``str``"""
    
    sorted_labels = sorted(labels.values(), key=lambda l: (
        len(l), [ord(c) for c in l][::-1]
    ))
    next_label = next_id(sorted_labels[-1] if labels else "@")
    a = sig[0] if polymer else (sig[0], "W") if non_polymer["is_water"] else sig
    if a in labels:
        label = labels[a]
    else:
        label = next_label
        labels[a] = next_label
    return label


def next_id(id):
    """Gets the next ID after a given ID. The order goes A to Z, then AA, BA
    etc. up to ZA, then AAA, BAA etc.
    
    :param str id: the current ID.
    :rtype: ``str``"""

    if all(char == "Z" for char in id):
        return "A" * (len(id) + 1)
    a_indices = []
    for i, char in enumerate(id):
        if char == "Z":
            a_indices.append(i)
        else:
            return "".join([
                "A" if n in a_indices else chr(ord(c) + 1) if n == i else c
                for n, c in enumerate(id)
            ])


def add_atom(mmcif, line, label, entity_id, seq_id, model_num):
    """Adds an atom the ``atom_site`` category.
    
    :param dict mmcif: the dictionary to update.
    :param str line: the ATOM or HETATM record.
    :param str label: the label_asym_id to use.
    :param str entity_id: the entity ID to use.
    :param str seq_id: the residue mumber to use.
    :param str model_num: the model number to use."""

    sig = (line[21], line[17:20].strip(), line[22:26].strip(), line[26].strip())
    mmcif["atom_site"].append({
        "group_PDB": line[:6].strip(), "id": line[6:11].strip(),
        "type_symbol": line[76:78].strip(), 
        "label_atom_id": line[12:16].strip(),
        "label_alt_id": line[16].strip() or ".", "label_comp_id": sig[1],
        "label_asym_id": label, "label_entity_id": entity_id,
        "label_seq_id": seq_id,
        "pdbx_PDB_ins_code": line[26].strip() or "?",
        "Cartn_x": line[30:38].strip(), "Cartn_y": line[38:46].strip(),
        "Cartn_z": line[46:54].strip(), "occupancy": line[54:60].strip(),
        "B_iso_or_equiv": line[60:66].strip(),
        "pdbx_formal_charge": line[78:80].strip() or "?",
        "auth_seq_id": line[22:26].strip(), "auth_comp_id": line[17:20].strip(),
        "auth_asym_id": sig[0], "auth_atom_id": line[12:16].strip(),
        "pdbx_PDB_model_num": model_num
    })


def build_atom_site_anisotrop(filestring, mmcif):
    """Creates the atom_site_anisotrop category in an mmCIF dictionary.
    
    :param str filestring: the contents of the .pdb file.
    :param dict mmcif: the dictionary to update."""

    anisou = re.findall(r"^ANISOU.+", filestring, re.M)
    if not anisou: return
    mmcif["atom_site_anisotrop"] = []
    atoms_by_id = {a["id"]: a for a in mmcif["atom_site"]}
    for a in anisou:
        atom_id = a[6:11].strip()
        atom = atoms_by_id.get(atom_id)
        if not atom: continue
        convert = lambda s: str(float(s) / 10000)
        mmcif["atom_site_anisotrop"].append({
            "id": atom_id, "type_symbol": atom["type_symbol"], 
            "pdbx_label_atom_id": atom["label_atom_id"], 
            "pdbx_label_alt_id": atom["label_alt_id"], 
            "pdbx_label_comp_id": atom["label_comp_id"], 
            "pdbx_label_asym_id": atom["label_asym_id"], 
            "pdbx_label_seq_id": atom["label_seq_id"], 
            "pdbx_PDB_ins_code": atom["pdbx_PDB_ins_code"], 
            "U[1][1]": convert(a[28:35].strip()),
            "U[2][2]": convert(a[35:42].strip()), 
            "U[3][3]": convert(a[42:49].strip()), 
            "U[1][2]": convert(a[49:56].strip()), 
            "U[1][3]": convert(a[56:63].strip()), 
            "U[2][3]": convert(a[63:70].strip()), 
            "pdbx_auth_seq_id": atom["auth_seq_id"], 
            "pdbx_auth_comp_id": atom["auth_comp_id"], 
            "pdbx_auth_asym_id": atom["auth_asym_id"], 
            "pdbx_auth_atom_id ": atom["auth_atom_id"], 
        })


def update_atom_ids(mmcif):
    """Ensures every atom is numbered sequentially, and the anisotrop category
    is updated likewise too.

    :param dict mmcif: the dictionary to update."""

    lookup = {}
    for i, atom in enumerate(mmcif.get("atom_site", []), start=1):
        lookup[atom["id"]] = str(i)
        atom["id"] = str(i)
    for aniso in mmcif.get("atom_site_anisotrop", []):
        aniso["id"] = lookup[aniso["id"]]


def build_struct_asym(mmcif):
    """Creates the struct_asym category in an mmCIF dictionary.
    
    :param str filestring: the contents of the .pdb file.
    :param dict mmcif: the dictionary to update."""

    lookup = {}
    atoms = mmcif.get("atom_site", [])
    if not atoms: return
    for atom in atoms:
        lookup[atom["label_asym_id"]] = atom["label_entity_id"]
    mmcif["struct_asym"] = [{
        "id": asym_id, "pdbx_blank_PDB_chainid_flag": "N",
        "pdbx_modified": "N", "entity_id": entity_id, "details": "?"
    } for asym_id, entity_id in lookup.items()]


def parse_annotation(filestring, mmcif):
    """Parses the annotation records.
    
    :param str filestring: the contents of the .pdb file.
    :param dict mmcif: the dictionary to update."""

    parse_helix(filestring, mmcif)
    parse_sheet(filestring, mmcif)
    parse_ssbond(filestring, mmcif)
    parse_link(filestring, mmcif)
    parse_cispep(filestring, mmcif)
    parse_site(filestring, mmcif)


def parse_helix(filestring, mmcif):
    """Parses HELIX records to produce the ``struct_conf`` category.
    
    :param str filestring: the contents of the .pdb file.
    :param dict mmcif: the dictionary to update."""

    lines = re.findall(r"^HELIX.+", filestring, re.M)
    if not lines: return
    mmcif["struct_conf"] = [{
        "conf_type_id": "HELX_P", "id": f"HELX_P{str(i)}", 
        "pdbx_PDB_helix_id": str(i), 
        "beg_label_comp_id": line[15:18].strip(), 
        "beg_label_asym_id": line[19].strip(), 
        "beg_label_seq_id": line[21:25].strip(), 
        "pdbx_beg_PDB_ins_code": line[25].strip() or "?", 
        "end_label_comp_id": line[27:30].strip(), 
        "end_label_asym_id": line[31].strip(), 
        "end_label_seq_id": line[33:37].strip(), 
        "pdbx_end_PDB_ins_code": line[37].strip() or "?", 
        "beg_auth_comp_id": line[15:18].strip(),
        "beg_auth_asym_id": line[19].strip(), 
        "beg_auth_seq_id": line[21:25].strip(), 
        "end_auth_comp_id": line[27:30].strip(), 
        "end_auth_asym_id": line[31].strip(), 
        "end_auth_seq_id": line[33:37].strip(), 
        "pdbx_PDB_helix_class": line[38:40].strip() or "?", 
        "details": line[40:70].strip() or "?",
        "pdbx_PDB_helix_length": line[71:76].strip(), 
    } for i, line in enumerate(lines, start=1)]


def parse_sheet(filestring, mmcif):
    """Parses SHEET records to produce the four beta sheet categories.
    
    :param str filestring: the contents of the .pdb file.
    :param dict mmcif: the dictionary to update."""

    lines = re.findall(r"^SHEET.+", filestring, re.M)
    if not lines: return
    mmcif["struct_sheet"] = []
    mmcif["struct_sheet_range"] = []
    for line in lines:
        strand_id, sheet_id = line[7:10].strip(), line[11:14].strip()
        strand_count = line[14:16].strip()
        if sheet_id not in [s["id"] for s in mmcif["struct_sheet"]]:
            mmcif["struct_sheet"].append({
                "id": sheet_id, "type": "?", "number_strands": strand_count,
                "details": "?"
            })
        mmcif["struct_sheet_range"].append({
            "sheet_id": sheet_id, "id": strand_id, 
            "beg_label_comp_id": line[17:20].strip(), 
            "beg_label_asym_id": line[21].strip(), 
            "beg_label_seq_id": line[22:26].strip(),
            "pdbx_beg_PDB_ins_code": line[26].strip() or "?", 
            "end_label_comp_id": line[28:31].strip(),
            "end_label_asym_id": line[32].strip(),
            "end_label_seq_id": line[33:37].strip(),
            "pdbx_end_PDB_ins_code": line[37].strip() or "?", 
            "beg_auth_comp_id": line[17:20].strip(),
            "beg_auth_asym_id": line[21].strip(),
            "beg_auth_seq_id": line[22:26].strip(),
            "end_auth_comp_id": line[28:31].strip(),
            "end_auth_asym_id": line[32].strip(),
            "end_auth_seq_id": line[33:37].strip(),
        })
    parse_inter_strand(lines, mmcif)


def parse_inter_strand(lines, mmcif):
    """Parses SHEET records to produce the two beta sheet categories which
    describe inter-strand relations.
    
    :param list lines: the SHEET lines.
    :param dict mmcif: the dictionary to update."""

    mmcif["struct_sheet_order"] = []
    mmcif["pdbx_struct_sheet_hbond"] = []
    for line in lines:
        line = line.ljust(80)
        strand_id, sheet_id = line[7:10].strip(), line[11:14].strip()
        if strand_id != "1":
            sense = "parallel" if line[38:40].strip() == "1" else "anti-parallel"
            mmcif["struct_sheet_order"].append({
                "sheet_id": sheet_id, "range_id_1": str(int(strand_id) - 1),
                "range_id_2": strand_id, "offset": "?", "sense": sense
            })
            mmcif["pdbx_struct_sheet_hbond"].append({
                "sheet_id": sheet_id,
                "range_id_1": str(int(strand_id) - 1),  "range_id_2": strand_id, 
                "range_1_label_atom_id": line[56:60].strip(), 
                "range_1_label_comp_id": line[60:63].strip(), 
                "range_1_label_asym_id": line[64].strip(), 
                "range_1_label_seq_id": line[65:69].strip(), 
                "range_1_PDB_ins_code": line[69].strip() or "?", 
                "range_1_auth_atom_id": line[56:60].strip(), 
                "range_1_auth_comp_id": line[60:63].strip(), 
                "range_1_auth_asym_id": line[64].strip(), 
                "range_1_auth_seq_id": line[65:69].strip(), 
                "range_2_label_atom_id": line[41:45].strip(), 
                "range_2_label_comp_id": line[45:48].strip(), 
                "range_2_label_asym_id": line[49].strip(), 
                "range_2_label_seq_id": line[33:37].strip(), 
                "range_2_PDB_ins_code": line[54].strip() or "?", 
                "range_2_auth_atom_id": line[41:45].strip(), 
                "range_2_auth_comp_id": line[45:48].strip(), 
                "range_2_auth_asym_id": line[49].strip(), 
                "range_2_auth_seq_id": line[33:37].strip(), 
            })


def parse_ssbond(filestring, mmcif):
    """Parses SSBOND records to produce the ``struct_conn`` category.
    
    :param str filestring: the contents of the .pdb file.
    :param dict mmcif: the dictionary to update."""

    lines = re.findall(r"^SSBOND.+", filestring, re.M)
    if not lines: return
    mmcif["struct_conn"] = []
    for n, line in enumerate(lines, start=1):
        mmcif["struct_conn"].append({
            "id": f"disulf{n}",  "conn_type_id": "disulf", 
            "pdbx_leaving_atom_flag": "?", "pdbx_PDB_id": "?", 
            "ptnr1_label_asym_id": line[15].strip(), 
            "ptnr1_label_comp_id": line[11:14].strip(), 
            "ptnr1_label_seq_id": line[17:21].strip(), 
            "ptnr1_label_atom_id": "?", "pdbx_ptnr1_label_alt_id": "?", 
            "pdbx_ptnr1_PDB_ins_code": line[21].strip() or "?", 
            "pdbx_ptnr1_standard_comp_id": "?", 
            "ptnr1_symmetry": line[59:65].strip(), 
            "ptnr2_label_asym_id": line[29].strip(), 
            "ptnr2_label_comp_id": line[25:28].strip(), 
            "ptnr2_label_seq_id": line[31:35].strip(), 
            "ptnr2_label_atom_id": "?", "pdbx_ptnr2_label_alt_id": "?", 
            "pdbx_ptnr2_PDB_ins_code": line[35].strip() or "?", 
            "ptnr1_auth_asym_id": line[15].strip(), 
            "ptnr1_auth_comp_id": line[11:14].strip(), 
            "ptnr1_auth_seq_id": line[17:21].strip(), 
            "ptnr2_auth_asym_id": line[29].strip(), 
            "ptnr2_auth_comp_id": line[25:28].strip(), 
            "ptnr2_auth_seq_id": line[31:35].strip(), 
            "ptnr2_symmetry": line[66:72].strip(), 
            "pdbx_ptnr3_label_atom_id": "?", "pdbx_ptnr3_label_seq_id": "?", 
            "pdbx_ptnr3_label_comp_id": "?", "pdbx_ptnr3_label_asym_id": "?", 
            "pdbx_ptnr3_label_alt_id": "?", "pdbx_ptnr3_PDB_ins_code": "?", 
            "details": "?", "pdbx_dist_value": line[73:78].strip(), 
            "pdbx_value_order": "?",
        })


def parse_link(filestring, mmcif):
    """Parses LINK records to produce or update the ``struct_conn`` category.
    
    :param str filestring: the contents of the .pdb file.
    :param dict mmcif: the dictionary to update."""

    lines = re.findall(r"^LINK .+", filestring, re.M)
    if not lines: return
    if "struct_conn" not in mmcif: mmcif["struct_conn"] = []
    for n, line in enumerate(lines, start=1):
        mmcif["struct_conn"].append({
            "id": f"covale{n}", "conn_type_id": "covale", 
            "pdbx_leaving_atom_flag": "?", "pdbx_PDB_id": "?", 
            "ptnr1_label_asym_id": line[21].strip(), 
            "ptnr1_label_comp_id": line[17:20].strip(), 
            "ptnr1_label_seq_id": line[22:26].strip(), 
            "ptnr1_label_atom_id": line[12:16].strip(), 
            "pdbx_ptnr1_label_alt_id": line[16].strip() or "?", 
            "pdbx_ptnr1_PDB_ins_code": line[26].strip() or "?", 
            "pdbx_ptnr1_standard_comp_id": "?", 
            "ptnr1_symmetry": line[59:65].strip(), 
            "ptnr2_label_asym_id": line[51].strip(), 
            "ptnr2_label_comp_id": line[47:50].strip(), 
            "ptnr2_label_seq_id": line[52:56].strip(), 
            "ptnr2_label_atom_id": line[42:46].strip(), 
            "pdbx_ptnr2_label_alt_id": line[46].strip() or "?", 
            "pdbx_ptnr2_PDB_ins_code": line[56].strip() or "?", 
            "ptnr1_auth_asym_id": line[21].strip(), 
            "ptnr1_auth_comp_id": line[17:20].strip(), 
            "ptnr1_auth_seq_id": line[22:26].strip(), 
            "ptnr2_auth_asym_id": line[51].strip(), 
            "ptnr2_auth_comp_id": line[47:50].strip(), 
            "ptnr2_auth_seq_id": line[52:56].strip(), 
            "ptnr2_symmetry": line[66:72].strip(), 
            "pdbx_ptnr3_label_atom_id": "?", "pdbx_ptnr3_label_seq_id": "?", 
            "pdbx_ptnr3_label_comp_id": "?", "pdbx_ptnr3_label_asym_id": "?", 
            "pdbx_ptnr3_label_alt_id": "?", "pdbx_ptnr3_PDB_ins_code": "?", 
            "details": "?", "pdbx_dist_value": line[73:78].strip(), 
            "pdbx_value_order": "?",
        })


def parse_site(filestring, mmcif):
    """Parses SITE records to produce the ``struct_site_gen`` category.
    
    :param str filestring: the contents of the .pdb file.
    :param dict mmcif: the dictionary to update."""

    lines = re.findall(r"^SITE.+", filestring, re.M)
    if not lines: return
    mmcif["struct_site_gen"] = []
    for line in lines:
        for n in range(4):
            start = 18 + (11 * n)
            section = line[start:start + 10].ljust(10)
            if not section.strip(): continue
            mmcif["struct_site_gen"].append({
                "id": str(len(mmcif["struct_site_gen"]) + 1), 
                "site_id": line[11:14].strip(), 
                "pdbx_num_res": line[15:17].strip(), 
                "label_comp_id": section[:3].strip(), 
                "label_asym_id": section[4].strip(), 
                "label_seq_id": section[5:9].strip(), 
                "pdbx_auth_ins_code": section[9].strip() or "?", 
                "auth_comp_id": section[:3].strip(), 
                "auth_asym_id": section[4].strip(), 
                "auth_seq_id": section[5:9].strip(), 
                "label_atom_id": ".",  "label_alt_id": "?", 
                "symmetry": "1_555", "details": "?",
            })


def parse_cispep(filestring, mmcif):
    """Parses CISPEP records to produce the ``struct_mon_prot_cis`` category.
    
    :param str filestring: the contents of the .pdb file.
    :param dict mmcif: the dictionary to update."""

    lines = re.findall(r"^CISPEP.+", filestring, re.M)
    if not lines: return
    mmcif["struct_mon_prot_cis"] = [{
        "pdbx_id": line[7:10].strip(), "label_comp_id": line[11:14].strip(),
        "label_seq_id": line[17:21].strip(), "label_asym_id": line[15].strip(),
        "label_alt_id": ".",
        "pdbx_PDB_ins_code": line[21].strip() or "?",
        "auth_comp_id": line[11:14].strip(), "auth_seq_id": line[17:21].strip(),
        "auth_asym_id": line[15].strip(),
        "pdbx_label_comp_id_2": line[25:28].strip(),
        "pdbx_label_seq_id_2": line[31:35].strip(),
        "pdbx_label_asym_id_2": line[29].strip(),
        "pdbx_PDB_ins_code_2": line[35].strip() or "?",
        "pdbx_auth_comp_id_2": line[25:28].strip(),
        "pdbx_auth_seq_id_2": line[31:35].strip(),
        "pdbx_auth_asym_id_2": line[29].strip(),
        "pdbx_PDB_model_num": str(int(line[43:46].strip() or 0) + 1),
        "pdbx_omega_angle": line[53:59].strip()
    } for line in lines]


def update_auth_and_label(mmcif):
    """When parsing the PDB, the chain IDs and sequence numbers will be
    converted to asym IDs and seq IDs for both auth and label variants, but
    really they only correspond to the auth variant. Once the atoms have been
    generated, we know how to fix the label values based on the auth values.

    :param dict mmcif: the dictionary to update."""

    update_auth_ids_in_pdbx_struct_assembly_gen(mmcif)
    auth_to_res = {}
    for a in mmcif["atom_site"]:
        auth_res = (a["auth_asym_id"], a["auth_seq_id"], a["pdbx_PDB_ins_code"])
        label_res = (a["label_asym_id"], a["label_seq_id"])
        auth_to_res[auth_res] = label_res
    for name, rows in mmcif.items():
        keys = list(rows[0].keys())
        templates = []
        if name == "atom_site": continue
        for key in keys:
            for label in ["label_asym_id", "label_seq_id"]:
                if label in key:
                    template = key.split(label)
                    if template in templates: continue
                    templates.append(template)
                    label_asym = key.replace("seq", "asym")
                    label_seq = key.replace("asym", "seq")
                    auth_asym = label_asym.replace("label", "auth")
                    auth_seq = label_seq.replace("label", "auth")
                    auth_ins = [k for k in keys if "ins" in k and
                        template[0] in k and template[1] in k][0]
                    if auth_ins not in keys or auth_seq not in keys or\
                        auth_asym not in keys: continue
                    for row in mmcif[name]:
                        sig = (row[auth_asym], row[auth_seq], row[auth_ins])
                        if sig in auth_to_res:
                            row[label_asym], row[label_seq] = auth_to_res[
                                (row[auth_asym], row[auth_seq], row[auth_ins])
                            ]


def update_auth_ids_in_pdbx_struct_assembly_gen(mmcif):
    """The assembly instructions need the label asym IDs of all molecules
    involved in a particular translation, but in PDB files the non-polymers are
    just treated as part of the associated polymer, so when generating that
    category, only auth asym IDs are present. This replaces them with a fuller
    list if label asym IDs, now that atoms have been generated.
    
    :param dict mmcif: the dictionary to update."""

    label_asym_to_auth_asym = {}
    for atom in mmcif.get("atom_site", []):
        label_asym_to_auth_asym[atom["label_asym_id"]] = atom["auth_asym_id"]
    auth_asym_to_label_asym = {a: [
        k for k, v in label_asym_to_auth_asym.items() if v == a
    ] for a in label_asym_to_auth_asym.values()}
    for row in mmcif.get("pdbx_struct_assembly_gen", []):
        auths = row["asym_id_list"].split(",")
        labels = [id for auth in auths for id in auth_asym_to_label_asym[auth]]
        row["asym_id_list"] = ",".join(sorted(
            labels, key=lambda l: list(label_asym_to_auth_asym.keys()).index(l)
        ))


def pdb_date_to_mmcif_date(date):
    """Turns a date in PDB format into one in mmCIF format. The datetime
    library is not used to avoid locale errors.
    
    :param str date: the date in PDB format.
    :rtype: ``str``"""

    if not date.strip(): return
    day, month, year = date.split("-")
    month = str(list(calendar.month_abbr).index(month.title())).zfill(2)
    if len(year) == 2:
        if int(year) > 50: year = "19" + year
        if int(year) <= 50: year = "20" + year
    return f"{year}-{month}-{day}"


def pdb_names_to_mmcif_names(lines):
    """Turns a set of names in PDB format into one in mmCIF format.
    
    :param list lines: the lines containing names
    :rtype: ``list``"""

    all_names = [name for l in lines for name in l.strip().split(",") if name]
    processed_names = []
    for name in all_names:
        if "." in name and "," not in name:
            names = [n.title() for n in name.split(".")]
            processed_names.append(f"{names[-1]}, {'.'.join(names[:-1])}")
        else: processed_names.append(name.title())
    return processed_names


def save_mmcif_dict(mmcif, path):
    """Saves an mmCIF dictionary to a .pdb file.

    :param dict mmcif: the dictionary to save.
    :param Path path: the location to save to."""

    lines = []
    lines += create_header_line(mmcif)
    lines += create_obslte_lines(mmcif)
    lines += create_title_lines(mmcif)
    lines += create_split_lines(mmcif)
    lines += create_caveat_lines(mmcif)
    lines += create_compnd_lines(mmcif)
    lines += create_source_lines(mmcif)
    lines += create_keywds_lines(mmcif)
    lines += create_expdta_lines(mmcif)
    lines += create_nummdl_lines(mmcif)
    lines += create_mdltyp_lines(mmcif)
    lines += create_author_lines(mmcif)
    lines += create_revdat_lines(mmcif)
    lines += create_sprsde_lines(mmcif)
    lines += create_jrnl_lines(mmcif)
    lines += create_remark_lines(mmcif)
    lines += create_dbref_lines(mmcif)
    lines += create_seqadv_lines(mmcif)
    lines += create_seqres_lines(mmcif)
    lines += create_modres_lines(mmcif)
    lines += create_het_lines(mmcif)
    lines += create_hetnam_lines(mmcif)
    lines += create_hetsyn_lines(mmcif)
    lines += create_formul_lines(mmcif)
    lines += create_helix_lines(mmcif)
    lines += create_sheet_lines(mmcif)
    lines += create_ssbond_lines(mmcif)
    lines += create_link_lines(mmcif)
    lines += create_cispep_lines(mmcif)
    lines += create_site_lines(mmcif)
    lines += create_cryst1_line(mmcif)
    lines += create_origxn_lines(mmcif)
    lines += create_scalen_lines(mmcif)
    lines += create_mtrixn_lines(mmcif)
    lines += create_atom_lines(mmcif)
    lines.append("END")
    with open(path, "w") as f:
        f.write("\n".join(lines))


def create_header_line(mmcif):
    """Creates the HEADER line from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""
    
    code = mmcif.get("entry", [{"id": ""}])[0]["id"]
    keyword = mmcif.get(
        "struct_keywords", [{"pdbx_keywords": ""}]
    )[0]["pdbx_keywords"][:40].upper()
    date = mmcif.get(
        "pdbx_database_status", [{"recvd_initial_deposition_date": ""}]
    )[0]["recvd_initial_deposition_date"]
    code = "    " if code == "?" else code[:4]
    keyword = (" " * 40) if keyword == "?" else keyword.ljust(40)
    date = (" " * 9) if date == "?" else create_pdb_date(date) or ""
    line = "HEADER" + (" " * 4) + keyword + date + (" " * 3) + code
    return [line] if line[6:].strip() else []


def create_obslte_lines(mmcif):
    """Creates the OBSLTE line from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    for row in mmcif.get("pdbx_database_PDB_obs_spr", []):
        if row["id"] == "OBSLTE":
            return ["OBSLTE     {:9} {:4}      {:4}".format(
                create_pdb_date(row["date"]),
                row["replace_pdb_id"], row["pdb_id"]
            )]
    return []


def create_title_lines(mmcif):
    """Creates the TITLE lines from a mmCIF dictionary.

    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    lines = []
    title = mmcif.get("struct", [{"title": "?"}])[0]["title"]
    if not title or title == "?": return []
    strings = split_lines(title.upper(), 60)
    for n, title in enumerate(strings, start=1):
        if n == 1:
            lines.append("TITLE     " + title.upper())
        else:
            lines.append("TITLE   " + str(n).rjust(2) + " " + title.upper())
    return lines


def create_split_lines(mmcif):
    """Creates the SPLIT lines from a mmCIF dictionary.

    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    lines = []
    codes = [s["db_id"][:4] for s in mmcif.get("pdbx_database_related", [])]
    line_count = math.ceil(len(codes) / 14)
    for n in range(line_count):
        lines.append("SPLIT   {:>2} {}".format(
            n + 1 if n else "", " ".join(codes[n * 14: (n + 1) * 14])
        ))
    return lines


def create_caveat_lines(mmcif):
    """Creates the CAVEAT lines from a mmCIF dictionary.

    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    lines = []
    if "database_PDB_caveat" not in mmcif: return []
    text = " ".join([r["text"] for r in mmcif["database_PDB_caveat"]]).upper()
    caveat_lines = split_lines(text, 60)
    for n, line in enumerate(caveat_lines):
        lines.append("CAVEAT  {:>2} {:4}    {}".format(
            "" if n == 0 else n + 1, mmcif["entry"][0]["id"], line
        ))
    return lines


def create_compnd_lines(mmcif):
    """Creates the COMPND lines from a mmCIF dictionary.

    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    lines = []
    for entity in mmcif["entity"]:
        if entity["type"] != "polymer": continue
        name_com = [e for e in mmcif.get("entity_name_com", [])
            if e["entity_id"] == entity["id"]]
        asym_lookup = {a["label_asym_id"]: a["auth_asym_id"]
            for a in mmcif["atom_site"]}
        asyms = [a["id"] for a in mmcif["struct_asym"]
            if a["entity_id"] == entity["id"]]
        mol = {
            "MOL_ID": entity["id"],
            "MOLECULE": entity["pdbx_description"].upper(),
            "CHAIN": ", ".join([asym_lookup[id] for id in asyms]),
            "SYNONYM": name_com[0]["name"].upper() if name_com else "?",
            "EC": entity["pdbx_ec"], "FRAGMENT": entity["pdbx_fragment"],
            "ENGINEERED": "YES" if entity["src_method"] == "man" else "?"
        }
        for key, value in mol.items():
            if value == "?": continue
            compnd = split_lines(f"{key}: {value};", 70)
            for line in compnd:
                if len(lines) == 0:
                    lines.append(f"COMPND    {line}")
                else:
                    lines.append(f"COMPND  {len(lines) + 1:>2} {line}")
    return lines


def create_source_lines(mmcif):
    """Creates the SOURCE lines from a mmCIF dictionary.

    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    lines = []
    for entity in mmcif["entity"]:
        if entity["type"] != "polymer": continue
        mol = {
            "MOL_ID": entity["id"],
            "SYNTHETIC": "YES" if entity["src_method"] == "syn" else "?"
        }
        if any(val != "?" for val in {
            k: v for k, v in mol.items() if k != "MOL_ID"
        }.values()):
            for key, value in mol.items():
                if value == "?": continue
                source = split_lines(f"{key}: {value};", 70)
                for line in source:
                    if len(lines) == 0:
                        lines.append(f"SOURCE    {line}")
                    else:
                        lines.append(f"SOURCE  {len(lines) + 1:>2} {line}")
    return lines


def create_keywds_lines(mmcif):
    """Creates the KEYWDS lines from a mmCIF dictionary.

    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    lines = []
    keywords = mmcif.get("struct_keywords", [{"text": "?"}])[0]["text"]
    if keywords == "?": return lines
    keywds = split_lines(keywords.upper(), 69)
    for line in keywds:
        if len(lines) == 0:
            lines.append(f"KEYWDS    {line}")
        else:
            lines.append(f"KEYWDS  {len(lines):>2} {line}")
    return lines


def create_expdta_lines(mmcif):
    """Creates the EXPDTA lines from a mmCIF dictionary.

    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    values = [r["method"] for r in mmcif.get("exptl", []) if r["method"] != "?"]
    if not values: return []
    text = "; ".join(values)
    expdta = split_lines(text.upper(), 69)
    lines = []
    for line in expdta:
        if len(lines) == 0:
            lines.append(f"EXPDTA    {line}")
        else:
            lines.append(f"EXPDTA  {len(lines):>2} {line}")
    return lines


def create_nummdl_lines(mmcif):
    """Creates the EXPDTA lines from a mmCIF dictionary.

    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    model_ids = set(a["pdbx_PDB_model_num"] for a in mmcif.get("atom_site", []))
    return [] if len(model_ids) == 1 else [f"NUMMDL    {len(model_ids)}"]


def create_mdltyp_lines(mmcif):
    """Creates the MDLTYP lines from a mmCIF dictionary.

    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    messages, lines = [], []
    message = mmcif.get(
        "struct", [{"pdbx_model_type_details": "?"}]
    )[0]["pdbx_model_type_details"]
    if message and message != "?": messages.append(message.upper())
    minimal_chains = []
    for row in  mmcif.get("pdbx_coordinate_model", []):
        if row["type"] not in minimal_chains:
            minimal_chains.append(row["type"])
    for minimal in minimal_chains:
        asyms = ", ".join([
            r["asym_id"] for r in mmcif.get("pdbx_coordinate_model", [])
                if r["type"] == minimal
        ])
        messages.append(f"{minimal.upper()}, CHAIN {asyms}")
    if not messages: return []
    mdltyp = split_lines("; ".join(messages), 60)
    for line in mdltyp:
        if len(lines) == 0:
            lines.append(f"MDLTYP    {line}")
        else:
            lines.append(f"MDLTYP  {len(lines) + 1:>2} {line}")
    return lines


def create_author_lines(mmcif):
    """Creates the AUTHOR lines from a mmCIF dictionary.

    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    lines = []
    names = [r["name"] for r in mmcif.get("audit_author", [])]
    if not names: return lines
    pdb_names = mmcif_names_to_pdb_names(names)
    author = split_lines(pdb_names, 65)
    for line in author:
        if len(lines) == 0:
            lines.append(f"AUTHOR    {line}")
        else:
            lines.append(f"AUTHOR  {len(lines) + 1:>2} {line}")
    return lines


def create_revdat_lines(mmcif):
    """Creates the REVDAT lines from a mmCIF dictionary.

    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    lines = []
    for row in mmcif.get("pdbx_audit_revision_history", []):
        lines.append("REVDAT {:>3}   {:9} {:4}".format(
            len(lines) + 1,
            create_pdb_date(row["revision_date"]),
            mmcif.get("entry", [{"id": ""}])[0]["id"]
        ))
    return lines[::-1]


def create_sprsde_lines(mmcif):
    """Creates the SPRSDE lines from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    for row in mmcif.get("pdbx_database_PDB_obs_spr", []):
        if row["id"] == "SPRSDE":
            return ["SPRSDE     {:9} {:4}      {:4}".format(
                create_pdb_date(row["date"]),
                row["replace_pdb_id"], row["pdb_id"]
            )]
    return []


def create_jrnl_lines(mmcif):
    """Creates the JRNL lines from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    lines = []
    lines += create_jrnl_auth_lines(mmcif)
    lines += create_jrnl_titl_lines(mmcif)
    lines += create_jrnl_edit_lines(mmcif)
    lines += create_jrnl_ref_lines(mmcif)
    lines += create_jrnl_publ_lines(mmcif)
    lines += create_jrnl_refn_lines(mmcif)
    lines += create_jrnl_pmid_lines(mmcif)
    lines += create_jrnl_doi_lines(mmcif)
    return lines


def create_jrnl_auth_lines(mmcif):
    """Creates the JRNL AUTH lines from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    lines = []
    names = [r["name"] for r in mmcif.get("citation_author", [])]
    if not names: return lines
    pdb_names = mmcif_names_to_pdb_names(names)
    author = split_lines(pdb_names, 60)
    for i, line in enumerate(author, start=1):
        lines.append(f"JRNL        AUTH{'' if i == 1 else i:>2} {line}")
    return lines


def create_jrnl_titl_lines(mmcif):
    """Creates the JRNL TITL lines from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    lines = []
    title = mmcif.get("citation", [{"title": "?"}])[0]["title"]
    if not title or title == "?": return lines
    title_lines = split_lines(title.upper(), 50)
    for i, line in enumerate(title_lines, start=1):
        lines.append(f"JRNL        TITL{'' if i == 1 else i:>2} {line}")
    return lines


def create_jrnl_edit_lines(mmcif):
    """Creates the JRNL EDIT lines from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    lines = []
    names = [r["name"] for r in mmcif.get("citation_editor", [])]
    if not names: return lines
    pdb_names = mmcif_names_to_pdb_names(names)
    author = split_lines(pdb_names, 60)
    for i, line in enumerate(author, start=1):
        lines.append(f"JRNL        EDIT{'' if i == 1 else i:>2} {line}")
    return lines


def create_jrnl_ref_lines(mmcif):
    """Creates the JRNL REF lines from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    lines = []
    citation = mmcif.get("citation", [{}])[0]
    journal = citation.get("journal_abbrev", "").replace("?", "")
    volume = citation.get("journal_volume", "").replace("?", "")
    page = citation.get("page_first", "").replace("?", "")
    year = citation.get("year", "").replace("?", "")
    if not (journal or volume or page or year): return lines
    line = "JRNL        REF {:2} {:28}  {:2}{:>4} {:>5} {:4}"
    pub_title_lines = split_lines(journal.upper(), 28)
    for i, title_line in enumerate(pub_title_lines, start=1):
        lines.append(line.format(
            "" if i == 1 else i, title_line,
            "V." if i == 1 and volume else "", volume.upper() if i == 1 else "",
            page if i == 1 else "", year if i == 1 else "",
        ).strip())
    return lines


def create_jrnl_publ_lines(mmcif):
    """Creates the JRNL PUBL lines from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    lines = []
    publ = mmcif.get("citation", [{"book_publisher": ""}])[0]["book_publisher"]
    if not publ or publ == "?": return lines
    publ_lines = split_lines(publ.upper(), 60)
    for i, line in enumerate(publ_lines, start=1):
        lines.append(f"JRNL        PUBL{'' if i == 1 else i:>2} {line}")
    return lines


def create_jrnl_refn_lines(mmcif):
    """Creates the JRNL REFN lines from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    ref = mmcif.get("citation", [{"journal_id_ISSN": ""}])[0]["journal_id_ISSN"]
    if not ref or ref == "?": return []
    return ["JRNL        REFN                   ISSN " + ref]


def create_jrnl_pmid_lines(mmcif):
    """Creates the JRNL PMID lines from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    pmid = mmcif.get(
        "citation", [{"pdbx_database_id_PubMed": ""}]
    )[0]["pdbx_database_id_PubMed"]
    if not pmid or pmid == "?": return []
    return ["JRNL        PMID   " + pmid]


def create_jrnl_doi_lines(mmcif):
    """Creates the JRNL DOI lines from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    doi = mmcif.get(
        "citation", [{"pdbx_database_id_DOI": ""}]
    )[0]["pdbx_database_id_DOI"]
    if not doi or doi == "?": return []
    return ["JRNL        DOI    " + doi]


def create_remark_lines(mmcif):
    """Creates the REMARK lines from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    lines = []
    lines += create_remark_2_lines(mmcif)
    lines += create_remark_3_lines(mmcif)
    lines += create_remark_350_lines(mmcif)
    lines += create_remark_465_lines(mmcif)
    lines += create_remark_800_lines(mmcif)
    return lines


def create_remark_2_lines(mmcif):
    """Creates the REMARK 2 lines from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    r = mmcif.get("reflns", [{"d_resolution_high": ""}])[0]["d_resolution_high"]
    try:
        return [
            "REMARK   2",
            "REMARK   2 RESOLUTION.    {:.2f} ANGSTROMS.".format(float(r))
        ]
    except ValueError: return []


def create_remark_3_lines(mmcif):
    """Creates the REMARK 3 lines from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    if "refine" not in mmcif: return []
    lines = []
    lines += create_remark_3_refinement_target_lines(mmcif)
    lines += create_remark_3_data_used_lines(mmcif)
    lines += create_remark_3_data_fit_lines(mmcif)
    lines += create_remark_3_bvalue_lines(mmcif)
    lines += create_remark_3_thermal_lines(mmcif)
    lines += create_remark_3_solvent_lines(mmcif)
    return lines


def create_remark_3_refinement_target_lines(mmcif):
    """Creates the refinement target REMARK 3 lines from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    target = mmcif["refine"][0].get("pdbx_stereochemistry_target_values", "?")
    target = "NULL" if target == "?" or not target else target
    return ["REMARK   3", "REMARK   3  REFINEMENT TARGET : " + target]


def create_remark_3_data_used_lines(mmcif):
    """Creates the data used REMARK 3 lines from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    refine = mmcif["refine"][0]
    get = lambda k: "NULL" if refine.get(k, "?") == "?" else refine[k]
    res_high = get("ls_d_res_high")
    res_low = get("ls_d_res_low")
    cutoff_high = get("pdbx_data_cutoff_high_absF")
    cutoff_low = get("pdbx_data_cutoff_low_absF")
    complete = get("ls_percent_reflns_obs")
    ref = get("ls_number_reflns_obs")
    return [
        "REMARK   3", "REMARK   3  DATA USED IN REFINEMENT.",                     
        "REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) : " + res_high,                     
        "REMARK   3   RESOLUTION RANGE LOW  (ANGSTROMS) : " + res_low,                  
        "REMARK   3   DATA CUTOFF            (SIGMA(F)) : " + cutoff_low,                  
        "REMARK   3   DATA CUTOFF HIGH         (ABS(F)) : " + cutoff_high,                  
        "REMARK   3   DATA CUTOFF LOW          (ABS(F)) : " + cutoff_low,                  
        "REMARK   3   COMPLETENESS (WORKING+TEST)   (%) : " + complete,                  
        "REMARK   3   NUMBER OF REFLECTIONS             : " + ref,
    ]


def create_remark_3_data_fit_lines(mmcif):
    """Creates the data fit REMARK 3 lines from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    refine = mmcif["refine"][0]
    get = lambda k: "NULL" if refine.get(k, "?") == "?" else refine[k]
    cv = get("pdbx_ls_cross_valid_method")
    select = get("pdbx_R_Free_selection_details")
    rwork = get("ls_R_factor_R_work")
    rfree = get("ls_R_factor_R_free")
    rfreepc = get("ls_percent_reflns_R_free")
    count = get("ls_number_reflns_R_free")
    error = get("ls_R_factor_R_free_error")
    return [
        "REMARK   3", "REMARK   3  FIT TO DATA USED IN REFINEMENT.",                     
        "REMARK   3   CROSS-VALIDATION METHOD          : " + cv,                     
        "REMARK   3   FREE R VALUE TEST SET SELECTION  : " + select,                     
        "REMARK   3   R VALUE            (WORKING SET) : " + rwork,                     
        "REMARK   3   FREE R VALUE                     : " + rfree,                     
        "REMARK   3   FREE R VALUE TEST SET SIZE   (%) : " + rfreepc,                     
        "REMARK   3   FREE R VALUE TEST SET COUNT      : " + count,                     
        "REMARK   3   ESTIMATED ERROR OF FREE R VALUE  : " + error,
    ]


def create_remark_3_bvalue_lines(mmcif):
    """Creates the bvalue REMARK 3 lines from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    refine = mmcif["refine"][0]
    get = lambda k: "NULL" if refine.get(k, "?") == "?" else refine[k]
    wilson = mmcif.get("reflns", [
        {"B_iso_Wilson_estimate": "?"}
    ])[0]["B_iso_Wilson_estimate"]
    if not wilson or wilson == "?": wilson = "NULL"
    mean = get("B_iso_mean")
    b11, b22 = get("aniso_B[1][1]"), get("aniso_B[2][2]")
    b33, b12 = get("aniso_B[3][3]"), get("aniso_B[1][2]")
    b13, b23 = get("aniso_B[1][3]"), get("aniso_B[2][3]")
    return [
        "REMARK   3", "REMARK   3  B VALUES.",                          
        "REMARK   3   FROM WILSON PLOT           (A**2) : " + wilson,                          
        "REMARK   3   MEAN B VALUE      (OVERALL, A**2) : " + mean,                          
        "REMARK   3   OVERALL ANISOTROPIC B VALUE.",                          
        "REMARK   3    B11 (A**2) : " + b11,                          
        "REMARK   3    B22 (A**2) : " + b22,                          
        "REMARK   3    B33 (A**2) : " + b33,                          
        "REMARK   3    B12 (A**2) : " + b12,                          
        "REMARK   3    B13 (A**2) : " + b13,                          
        "REMARK   3    B23 (A**2) : " + b23,
    ]


def create_remark_3_thermal_lines(mmcif):
    """Creates the thermal REMARK 3 lines from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    model = mmcif["refine"][0].get("pdbx_isotropic_thermal_model")
    if not model or model == "?": model = "NULL"
    return ["REMARK   3", "REMARK   3  ISOTROPIC THERMAL MODEL : " + model]


def create_remark_3_solvent_lines(mmcif):
    """Creates the solvent REMARK 3 lines from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    refine = mmcif["refine"][0]
    get = lambda k: "NULL" if refine.get(k, "?") == "?" else refine[k]
    method = get("solvent_model_details")
    ksol = get("solvent_model_param_ksol")
    bsol = get("solvent_model_param_bsol")
    return [
        "REMARK   3",
        "REMARK   3  BULK SOLVENT MODELING.",
        "REMARK   3   METHOD USED : " + method,
        "REMARK   3   KSOL        : " + ksol,
        "REMARK   3   BSOL        : " + bsol,
    ]


def create_remark_350_lines(mmcif):
    """Creates the REMARK 350 lines from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    lines = []
    assemblies = mmcif.get("pdbx_struct_assembly", [])
    if not assemblies: return []
    asym_lookup = {
        a["label_asym_id"]: a["auth_asym_id"] for a in mmcif["atom_site"]
    }
    r3 = "REMARK 350"
    lines = [
        r3, f"{r3} COORDINATES FOR A COMPLETE MULTIMER REPRESENTING THE KNOWN",
        f"{r3} BIOLOGICALLY SIGNIFICANT OLIGOMERIZATION STATE OF THE",
        f"{r3} MOLECULE CAN BE GENERATED BY APPLYING BIOMT TRANSFORMATIONS",
        f"{r3} GIVEN BELOW.  BOTH NON-CRYSTALLOGRAPHIC AND",
        f"{r3} CRYSTALLOGRAPHIC OPERATIONS ARE GIVEN."
    ]
    for assembly in mmcif.get("pdbx_struct_assembly", []):
        lines += create_single_assembly_lines(assembly, mmcif, asym_lookup)
    return lines


def create_single_assembly_lines(assembly, mmcif, asym_lookup):
    """Creates the REMARK 350 lines for a single assembly, from an assembly row
    dictionary, the full mmCIF dictionary, and a mapping of label_asym_ids to
    auth_asym_ids.
    
    :param dict assembly: the assembly dictionary to parse.
    :param dict mmcif: the full dictionary to parse.
    :param dict asym_lookup: mapping of label_asym_ids to auth_asym_ids.
    :rtype: ``list``"""

    lines = []
    lines.append("REMARK 350")
    lines.append("REMARK 350 BIOMOLECULE: " + str(assembly["id"]))
    if assembly.get("method_details", "?") != "?":
        lines.append("REMARK 350 SOFTWARE USED: " + assembly["method_details"])
    properties = [
        ["ABSA (A^2)", "TOTAL BURIED SURFACE AREA", "ANGSTROM**2"],
        ["SSA (A^2)", "SURFACE AREA OF THE COMPLEX", "ANGSTROM**2"],
        ["MORE", "CHANGE IN SOLVENT FREE ENERGY", "KCAL/MOL"],
    ]
    for type, text, units in properties:
        for row in mmcif.get("pdbx_struct_assembly_prop", []):
            if row["biol_id"] == assembly["id"] and row["type"] == type:
                lines.append(f"REMARK 350 {text}: {row['value']} {units}")
    for row in mmcif.get("pdbx_struct_assembly_gen", []):
        if row["assembly_id"] == assembly["id"]:
            lines += create_assembly_gen_lines(row, mmcif, asym_lookup)
    return lines


def create_assembly_gen_lines(gen, mmcif, asym_lookup):
    lines = []
    asym_ids = gen["asym_id_list"].split(",")
    chain_ids = sorted(set(asym_lookup[id] for id in asym_ids))
    chain_lines = split_lines(", ".join(chain_ids), 27)
    for i, line in enumerate(chain_lines):
        if i == 0:
            lines.append(f"REMARK 350 APPLY THE FOLLOWING TO CHAINS: {line},")
        else:
            lines.append(f"REMARK 350                    AND CHAINS: {line},")
    line = "REMARK 350   BIOMT{:<2} {:>2} {} {} {}       {}"
    for op_num, operation in enumerate(get_operations(mmcif, gen), start=1):
        for row_num, row in enumerate(operation[:3], start=1):
            values = [
                ("-" if v < 0 else " ") + "{0:.10f}".format(abs(v))[:8]
                for v in row
            ]
            lines.append(line.format(row_num, op_num, *values))
    return lines


def create_remark_465_lines(mmcif):
    """Creates the REMARK 465 lines from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    residues = mmcif.get("pdbx_unobs_or_zero_occ_residues", [])
    if not residues: return []
    lines = [
        "REMARK 465",
        "REMARK 465 MISSING RESIDUES",
        "REMARK 465 THE FOLLOWING RESIDUES WERE NOT LOCATED IN THE",
        "REMARK 465 EXPERIMENT. (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN",
        "REMARK 465 IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)",
        "REMARK 465",
        "REMARK 465   M RES C SSSEQI"
    ]
    for residue in residues:
        lines.append("REMARK 465     {:3} {:1} {:>5}{:1}".format(
            residue["auth_comp_id"], residue["auth_asym_id"],
            residue["auth_seq_id"], residue["PDB_ins_code"].replace("?", ""),
        ))
    return lines


def create_remark_800_lines(mmcif):
    """Creates the REMARK 800 lines from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    sites = mmcif.get("struct_site", [])
    if not sites: return []
    lines = ["REMARK 800", "REMARK 800 SITE"]
    for site in sites:
        id = site["id"]
        evidence = site["pdbx_evidence_code"]
        details = site["details"]
        if id and id != "?":
            lines.append(f"REMARK 800 SITE_IDENTIFIER: {id.upper()}")
        else:
            continue
        if evidence and evidence != "?":
            lines.append(f"REMARK 800 EVIDENCE_CODE: {evidence.upper()}")
        if details and details != "?":
            lines.append(f"REMARK 800 SITE_DESCRIPTION: {details.upper()}")
    return lines


def create_dbref_lines(mmcif):
    """Creates the DBREF (including DBREF1 and DBREF2) lines from a mmCIF
    dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    lines = []
    code = mmcif.get("entry", [{"id": ""}])[0]["id"]
    line = "DBREF  {:4} {:1} {:>4}{:1} {:>4}{:1} {:6} {:8} {:12} {:>5}{:1} {:>5}{:1}"
    for seq in mmcif.get("struct_ref_seq", []):
        refs = [r for r in mmcif.get("struct_ref", []) if r["id"] == seq["ref_id"]]
        if not refs: continue
        ref = refs[0]
        multi = len(ref["pdbx_db_accession"]) > 8 or len(ref["db_code"]) > 12
        if multi:
            lines += create_dbrefn_lines(code, ref, seq)
        else:
            lines.append(line.format(
                code, seq["pdbx_strand_id"], seq["pdbx_auth_seq_align_beg"],
                seq["pdbx_seq_align_beg_ins_code"].replace("?", ""),
                seq["pdbx_auth_seq_align_end"],
                seq["pdbx_seq_align_end_ins_code"].replace("?", ""),
                ref["db_name"], ref["pdbx_db_accession"], ref["db_code"],
                seq["db_align_beg"],
                seq["pdbx_db_align_beg_ins_code"].replace("?", ""),
                seq["db_align_end"],
                seq["pdbx_db_align_end_ins_code"].replace("?", ""),
            ))
    return lines


def create_dbrefn_lines(code, ref, seq):
    """Creates the DBREF DBREF1 and DBREF2 lines.
    
    :param str code: the ID of the mmCIF.
    :param dict ref: a struct_ref row.
    :param dict mmcif: a struct_ref_seq row.
    :rtype: ``list``"""

    line1 = "DBREF1 {:4} {:1} {:>4}{:1} {:>4}{:1} {:6}               {:10}"
    line2 = "DBREF2 {:4} {:1}     {:22}     {:>10}{:1} {:>10}{:1}"
    lines = []
    lines.append(line1.format(
        code, seq["pdbx_strand_id"], seq["pdbx_auth_seq_align_beg"],
        seq["pdbx_seq_align_beg_ins_code"].replace("?", ""),
        seq["pdbx_auth_seq_align_end"],
        seq["pdbx_seq_align_end_ins_code"].replace("?", ""),
        ref["db_name"], ref["db_code"]
    ))
    lines.append(line2.format(
        code, seq["pdbx_strand_id"], ref["pdbx_db_accession"],
        seq["db_align_beg"], seq["pdbx_db_align_beg_ins_code"].replace("?", ""),
        seq["db_align_end"], seq["pdbx_db_align_end_ins_code"].replace("?", ""),
    ))
    return lines


def create_seqadv_lines(mmcif):
    """Creates the SEQADV lines from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    lines = []
    line = "SEQADV {:4} {:>3} {:1} {:>4}{:1} {:4} {:9} {:>3} {:>5} {:20}"
    for seqadv in mmcif.get("struct_ref_seq_dif", []):
        lines.append(line.format(
            seqadv["pdbx_pdb_id_code"], seqadv["mon_id"],
            seqadv["pdbx_pdb_strand_id"], seqadv["pdbx_auth_seq_num"],
            seqadv["pdbx_pdb_ins_code"].replace("?", ""), 
            seqadv["pdbx_seq_db_name"], seqadv["pdbx_seq_db_accession_code"],
            seqadv["db_mon_id"].replace("?", ""),
            seqadv["pdbx_seq_db_seq_num"].replace("?", ""),
            seqadv["details"].replace("?", ""),
        ).strip())
    return lines


def create_seqres_lines(mmcif):
    """Creates the SEQRES lines from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    lines = []
    line = "SEQRES {:3} {:1} {:4} " + (" {:>3}" * 13)
    poly_seq = mmcif.get("entity_poly_seq", [])
    for entity in mmcif.get("entity_poly", []):
        chain_ids = entity["pdbx_strand_id"].split(",")
        res = [r for r in poly_seq if r["entity_id"] == entity["entity_id"]]
        line_count = math.ceil(len(res) / 13)
        for chain_id in chain_ids:
            for n in range(line_count):
                lineres = res[n * 13:(n + 1) * 13]
                lines.append(line.format(n + 1, chain_id, len(res), *(
                    [r["mon_id"] for r in lineres] + [""] * (13 - len(lineres))
                )).strip())
    return lines


def create_modres_lines(mmcif):
    """Creates the MODRES lines from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    lines = []
    code = mmcif.get("entry", [{"id": ""}])[0]["id"] or ""
    line = "MODRES {:4} {:>3} {:1} {:>4}{:1} {:>3}  {:40}"
    for modres in mmcif.get("pdbx_struct_mod_residue", []):
        lines.append(line.format(
            code, modres["auth_comp_id"], modres["auth_asym_id"],
            modres["auth_seq_id"], modres["PDB_ins_code"].replace("?", ""),
            modres["parent_comp_id"], modres["details"].replace("?", ""),
        ).strip())
    return lines
    

def create_het_lines(mmcif):
    """Creates the HET lines from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    lines = []
    sig_counts = get_sig_counts(mmcif)
    for s, count in sorted(sig_counts.items(), key=lambda s: s[0][3]):
        lines.append(f"HET    {s[3]:3}  {s[0]:1}{s[1]:>4}{s[2]:1}  {count:>5}")
    return lines


def create_hetnam_lines(mmcif):
    """Creates the HETNAM lines from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    lines = []
    if "chem_comp" not in mmcif: return []
    chem_comp = [c for c in mmcif["chem_comp"] if c["mon_nstd_flag"] != "y"]
    lookup = {e["id"]: e["type"] for e in mmcif["entity"]}
    lookup = {s["comp_id"]: lookup[s["entity_id"]] for s in
        mmcif.get("pdbx_entity_nonpoly", [])}
    for chem in chem_comp:
        if chem["name"] != "?" and lookup.get(chem["id"]) != "water":
            strings = split_lines(chem["name"], 55)
            for n, hetnam in enumerate(strings, start=1):
                if n == 1:
                    lines.append(f"HETNAM     {chem['id']:3} " + hetnam)
                else:
                    lines.append(f"HETNAM  {n:2} {chem['id']:3}  " + hetnam)
    return lines


def create_hetsyn_lines(mmcif):
    """Creates the HETSYN lines from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    lines = []
    if "chem_comp" not in mmcif: return []
    chem_comp = [c for c in mmcif["chem_comp"] if c["mon_nstd_flag"] != "y"]
    lookup = {e["id"]: e["type"] for e in mmcif["entity"]}
    lookup = {s["comp_id"]: lookup[s["entity_id"]] for s in
        mmcif.get("pdbx_entity_nonpoly", [])}
    for chem in chem_comp:
        if chem["pdbx_synonyms"] != "?" and lookup.get(chem["id"]) != "water":
            strings = split_lines(chem["pdbx_synonyms"].replace(", ", "; "), 55)
            for n, hetsyn in enumerate(strings, start=1):
                if n == 1:
                    lines.append(f"HETSYN     {chem['id']:3} " + hetsyn)
                else:
                    lines.append(f"HETSYN  {n:2} {chem['id']:3}  " + hetsyn)
    return lines


def create_formul_lines(mmcif):
    """Creates the FORMUL lines from a mmCIF dictionary. A chemical component is
    given a FORMUL line if it isn't a MODRES, and has a formula.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    lines = []
    chem_comp = [c for c in mmcif["chem_comp"] if c["mon_nstd_flag"] != "y"]
    lookup = {e["id"]: e["type"] for e in mmcif["entity"]}
    lookup = {s["comp_id"]: lookup[s["entity_id"]] for s in
        mmcif.get("pdbx_entity_nonpoly", [])}
    name_to_entity_id = {s["comp_id"]: s["entity_id"] for s in
        mmcif.get("pdbx_entity_nonpoly", [])}
    entity_ids = sorted(set(s["entity_id"] for s in mmcif["struct_asym"]))
    sig_counts = get_sig_counts(mmcif, include_water=True, representative=True)
    id_counts = {sig[3]: sum(
        v for k, v in sig_counts.items() if k[3] == sig[3]
    ) for sig in sig_counts}
    chem_comp = sorted(chem_comp, key=lambda c: \
        list(id_counts.keys()).index(c["id"]) if c["id"] in id_counts else 0)
    for chem in chem_comp:
        if chem["formula"] == "?": continue
        entity_id = name_to_entity_id.get(chem["id"])
        num = (entity_ids.index(entity_id) + 1) if entity_id else ""
        char = "*" if lookup.get(chem["id"]) == "water" else " "
        formula = chem["formula"]
        if id_counts.get(chem["id"], 0) > 1:
            formula = f"{id_counts[chem['id']]}({chem['formula']})"
        strings = split_lines(formula, 50)
        for n, formul in enumerate(strings, start=1):
            start = f"FORMUL  {num:2}  {chem['id']:3}"
            if n == 1:
                lines.append(f"{start}   {char}{formul}")
            else:
                lines.append(f"{start} {n:2}{char}{formul}")
    return lines


def get_sig_counts(mmcif, include_water=False, representative=False):
    """Gets all distinct ligands from the atoms, along with how many atoms they
    have.
    
    :param dict mmcif: the dictionary to parse.
    :param bool include_water: whether waters count. 
    :param bool representative: will not count atoms if true. 
    :rtype: ``dict``"""

    sigs = []
    names = [c["id"] for c in mmcif["chem_comp"] if c["mon_nstd_flag"] != "y"]
    lookup = {e["id"]: e["type"] for e in mmcif["entity"]}
    for line in mmcif["atom_site"]:
        sig = (
            line["auth_asym_id"], line["auth_seq_id"],
            line["pdbx_PDB_ins_code"].replace("?", ""), line["auth_comp_id"]
        )
        if not include_water and lookup[line["label_entity_id"]] == "water":
            continue
        if representative and sig in sigs: continue
        if sig[3] in names: sigs.append(sig)
    return Counter(sigs)


def create_helix_lines(mmcif):
    """Creates the HELIX lines from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    line = "HELIX  {:>3} {:>3} {:3} {:1} {:>4}{:1} {:3} " + \
    "{:1} {:>4}{:1}{:>2}{:30} {:>5}"
    return [line.format(
        helix["pdbx_PDB_helix_id"], helix["pdbx_PDB_helix_id"],
        helix["beg_auth_comp_id"], helix["beg_auth_asym_id"],
        helix["beg_auth_seq_id"],
        helix["pdbx_beg_PDB_ins_code"].replace("?", ""),
        helix["end_auth_comp_id"], helix["end_auth_asym_id"],
        helix["end_auth_seq_id"],
        helix["pdbx_end_PDB_ins_code"].replace("?", ""),
        helix["pdbx_PDB_helix_class"], "", helix["pdbx_PDB_helix_length"]
    ) for helix in mmcif.get("struct_conf", [])]


def create_sheet_lines(mmcif):
    """Creates the SHEET lines from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    lines = []
    line = "SHEET  {:>3} {:>3}{:>2} {:>3} {:1}{:>4}{:1} {:>3} {:1}{:>4}{:1}" +\
    "{:>2} {:<4}{:>3} {:1}{:>4}{:1} {:<4}{:>3} {:1}{:>4}{:1}"
    for strand in mmcif.get("struct_sheet_range", []):
        order = [r for r in mmcif["struct_sheet_order"] if
            r["sheet_id"] == strand["sheet_id"] and r["range_id_2"] == strand["id"]]
        order = order[0] if order else None
        hbond = [r for r in mmcif["pdbx_struct_sheet_hbond"] if
            r["sheet_id"] == strand["sheet_id"] and r["range_id_2"] == strand["id"]]
        hbond = hbond[0] if hbond else None
        count = len([r for r in mmcif["struct_sheet_range"] if
            r["sheet_id"] == strand["sheet_id"]])
        lines.append(line.format(
            strand["id"], strand["sheet_id"], str(count),
            strand["beg_auth_comp_id"], strand["beg_auth_asym_id"],
            strand["beg_auth_seq_id"],
            strand["pdbx_beg_PDB_ins_code"].replace("?", ""),
            strand["end_auth_comp_id"], strand["end_auth_asym_id"],
            strand["end_auth_seq_id"],
            strand["pdbx_end_PDB_ins_code"].replace("?", ""),
            0 if not order else 1 if order["sense"] == "parallel" else -1,
            f"{'' if len(hbond['range_2_auth_atom_id']) == 4 else ' '}" + \
            f"{hbond['range_2_auth_atom_id']}" if hbond else "",
            hbond["range_2_auth_comp_id"] if hbond else "",
            hbond["range_2_auth_asym_id"] if hbond else "",
            hbond["range_2_auth_seq_id"] if hbond else "",
            hbond["range_2_PDB_ins_code"].replace("?", "") if hbond else "",
            f"{'' if len(hbond['range_1_auth_atom_id']) == 4 else ' '}" + \
            f"{hbond['range_1_auth_atom_id']}" if hbond else "",
            hbond["range_1_auth_comp_id"] if hbond else "",
            hbond["range_1_auth_asym_id"] if hbond else "",
            hbond["range_1_auth_seq_id"] if hbond else "",
            hbond["range_1_PDB_ins_code"].replace("?", "") if hbond else "",
        ).strip())
    return lines


def create_ssbond_lines(mmcif):
    """Creates the SSBOND lines from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    line = "SSBOND {:>3} {:>3} {:1} {:>4}{:1}   {:>3} {:1} {:>4}{:1}" +\
    "                       {:>6} {:>6} {:>5}"
    lines = []
    for ssbond in mmcif.get("struct_conn", []):
        if ssbond["conn_type_id"] != "disulf": continue
        lines.append(line.format(
            str(len(lines) + 1), ssbond["ptnr1_auth_comp_id"],
            ssbond["ptnr1_auth_asym_id"], ssbond["ptnr1_auth_seq_id"],
            ssbond["pdbx_ptnr1_PDB_ins_code"].replace("?", ""),
            ssbond["ptnr2_auth_comp_id"], ssbond["ptnr2_auth_asym_id"],
            ssbond["ptnr2_auth_seq_id"],
            ssbond["pdbx_ptnr2_PDB_ins_code"].replace("?", ""),
            ssbond["ptnr1_symmetry"], ssbond["ptnr2_symmetry"],
            ssbond["pdbx_dist_value"],
        ))
    return lines


def create_link_lines(mmcif):
    """Creates the LINK lines from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    line = "LINK        {:>4}{:1}{:>3} {:1}{:>4}{:1}               " + \
    "{:>4}{:1}{:>3} {:1}{:>4}{:1}  {:>6} {:>6} {:>5}"
    lines = []
    for link in mmcif.get("struct_conn", []):
        if link["conn_type_id"] != "covale": continue
        lines.append(line.format(
            link["ptnr1_label_atom_id"],
            link["pdbx_ptnr1_label_alt_id"].replace("?", ""),
            link["ptnr1_auth_comp_id"], link["ptnr1_auth_asym_id"],
            link["ptnr1_auth_seq_id"],
            link["pdbx_ptnr1_PDB_ins_code"].replace("?", ""),
            link["ptnr2_label_atom_id"],
            link["pdbx_ptnr2_label_alt_id"].replace("?", ""),
            link["ptnr2_auth_comp_id"], link["ptnr2_auth_asym_id"],
            link["ptnr2_auth_seq_id"],
            link["pdbx_ptnr2_PDB_ins_code"].replace("?", ""),
            link["ptnr1_symmetry"], link["ptnr2_symmetry"],
            link["pdbx_dist_value"],
        ))
    return lines
    

def create_cispep_lines(mmcif):
    """Creates the CISPEP lines from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    line = "CISPEP {:>3} {:>3} {:1} {:>4}{:1}   {:>3} {:1} " +\
    "{:>4}{:1}       {:>3}       {:>6}"
    return [line.format(
        cispep["pdbx_id"], cispep["auth_comp_id"], cispep["auth_asym_id"],
        cispep["auth_seq_id"], cispep["pdbx_PDB_ins_code"].replace("?", ""),
        cispep["pdbx_auth_comp_id_2"], cispep["pdbx_auth_asym_id_2"],
        cispep["pdbx_auth_seq_id_2"],
        cispep["pdbx_PDB_ins_code_2"].replace("?", ""),
        str(int(cispep["pdbx_PDB_model_num"]) - 1), cispep["pdbx_omega_angle"],
    ) for cispep in mmcif.get("struct_mon_prot_cis", [])]


def create_site_lines(mmcif):
    """Creates the SITE lines from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    line = "SITE   {:>3} {:>3} {:>2} {:3} {:1}{:>4}{:1} {:3} {:1}{:>4}{:1} " +\
    "{:3} {:1}{:>4}{:1} {:3} {:1}{:>4}{:1}"
    lines, site_ids = [], []
    for row in mmcif.get("struct_site_gen", []):
        if row["site_id"] not in site_ids: site_ids.append(row["site_id"])
    for site_id in site_ids:
        rows = [r for r in mmcif["struct_site_gen"] if r["site_id"] == site_id]
        row_count = math.ceil(len(rows) / 4)
        for line_num in range(row_count):
            line_rows = rows[line_num * 4:line_num * 4 + 4]
            values = [line_num + 1, site_id, line_rows[0]["pdbx_num_res"]]
            for row in line_rows:
                values += [
                    row["auth_comp_id"], row["auth_asym_id"], row["auth_seq_id"],
                    row["pdbx_auth_ins_code"].replace("?", ""),
                ]
            values += [""] * (19 - len(values))
            lines.append(line.format(*values).strip())
    return lines


def create_cryst1_line(mmcif):
    """Creates the CRYST1 lines from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    if "cell" not in mmcif: return []
    line = "CRYST1{:>9}{:>9}{:>9}{:>7}{:>7}{:>7} {:<11}{:>4}"
    line = line.format(
        mmcif["cell"][0].get("length_a", "?"),
        mmcif["cell"][0].get("length_b", "?"),
        mmcif["cell"][0].get("length_c", "?"),
        mmcif["cell"][0].get("angle_alpha", "?"),
        mmcif["cell"][0].get("angle_beta", "?"),
        mmcif["cell"][0].get("angle_gamma", "?"),
        mmcif.get("symmetry", [{}])[0].get("space_group_name_H-M", "?"),
        mmcif["cell"][0].get("Z_pdb", "?"),
    ).replace("?", " ")
    return [line.strip()] if line[6:].strip() else []


def create_origxn_lines(mmcif):
    """Creates the ORIGX lines from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    if "database_PDB_matrix" not in mmcif: return []
    lines = []
    line = "ORIGX{:1}    {:>10}{:>10}{:>10}     {:>10}"
    for n in range(1, 4):
        lines.append(line.format(
            n, mmcif["database_PDB_matrix"][0][f"origx[{n}][1]"],
            mmcif["database_PDB_matrix"][0][f"origx[{n}][2]"],
            mmcif["database_PDB_matrix"][0][f"origx[{n}][3]"],
            mmcif["database_PDB_matrix"][0][f"origx_vector[{n}]"],
        ).replace("?", " "))
    return lines


def create_scalen_lines(mmcif):
    """Creates the SCALE lines from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    if "atom_sites" not in mmcif: return []
    lines = []
    line = "SCALE{:1}    {:>10}{:>10}{:>10}     {:>10}"
    for n in range(1, 4):
        lines.append(line.format(
            n, mmcif["atom_sites"][0][f"fract_transf_matrix[{n}][1]"],
            mmcif["atom_sites"][0][f"fract_transf_matrix[{n}][2]"],
            mmcif["atom_sites"][0][f"fract_transf_matrix[{n}][3]"],
            mmcif["atom_sites"][0][f"fract_transf_vector[{n}]"],
        ).replace("?", " "))
    return lines


def create_mtrixn_lines(mmcif):
    """Creates the MTRIX lines from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    lines = []
    line = "MTRIX{:1} {:3}{:>10}{:>10}{:>10}     {:>10}    {:1}"
    for i, row in enumerate(mmcif.get("struct_ncs_oper", []), start=1):
        for n in range(1, 4):
            lines.append(line.format(
                n, i, row[f"matrix[{n}][1]"], row[f"matrix[{n}][2]"],
                row[f"matrix[{n}][3]"], row[f"vector[{n}]"],
                "1" if row["code"] == "given" else " "
            ).replace("?", " "))
    return lines


def create_atom_lines(mmcif):
    """Creates the ATOM, HETATM, MODEL, ENDMDL and TER lines from a mmCIF
    dictionary.
    
    :param dict mmcif: the dictionary to parse.
    :rtype: ``list``"""

    if not mmcif.get("atom_site", []): return []
    lines, model_num = [], 0
    atom_id, asym_id, entity_id = 1, mmcif["atom_site"][0]["label_asym_id"], ""
    aniso_lookup = {a["id"]: a for a in mmcif.get("atom_site_anisotrop", [])}
    model_nums = set(a["pdbx_PDB_model_num"] for a in  mmcif["atom_site"])
    lookup = {e["id"]: e["type"] for e in mmcif["entity"]}
    add_ter = lambda: lines.append(f"TER   {atom_id:>5}      {lines[-1][17:26]}")
    for atom in mmcif["atom_site"]:
        if int(atom["pdbx_PDB_model_num"]) > model_num and len(model_nums) > 1:
            if model_num != 0:
                if lookup[entity_id] == "polymer": add_ter()
                lines.append("ENDMDL")
            model_num += 1
            lines.append(f"MODEL     {model_num:>4}")
            atom_id = 1
        if atom["label_asym_id"] != asym_id and lookup[entity_id] == "polymer":
            add_ter()
            atom_id += 1
        asym_id, entity_id = atom["label_asym_id"], atom["label_entity_id"]
        lines.append(create_atom_line(atom, atom_id))
        line = create_aniso_line(atom, aniso_lookup.get(atom["id"]), atom_id)
        if line: lines.append(line)
        is_last = atom == mmcif["atom_site"][-1]
        atom_id += 1
        if is_last and lookup[entity_id] == "polymer": add_ter()
    if len(model_nums) > 1:
        if lookup[entity_id] == "polymer" and not lines[-1].startswith("TER"):
            add_ter()
        lines.append("ENDMDL")
    return lines


def create_atom_line(atom, atom_id):
    """Creates an ATOM or HETATM line from a mmCIF atom_site dictionary.
    
    :param dict atom: the atom to parse.
    :param int atom_id: the atom serial number.
    :rtype: ``str"""

    line = "{:6}{:>5} {:<4} {:3} {:1}{:>4}{:1}   "
    line += "{:>8}{:>8}{:>8}  1.00{:>6}          {:>2}{:2}"
    line = line.format(
        atom["group_PDB"], atom_id,
        f"{'' if len(atom['label_atom_id']) == 4 else ' '}{atom['label_atom_id']}",
        atom["auth_comp_id"], atom["auth_asym_id"],
        atom["auth_seq_id"], atom["pdbx_PDB_ins_code"],
        "{:.3f}".format(float(atom["Cartn_x"])),
        "{:.3f}".format(float(atom["Cartn_y"])),
        "{:.3f}".format(float(atom["Cartn_z"])),
        "{:.2f}".format(float(atom["B_iso_or_equiv"])),
        atom["type_symbol"] or "", atom["pdbx_formal_charge"][::-1],
    )
    return line.replace("?", " ")


def create_aniso_line(atom, aniso, atom_id):
    """Creates an ANISOU line from a mmCIF atom_site dictionary.
    
    :param dict atom: the atom to parse.
    :param dict aniso: the aniso information.
    :param int atom_id: the atom serial number.
    :rtype: ``str"""

    if not aniso: return
    line = "ANISOU{:5} {:4} {:3} {:1}{:>4}{:1} "
    line += "{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}      {:>2}{:2}"
    convert = lambda s: int(float(s) * 10000)
    line = line.format(
        atom_id, f"{'' if len(atom['label_atom_id']) == 4 else ' '}{atom['label_atom_id']}",
        atom["auth_comp_id"],
        atom["auth_asym_id"], atom["auth_seq_id"], atom["pdbx_PDB_ins_code"],
        convert(aniso["U[1][1]"]), convert(aniso["U[2][2]"]),
        convert(aniso["U[3][3]"]), convert(aniso["U[1][2]"]),
        convert(aniso["U[1][3]"]), convert(aniso["U[2][3]"]),
        atom["type_symbol"] or "", atom["pdbx_formal_charge"][::-1],
    )
    return line.replace("?", " ")


def create_pdb_date(date):
    """Takes a mmCIF formatted date and returns a PDB formatted date.
    
    :param str date: the date to parse.
    :rtype: ``str``"""

    if not date.strip(): return
    year, month, day = date.split("-")
    month = list(calendar.month_abbr)[int(month)].upper()
    return f"{day}-{month}-{year[2:]}"


def split_lines(string, length):
    """Splits a string into multiple lines based on a maximum length. The break
    points will be picked based on the location of spaces and commas.
    
    :param str string: the string to split.
    :param int length: the maximum length of each line.
    :rtype: ``list``"""

    strings = []
    has_spaces = " " in string
    while string:
        if len(string) <= length:
            strings.append(string)
            break
        first = string[:length]
        last_space = first[::-1].find(" " if has_spaces else ",")
        if last_space <= length and last_space != -1:
            first = string[:length - last_space]
            string = string[length - last_space:].lstrip()
        else:
            string = string[length:]
        strings.append(first.rstrip())
    return strings


def mmcif_names_to_pdb_names(names):
    """Converts a comma separated list of names from mmCIF format to a list of
    names in PDB format as a single string.
    
    :param str names: the names to convert.
    :rtype: ``str``"""
    
    pdb_names = []
    for name in names:
        surname, initials = name.split(", ")
        if initials.endswith("."): initials = initials[:-1]
        pdb_names.append(f"{initials}.{surname}".upper())
    return ",".join(pdb_names)