import re
import math
import calendar
import itertools
from collections import Counter
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

    mmcif["exptl"] = [{"entry_id": mmcif["entry"][0]["id"], "method": "?"}]
    exptl = re.findall(f"^EXPDTA.+", filestring, re.M)
    if not exptl: return
    mmcif["exptl"][0]["method"] = " ".join([l[10:79].strip() for l in exptl])


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
    authors = process_names(author_lines)
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
        names = process_names(lines)
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
        (r"DATA CUTOFF HIGH.+?\(ABS\(F\)\).+?\:(.+)", "pdbx_data_cutoff_high_rms_absF"),
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
    oper = "pdbx_struct_oper_list"
    for tf in gen["transformations"]:
        op = {
            "type": "?", "pdbx_struct_oper_list.type": "?",
            "pdbx_struct_oper_list.name": "?",
            "pdbx_struct_oper_list.symmetry_operation": "x,y,z",
        }
        for i in range(3):
            for j in range(3):
                op[f"{oper}.matrix[{i + 1}][{j + 1}]"] = str(tf["matrix"][i][j])
            op[f"{oper}.vector[{i + 1}]"] = str(tf["vector"][i])
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
        "asym_id": ",".join(gen["chains"])
    })
    

def parse_remark_465(filestring, mmcif):
    """Parses REMARK 465 to produce the ``pdbx_unobs_or_zero_occ_residues``
    category.
    
    :param str filestring: the contents of the .pdb file.
    :param dict mmcif: the dictionary to update."""

    records = re.findall(r"^REMARK 465 .+", filestring, re.M)
    if not records[6:]: return
    mmcif["pdbx_unobs_or_zero_occ_residues"] = [{
        "id": str(i), "PDB_model_num": "1", "polymer_flag": "Y",
        "occupancy_flag": "1", "auth_asym_id": rec[19].strip() or "?",
        "auth_comp_id": rec[15:18].strip() or "?",
        "auth_seq_id": rec[21:26].strip() or "?",
        "PDB_ins_code": (rec[26].strip() if len(rec) > 26 else "") or "?",
        "label_asym_id": rec[19].strip() or "?",
        "label_comp_id": rec[15:18].strip() or "?",
        "label_seq_id": rec[21:26].strip() or "?"
    } for i, rec in enumerate(records[7:], start=1)]


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
    information about the entity."""

    polymer_entities = {}
    for record_name in ["COMPND", "SOURCE"]:
        records = re.findall(rf"^{record_name}.+", filestring, re.M)
        molecules, molecule = [], ""
        for record in records:
            if "MOL_ID" in record:
                if molecule: molecules.append(molecule)
                molecule = ""
            molecule += record[10:].strip() + " "
        if molecule: molecules.append(molecule)
        molecules = [parse_entity_string(mol) for mol in molecules]
        for molecule in molecules:
            molecule["id"] = molecule["MOL_ID"]
            del molecule["MOL_ID"]
            if molecule["id"] not in polymer_entities:
                polymer_entities[molecule["id"]] = {}
            polymer_entities[molecule["id"]].update(molecule)
    polymer_entities = {k: v for k, v in polymer_entities.items() if "CHAIN" in v}
    return polymer_entities


def parse_entity_string(s):
    """Takes a molecule string from a COMPND or SOURCE record and parses it into
    a dict, converting values where appropriate.

    :param str s: the string to convert.
    :rtype: ``dict``"""

    molecule = {"molecules": {}}
    entries = [entry.strip() for entry in s.split(";")]
    for entry in entries:
        key =  entry.split(":")[0].strip()
        value = ":".join([s.strip() for s in entry.split(":")[1:]])
        if value == "YES":value = True
        if value == "NO": value = False
        if key == "CHAIN": value = tuple([
            c.strip() for c in value.split(",")
        ])
        molecule[key] = value
    return molecule


def parse_dbref(filestring):
    """Parses the DBREF records to get external DB info for each polymer
    molecule. A dictionary of polymers is started, with each being given a
    dbrefs list here."""

    lines = re.findall(r"^DBREF.+", filestring, re.M)
    polymers = {}
    for line in lines:
        chain_id = line[12]
        if chain_id not in polymers: polymers[chain_id] = {"dbrefs": []}
        polymers[chain_id]["dbrefs"].append({
            "start": line[14:18].strip(),
            "start_insert": line[18:19].strip(),
            "end": line[20:24].strip(),
            "end_insert": line[24:25].strip(),
            "database": line[26:32].strip(),
            "accession": line[33:41].strip(),
            "id": line[42:54].strip(),
            "db_start": line[55:60].strip(),
            "db_start_insert": line[60:61].strip(),
            "db_end": line[62:67].strip(),
            "db_end_insert": line[67:68].strip(),
        })
    return polymers


def parse_seqadv(filestring, polymers):
    """Parses the SEQADV records and updates a polymers dictionary with any
    sequence modifications for each polymer, in a 'differences' list."""

    lines = re.findall(r"^SEQADV.+", filestring, re.M)
    for line in lines:
        chain_id = line[16]
        if chain_id not in polymers: polymers[chain_id] = {}
        if "differences" not in polymers[chain_id]:
            polymers[chain_id]["differences"] = []
        polymers[chain_id]["differences"].append({
            "name": line[12:15].strip(),
            "number": line[18:22].strip(),
            "insert": line[22].strip(),
            "database": line[24:28].strip(),
            "accession": line[29:38].strip(),
            "db_name": line[39:42].strip(),
            "db_number": line[43:48].strip(),
            "comment": line[49:70].strip(),
        })


def parse_seqres(filestring, polymers):
    """Parses the SEQRES records to get the list of residue names for each
    polymer molecule in a polymers dictionary, adding them to a 'residues'
    list."""

    lines = re.findall(r"^SEQRES.+", filestring, re.M)
    for line in lines:
        chain_id = line[11]
        if chain_id not in polymers: polymers[chain_id] = {}
        if "residues" not in polymers[chain_id]:
            polymers[chain_id]["residues"] = []
        polymers[chain_id]["residues"] += line[19:70].strip().split()


def parse_modres(filestring, polymers):
    lines = re.findall(r"^MODRES.+", filestring, re.M)
    for line in lines:
        chain_id = line[16]
        if chain_id not in polymers: polymers[chain_id] = {}
        if "modified" not in polymers[chain_id]:
            polymers[chain_id]["modified"] = []
        polymers[chain_id]["modified"].append({
            "name": line[12:15].strip(),
            "number": line[18:22].strip(),
            "insert": line[22].strip(),
            "standard_name": line[24:27].strip(),
            "comment": line[29:70].strip(),
        })


def parse_het(filestring):
    non_polymer_entities, non_polymers = {}, {}
    lines = re.findall(r"^HET .+", filestring, re.M)
    for line in lines:
        name = line[7:10].strip()
        non_polymer_entities[name] = {}
        sig = (line[12], name, line[13:17].strip(), line[17].strip())
        non_polymers[sig] = {}
    return non_polymer_entities, non_polymers


def parse_hetnam(filestring, non_polymer_entities):
    lines = re.findall(r"^HETNAM.+", filestring, re.M)
    names = [l[11:14].strip() for l in lines]
    for name in names:
        if name not in non_polymer_entities: non_polymer_entities[name] = {}
        non_polymer_entities[name]["name"] = " ".join([
            l[15:70].strip() for l in lines if l[11:14].strip() == name
        ])


def parse_hetsyn(filestring, non_polymer_entities):
    lines = re.findall(r"^HETSYN.+", filestring, re.M)
    names = [l[11:14].strip() for l in lines]
    for name in names:
        if name not in non_polymer_entities:
            non_polymer_entities[name] = {"name": ""}
        non_polymer_entities[name]["synonyms"] = " ".join([
            l[15:70].strip() for l in lines if l[11:14].strip() == name
        ]).split(";")


def parse_formul(filestring, non_polymer_entities):
    lines = re.findall(r"^FORMUL.+", filestring, re.M)
    names = [l[12:15].strip() for l in lines]
    for name in names:
        if name not in non_polymer_entities:
            non_polymer_entities[name] = {"name": "", "synonyms": []}
        match = [l for l in lines if l[12:15].strip() == name]
        non_polymer_entities[name]["formula"] = "".join(
            [l[19:70].strip() for l in match]
        )
        non_polymer_entities[name]["is_water"] = match[0][18] == "*"


def update_entities_from_atoms(filestring, polymer_entities, polymers, non_polymer_entities, non_polymers):
    lines = re.findall(r"^ATOM.+|^HETATM.+|^TER", filestring, re.M)
    sig, sigs, chain_id = None, [], ""
    for line in lines:
        if line == "TER":
            # All residues so far have been for a polymer
            if chain_id not in polymers:
                polymers[chain_id] = {}
            polymers[chain_id]["observed_residues"] = [r[1] for r in sigs]
            for entity in polymer_entities.values():
                if chain_id in entity["CHAIN"]: break
            else:
                new_entity_id = int(list(polymer_entities.keys() or [0])[-1]) + 1
                polymer_entities[str(new_entity_id)] = {"CHAIN": (chain_id,)}
            for sig in sigs:
                if sig in non_polymers: del non_polymers[sig]
            sig, sigs = None, []
        else:
            # Make a note of this residue
            line_sig = residue_sig(line)
            if line_sig != sig: sigs.append(line_sig)
            sig = line_sig
            chain_id = line[21]
    # Any residues left are non-polymers
    for order, sig in enumerate(sigs):
        if sig[1] not in non_polymer_entities:
            non_polymer_entities[sig[1]] = {"order": order}
        else:
             non_polymer_entities[sig[1]]["order"] = order
        if sig not in non_polymers: non_polymers[sig] = {}


def finalize_entities(polymer_entities, non_polymer_entities):
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


def finalize_polymers(polymers):
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
    

def add_molecules_to_entities(polymer_entities, polymers, non_polymer_entities, non_polymers):
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
    '''non_pol_without_molecules = []
    for id, entity in non_polymer_entities.items():
        if not entity["molecules"]: non_pol_without_molecules.append(id)
    for id in non_pol_without_molecules: del non_polymer_entities[id]'''


def build_entity_category(polymer_entities, non_polymer_entities, mmcif):
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
        mmcif["entity"][-1]["type"] = "polymer"
        mmcif["entity"][-1]["pdbx_description"] = entity.get("MOLECULE", "?")
        mmcif["entity"][-1]["pdbx_number_of_molecules"] = str(len(entity["CHAIN"]))
        mmcif["entity"][-1]["pdbx_ec"] = entity.get("EC", "?")
        mmcif["entity"][-1]["pdbx_mutation"] = "?"
        mmcif["entity"][-1]["pdbx_fragment"] = entity.get("FRAGMENT", "?")
        mmcif["entity"][-1]["details"] = entity.get("OTHER_DETAILS", "?")
        mmcif["entity"][-1]["src_method"] = "syn" if entity.get("SYNTHETIC") \
            else "man" if entity.get("ENGINEERED") else "nat"
    for entity in non_polymer_entities.values():
        if not entity["molecules"]: continue
        entity_id = str(len(mmcif["entity"]) + 1)
        mmcif["entity"].append({**entity_template})
        mmcif["entity"][-1]["id"] = entity["id"] = entity_id
        mmcif["entity"][-1]["type"] = "water" if entity["is_water"] else "non-polymer"
        mmcif["entity"][-1]["pdbx_description"] = entity.get("name") or "?"
        if mmcif["entity"][-1]["pdbx_description"] == "?" and mmcif["entity"][-1]["type"] == "water":
            mmcif["entity"][-1]["pdbx_description"] = "water"
        mmcif["entity"][-1]["pdbx_number_of_molecules"] = str(len(entity["molecules"]))
        mmcif["entity"][-1]["src_method"] = "nat" if entity["is_water"] else "syn"
    if any(e.get("SYNONYM") for e in polymer_entities.values()):
        mmcif["entity_name_com"] = [{
            "entity_id": entity["id"], "name": entity.get("SYNONYM", "?")
        } for entity in polymer_entities.values()]


def build_entity_poly(polymer_entities, mmcif):
    mmcif["entity_poly"] = []
    for entity in polymer_entities.values():
        if not entity["molecules"]: continue
        polymer = list(entity["molecules"].values())[0]
        nucleocount = len([r for r in polymer["residues"] if r in "ACGTU"])
        is_nucleotide = nucleocount / len(polymer["residues"]) > 0.75
        sequence, can_sequence = [], []
        modified_lookup = {m["name"]: m["standard_name"] for m in polymer["modified"]}
        for res in polymer["residues"]:
            if res in modified_lookup:
                sequence.append(f"({res})")
                can_sequence.append(modified_lookup[res])
            else:
                sequence.append(CODES.get(res, f"({res})"))
                can_sequence.append(CODES.get(res, "X"))
        mmcif["entity_poly"].append({
            "entity_id": entity["id"],
            "type": "polyribonucleotide" if is_nucleotide else "polypeptide(L)",
            "nstd_linkage": "no", "nstd_monomer": "yes" if polymer["modified"] else "no",
            "pdbx_seq_one_letter_code": "".join(sequence),
            "pdbx_seq_one_letter_code_can": "".join(can_sequence),
            "pdbx_strand_id": ",".join(entity["CHAIN"]),
            "pdbx_target_identifier": "?"
        })
    if not mmcif["entity_poly"]: del mmcif["entity_poly"]


def build_entity_poly_seq(polymer_entities, mmcif):
    mmcif["entity_poly_seq"] = []
    for entity in polymer_entities.values():
        if not entity["molecules"]: continue
        polymer = list(entity["molecules"].values())[0]
        for i, residue in enumerate(polymer["residues"], start=1):
            mmcif["entity_poly_seq"].append({
                "entity_id": entity["id"],
                "num": str(i),
                "mon_id": residue,
                "hetero": "n"
            })
    if not mmcif["entity_poly_seq"]: del mmcif["entity_poly_seq"]


def build_struct_ref(polymer_entities, mmcif):
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
                "db_name": dbref["database"], "db_code": dbref["id"],
                "entity_id": entity["id"], "pdbx_seq_one_letter_code": "?",
                "pdbx_align_begin": dbref["start"],
                "pdbx_db_accession": dbref["accession"], "pdbx_db_isoform": "?"
            })
    if not mmcif["struct_ref"]: del mmcif["struct_ref"]


def build_struct_ref_seq(polymer_entities, mmcif):
    mmcif["struct_ref_seq"] = []
    for entity in polymer_entities.values():
        for molecule in entity["molecules"].items():
            chain_id, molecule = molecule
            for dbref in molecule["dbrefs"]:
                mmcif["struct_ref_seq"].append({
                    "align_id": str(len(mmcif["struct_ref_seq"]) + 1), "ref_id": "1",
                    "pdbx_PDB_id_code": mmcif.get("entry", [{"id": ""}])[0]["id"] or "?",
                    "pdbx_strand_id": chain_id,
                    "seq_align_beg": "?",
                    "pdbx_seq_align_beg_ins_code": dbref["start_insert"] or "?",
                    "seq_align_end": "?",
                    "pdbx_seq_align_end_ins_code": dbref["end_insert"] or "?",
                    "pdbx_db_accession": dbref["accession"],
                    "db_align_beg": dbref["db_start"],
                    "pdbx_db_align_beg_ins_code": dbref["db_start_insert"] or "?",
                    "db_align_end": dbref["db_end"],
                    "pdbx_db_align_end_ins_code": dbref["db_end_insert"] or "?",
                    "pdbx_auth_seq_align_beg": dbref["start"],
                    "pdbx_auth_seq_align_end": dbref["end"]
                })
    if not mmcif["struct_ref_seq"]: del mmcif["struct_ref_seq"]


def build_struct_ref_seq_dif(polymer_entities, mmcif):
    mmcif["struct_ref_seq_dif"] = []
    for entity in polymer_entities.values():
        for mol_id, molecule in entity["molecules"].items():
            for diff in molecule["differences"]:
                for align in mmcif["struct_ref_seq"]:
                    if align["pdbx_strand_id"] == mol_id:
                        align_id = align["align_id"]
                        break
                else:
                    align_id = "?"
                mmcif["struct_ref_seq_dif"].append({
                    "align_id": align_id,
                    "pdbx_pdb_id_code":  mmcif.get("entry", [{"id": ""}])[0]["id"] or "?",
                    "mon_id": diff["name"] or "?",
                    "pdbx_pdb_strand_id": mol_id,
                    "seq_num": "?",
                    "pdbx_pdb_ins_code": diff["insert"] or "?",
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
    mmcif["pdbx_struct_mod_residue"] = []
    for entity in polymer_entities.values():
        for mol_id, molecule in entity["molecules"].items():
            for mod in molecule["modified"]:
                mmcif["pdbx_struct_mod_residue"].append({
                    "id": str(len(mmcif["pdbx_struct_mod_residue"]) + 1),
                    "label_asym_id": "?",
                    "label_comp_id": mod["name"],
                    "label_seq_id": "?",
                    "auth_asym_id": mol_id,
                    "auth_comp_id": mod["name"],
                    "auth_seq_id": mod["number"],
                    "PDB_ins_code": mod["insert"] or "?",
                    "parent_comp_id": mod["standard_name"],
                    "details": mod["comment"] or "?",
                })
    if not mmcif["pdbx_struct_mod_residue"]: del mmcif["pdbx_struct_mod_residue"]



def build_pdbx_entity_nonpoly(non_polymers, mmcif):
    mmcif["pdbx_entity_nonpoly"] = []
    for name, entity in non_polymers.items():
        if not entity["molecules"]: continue
        entity_row = [e for e in mmcif["entity"] if e["id"] == entity["id"]][0]
        mmcif["pdbx_entity_nonpoly"].append({
            "entity_id": entity["id"],
            "name": entity_row["pdbx_description"],
            "comp_id": name,
        })
    if not mmcif["pdbx_entity_nonpoly"]: del mmcif["pdbx_entity_nonpoly"]


def build_chem_comp(non_polymer_entities, mmcif):
    mmcif["chem_comp"] = []
    residues = {r["mon_id"] for r in mmcif.get("entity_poly_seq", [])}
    formulae_lookup = {}
    for name, entity in non_polymer_entities.items():
        if not entity["molecules"]:
            formulae_lookup[name] = entity["formula"]
            continue
        entity_row = [e for e in mmcif["entity"] if e["id"] == entity["id"]][0]
        formula = entity["formula"] or "?"
        if re.match(r"^\d+\(", formula): formula = formula[formula.find("(") + 1:-1]
        weight = "?"
        mmcif["chem_comp"].append({
            "id": name, "type": "non-polymer",
            "mon_nstd_flag": ".", "name": entity_row["pdbx_description"].upper(),
            "pdbx_synonyms": ",".join(entity.get("synonyms", [])) or "?",
            "formula": formula,
            "formula_weight": f"{formula_to_weight(formula):.3f}"
        })
    for res in residues:
        weight = RESIDUE_MASSES.get(res)
        weight = f"{weight:.3f}" if weight else "?"
        mmcif["chem_comp"].append({
            "id": res, "type": "L-peptide linking",
            "mon_nstd_flag": "y" if res in FORMULAE else "n",
            "name": FULL_NAMES.get(res, "?").upper(),
            "pdbx_synonyms": "?", "formula": FORMULAE.get(res, formulae_lookup.get(res, "?")),
            "formula_weight": weight
        })
    mmcif["chem_comp"].sort(key=lambda c: c["id"])


def formula_to_weight(formula):
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
    lines = re.findall(r"^ATOM.+|^HETATM.+", filestring, re.M)
    elements = sorted(set(line[76:78].strip() or "?" for line in lines))
    mmcif["atom_type"] = [{"symbol": element} for element in elements]


def build_atom_site(filestring, polymer_entities, non_polymer_entities, mmcif):
    lines = re.findall(r"^ATOM.+|^HETATM.+|^ENDMDL", filestring, re.M)
    mmcif["atom_site"] = []
    labels = {}
    model_index, residue_index, current_sig, current_chain_id = 0, -1, None, None
    for line in lines:
        if line == "ENDMDL":
            model_index += 1
            residue_index = -1
            current_chain_id = None
        else:
            sig = residue_sig(line)
            if sig != current_sig: residue_index += 1
            current_sig = sig
            if sig[0] != current_chain_id: residue_index = 0
            current_chain_id = sig[0]
            polymer_entity, polymer, non_polymer_entity = get_atom_entity(
                sig, polymer_entities, non_polymer_entities
            )
            label = get_atom_label(sig, polymer, non_polymer_entity, labels)
            mmcif["atom_site"].append({
                "group_PDB": line[:6].strip(),
                "id": line[6:11].strip(),
                "type_symbol": line[76:78].strip(),
                "label_atom_id": line[12:16].strip(),
                "label_alt_id": line[16].strip() or ".",
                "label_comp_id": sig[1],
                "label_asym_id": label,
                "label_entity_id": (polymer_entity or non_polymer_entity)["id"],
                "label_seq_id": str(polymer["alignment"][residue_index] + 1) if polymer else ".",
                "pdbx_PDB_ins_code": line[26].strip() or "?",
                "Cartn_x": line[30:38].strip(),
                "Cartn_y": line[38:46].strip(),
                "Cartn_z": line[46:54].strip(),
                "occupancy": line[54:60].strip(),
                "B_iso_or_equiv": line[60:66].strip(),
                "pdbx_formal_charge": line[78:80].strip() or "?",
                "auth_seq_id": line[22:26].strip(),
                "auth_comp_id": line[17:20].strip(),
                "auth_asym_id": sig[0],
                "auth_atom_id": line[12:16].strip(),
                "pdbx_PDB_model_num": str(model_index + 1)
            })


def get_atom_entity(sig, polymer_entities, non_polymer_entities):
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
    sorted_labels = sorted(labels.values(), key=lambda l: (len(l), [ord(c) for c in l][::-1]))
    next_label = next_id(sorted_labels[-1] if labels else "@")
    if polymer:
        if sig[0] in labels:
            label = labels[sig[0]]
        else:
            label = next_label
            labels[sig[0]] = next_label
    elif not non_polymer["is_water"]:
        if sig in labels:
            label = labels[sig]
        else:
            label = next_label
            labels[sig] = next_label
    else:
        chain_id = (sig[0], "W")
        if chain_id in labels:
            label = labels[chain_id]
        else:
            label = next_label
            labels[chain_id] = next_label
    return label


def next_id(id):
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


def build_atom_site_anisotrop(filestring, mmcif):
    anisou = re.findall(r"^ANISOU.+", filestring, re.M)
    mmcif["atom_site_anisotrop"] = []
    atoms_by_id = {a["id"]: a for a in mmcif["atom_site"]}
    for a in anisou:
        atom_id = a[6:11].strip()
        atom = atoms_by_id[atom_id]
        convert = lambda s: str(float(s) / 10000)
        mmcif["atom_site_anisotrop"].append({
            "id": atom_id, 
            "type_symbol": atom["type_symbol"], 
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
    if not mmcif["atom_site_anisotrop"]: del mmcif["atom_site_anisotrop"]


def update_atom_ids(mmcif):
    lookup = {}
    for i, atom in enumerate(mmcif["atom_site"], start=1):
        lookup[atom["id"]] = str(i)
        atom["id"] = str(i)
    for aniso in mmcif.get("atom_site_anisotrop", []):
        aniso["id"] = lookup[aniso["id"]]


def build_struct_asym(mmcif):
    lookup = {}
    for atom in mmcif["atom_site"]:
        lookup[atom["label_asym_id"]] = atom["label_entity_id"]
    mmcif["struct_asym"] = [{
        "id": asym_id, "pdbx_blank_PDB_chainid_flag": "N",
        "pdbx_modified": "N", "entity_id": entity_id, "details": "?"
    } for asym_id, entity_id in lookup.items()]


def parse_annotation(filestring, mmcif):
    parse_helix(filestring, mmcif)
    parse_sheet(filestring, mmcif)
    parse_ssbond(filestring, mmcif)
    parse_link(filestring, mmcif)
    parse_cispep(filestring, mmcif)
    parse_site(filestring, mmcif)


def parse_helix(filestring, mmcif):
    lines = re.findall(r"^HELIX.+", filestring, re.M)
    if not lines: return
    mmcif["struct_conf"] = [{
        "conf_type_id": "HELX_P", 
        "id": f"HELX_P{str(i)}", 
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
        "pdbx_PDB_helix_class": line[38:40].strip(), 
        "details": "?", 
        "pdbx_PDB_helix_length": line[71:76].strip(), 
    } for i, line in enumerate(lines, start=1)]


def parse_sheet(filestring, mmcif):
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
    mmcif["struct_sheet_order"] = []
    mmcif["pdbx_struct_sheet_hbond"] = []
    for line in lines:
        line = line.ljust(80)
        strand_id, sheet_id = line[7:10].strip(), line[11:14].strip()
        if strand_id != "1":
            mmcif["struct_sheet_order"].append({
                "sheet_id": sheet_id,
                "range_id_1": str(int(strand_id) - 1), "range_id_2": strand_id,
                "offset": "?",
                "sense": "parallel" if line[38:40].strip() == "1" else "anti-parallel"
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
    lines = re.findall(r"^SSBOND.+", filestring, re.M)
    if not lines: return
    mmcif["struct_conn"] = []
    for n, line in enumerate(lines, start=1):
        mmcif["struct_conn"].append({
            "id": f"disulf{n}", 
            "conn_type_id": "disulf", 
            "pdbx_leaving_atom_flag": "?", 
            "pdbx_PDB_id": "?", 
            "ptnr1_label_asym_id": line[15].strip(), 
            "ptnr1_label_comp_id": line[11:14].strip(), 
            "ptnr1_label_seq_id": line[17:21].strip(), 
            "ptnr1_label_atom_id": "?", 
            "pdbx_ptnr1_label_alt_id": "?", 
            "pdbx_ptnr1_PDB_ins_code": line[21].strip() or "?", 
            "pdbx_ptnr1_standard_comp_id": "?", 
            "ptnr1_symmetry": line[59:65].strip(), 
            "ptnr2_label_asym_id": line[29].strip(), 
            "ptnr2_label_comp_id": line[25:28].strip(), 
            "ptnr2_label_seq_id": line[31:35].strip(), 
            "ptnr2_label_atom_id": "?", 
            "pdbx_ptnr2_label_alt_id": "?", 
            "pdbx_ptnr2_PDB_ins_code": line[36].strip() or "?", 
            "ptnr1_auth_asym_id": line[15].strip(), 
            "ptnr1_auth_comp_id": line[11:14].strip(), 
            "ptnr1_auth_seq_id": line[17:21].strip(), 
            "ptnr2_auth_asym_id": line[29].strip(), 
            "ptnr2_auth_comp_id": line[25:28].strip(), 
            "ptnr2_auth_seq_id": line[31:35].strip(), 
            "ptnr2_symmetry": line[66:72].strip(), 
            "pdbx_ptnr3_label_atom_id": "?", 
            "pdbx_ptnr3_label_seq_id": "?", 
            "pdbx_ptnr3_label_comp_id": "?", 
            "pdbx_ptnr3_label_asym_id": "?", 
            "pdbx_ptnr3_label_alt_id": "?", 
            "pdbx_ptnr3_PDB_ins_code": "?", 
            "details": "?", 
            "pdbx_dist_value": line[73:78].strip(), 
            "pdbx_value_order": "?",
        })


def parse_link(filestring, mmcif):
    lines = re.findall(r"^LINK.+", filestring, re.M)
    if not lines: return
    if "struct_conn" not in mmcif: mmcif["struct_conn"] = []
    for n, line in enumerate(lines, start=1):
        mmcif["struct_conn"].append({
            "id": f"covale{n}", 
            "conn_type_id": "covale", 
            "pdbx_leaving_atom_flag": "?", 
            "pdbx_PDB_id": "?", 
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
            "pdbx_ptnr3_label_atom_id": "?", 
            "pdbx_ptnr3_label_seq_id": "?", 
            "pdbx_ptnr3_label_comp_id": "?", 
            "pdbx_ptnr3_label_asym_id": "?", 
            "pdbx_ptnr3_label_alt_id": "?", 
            "pdbx_ptnr3_PDB_ins_code": "?", 
            "details": "?", 
            "pdbx_dist_value": line[73:78].strip(), 
            "pdbx_value_order": "?",
        })


def parse_site(filestring, mmcif):
    lines = re.findall(r"^SITE.+", filestring, re.M)
    if not lines: return
    mmcif["struct_site_gen"] = []
    for line in lines:
        # TODO: UPDATE struct_site table too
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
                "label_atom_id": ".", 
                "label_alt_id": "?", 
                "symmetry": "1_555", 
                "details": "?",
            })


def parse_cispep(filestring, mmcif):
    lines = re.findall(r"^CISPEP.+", filestring, re.M)
    if not lines: return
    mmcif["struct_mon_prot_cis"] = [{
        "pdbx_id": line[7:10].strip(),
        "label_comp_id": line[11:14].strip(),
        "label_seq_id": line[17:21].strip(),
        "label_asym_id": line[15].strip(),
        "label_alt_id": ".",
        "pdbx_PDB_ins_code": line[21].strip() or "?",
        "auth_comp_id": line[11:14].strip(),
        "auth_seq_id": line[17:21].strip(),
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
    auth_asym_to_label_asym = {}
    label_asym_to_auth_asym = {}
    auth_res_to_label_res = {}
    for atom in mmcif["atom_site"]:
        label_asym_to_auth_asym[atom["label_asym_id"]] = atom["auth_asym_id"]
        auth_res = (atom["auth_asym_id"], atom["auth_seq_id"], atom["pdbx_PDB_ins_code"])
        label_res = (atom["label_asym_id"], atom["label_seq_id"])
        auth_res_to_label_res[auth_res] = label_res
    auth_asym_to_label_asym = {a: [
        k for k, v in label_asym_to_auth_asym.items() if v == a
    ] for a in label_asym_to_auth_asym.values()}
    for row in mmcif.get("pdbx_struct_assembly_gen", []):
        auths = row["asym_id"].split(",")
        labels = [id for auth in auths for id in auth_asym_to_label_asym[auth]]
        row["asym_id"] = ",".join(sorted(labels, key=lambda l: list(label_asym_to_auth_asym.keys()).index(l)))
    
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
                    label_asym_key = key.replace("seq", "asym")
                    label_seq_key = key.replace("asym", "seq")
                    auth_asym_key = label_asym_key.replace("label", "auth")
                    auth_seq_key = label_seq_key.replace("label", "auth")
                    auth_ins_key = [k for k in keys if "ins" in k and template[0] in k and template[1] in k][0]
                    if auth_ins_key not in keys or auth_seq_key not in keys or auth_asym_key not in keys: continue

                    for row in mmcif[name]:
                        sig = (row[auth_asym_key], row[auth_seq_key], row[auth_ins_key])
                        if sig in auth_res_to_label_res:
                            row[label_asym_key], row[label_seq_key] = auth_res_to_label_res[
                                (row[auth_asym_key], row[auth_seq_key], row[auth_ins_key])
                            ]


def pdb_date_to_mmcif_date(date):
    if not date.strip(): return
    day, month, year = date.split("-")
    month = str(list(calendar.month_abbr).index(month.title())).zfill(2)
    if len(year) == 2:
        if int(year) > 50: year = "19" + year
        if int(year) <= 50: year = "20" + year
    return f"{year}-{month}-{day}"


def process_names(lines):
    all_names = [name for line in lines for name in line.split(",") if name]
    processed_names = []
    for name in all_names:
        if "." in name and "," not in name:
            names = [n.title() for n in name.split(".")]
            processed_names.append(f"{names[-1]}, {'.'.join(names[:-1])}")
        else: processed_names.append(name.title())
    return processed_names


def residue_sig(atom):
    return (atom[21], atom[17:20].strip(), atom[22:26].strip(), atom[26].strip())


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
    lines += create_origix_lines(mmcif)
    lines += create_scalen_lines(mmcif)
    lines += create_mtrixn_lines(mmcif)
    lines += create_atom_lines(mmcif)
    lines.append("END")
    with open(path, "w") as f:
        f.write("\n".join(lines))


def create_header_line(mmcif):
    """Creates the HEADER line from a mmCIF dictionary.
    
    :param dict mmcif: the dictionary to update.
    :rtye: ``list``"""
    
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
    
    :param dict mmcif: the dictionary to update.
    :rtye: ``list``"""

    for row in mmcif.get("pdbx_database_PDB_obs_spr", []):
        if row["id"] == "OBSLTE":
            return ["OBSLTE     {:9} {:4}      {:4}".format(
                create_pdb_date(row["date"]),
                row["replace_pdb_id"], row["pdb_id"]
            )]
    return []


def create_title_lines(mmcif):
    """Creates the TITLE lines from a mmCIF dictionary.

    :param dict mmcif: the dictionary to update.
    :rtye: ``list``"""

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

    :param dict mmcif: the dictionary to update.
    :rtye: ``list``"""

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

    :param dict mmcif: the dictionary to update.
    :rtye: ``list``"""

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

    :param dict mmcif: the dictionary to update.
    :rtye: ``list``"""

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
            "EC": entity["pdbx_ec"],
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

    :param dict mmcif: the dictionary to update.
    :rtye: ``list``"""

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

    :param dict mmcif: the dictionary to update.
    :rtye: ``list``"""

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


def create_seqadv_lines(mmcif):
    lines = []
    line = "SEQADV {:4} {:>3} {:1} {:>4}{:1} {:4} {:9} {:>3} {:>5} {:20}"
    for seqadv in mmcif.get("struct_ref_seq_dif", []):
        lines.append(line.format(
            seqadv["pdbx_pdb_id_code"], seqadv["mon_id"],
            seqadv["pdbx_pdb_strand_id"], seqadv["pdbx_auth_seq_num"],
            seqadv["pdbx_pdb_ins_code"].replace("?", ""), seqadv["pdbx_seq_db_name"],
            seqadv["pdbx_seq_db_accession_code"], seqadv["db_mon_id"].replace("?", ""),
            seqadv["pdbx_seq_db_seq_num"].replace("?", ""), seqadv["details"].replace("?", ""),
        ))
    return lines


def create_seqres_lines(mmcif):
    lines = []
    line = "SEQRES {:3} {:1} {:4}  {:>3} {:>3} {:>3} {:>3} {:>3} {:>3} {:>3} {:>3} {:>3} {:>3} {:>3} {:>3} {:>3}"
    for entity in mmcif.get("entity_poly", []):
        chain_ids = entity["pdbx_strand_id"].split(",")
        residues = [r for r in mmcif.get("entity_poly_seq", []) if r["entity_id"] == entity["entity_id"]]
        line_count = math.ceil(len(residues) / 13)
        for chain_id in chain_ids:
            for n in range(line_count):
                line_res = residues[n * 13:(n + 1) * 13]
                names = [r["mon_id"] for r in line_res] + [""] * (13 - len(line_res))
                lines.append(line.format(n + 1, chain_id, len(residues), *names))
    return lines


def create_modres_lines(mmcif):
    lines = []
    line = "MODRES {:3} {:>3} {:1} {:>4}{:1} {:>3}  {:40}"
    for modres in mmcif.get("pdbx_struct_mod_residue", []):
        lines.append(line.format(
            mmcif["entry"][0]["id"], modres["auth_comp_id"], modres["auth_asym_id"],
            modres["auth_seq_id"], modres["PDB_ins_code"].replace("?", ""),
            modres["parent_comp_id"], modres["details"].replace("?", ""),
        ))
    return lines
    

def create_het_lines(mmcif):
    lines = []
    sig_counts = get_sig_counts(mmcif)
    for sig, count in sorted(sig_counts.items(), key=lambda s: s[0][3]):
        lines.append(f"HET    {sig[3]:3}  {sig[0]:1}{sig[1]:>4}{sig[2]:1}  {count:>5}")
    return lines


def create_hetnam_lines(mmcif):
    lines = []
    chem_comp = [c for c in mmcif["chem_comp"] if c["mon_nstd_flag"] != "y"]
    lookup = {e["id"]: e["type"] for e in mmcif["entity"]}
    lookup = {s["comp_id"]: lookup[s["entity_id"]] for s in mmcif.get("pdbx_entity_nonpoly", [])}
    for chem in chem_comp:
        if chem["name"] != "?" and lookup[chem["id"]] != "water":
            strings = split_lines(chem["name"], 55)
            for n, hetnam in enumerate(strings, start=1):
                if n == 1:
                    lines.append(f"HETNAM     {chem['id']:3} " + hetnam)
                else:
                    lines.append(f"HETNAM  {n:2} {chem['id']:3}  " + hetnam)
    return lines


def create_hetsyn_lines(mmcif):
    lines = []
    chem_comp = [c for c in mmcif["chem_comp"] if c["mon_nstd_flag"] != "y"]
    lookup = {e["id"]: e["type"] for e in mmcif["entity"]}
    lookup = {s["comp_id"]: lookup[s["entity_id"]] for s in mmcif.get("pdbx_entity_nonpoly", [])}
    for chem in chem_comp:
        if chem["pdbx_synonyms"] != "?" and lookup[chem["id"]] != "water":
            strings = split_lines(chem["pdbx_synonyms"].replace(", ", "; "), 55)
            for n, hetsyn in enumerate(strings, start=1):
                if n == 1:
                    lines.append(f"HETSYN     {chem['id']:3} " + hetsyn)
                else:
                    lines.append(f"HETSYN  {n:2} {chem['id']:3}  " + hetsyn)
    return lines


def create_formul_lines(mmcif):
    lines = []
    chem_comp = [c for c in mmcif["chem_comp"] if c["mon_nstd_flag"] != "y"]
    lookup = {e["id"]: e["type"] for e in mmcif["entity"]}
    lookup = {s["comp_id"]: lookup[s["entity_id"]] for s in mmcif.get("pdbx_entity_nonpoly", [])}
    name_to_entity_id = {s["comp_id"]: s["entity_id"] for s in mmcif.get("pdbx_entity_nonpoly", [])}
    entity_ids = sorted(s["entity_id"] for s in mmcif["struct_asym"])
    sig_counts = get_sig_counts(mmcif, include_water=True, representative=True)
    name_counts = {sig[3]: sum(
        v for k, v in sig_counts.items() if k[3] == sig[3]
    ) for sig in sig_counts}
    for chem in sorted(chem_comp, key=lambda c: list(name_counts.keys()).index(c["id"]) if c["id"] in name_counts else 0):
        if chem["formula"] != "?":
            entity_id = name_to_entity_id.get(chem["id"])
            if not entity_id: continue
            number = entity_ids.index(entity_id) + 1
            char = "*" if lookup[chem["id"]] == "water" else " "
            formula = f"{name_counts[chem['id']]}({chem['formula']})"
            strings = split_lines(formula, 50)
            for n, formul in enumerate(strings, start=1):
                if n == 1:
                    lines.append(f"FORMUL  {number:2}  {chem['id']:3}   {char}{formul}")
                else:
                    lines.append(f"FORMUL  {number:2}  {chem['id']:3} {n:2}{char}{formul}")
    return lines


def get_sig_counts(mmcif, include_water=False, representative=False):
    sigs = []
    names = [c["id"] for c in mmcif["chem_comp"]  if c["mon_nstd_flag"] != "y"]
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
    line = "HELIX  {:>3} {:>3} {:3} {:1} {:>4}{:1} {:3} {:1} {:>4}{:1}{:>2}{:30} {:>5}"
    return [line.format(
        helix["pdbx_PDB_helix_id"],
        helix["pdbx_PDB_helix_id"],
        helix["beg_auth_comp_id"],
        helix["beg_auth_asym_id"],
        helix["beg_auth_seq_id"],
        helix["pdbx_beg_PDB_ins_code"].replace("?", ""),
        helix["end_auth_comp_id"],
        helix["end_auth_asym_id"],
        helix["end_auth_seq_id"],
        helix["pdbx_end_PDB_ins_code"].replace("?", ""),
        helix["pdbx_PDB_helix_class"],
        "", helix["pdbx_PDB_helix_length"]
    ) for helix in mmcif.get("struct_conf", [])]


def create_sheet_lines(mmcif):
    lines = []
    line = "SHEET  {:>3} {:>3}{:>2} {:>3} {:1}{:>4}{:1} {:>3} {:1}{:>4}{:1}{:>2} {:<4}{:>3} {:1}{:>4}{:1} {:<4}{:>3} {:1}{:>4}{:1}"
    for strand in mmcif.get("struct_sheet_range", []):
        order = [r for r in mmcif["struct_sheet_order"] if r["sheet_id"] == strand["sheet_id"] and r["range_id_2"] == strand["id"]]
        order = order[0] if order else None
        hbond = [r for r in mmcif["pdbx_struct_sheet_hbond"] if r["sheet_id"] == strand["sheet_id"] and r["range_id_2"] == strand["id"]]
        hbond = hbond[0] if hbond else None
        lines.append(line.format(
            strand["id"],
            strand["sheet_id"],
            str(len([r for r in mmcif["struct_sheet_range"] if r["sheet_id"] == strand["sheet_id"]])),
            strand["beg_auth_comp_id"],
            strand["beg_auth_asym_id"],
            strand["beg_auth_seq_id"],
            strand["pdbx_beg_PDB_ins_code"].replace("?", ""),
            strand["end_auth_comp_id"],
            strand["end_auth_asym_id"],
            strand["end_auth_seq_id"],
            strand["pdbx_end_PDB_ins_code"].replace("?", ""),
            0 if not order else 1 if order["sense"] == "parallel" else -1,
            f"{'' if len(hbond['range_2_auth_atom_id']) == 4 else ' '}{hbond['range_2_auth_atom_id']}" if hbond else "",
            hbond["range_2_auth_comp_id"] if hbond else "",
            hbond["range_2_auth_asym_id"] if hbond else "",
            hbond["range_2_auth_seq_id"] if hbond else "",
            hbond["range_2_PDB_ins_code"].replace("?", "") if hbond else "",
            f"{'' if len(hbond['range_1_auth_atom_id']) == 4 else ' '}{hbond['range_1_auth_atom_id']}" if hbond else "",
            hbond["range_1_auth_comp_id"] if hbond else "",
            hbond["range_1_auth_asym_id"] if hbond else "",
            hbond["range_1_auth_seq_id"] if hbond else "",
            hbond["range_1_PDB_ins_code"].replace("?", "") if hbond else "",
        ))
    return lines


def create_ssbond_lines(mmcif):
    line = "SSBOND {:>3} {:>3} {:1} {:>4}{:1}   {:>3} {:1} {:>4}{:1}                       {:>6} {:>6} {:>5}"
    lines = []
    for ssbond in mmcif.get("struct_conn", []):
        if ssbond["conn_type_id"] != "disulf": continue
        lines.append(line.format(
            str(len(lines) + 1),
            ssbond["ptnr1_auth_comp_id"],
            ssbond["ptnr1_auth_asym_id"],
            ssbond["ptnr1_auth_seq_id"],
            ssbond["pdbx_ptnr1_PDB_ins_code"].replace("?", ""),
            ssbond["ptnr2_auth_comp_id"],
            ssbond["ptnr2_auth_asym_id"],
            ssbond["ptnr2_auth_seq_id"],
            ssbond["pdbx_ptnr2_PDB_ins_code"].replace("?", ""),
            ssbond["ptnr1_symmetry"],
            ssbond["ptnr2_symmetry"],
            ssbond["pdbx_dist_value"],
        ))
    return lines


def create_link_lines(mmcif):
    line = "LINK        {:>4}{:1}{:>3} {:1}{:>4}{:1}               {:>4}{:1}{:>3} {:1}{:>4}{:1}  {:>6} {:>6} {:>5}"
    lines = []
    for link in mmcif.get("struct_conn", []):
        if link["conn_type_id"] != "covale": continue
        lines.append(line.format(
            link["ptnr1_label_atom_id"],
            link["pdbx_ptnr1_label_alt_id"].replace("?", ""),
            link["ptnr1_auth_comp_id"],
            link["ptnr1_auth_asym_id"],
            link["ptnr1_auth_seq_id"],
            link["pdbx_ptnr1_PDB_ins_code"].replace("?", ""),
            link["ptnr2_label_atom_id"],
            link["pdbx_ptnr2_label_alt_id"].replace("?", ""),
            link["ptnr2_auth_comp_id"],
            link["ptnr2_auth_asym_id"],
            link["ptnr2_auth_seq_id"],
            link["pdbx_ptnr2_PDB_ins_code"].replace("?", ""),
            link["ptnr1_symmetry"],
            link["ptnr2_symmetry"],
            link["pdbx_dist_value"],
        ))
    return lines
    

def create_cispep_lines(mmcif):
    line = "CISPEP {:>3} {:>3} {:1} {:>4}{:1}   {:>3} {:1} {:>4}{:1}       {:>3}       {:>6}"
    return [line.format(
        cispep["pdbx_id"],
        cispep["auth_comp_id"],
        cispep["auth_asym_id"],
        cispep["auth_seq_id"],
        cispep["pdbx_PDB_ins_code"].replace("?", ""),
        cispep["pdbx_auth_comp_id_2"],
        cispep["pdbx_auth_asym_id_2"],
        cispep["pdbx_auth_seq_id_2"],
        cispep["pdbx_PDB_ins_code_2"].replace("?", ""),
        str(int(cispep["pdbx_PDB_model_num"]) - 1),
        cispep["pdbx_omega_angle"],
    ) for cispep in mmcif.get("struct_mon_prot_cis", [])]


def create_site_lines(mmcif):
    line = "SITE   {:>3} {:>3} {:>2} {:3} {:1}{:>4}{:1} {:3} {:1}{:>4}{:1} {:3} {:1}{:>4}{:1} {:3} {:1}{:>4}{:1}"
    lines = []
    site_ids = []
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
                    row["auth_comp_id"],
                    row["auth_asym_id"],
                    row["auth_seq_id"],
                    row["pdbx_auth_ins_code"].replace("?", ""),
                ]
            values += [""] * (19 - len(values))
            lines.append(line.format(*values))
    return lines


def create_cryst1_line(mmcif):
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
    return [line]


def create_origix_lines(mmcif):
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
    lines = []
    line = "MTRIX{:1} {:3}{:>10}{:>10}{:>10}     {:>10}    {:1}"
    for i, row in enumerate(mmcif.get("struct_ncs_oper", []), start=1):
        for n in range(1, 4):
            lines.append(line.format(
                n, i, row[f"matrix[{n}][1]"],
                row[f"matrix[{n}][2]"],
                row[f"matrix[{n}][3]"],
                row[f"vector[{n}]"],
                "1" if row["code"] == "given" else " "
            ).replace("?", " "))
    return lines


def create_atom_lines(mmcif):
    if not mmcif.get("atom_site", []): return []
    lines = []
    model_num = 0
    atom_id, asym_id, entity_id = 1, mmcif["atom_site"][0]["label_asym_id"], ""
    aniso_lookup = {a["id"]: a for a in mmcif.get("atom_site_anisotrop", [])}
    model_nums = set(a["pdbx_PDB_model_num"] for a in  mmcif["atom_site"])
    lookup = {e["id"]: e["type"] for e in mmcif["entity"]}
    for atom in mmcif["atom_site"]:
        if int(atom["pdbx_PDB_model_num"]) > model_num and len(model_nums) > 1:
            if model_num != 0: lines.append("ENDMDL")
            model_num += 1
            lines.append(f"MODEL     {model_num:>4}")
            atom_id = 1
        if atom["label_asym_id"] != asym_id and lookup[entity_id] == "polymer":
            lines.append(f"TER   {atom_id:>5}      {lines[-1][17:26]}")
            atom_id += 1
        asym_id = atom["label_asym_id"]
        entity_id = atom["label_entity_id"]
        lines.append(create_atom_line(atom, atom_id))
        line = create_aniso_line(atom, aniso_lookup.get(atom["id"]), atom_id)
        if line: lines.append(line)
        atom_id += 1
    return lines


def create_atom_line(atom, atom_id):
    line = "{:6}{:>5} {:<4} {:3} {:1}{:>4}{:1}   "
    line += "{:>8}{:>8}{:>8}  1.00{:>6}          {:>2}{:2}"
    line = line.format(
        atom["group_PDB"],
        atom_id,
        f"{'' if len(atom['label_atom_id']) == 4 else ' '}{atom['label_atom_id']}",
        atom["auth_comp_id"],
        atom["auth_asym_id"], atom["auth_seq_id"], atom["pdbx_PDB_ins_code"],
        "{:.3f}".format(float(atom["Cartn_x"])),
        "{:.3f}".format(float(atom["Cartn_y"])),
        "{:.3f}".format(float(atom["Cartn_z"])),
        "{:.2f}".format(float(atom["B_iso_or_equiv"])),
        atom["type_symbol"] or "", atom["pdbx_formal_charge"][::-1],
    )
    return line.replace("?", " ")


def create_aniso_line(atom, aniso, atom_id):
    if not aniso: return
    line = "ANISOU{:5} {:4} {:3} {:1}{:>4}{:1} "
    line += "{:>7}{:>7}{:>7}{:>7}{:>7}{:>7}      {:>2}{:2}"
    line = line.format(
        atom_id, f"{'' if len(atom['label_atom_id']) == 4 else ' '}{atom['label_atom_id']}",
        atom["auth_comp_id"],
        atom["auth_asym_id"], atom["auth_seq_id"], atom["pdbx_PDB_ins_code"],
        int(float(aniso["U[1][1]"]) * 10000),
        int(float(aniso["U[2][2]"]) * 10000),
        int(float(aniso["U[3][3]"]) * 10000),
        int(float(aniso["U[1][2]"]) * 10000),
        int(float(aniso["U[1][3]"]) * 10000),
        int(float(aniso["U[2][3]"]) * 10000),
        atom["type_symbol"] or "",
        atom["pdbx_formal_charge"][::-1],
    )
    return line.replace("?", " ")


def create_pdb_date(date):
    if not date.strip(): return
    year, month, day = date.split("-")
    month = list(calendar.month_abbr)[int(month)].upper()
    return f"{day}-{month}-{year[2:]}"


def split_lines(string, length):
    if len(string) <= length: return [string]
    strings = []
    while string:
        first = string[:length]
        last_space = first[::-1].find(" ")
        first = string[:length - last_space]
        string = string[length - last_space:].lstrip()
        strings.append(first)
    return strings