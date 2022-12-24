"""Functions for handling the mmCIF file type."""

from collections import deque

def mmcif_string_to_mmcif_dict(filestring):
    """Takes the contents of a .cif file and returns a fully parsed dictionary
    representation of it.
    
    :param str filestring: the contents of the .cif file.
    :rtype: ``dict``"""

    mmcif_dict = {}
    sections = get_sections_from_filestring(filestring.strip())
    for section in sections:
        if section[0] == "loop_":
            add_loop_section(section, mmcif_dict)
        else:
            add_non_loop_section(section, mmcif_dict)
    return mmcif_dict


def get_sections_from_filestring(filestring):
    """Returns a list of 'sections' in a .cif file, where each section is a list
    of lines for a particular category.
    
    :param str filestring: the contents of the .cif file.
    :rtype: ``list``"""

    sections = filestring.split("\n#")[1:-1]
    sections = [[
        line.rstrip() for line in s.strip().split("\n")
    ] for s in sections]
    return sections


def add_non_loop_section(section, mmcif_dict):
    """Takes a non-loop section and updates an mmCIF dictionary with its
    contents. The category will be a single item list.
    
    :param list section: the lines of the section.
    :param dict mmcif_dict: the dictionary to update."""

    row = {}
    section = deque(section)
    isquote = lambda s: any(s[0] == c and s[-1] == c for c in "'\"")
    while section:
        line = section[0]
        category_name, key = line.split()[0].split(".")
        if len(section) > 2 and section[1].startswith(";"):
            section.popleft()
            value = get_semicolon_value(section)
        elif len(section) > 1 and " " not in line:
            value = section[1].lstrip()
            section.popleft()
        else:
            value = line[len(category_name) + len(key) + 2:].strip()
        if isquote(value): value = value[1:-1]
        row[key] = value
        section.popleft()
    mmcif_dict[category_name[1:]] = [row]


def add_loop_section(section, mmcif_dict):
    """Takes a loop section and updates an mmCIF dictionary with its contents.
    The category will be a multi item list.
    
    :param list section: the lines of the section.
    :param dict mmcif_dict: the dictionary to update."""

    keys, values = [], []
    section = deque(section[1:])
    while section:
        line = section[0]
        if line[0] == "_":
            category_name, key = line.strip().split(".")
            keys.append(key)
        elif line[0] == ";":
            values.append(get_semicolon_value(section))
        else:
            values += get_line_values(line)
        section.popleft()
    row_count = len(values) // len(keys)
    mmcif_dict[category_name[1:]] = [{
        key: values[n * len(keys) + i] for i, key in enumerate(keys)
    } for n in range(row_count)]


def get_line_values(line):
    """Takes an mmCIF lines representing a category table row and breaks into
    individual cells. Quote strings are handled.
    
    :param str line: the line to break up.
    :rtype: ``str``"""

    if "'" not in line and '"' not in line: return line.split()
    values = []
    while line.strip():
        value = line.split()[0]
        if value[0] in "'\"":
            end = line[1:].find(value[0])
            values.append(line[1:end + 1])
            line = line[len(values[-1]) + 2:].lstrip()
        else:
            values.append(value)
            line = line[len(value):].lstrip()
    return values


def get_semicolon_value(lines):
    """Takes a list of lines (as deque) where the first line is a semicolon
    value, and returns the full semi-colon value over however many lines it is
    over, and pops lines off the deque so that the first line is the single
    semi-colon at the end.
    
    :param deque lines: a list of lines where the first line starts with `;`.
    :rtype: ``str``"""
    
    values = []
    values.append(lines[0])
    lines.popleft()
    while lines[0] != ";":
        values.append(lines[0])
        lines.popleft()
    join = "\n" if any(" " in val for val in values) else ""
    return join.join(values)[1:].strip()


def save_mmcif_dict(mmcif_dict, path):
    """Saves an mmCIF dictionary to a .cif file.

    :param dict mmcif_dict: the dictionary to save.
    :param Path path: the location to save to."""

    code = mmcif_dict.get("entry", [{"id": ""}])[0]["id"]
    lines = [f"data_{code}" if code else "data_XXXX"]
    for category_name, rows in mmcif_dict.items():
        if len(rows): lines.append("#")
        if len(rows) == 1:
            lines += get_non_loop_lines(category_name, rows)
        else:
            lines += get_loop_lines(category_name, rows)
    lines.append("#")
    with open(path, "w") as f:
        f.write("\n".join(lines))


def get_non_loop_lines(name, rows):
    """Gets the mmCIF lines for a single row section. Values will be enclosed in
    quotes where required, and put onto separate semi-colon lines where
    required.
    
    :param str name: the name of the category.
    :param list rows: the mmCIF rows to process.
    :rtype: ``list``"""

    lines = []
    keys = list(rows[0].keys())
    max_length = max(len(k) for k in keys)
    for key in keys:
        value = format_value(rows[0][key])
        if isinstance(value, list):
            lines.append(f"_{name}.{key.ljust(max_length)}")
            lines += value
        else:
            lines.append(f"_{name}.{key.ljust(max_length)} {value}")
    return lines


def get_loop_lines(name, rows):
    """Gets the mmCIF lines for a multi-row section. Values will be enclosed in
    quotes where required, and put onto separate semi-colon lines where
    required.
    
    :param str name: the name of the category.
    :param list rows: the mmCIF rows to process.
    :rtype: ``list``"""

    if len(rows) == 0: return []
    keys = rows[0].keys()
    lines = ["loop_", *[f"_{name}.{key}" for key in keys]]
    grid = [[format_value(row[key]) for key in keys] for row in rows]
    max_lengths = [max([
        len(row[col_num]) for row in grid if not isinstance(row[col_num], list)
    ] or [0]) for col_num in range(len(keys))]
    for row in grid:
        line = []
        for i, cell in enumerate(row):
            if isinstance(cell, list):
                if line: lines.append(" ".join(line))
                lines += cell
                line = []
            else:
                line.append(cell.ljust(max_lengths[i]))
        if line: lines.append(" ".join(line))
    return lines


def format_value(s):
    """Formats a mmCIF value for encoding in a filestring. If the value contains
    spaces or double quotes, it will be enclosed in single quotes. If it
    contains double quotes it will be enclosed in single quotes. And if it
    contains line breaks, or two different kinds of quotes, a list of lines will
    be returned which should be inserted to represent that value.
    
    :param str s: the string to encode.
    :rtype: ``str``"""

    if "\n" in s or ("'" in s and '"' in s):
        lines = s.split("\n")
        lines[0] = f";{lines[0]}"
        lines.append(";")
        return lines
    elif "'" in s:
        return f'"{s}"'
    elif '"' in s or " " in s:
        return f"'{s}'"
    else: return s


CATEGORIES = {
    "reflns": [
        "entry_id",
        "observed_criterion_sigma_I",
        "observed_criterion_sigma_F",
        "d_resolution_low",
        "d_resolution_high",
        "number_obs",
        "number_all",
        "percent_possible_obs",
        "pdbx_Rmerge_I_obs",
        "pdbx_Rsym_value",
        "pdbx_netI_over_sigmaI",
        "B_iso_Wilson_estimate",
        "pdbx_redundancy",
        "R_free_details",
        "limit_h_max",
        "limit_h_min",
        "limit_k_max",
        "limit_k_min",
        "limit_l_max",
        "limit_l_min",
        "observed_criterion_F_max",
        "observed_criterion_F_min",
        "pdbx_chi_squared",
        "pdbx_scaling_rejects",
        "pdbx_diffrn_id",
        "pdbx_ordinal",
        "d_res_high",
        "d_res_low",
        "percent_possible_all",
        "Rmerge_I_obs",
        "pdbx_Rsym_value",
        "meanI_over_sigI_obs",
        "pdbx_redundancy",
        "percent_possible_obs",
        "number_unique_all",
        "number_measured_all",
        "number_measured_obs",
        "number_unique_obs",
        "pdbx_chi_squared",
        "pdbx_diffrn_id",
        "pdbx_ordinal"
    ],
    "reflns_shell": [
        "d_res_high",
        "d_res_low",
        "percent_possible_all",
        "Rmerge_I_obs",
        "pdbx_Rsym_value",
        "meanI_over_sigI_obs",
        "pdbx_redundancy",
        "percent_possible_obs",
        "number_unique_all",
        "number_measured_all",
        "number_measured_obs",
        "number_unique_obs",
        "pdbx_chi_squared",
        "pdbx_diffrn_id",
        "pdbx_ordinal"
    ],
    "refine": [
        "entry_id",
        "ls_number_reflns_obs",
        "ls_number_reflns_all",
        "pdbx_ls_sigma_I",
        "pdbx_ls_sigma_F",
        "pdbx_data_cutoff_high_absF",
        "pdbx_data_cutoff_low_absF",
        "ls_d_res_low",
        "ls_d_res_high",
        "ls_percent_reflns_obs",
        "ls_R_factor_obs",
        "ls_R_factor_all",
        "ls_R_factor_R_work",
        "ls_R_factor_R_free",
        "ls_R_factor_R_free_error",
        "ls_R_factor_R_free_error_details",
        "ls_percent_reflns_R_free",
        "ls_number_reflns_R_free",
        "ls_number_parameters",
        "ls_number_restraints",
        "occupancy_min",
        "occupancy_max",
        "correlation_coeff_Fo_to_Fc",
        "correlation_coeff_Fo_to_Fc_free",
        "B_iso_mean",
        "aniso_B[1][1]",
        "aniso_B[2][2]",
        "aniso_B[3][3]",
        "aniso_B[1][2]",
        "aniso_B[1][3]",
        "aniso_B[2][3]",
        "solvent_model_details",
        "solvent_model_param_ksol",
        "solvent_model_param_bsol",
        "pdbx_solvent_vdw_probe_radii",
        "pdbx_solvent_ion_probe_radii",
        "pdbx_solvent_shrinkage_radii",
        "pdbx_ls_cross_valid_method",
        "details",
        "pdbx_starting_model",
        "pdbx_method_to_determine_struct",
        "pdbx_isotropic_thermal_model",
        "pdbx_stereochemistry_target_values",
        "pdbx_stereochem_target_val_spec_case",
        "pdbx_R_Free_selection_details",
        "pdbx_overall_ESU_R_Free",
        "overall_SU_B",
        "ls_redundancy_reflns_obs",
        "B_iso_min",
        "B_iso_max",
        "overall_SU_R_Cruickshank_DPI",
        "overall_SU_R_free",
        "overall_SU_ML",
        "pdbx_overall_ESU_R",
        "pdbx_data_cutoff_high_rms_absF",
        "pdbx_refine_id",
        "pdbx_overall_phase_error",
        "ls_wR_factor_R_free",
        "ls_wR_factor_R_work",
        "overall_FOM_free_R_set",
        "overall_FOM_work_R_set",
        "pdbx_diffrn_id",
        "pdbx_TLS_residual_ADP_flag",
        "pdbx_overall_SU_R_free_Cruickshank_DPI",
        "pdbx_overall_SU_R_Blow_DPI",
        "pdbx_overall_SU_R_free_Blow_DPI",
        "entry_id",
        "Luzzati_coordinate_error_obs",
        "Luzzati_sigma_a_obs",
        "Luzzati_d_res_low_obs",
        "Luzzati_coordinate_error_free",
        "Luzzati_sigma_a_free",
        "Luzzati_d_res_low_free",
        "number_disordered_residues",
        "occupancy_sum_hydrogen",
        "occupancy_sum_non_hydrogen",
        "pdbx_Luzzati_d_res_high_obs",
        "pdbx_refine_id",
        "pdbx_refine_id",
        "cycle_id",
        "pdbx_number_atoms_protein",
        "pdbx_number_atoms_nucleic_acid",
        "pdbx_number_atoms_ligand",
        "number_atoms_solvent",
        "number_atoms_total",
        "d_res_high",
        "d_res_low",
        "type",
        "dev_ideal",
        "dev_ideal_target",
        "weight",
        "number",
        "pdbx_refine_id",
        "pdbx_restraint_function",
        "pdbx_total_number_of_bins_used",
        "d_res_high",
        "d_res_low",
        "number_reflns_R_work",
        "R_factor_R_work",
        "percent_reflns_obs",
        "R_factor_R_free",
        "R_factor_R_free_error",
        "percent_reflns_R_free",
        "number_reflns_R_free",
        "number_reflns_obs",
        "redundancy_reflns_obs",
        "number_reflns_all",
        "pdbx_refine_id",
        "R_factor_all"
    ],
    "refine_analyze": [
        "entry_id",
        "Luzzati_coordinate_error_obs",
        "Luzzati_sigma_a_obs",
        "Luzzati_d_res_low_obs",
        "Luzzati_coordinate_error_free",
        "Luzzati_sigma_a_free",
        "Luzzati_d_res_low_free",
        "number_disordered_residues",
        "occupancy_sum_hydrogen",
        "occupancy_sum_non_hydrogen",
        "pdbx_Luzzati_d_res_high_obs",
        "pdbx_refine_id"
    ],
    "refine_hist": [
        "pdbx_refine_id",
        "cycle_id",
        "pdbx_number_atoms_protein",
        "pdbx_number_atoms_nucleic_acid",
        "pdbx_number_atoms_ligand",
        "number_atoms_solvent",
        "number_atoms_total",
        "d_res_high",
        "d_res_low"
    ],
    "refine_ls_restr": [
        "type",
        "dev_ideal",
        "dev_ideal_target",
        "weight",
        "number",
        "pdbx_refine_id",
        "pdbx_restraint_function"
    ],
    "refine_ls_shell": [
        "pdbx_total_number_of_bins_used",
        "d_res_high",
        "d_res_low",
        "number_reflns_R_work",
        "R_factor_R_work",
        "percent_reflns_obs",
        "R_factor_R_free",
        "R_factor_R_free_error",
        "percent_reflns_R_free",
        "number_reflns_R_free",
        "number_reflns_obs",
        "redundancy_reflns_obs",
        "number_reflns_all",
        "pdbx_refine_id",
        "R_factor_all"
    ]
}