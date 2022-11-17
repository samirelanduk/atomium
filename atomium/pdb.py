import re
import calendar

def pdb_string_to_mmcif_dict(filestring):
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
    parse_scalen(filestring, mmcif)
    parse_mtrixn(filestring, mmcif)
    return mmcif


def parse_header(filestring, mmcif):
    header = re.search(f"HEADER.+", filestring, re.M).group(0)
    code = (header[62:66].strip() or "1XXX") if header else "1XXX"
    if set(code) == {"-"}: code = "1XXX"
    mmcif["entry"] = [{"id": code}]
    date = header[50:59].strip()
    mmcif["pdbx_database_status"] = [{
        "status_code": "REL", "entry_id": "1LOL",
        "recvd_initial_deposition_date": pdb_date_to_mmcif_date(date)
    }]
    keyword = (header[10:50].strip() or "?") if header else "?"
    if set(keyword) == {"-"} or keyword in ["NULL", "NONE"]: keyword = "?"
    mmcif["struct_keywords"] = [{
        "entry_id": code, "pdbx_keywords": keyword, "text": "?"
    }]


def parse_obslte(filestring, mmcif):
    obslte_lines = re.findall(r"^OBSLTE.+", filestring, re.M)
    if obslte_lines:
        codes = [code for l in obslte_lines for code in l[31:].strip().split()]
        mmcif["pdbx_database_PDB_obs_spr"] = [{
            "id": "OBSLTE", "details": "?",
            "date": pdb_date_to_mmcif_date(obslte_lines[0][11:20].strip()),
            "pdb_id": " ".join(codes),
            "replace_pdb_id": obslte_lines[0][21:25].strip()
        }]


def parse_title(filestring, mmcif):    
    mmcif["struct"] = [{"entry_id": mmcif["entry"][0]["id"], **{
        k: "?" for k in [
            "title", "pdbx_descriptor", "pdbx_model_details", "pdbx_CASP_flag",
            "pdbx_model_type_details"
        ]
    }}]
    title_lines = re.findall(r"^TITLE.+", filestring, re.M)
    title = " ".join([l[10:80] for l in title_lines]).strip()
    mmcif["struct"][0]["title"] = " ".join(title.split())


def parse_split(filestring, mmcif):
    split_lines = re.findall(r"^SPLIT.+", filestring, re.M)
    if split_lines:
        codes = [code for l in split_lines for code in l[11:].strip().split()]
        mmcif["pdbx_database_related"] = [{
            "db_name": "PDB", "db_id": code, "content_type": "split",
            "details": f"Split {n}"
        } for n, code in enumerate(codes, start=1)]


def parse_caveat(filestring, mmcif):    
    caveat_lines = re.findall(r"^CAVEAT.+", filestring, re.M)
    caveat = " ".join([l[19:80] for l in caveat_lines]).strip()
    mmcif["database_PDB_caveat"] = [{"text": " ".join(caveat.split())}]


def parse_keywds(filestring, mmcif):
    lines = re.findall(r"^KEYWDS.+", filestring, re.M)
    if not lines: return
    keywords = " ".join([l[10:].strip() for l in lines]).strip()
    mmcif["struct_keywords"][0]["text"] = keywords


def parse_expdta(filestring, mmcif):
    mmcif["exptl"] = [{"entry_id": mmcif["entry"][0]["id"], "method": "?"}]
    exptl = re.findall(f"^EXPDTA.+", filestring, re.M)
    if not exptl: return
    mmcif["exptl"][0]["method"] = " ".join([l[10:79].strip() for l in exptl])


def parse_mdltyp(filestring, mmcif):
    mdltyp = re.findall(f"^MDLTYP.+", filestring, re.M)
    if not mdltyp: return
    text = " ".join([l[10:80].strip() for l in mdltyp]).strip()
    features = text.split(";")
    sections = []
    for feature in features:
        match = re.match(r"(.+?), CHAIN (.+)", feature)
        if match:
            if "pdbx_coordinate_model" not in mmcif:
                mmcif["pdbx_coordinate_model"] = []
            for chain in match[2].split(","):
                mmcif["pdbx_coordinate_model"].append({
                    "asym_id": chain.strip(), "type": match[1].strip()
                })
        else:
            sections.append(feature.strip())
    if sections:
        mmcif["struct"][0]["pdbx_model_type_details"] = " ; ".join(sections)


def parse_author(filestring, mmcif):
    lines = re.findall(r"^AUTHOR.+", filestring, re.M)
    if not lines: return
    author_lines = [line[10:].strip() for line in lines]
    authors = process_names(author_lines)
    mmcif["audit_author"] = [{
        "name": author, "pdbx_ordinal": str(num)
    }  for num, author in enumerate(authors, start=1)]


def parse_revdat(filestring, mmcif):
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
    sprsde_lines = re.findall(r"^SPRSDE.+", filestring, re.M)
    if sprsde_lines:
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
    title_lines = [l for l in journal_lines if l[12:16] == "TITL"]
    if not title_lines: return
    title = " ".join([l[19:80] for l in title_lines])
    mmcif["citation"][0]["title"] = " ".join(title.split())


def parse_journal_references(journal_lines, mmcif):
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
        mmcif["citation"][0]["book_publisher"] = lines[0][19:80].strip() or "?"


def parse_journal_ids(journal_lines, mmcif):
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
    parse_remark_2(filestring, mmcif)
    parse_remark_3(filestring, mmcif)
    parse_remark_200(filestring, mmcif)
    parse_remark_465(filestring, mmcif)
    parse_remark_470(filestring, mmcif)
    parse_remark_480(filestring, mmcif)
    parse_remark_800(filestring, mmcif)


def parse_remark_2(filestring, mmcif):
    reflns = [
        "B_iso_Wilson_estimate", "entry_id", "data_reduction_details",
        "data_reduction_method", "d_resolution_high", "d_resolution_low",
        "details", "limit_h_max", "limit_h_min", "limit_k_max", "limit_k_min",
        "limit_l_max", "limit_l_min", "number_all", "number_obs",
        "observed_criterion", "observed_criterion_F_max",
        "observed_criterion_F_min", "observed_criterion_I_max",
        "observed_criterion_I_min", "observed_criterion_sigma_F",
        "observed_criterion_sigma_I", "percent_possible_obs", "R_free_details",
        "Rmerge_F_all", "Rmerge_F_obs", "Friedel_coverage", "number_gt",
        "threshold_expression", "pdbx_redundancy", "pdbx_Rmerge_I_obs",
        "pdbx_Rmerge_I_all", "pdbx_Rsym_value", "pdbx_netI_over_av_sigmaI",
        "pdbx_netI_over_sigmaI", "pdbx_res_netI_over_av_sigmaI_2",
        "pdbx_res_netI_over_sigmaI_2", "pdbx_chi_squared",
        "pdbx_scaling_rejects", "pdbx_d_res_high_opt", "pdbx_d_res_low_opt",
        "pdbx_d_res_opt_method", "phase_calculation_details", "pdbx_Rrim_I_all",
        "pdbx_Rpim_I_all", "pdbx_d_opt", "pdbx_number_measured_all",
        "pdbx_diffrn_id", "pdbx_ordinal", "pdbx_CC_half", "pdbx_R_split"
    ]
    mmcif["reflns"] = [{
        **{key: "?" for key in reflns}, "entry_id": mmcif["entry"][0]["id"]
    }]
    records = re.findall(r"^REMARK   2 .+", filestring, re.M)
    for r in records:
        if "RESOLUTION" in r:
            mmcif["reflns"][0]["d_resolution_high"] = r[10:].strip().split()[1]
            break


def parse_remark_3(filestring, mmcif):
    records = re.findall(r"^REMARK   3 .+", filestring, re.M)
    if not records: return
    string = "\n".join(records)
    for regex, key in [
        [r"RESOLUTION RANGE LOW.+?\:(.+)", "d_resolution_low"],
        [r"RESOLUTION RANGE HIGH.+?\:(.+)", "d_resolution_high"],
        [r"R VALUE[ ]+?\(WORKING SET\) \: (.+)", "ls_R_factor_R_work"],
        [r"FREE R VALUE[ ]+?\:(.+)", "ls_R_factor_R_free"],
        [
            r"FREE R VALUE TEST SET SIZE[ ]+?\(%\)[ ]+?\:(.+)",
            "ls_percent_reflns_R_free"
        ],
        [r"FREE R VALUE TEST SET COUNT.+?\:(.+)", "ls_number_reflns_R_free"],
    ]:
        value = re.search(regex, string)
        mmcif["reflns"][0][key] = value[1].strip() if value else "?"


def parse_remark_200(filestring, mmcif):
    records = re.findall(r"^REMARK 200 .+", filestring, re.M)
    if not records: return
    string = "\n".join(records)
    temp = re.search(r"TEMPERATURE.+?\: (\d+)", string)
    detector = re.search(r"DETECTOR TYPE[ ]+\:[ ]+(.+)", string)
    type_ = re.search(r"DETECTOR MANUFACTURER.+?\:[ ]+(.+)", string)
    collection = re.search(r"DATE OF DATA COLLECTION.+?\:[ ]+(.+)", string)
    details = re.search(r"OPTICS.+?\:(.+?)  ", string)
    mmcif["diffrn"] = [{
        "id": "1", "ambient_temp": temp[1].strip() if temp else "?",
        "ambient_temp_details": "?", "crystal_id": "1"
    }]
    mmcif["diffrn_detector"] = [{
        "diffrn_id": "1",
        "detector": detector[1].strip() if detector else "?",
        "type": type_[1].strip() if type_ else "?",
        "pdbx_collection_date": pdb_date_to_mmcif_date(collection[1].strip())\
            if collection and collection[1].strip() != "NULL" else "?",
        "details": details[1].strip() if details else "?"
    }]


def parse_remark_465(filestring, mmcif):
    records = re.findall(r"^REMARK 465 .+", filestring, re.M)
    if not records: return
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


def parse_remark_470(filestring, mmcif):
    records = re.findall(r"^REMARK 470 .+", filestring, re.M)
    if not records[6:]: return
    mmcif["pdbx_unobs_or_zero_occ_atoms"] = []
    for rec in records[6:]:
        for atom in rec[27:].strip().split():
            mmcif["pdbx_unobs_or_zero_occ_atoms"].append({
                "id": str(len(mmcif["pdbx_unobs_or_zero_occ_atoms"]) + 1),
                "PDB_model_num": "1", "polymer_flag": "Y", "occupancy_flag": "0",
                "auth_asym_id": rec[19].strip() or "?",
                "auth_comp_id": rec[15:18].strip() or "?",
                "auth_seq_id": rec[21:26].strip() or "?",
                "PDB_ins_code": rec[26].strip() or "?",
                "auth_atom_id": atom, "label_alt_id": "?",
                "label_asym_id": rec[19].strip() or "?",
                "label_comp_id": rec[15:18].strip() or "?",
                "label_seq_id": rec[21:26].strip() or "?", "label_atom_id": atom
            })


def parse_remark_480(filestring, mmcif):
    records = re.findall(r"^REMARK 480 .+", filestring, re.M)
    if not records[8:]: return
    if "pdbx_unobs_or_zero_occ_atoms" not in mmcif:
        mmcif["pdbx_unobs_or_zero_occ_atoms"] = []
    for rec in records[8:]:
        for atom in rec[27:].strip().split():
            mmcif["pdbx_unobs_or_zero_occ_atoms"].append({
                "id": str(len(mmcif["pdbx_unobs_or_zero_occ_atoms"]) + 1),
                "PDB_model_num": "1", "polymer_flag": "Y", "occupancy_flag": "0",
                "auth_asym_id": rec[19].strip() or "?",
                "auth_comp_id": rec[15:18].strip() or "?",
                "auth_seq_id": rec[21:26].strip() or "?",
                "PDB_ins_code": rec[26].strip() or "?",
                "auth_atom_id": atom, "label_alt_id": "?",
                "label_asym_id": rec[19].strip() or "?",
                "label_comp_id": rec[15:18].strip() or "?",
                "label_seq_id": rec[21:26].strip() or "?", "label_atom_id": atom
            })
    mmcif["pdbx_unobs_or_zero_occ_atoms"].sort(key=lambda r: int(
        "".join([char for char in r["auth_seq_id"] if char.isdigit()])
    ))
    for i, row in enumerate(mmcif["pdbx_unobs_or_zero_occ_atoms"], start=1):
        row["id"] = str(i)


def parse_remark_800(filestring, mmcif):
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
    cryst1 = re.search(r"CRYST1.+", filestring, re.M)
    mmcif["cell"] = [{**{k: "?" for k in [
        "entry_id", "length_a", "length_b", "length_c", "angle_alpha",
        "angle_beta", "angle_gamma", "Z_pdb", "pdbx_unique_axis"
    ]}, "entry_id": mmcif["entry"][0]["id"], }]
    mmcif["symmetry"] = [{
        **{k: "?" for k in [
            "entry_id", "space_group_name_H-M",
            "pdbx_full_space_group_name_H-M", "cell_setting", "Int_Tables_number"
        ]},
        "entry_id": mmcif["entry"][0]["id"],
    }]
    if cryst1:
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
    mmcif["database_PDB_matrix"] = [{"entry_id": mmcif["entry"][0]["id"], **{
        k: "?" for k in [
        "origx[1][1]", "origx[2][1]", "origx[3][1]", "origx[1][2]", 
        "origx[2][2]", "origx[3][2]", "origx[1][3]", "origx[2][3]", 
        "origx[3][3]", "origx_vector[1]", "origx_vector[2]", "origx_vector[3]", 
    ]}}]
    for rec in re.findall(r"ORIGX.+", filestring, re.M):
        n = rec[5]
        mmcif["database_PDB_matrix"][0][f"origx[{n}][1]"] = rec[10:20].strip() or "?"
        mmcif["database_PDB_matrix"][0][f"origx[{n}][2]"] = rec[20:30].strip() or "?"
        mmcif["database_PDB_matrix"][0][f"origx[{n}][3]"] = rec[30:40].strip() or "?"
        mmcif["database_PDB_matrix"][0][f"origx_vector[{n}]"] = rec[45:55].strip() or "?"


def parse_scalen(filestring, mmcif):
    mmcif["atom_sites"] = [{"entry_id": mmcif["entry"][0]["id"], **{
        k: "?" for k in [
        "fract_transf_matrix[1][1]", "fract_transf_matrix[2][1]",
        "fract_transf_matrix[3][1]", "fract_transf_matrix[1][2]", 
        "fract_transf_matrix[2][2]", "fract_transf_matrix[3][2]",
        "fract_transf_matrix[1][3]", "fract_transf_matrix[2][3]", 
        "fract_transf_matrix[3][3]", "fract_transf_vector[1]",
        "fract_transf_vector[2]", "fract_transf_vector[3]", 
    ]}}]
    for rec in re.findall(r"SCALE.+", filestring, re.M):
        n = rec[5]
        if not n.strip(): continue
        mmcif["atom_sites"][0][f"fract_transf_matrix[{n}][1]"] = rec[10:20].strip() or "?"
        mmcif["atom_sites"][0][f"fract_transf_matrix[{n}][2]"] = rec[20:30].strip() or "?"
        mmcif["atom_sites"][0][f"fract_transf_matrix[{n}][3]"] = rec[30:40].strip() or "?"
        mmcif["atom_sites"][0][f"fract_transf_vector[{n}]"] = rec[45:55].strip() or "?"


def parse_mtrixn(filestring, mmcif):
    lines = re.findall(r"MTRIX\d.+", filestring, re.M)
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

    for rec in re.findall(r"MTRIX.+", filestring, re.M):
        pass


    """ mmcif["struct_ncs_oper"] = [{"entry_id": mmcif["entry"][0]["id"], **{
        k: "?" for k in [
        "matrix[1][1]", "matrix[2][1]", "matrix[3][1]", "matrix[1][2]", 
        "matrix[2][2]", "matrix[3][2]", "matrix[1][3]", "matrix[2][3]", 
        "matrix[3][3]", "vector[1]", "vector[2]", "vector[3]", 
    ]}}] """


def pdb_date_to_mmcif_date(date):
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