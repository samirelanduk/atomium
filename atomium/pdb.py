import re
import calendar

def pdb_string_to_mmcif_dict(filestring):
    mmcif = {}
    parse_header(filestring, mmcif)
    parse_obslte(filestring, mmcif)
    parse_title(filestring, mmcif)
    parse_split(filestring, mmcif)
    parse_sprsde(filestring, mmcif)
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
    if title and set(title) != {"-"} and title not in ["NULL", "NONE"]:
        mmcif["struct"][0]["title"] = " ".join(title.split())


def parse_split(filestring, mmcif):
    split_lines = re.findall(r"^SPLIT.+", filestring, re.M)
    if split_lines:
        codes = [code for l in split_lines for code in l[11:].strip().split()]
        mmcif["pdbx_database_related"] = [{
            "db_name": "PDB", "db_id": code, "content_type": "split",
            "details": f"Split {n}"
        } for n, code in enumerate(codes, start=1)]


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
        



def pdb_date_to_mmcif_date(date):
    day, month, year = date.split("-")
    month = str(list(calendar.month_abbr).index(month.title())).zfill(2)
    if len(year) == 2:
        if int(year) > 50: year = "19" + year
        if int(year) <= 50: year = "20" + year
    return f"{year}-{month}-{day}"