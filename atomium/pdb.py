import re
import calendar

def pdb_string_to_mmcif_dict(filestring):
    mmcif = {}
    parse_header(filestring, mmcif)
    return mmcif


def parse_header(filestring, mmcif):
    header = re.search(f"HEADER.+", filestring, re.M).group(0)
    code = (header[62:66].strip() or "1XXX") if header else "1XXX"
    if set(code) == {"-"}: code = "1XXX"
    mmcif["entry"] = [{"id": code}]
    date = header[50:59].strip()
    day, month, year = date.split("-")
    month = str(list(calendar.month_abbr).index(month.title())).zfill(2)
    if len(year) == 2:
        if int(year) > 50: year = "19" + year
        if int(year) <= 50: year = "20" + year
    mmcif["pdbx_database_status"] = [{
        "status_code": "REL", "entry_id": "1LOL",
        "recvd_initial_deposition_date": f"{year}-{month}-{day}"
    }]
    keyword = (header[10:50].strip() or "?") if header else "?"
    if set(keyword) == {"-"} or keyword in ["NULL", "NONE"]: keyword = "?"
    mmcif["struct_keywords"] = [{
        "entry_id": code, "pdbx_keywords": keyword, "text": "?"
    }]