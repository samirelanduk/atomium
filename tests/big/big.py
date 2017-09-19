import requests
import xml.etree.ElementTree as ET

def get_all_codes():
    response = requests.get("https://www.rcsb.org/pdb/rest/getCurrent")
    root = ET.fromstring(response.text)
    codes = []
    for child in root:
        codes.append(child.attrib["structureId"])
    return codes
