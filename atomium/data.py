CODES = {
    "VAL": "V", "ILE": "I", "LEU": "L", "GLU": "E", "GLN": "Q", "ASP": "D",
    "ASN": "N", "HIS": "H", "TRP": "W", "PHE": "F", "TYR": "Y", "ARG": "R",
    "LYS": "K", "SER": "S", "THR": "T", "MET": "M", "ALA": "A", "GLY": "G",
    "PRO": "P", "CYS": "C", "HIP": "H", "HIE": "H",
    "DA": "A", "DG": "G", "DC": "C", "DT": "T", "A": "A", "G": "G", "C": "C",
    "U": "U"
}

FULL_NAMES = {
    "GLY": "glycine", "ALA": "alanine", "VAL": "valine", "LEU": "leucine",
    "ILE": "isoleucine", "MET": "methionine", "PHE": "phenylalanine",
    "TRP": "tryptophan", "PRO": "proline", "SER": "serine", "THR": "threonine",
    "CYS": "cysteine", "TYR": "tyrosine", "ASN": "asparagine",
    "GLN": "glutamine", "ASP": "aspartic acid", "GLU": "glutamic acid",
    "LYS": "lysine", "ARG": "arginine", "HIS": "histidine",
    "A": "ADENOSINE-5'-MONOPHOSPHATE", "C": "CYTIDINE-5'-MONOPHOSPHATE",
    "G": "GUANOSINE-5'-MONOPHOSPHATE", "U": "URIDINE-5'-MONOPHOSPHATE",
    "DA": "2'-DEOXYADENOSINE-5'-MONOPHOSPHATE",
    "DC": "2'-DEOXYCYTIDINE-5'-MONOPHOSPHATE",
    "DG": "2'-DEOXYGUANOSINE-5'-MONOPHOSPHATE",
    "DT": "THYMIDINE-5'-MONOPHOSPHATE"
}

FORMULAE = {
    "GLY": "C2 H5 N O2", "ALA": "C3 H7 N O2", "VAL": "C5 H11 N O2",
    "LEU": "C6 H13 N O2", "ILE": "C6 H13 N O2", "MET": "C5 H11 N O2 S",
    "PHE": "C9 H11 N O2", "TRP": "C11 H12 N2 O2", "PRO": "C5 H9 N O2",
    "SER": "C3 H7 N O3", "THR": "C4 H9 N O3", "CYS": "C3 H7 N O2 S",
    "TYR": "C9 H11 N O3", "ASN": "C4 H8 N2 O3", "GLN": "C5 H10 N2 O3",
    "ASP": "C4 H7 N O4", "GLU": "C5 H9 N O4", "LYS": "C6 H15 N2 O2 1",
    "ARG": "C6 H15 N4 O2 1", "HIS": "C6 H10 N3 O2 1", "A": "C10 H14 N5 O7 P",
    "C": "C9 H14 N3 O8 P", "G": "C10 H14 N5 O8 P", "U": "C9 H13 N2 O9 P",
    "DA": "C10 H14 N5 O6 P", "DC": "C9 H14 N3 O7 P", "DG": "C10 H14 N5 O7 P",
    "DT": "C10 H15 N2 O8 P"
}

RESIDUE_MASSES = {
    "GLY": 75.067, "ALA": 89.093, "VAL": 117.146, "LEU": 131.173,
    "ILE": 131.173, "MET": 149.211, "PHE": 165.189, "TRP": 204.225,
    "PRO": 115.130, "SER": 115.130, "THR": 119.119, "CYS": 121.158,
    "TYR": 181.189, "ASN": 132.118, "GLN": 146.144, "ASP": 133.103,
    "GLU": 147.129, "LYS": 147.195, "ARG": 175.209, "HIS": 156.162,
    "A": 347.22, "C": 323.20, "G": 363.22, "U": 324.18,
    "DA": 331.22, "DC": 307.20, "DG": 347.22, "DT":322.21
}

WATER_NAMES = ["HOH", "WAT"]

PERIODIC_TABLE = {
    "H": 1.0079, "HE": 4.0026, "LI": 6.941, "BE": 9.0122, "B": 10.811,
    "C": 12.0107, "N": 14.0067, "O": 15.9994, "F": 18.9984, "NE": 20.1797,
    "NA": 22.9897, "MG": 24.305, "AL": 26.9815, "SI": 28.0855, "P": 30.9738,
    "S": 32.065, "CL": 35.453, "K": 39.0983, "AR": 39.948, "CA": 40.078,
    "SC": 44.9559, "TI": 47.867, "V": 50.9415, "CR": 51.9961, "MN": 54.938,
    "FE": 55.845, "NI": 58.6934, "CO": 58.9332, "CU": 63.546, "ZN": 65.39,
    "GA": 69.723, "GE": 72.64, "AS": 74.9216, "SE": 78.96, "BR": 79.904,
    "KR": 83.8, "RB": 85.4678, "SR": 87.62, "Y": 88.9059, "ZR": 91.224,
    "NB": 92.9064, "MO": 95.94, "TC": 98, "RU": 101.07, "RH": 102.9055,
    "PD": 106.42, "AG": 107.8682, "CD": 112.411, "IN": 114.818, "SN": 118.71,
    "SB": 121.76, "I": 126.9045, "TE": 127.6, "XE": 131.293, "CS": 132.9055,
    "BA": 137.327, "LA": 138.9055, "CE": 140.116, "PR": 140.9077, "ND": 144.24,
    "PM": 145, "SM": 150.36, "EU": 151.964, "GD": 157.25, "TB": 158.9253,
    "DY": 162.5, "HO": 164.9303, "ER": 167.259, "TM": 168.9342, "YB": 173.04,
    "LU": 174.967, "HF": 178.49, "TA": 180.9479, "W": 183.84, "RE": 186.207,
    "OS": 190.23, "IR": 192.217, "PT": 195.078, "AU": 196.9665, "HG": 200.59,
    "TL": 204.3833, "PB": 207.2, "BI": 208.9804, "PO": 209, "AT": 210,
    "RN": 222, "FR": 223, "RA": 226, "AC": 227, "PA": 231.0359, "TH": 232.0381,
    "NP": 237, "U": 238.0289, "AM": 243, "PU": 244, "CM": 247, "BK": 247,
    "CF": 251, "ES": 252, "FM": 257, "MD": 258, "NO": 259, "RF": 261, "LR": 262,
    "DB": 262, "BH": 264, "SG": 266, "MT": 268, "RG": 272, "HS": 277
}

COVALENT_RADII = {
    "H": 0.31, "HE": 0.28, "LI": 1.28, "BE": 0.96, "B": 0.85, "C": 0.76,
    "N": 0.71, "O": 0.66, "F": 0.57, "NE": 0.58, "NA": 1.66, "MG": 1.41,
    "AL": 1.21, "SI": 1.11, "P": 1.07, "S": 1.05, "CL": 1.02, "AR": 1.06,
    "K": 2.03, "CA": 1.76, "SC": 1.7, "TI": 1.6, "V": 1.53, "CR": 1.39,
    "MN": 1.39, "FE": 1.32, "CO": 1.26, "NI": 1.24, "CU": 1.32, "ZN": 1.22,
    "GA": 1.22, "GE": 1.2, "AS": 1.19, "SE": 1.2, "BR": 1.2, "KR": 1.16,
    "RB": 2.2, "SR": 1.95, "Y": 1.9, "ZR": 1.75, "NB": 1.64, "MO": 1.54,
    "TC": 1.47, "RU": 1.46, "RH": 1.42, "PD": 1.39, "AG": 1.45, "CD": 1.44,
    "IN": 1.42, "SN": 1.39, "SB": 1.39, "TE": 1.38, "I": 1.39, "XE": 1.4,
    "CS": 2.44, "BA": 2.15, "LA": 2.07, "CE": 2.04, "PR": 2.03, "ND": 2.01,
    "PM": 1.99, "SM": 1.98, "EU": 1.98, "GD": 1.96, "TB": 1.94, "DY": 1.92,
    "HO": 1.92, "ER": 1.89, "TM": 1.9, "YB": 1.87, "LU": 1.87, "HF": 1.75,
    "TA": 1.7, "W": 1.62, "RE": 1.51, "OS": 1.44, "IR": 1.41, "PT": 1.36,
    "AU": 1.36, "HG": 1.32, "TL": 1.45, "PB": 1.46, "BI": 1.48, "PO": 1.4,
    "AT": 1.5, "RN": 1.5, "FR": 2.6, "RA": 2.21, "AC": 2.15, "TH": 2.06,
    "PA": 2.0, "U": 1.96, "NP": 1.9, "PU": 1.87, "AM": 1.8, "CM": 1.69
}

METALS = [
    "LI", "BE", "NA", "MG", "AL", "K", "CA", "SC", "TI", "V", "CR", "MN", "FE",
    "CO", "NI", "CU", "ZN", "HA", "RB", "SR", "Y", "ZR", "NB", "MO", "TC", "RU",
    "RH", "PD", "AG", "CD", "IN", "SN", "CS", "BA", "LA", "CE", "PR", "ND",
    "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU", "HF",
    "TA", "W", "RE", "OS", "IR", "PT", "AU", "HG", "TL", "PB", "BI", "PO", "FR",
    "RA", "AC", "TH", "PA", "U", "NP", "PU", "AM", "CM", "BK", "CF", "ES", "FM",
    "MD", "NO", "LR", "RF", "DB", "SG", "BH", "HS", "MT", "DS", "RG", "CN",
    "UUT", "FL", "LV"
]