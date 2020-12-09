"""Contains logic for turning data dictionaies into a parsed Python objects."""

from .structures import *

class File:
    """When a file is parsed, the result is a ``File``. It contains the
    structure of interest, as well as meta information.

    :param str filetype: the type of file that was parsed to make this."""

    def __init__(self, filetype):
        self._filetype = filetype
        self._models = []


    def __repr__(self):
        return "<{}.{} File>".format(self._code or "", self._filetype)


    @property
    def filetype(self):
        """The filetype that this File was created from, such as .pdb or
        .cif.

        :rtype: ``str``"""

        return self._filetype


    @property
    def code(self):
        """The unique database identifer for this structure.

        :rtype: ``str``"""

        return self._code


    @property
    def title(self):
        """The structure's text description.

        :rtype: ``str``"""

        return self._title


    @property
    def deposition_date(self):
        """The date the structure was submitted for publication.

        :rtype: ``datetime.date``"""

        return self._deposition_date


    @property
    def classification(self):
        """The structure's formal classification.

        :rtype: ``str``"""

        return self._classification


    @property
    def keywords(self):
        """The structure's keyword descriptors.

        :rtype: ``list``"""

        return self._keywords


    @property
    def authors(self):
        """The structure's authors.

        :rtype: ``list``"""

        return self._authors


    @property
    def technique(self):
        """The structure's experimental technique.

        :rtype: ``str``"""

        return self._technique


    @property
    def source_organism(self):
        """The structure's original organism.

        :rtype: ``float``"""

        return self._source_organism


    @property
    def expression_system(self):
        """The organism the structure was expressed in.

        :rtype: ``float``"""

        return self._expression_system


    @property
    def missing_residues(self):
        """The residues that should be in the model but aren't.

        :rtype: ``list``"""

        return self._missing_residues


    @property
    def resolution(self):
        """The structure's resolution.

        :rtype: ``float``"""

        return self._resolution


    @property
    def rvalue(self):
        """The structure's R-value.

        :rtype: ``float``"""

        return self._rvalue


    @property
    def rfree(self):
        """The structure's R-free value.

        :rtype: ``float``"""

        return self._rfree


    @property
    def assemblies(self):
        """The structure's biological assembly instructions.

        :rtype: ``list``"""

        return self._assemblies


    @property
    def models(self):
        """The structure's models.

        :rtype: ``list``"""

        return self._models


    @property
    def model(self):
        """The structure's first model (and only model if it has only one).

        :rtype: ``Model``"""

        return self._models[0]


    def generate_assembly(self, id):
        """Generates a new model from the existing model using one of the file's
        set of assembly instructions (for which you provide the ID).

        For example:

            >>> import atomium
            >>> pdb = atomium.fetch('1xda')
            >>> pdb.model
            <Model (8 chains, 16 ligands)>
            >>> pdb.generate_assembly(1)
            <Model (2 chains, 4 ligands)>
            >>> pdb.generate_assembly(5)
            <Model (12 chains, 24 ligands)>

        :param int id: the ID of the assembly to generate.
        :rtype: ``Model``"""
        
        m = self._models[0]
        for assembly in self._assemblies:
            if assembly["id"] == id: break
        else:
            raise ValueError(f"No assembly with ID {id}")
        all_structures = []
        for t in assembly["transformations"]:
            structures = {}
            for chain_id in t["chains"]:
                for obj in list(m.chains()) + list(m.ligands() | m.waters()):
                    if obj._internal_id == chain_id:
                        copy = obj.copy()
                        if isinstance(copy, Ligand):
                            copy._chain = structures.get(obj.chain)
                        structures[obj] = copy
            atoms = set()
            for s in structures.values(): atoms.update(s.atoms())
            Atom.transform_atoms(t["matrix"], *atoms)
            Atom.translate_atoms(t["vector"], *atoms)
            all_structures += structures.values()
        return Model(*all_structures)


def data_dict_to_file(data_dict, filetype):
    """Turns an atomium data dictionary into a :py:class:`.File`.

    :param dict data_dict: the data dictionary to parse.
    :param str filetype: the file type that is being converted.
    :rtype: ``File``"""

    f = File(filetype)
    for key in data_dict.keys():
        if key != "models":
            for subkey, value in data_dict[key].items():
                setattr(f, "_" + subkey, value)
    f._models = [model_dict_to_model(m) for m in data_dict["models"]]
    return f


def model_dict_to_model(model_dict):
    """Takes a model dictionary and turns it into a fully processed
    :py:class:`.Model` object.

    :param dict model_dict: the model dictionary.
    :rtype: ``Model``"""


    chains = create_chains(model_dict)
    ligands = create_ligands(model_dict, chains)
    waters = create_ligands(model_dict, chains, water=True)
    model = Model(*(chains + ligands + waters))
    return model


def create_chains(model_dict):
    """Creates a list of :py:class:`.Chain` objects from a model dictionary.

    :param dict model_dict: the model dictionary.
    :rtype: ``list``"""

    chains = []
    for chain_id, chain in model_dict["polymer"].items():
        res = [create_het(r, i) for i, r in sorted(
         chain["residues"].items(), key=lambda x: x[1]["number"]
        )]
        res_by_id = {r.id: r for r in res}
        for res1, res2 in zip(res[:-1], res[1:]):
            res1._next, res2._previous = res2, res1
        chains.append(Chain(
         *res, id=chain_id,
         helices=[[res_by_id[r] for r in h] for h in chain["helices"]],
         strands=[[res_by_id[r] for r in s] for s in chain["strands"]],
         internal_id=chain["internal_id"], sequence=chain["sequence"]
        ))
    return chains


def create_ligands(model_dict, chains, water=False):
    """Creates a list of :py:class:`.Ligand` objects from a model dictionary.

    :param dict model_dict: the model dictionary.
    :param list chains: a list of :py:class:`.Chain` objects to assign by ID.
    :param bool water: if `True``, water ligands will be made.
    :rtype: ``list``"""

    ligands = []
    for lig_id, lig in model_dict["water" if water else "non-polymer"].items():
        chain = None
        for c in chains:
            if c._id == lig["polymer"]:
                chain = c
                break
        ligands.append(
         create_het(lig, lig_id, ligand=True, chain=chain, water=water)
        )
    return ligands


def create_het(d, id, ligand=False, chain=None, water=False):
    """Creates a :py:class:`.Residue` or :py:class:`.Ligand` from some
    atom-containing dictionary.

    If there is multiple occupancy, only one position will be used.

    :param dict d: the dictionary to parse.
    :param str id: the ID of the structure to make.
    :param bool ligand: if ``True`` a ligand will be made, not a residue.
    :param Chain chain: the :py:class:`.Chain` to assign if a ligand.
    :param bool water: if ``True``, the ligand will be a water ligand.
    :rtype: ``Residue`` or ``Ligand``"""

    alt_loc = None
    if any([atom["occupancy"] < 1 for atom in d["atoms"].values()]):
        if any([atom["alt_loc"] for atom in d["atoms"].values()]):
            alt_loc = sorted([atom["alt_loc"] for atom in d["atoms"].values()
             if atom["alt_loc"]])[0]
    atoms = [atom_dict_to_atom(a, i) for i, a in d["atoms"].items()
     if a["occupancy"] == 1 or a["alt_loc"] is None or a["alt_loc"] == alt_loc]
    if ligand:
        return Ligand(*atoms, id=id, name=d["name"], chain=chain,
         internal_id=d["internal_id"], water=water, full_name=d["full_name"])
    else:
        return Residue(*atoms, id=id, name=d["name"], full_name=d["full_name"])


def atom_dict_to_atom(d, atom_id):
    """Creates an :py:class:`.Atom` from an atom dictionary.

    :param dict d: the atom dictionary.
    :param int id: the atom's ID.
    :rtype: ``Atom``"""

    return Atom(
     d["element"], d["x"], d["y"], d["z"], atom_id,
     d["name"], d["charge"], d["bvalue"], d["anisotropy"], d["is_hetatm"]
    )



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
 "TL": 204.3833, "PB": 207.2, "BI": 208.9804, "PO": 209, "AT": 210, "RN": 222,
 "FR": 223, "RA": 226, "AC": 227, "PA": 231.0359, "TH": 232.0381, "NP": 237,
 "U": 238.0289, "AM": 243, "PU": 244, "CM": 247, "BK": 247, "CF": 251,
 "ES": 252, "FM": 257, "MD": 258, "NO": 259, "RF": 261, "LR": 262, "DB": 262,
 "BH": 264, "SG": 266, "MT": 268, "RG": 272, "HS": 277
}

COVALENT_RADII = {
 "H": 0.31, "HE": 0.28, "LI": 1.28, "BE": 0.96, "B": 0.85, "C": 0.76,
 "N": 0.71, "O": 0.66, "F": 0.57, "NE": 0.58, "NA": 1.66, "MG": 1.41,
 "AL": 1.21, "SI": 1.11, "P": 1.07, "S": 1.05, "CL": 1.02, "AR": 1.06,
 "K": 2.03, "CA": 1.76, "SC": 1.7, "TI": 1.6, "V": 1.53, "CR": 1.39, "MN": 1.39,
 "FE": 1.32, "CO": 1.26, "NI": 1.24, "CU": 1.32, "ZN": 1.22, "GA": 1.22,
 "GE": 1.2, "AS": 1.19, "SE": 1.2, "BR": 1.2, "KR": 1.16, "RB": 2.2, "SR": 1.95,
 "Y": 1.9, "ZR": 1.75, "NB": 1.64, "MO": 1.54, "TC": 1.47, "RU": 1.46,
 "RH": 1.42, "PD": 1.39, "AG": 1.45, "CD": 1.44, "IN": 1.42, "SN": 1.39,
 "SB": 1.39, "TE": 1.38, "I": 1.39, "XE": 1.4, "CS": 2.44, "BA": 2.15,
 "LA": 2.07, "CE": 2.04, "PR": 2.03, "ND": 2.01, "PM": 1.99, "SM": 1.98,
 "EU": 1.98, "GD": 1.96, "TB": 1.94, "DY": 1.92, "HO": 1.92, "ER": 1.89,
 "TM": 1.9, "YB": 1.87, "LU": 1.87, "HF": 1.75, "TA": 1.7, "W": 1.62,
 "RE": 1.51, "OS": 1.44, "IR": 1.41, "PT": 1.36, "AU": 1.36, "HG": 1.32,
 "TL": 1.45, "PB": 1.46, "BI": 1.48, "PO": 1.4, "AT": 1.5, "RN": 1.5, "FR": 2.6,
 "RA": 2.21, "AC": 2.15, "TH": 2.06, "PA": 2.0, "U": 1.96, "NP": 1.9,
 "PU": 1.87, "AM": 1.8, "CM": 1.69
}

METALS = [
 "LI", "BE", "NA", "MG", "AL", "K", "CA", "SC", "TI", "V", "CR", "MN", "FE",
 "CO", "NI", "CU", "ZN", "HA", "RB", "SR", "Y", "ZR", "NB", "MO", "TC", "RU",
 "RH", "PD", "AG", "CD", "IN", "SN", "CS", "BA", "LA", "CE", "PR", "ND", "PM",
 "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU", "HF", "TA", "W",
 "RE", "OS", "IR", "PT", "AU", "HG", "TL", "PB", "BI", "PO", "FR", "RA", "AC",
 "TH", "PA", "U", "NP", "PU", "AM", "CM", "BK", "CF", "ES", "FM", "MD", "NO",
 "LR", "RF", "DB", "SG", "BH", "HS", "MT", "DS", "RG", "CN", "UUT", "FL", "LV"
]

FULL_NAMES = {
 "GLY": "glycine", "ALA": "alanine", "VAL": "valine", "LEU": "leucine",
 "ILE": "isoleucine", "MET": "methionine", "PHE": "phenylalanine",
 "TRP": "tryptophan", "PRO": "proline", "SER": "serine", "THR": "threonine",
 "CYS": "cysteine", "TYR": "tyrosine", "ASN": "asparagine", "GLN": "glutamine",
 "ASP": "aspartic acid", "GLU": "glutamic acid", "LYS": "lysine",
 "ARG": "arginine", "HIS": "histidine", "HOH": "water"
}

CODES = {
 "VAL": "V", "ILE": "I", "LEU": "L", "GLU": "E", "GLN": "Q", "ASP": "D",
 "ASN": "N", "HIS": "H", "TRP": "W", "PHE": "F", "TYR": "Y", "ARG": "R",
 "LYS": "K", "SER": "S", "THR": "T", "MET": "M", "ALA": "A", "GLY": "G",
 "PRO": "P", "CYS": "C", "HIP": "H", "HIE": "H",
 "DA": "A", "DG": "G", "DC": "C", "DT": "T", "A": "A", "G": "G", "C": "C",
 "U": "U"
}
