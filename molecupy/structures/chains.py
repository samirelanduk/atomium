from .molecules import AtomicStructure, Residue, SmallMolecule
from . import matrix
from ..exceptions import NoResiduesError

class ResiduicStructure(AtomicStructure):

    def __init__(self, *residues):
        if len(residues) == 0:
            raise NoResiduesError("Cannot make a ResiduicStructure with no residues")
        for residue in residues:
            if not isinstance(residue, Residue):
                raise TypeError(
                 "Can only make ResiduicStructure with Residues, not '%s'" % str(residue)
                )
        self._residues = set(residues)


    def __repr__(self):
        return "<%s (%i residues)>" % (self.__class__.__name__, len(self._residues))


    def __getattr__(self, attribute):
        if attribute == "_atoms":
            atoms = set()
            for residue in self.residues():
                atoms.update(residue.atoms())
            return atoms
        else:
            return self.__getattribute__(attribute)


    def residues(self, include_missing=True):
        if include_missing:
            return set(self._residues)
        else:
            return set([res for res in self._residues if not res.is_missing()])


    def add_residue(self, residue):
        if not isinstance(residue, Residue):
            raise TypeError(
             "Can only add Residues to ResiduicStructures, not '%s'" % str(residue)
            )
        self._residues.add(residue)


    def remove_residue(self, residue):
        self._residues.remove(residue)


    def get_residue_by_id(self, residue_id):
        if not isinstance(residue_id, str):
            raise TypeError(
             "Residue ID search must be by str, not '%s'" % str(residue_id)
            )
        for residue in self.residues():
            if residue.residue_id() == residue_id:
                return residue


    def get_residues_by_name(self, residue_name, include_missing=True):
        if not isinstance(residue_name, str):
            raise TypeError(
             "Residue name search must be by str, not '%s'" % str(residue_name)
            )
        return set([
         residue for residue in self.residues(include_missing=include_missing)
          if residue.residue_name() == residue_name
        ])


    def get_residue_by_name(self, residue_name, include_missing=True):
        if not isinstance(residue_name, str):
            raise TypeError(
             "Residue name search must be by str, not '%s'" % str(residue_name)
            )
        for residue in self.residues(include_missing=include_missing):
            if residue.residue_name() == residue_name:
                return residue



class ResiduicSequence(ResiduicStructure):

    def __init__(self, *residues):
        ResiduicStructure.__init__(self, *residues)
        self._residues = list(residues)


    def residues(self, include_missing=True):
        if include_missing:
            return list(self._residues)
        else:
            return [res for res in self._residues if not res.is_missing()]


    def add_residue(self, residue):
        if not isinstance(residue, Residue):
            raise TypeError(
             "Can only add Residues to ResiduicSequences, not '%s'" % str(residue)
            )
        self._residues.append(residue)


    def sequence_string(self, include_missing=True):
        return "".join([RESIDUES.get(
         res.residue_name().upper(), "X"
        ) for res in self.residues(include_missing=include_missing)])



class Chain(ResiduicSequence):

    def __init__(self, chain_id, *residues):
        if not isinstance(chain_id, str):
            raise TypeError("'%s' is not a valid chain_id" % str(chain_id))
        self._chain_id = chain_id
        ResiduicSequence.__init__(self, *residues)
        for residue in self._residues:
            residue._chain = self


    def __repr__(self):
        return "<Chain %s (%i residues)>" % (self._chain_id, len(self._residues))


    def chain_id(self):
        return self._chain_id


    def add_residue(self, residue):
        ResiduicSequence.add_residue(self, residue)
        residue._chain = self


    def remove_residue(self, residue):
        ResiduicSequence.remove_residue(self, residue)
        residue._chain = None


    generate_residue_distance_matrix = matrix.generate_residue_distance_matrix



class BindSite(ResiduicStructure):

    def __init__(self, site_id, *residues):
        if not isinstance(site_id, str):
            raise TypeError("'%s' is not a valid site_id" % str(site_id))
        self._site_id = site_id
        self._ligand = None
        ResiduicStructure.__init__(self, *residues)


    def __repr__(self):
        return "<BindSite %s (%i residues)>" % (self._site_id, len(self._residues))


    def site_id(self):
        return self._site_id


    def ligand(self, ligand=None):
        if ligand is None:
            return self._ligand
        else:
            if not isinstance(ligand, SmallMolecule):
                raise TypeError(
                 "'%s' is not a valid ligand" % str(ligand)
                )
            self._ligand = ligand
            ligand._bind_site = self


    def continuous_sequence(self):
        if len(set([res.chain() for res in self.residues() if res.chain()])) == 1:
            chain = list(self.residues())[0].chain()
            min_index = min([chain.residues().index(res) for res in self.residues()])
            max_index = max([chain.residues().index(res) for res in self.residues()])
            return ResiduicSequence(*chain.residues()[min_index:max_index])
        else:
            return None



RESIDUES = {
 "GLY": "G", "ALA": "A", "LEU": "L", "MET": "M", "PHE": "F",
 "TRP": "W", "LYS": "K", "GLN": "Q", "GLU": "E", "SER": "S",
 "PRO": "P", "VAL": "V", "ILE": "I", "CYS": "C", "TYR": "Y",
 "HIS": "H", "ARG": "R", "ASN": "N", "ASP": "D", "THR": "T"
}
