from unittest import TestCase
from molecupy import exceptions
from molecupy.structures import PdbSmallMolecule, PdbAtom, AtomicStructure
from molecupy.structures import PdbResidue, PdbChain, PdbSite, PdbModel

class SmallMoleculeTest(TestCase):

    def setUp(self):
        self.atom1 = PdbAtom(1.0, 1.0, 1.0, "H", atom_id=1, atom_name="H1")
        self.atom2 = PdbAtom(1.0, 1.0, 2.0, "C", atom_id=2, atom_name="CA")


    def check_valid_small_molecule(self, small_molecule):
        self.assertIsInstance(small_molecule, PdbSmallMolecule)
        self.assertIsInstance(small_molecule, AtomicStructure)
        for atom in small_molecule.atoms:
            self.assertEqual(atom.molecule, small_molecule)
        self.assertIsInstance(small_molecule.molecule_id, str)
        self.assertIsInstance(small_molecule.molecule_name, str)
        small_molecule.model
        self.assertRegex(str(small_molecule), r"<SmallMolecule \((.+)\)>")



class SmallMoleculeCreationTest(SmallMoleculeTest):

    def test_can_create_small_molecule(self):
        small_molecule = PdbSmallMolecule("A1", "HET", self.atom1, self.atom2)
        self.check_valid_small_molecule(small_molecule)


    def test_molecule_id_must_be_str(self):
        with self.assertRaises(TypeError):
            small_molecule = PdbSmallMolecule(1.1, "HET", self.atom1, self.atom2)
        with self.assertRaises(TypeError):
            small_molecule = PdbSmallMolecule(1, "HET", self.atom1, self.atom2)


    def test_molecule_name_must_be_str(self):
        with self.assertRaises(TypeError):
            small_molecule = PdbSmallMolecule("A1", 1, self.atom1, self.atom2)



class BindingSiteDetectionTests(SmallMoleculeTest):

    def setUp(self):
        SmallMoleculeTest.setUp(self)
        atom3 = PdbAtom(1.0, 1.0, 3.0, "C", 3, "C")
        self.residue = PdbResidue("A1", "VAL", atom3)
        self.chain = PdbChain("A", self.residue)
        self.site = PdbSite("S1", self.residue)
        self.small_molecule = PdbSmallMolecule("A1", "HET", self.atom1, self.atom2)
        self.site.ligand = self.small_molecule
        self.model = PdbModel()
        self.model.add_small_molecule(self.small_molecule)
        self.model.add_chain(self.chain)
        self.model.add_site(self.site)

    def test_small_molecule_can_get_annotated_site(self):
        self.assertEqual(
         self.small_molecule.get_binding_site(),
         self.site
        )
