import math
from tests.integration.base import IntegratedTest
from atomium.structures import Model, Atom, Residue, Chain, Molecule, Complex
import atomium

class CreationTests(IntegratedTest):
    """Tests relating to creating and modifying structures."""

    def setUp(self):
        IntegratedTest.setUp(self)
        self.atoms = [
         Atom("N", id=1, name="N", bfactor=1.5),
         Atom("C", 1, 0, 0, id=2, name="CA", charge=0.9, bfactor=1.5),
         Atom("C", 2, 0, 0, id=3, name="C", charge=-0.3, bfactor=1.5),
         Atom("O", 0, 1, 0, id=4, name="OA", bfactor=1.5),
         Atom("C", 1, 1, 0, id=5, name="CB", charge=0.5, bfactor=1.5),
         Atom("O", 2, 1, 0, id=6, name="OB", bfactor=1.5),
         Atom("S", 0, 2, 0, id=7, name="S1", bfactor=1.5),
         Atom("C", 1, 2, 0, id=8, name="CC", bfactor=1.5),
         Atom("S", 2, 2, 0, id=9, name="S2", bfactor=1.5, anisotropy=[0.03, 0.4, 1.01, 0.2, 0.3, 1]),
         Atom("N", 0, 0, 1, id=10, name="N", bfactor=0.5),
         Atom("C", 1, 0, 1, id=11, name="CA", charge=0.5, bfactor=0.5),
         Atom("C", 2, 0, 1, id=12, name="C", bfactor=0.5),
         Atom("O", 0, 1, 1, id=13, name="OA", bfactor=0.5),
         Atom("C", 1, 1, 1, id=14, name="CB", charge=-2, bfactor=0.5),
         Atom("O", 2, 1, 1, id=15, name="OB", bfactor=0.5),
         Atom("P", 0, 2, 1, id=16, name="P1", bfactor=0.5),
         Atom("H", 1, 2, 1, id=17, name="CC", charge=0.5, bfactor=0.5),
         Atom("H", 2, 2, 1, id=18, name="P2", bfactor=0.5),
         Atom("N", 0, 0, 2, id=19, name="N", bfactor=1.5),
         Atom("C", 1, 0, 2, id=20, name="CA", bfactor=1.5),
         Atom("C", 2, 0, 2, id=21, name="CO", bfactor=1.5),
         Atom("O", 0, 1, 2, id=22, name="OA", bfactor=1.5),
         Atom("C", 1, 1, 2, id=23, name="CB", charge=0.5, bfactor=1.5),
         Atom("O", 2, 1, 2, id=24, name="OB", bfactor=1.5),
         Atom("Fe", 0, 2, 2, name="F1", bfactor=1.5),
         Atom("C", 1, 2, 2,  name="CC", bfactor=1.5),
         Atom("Fe", 2, 2, 2, name="F2", bfactor=1.5),
        ]


    def test_atom_creation(self):
        """Manipulation of pure atoms."""

        # Basic properties
        self.assertEqual(self.atoms[-1].element, "Fe")
        self.assertEqual(self.atoms[-1].name, "F2")
        self.assertEqual(self.atoms[0].id, 1)
        self.assertEqual(self.atoms[-3].id, 0)
        self.assertEqual(self.atoms[-2].id, 0)
        self.assertEqual(self.atoms[-1].id, 0)
        self.assertEqual(self.atoms[13].charge, -2)
        self.assertEqual(self.atoms[13].location, (1, 1, 1))
        self.assertEqual(self.atoms[0].location, (0, 0, 0))
        self.assertEqual(self.atoms[13].bfactor, 0.5)
        self.assertEqual(self.atoms[13].anisotropy, [0] * 6)
        self.assertEqual(self.atoms[8].anisotropy, [0.03, 0.4, 1.01, 0.2, 0.3, 1])
        self.atoms[0].x = 4
        self.assertEqual(self.atoms[0].location, (4, 0, 0))
        self.atoms[0].x -= 4
        self.assertEqual(self.atoms[0].location, (0, 0, 0))

        # Transformations
        self.atoms[0].translate(0.5, -1.9, 12)
        self.assertEqual(self.atoms[0].location, (0.5, -1.9, 12))
        self.atoms[0].translate(0.4, -2.8, 9)
        self.assertEqual(self.atoms[0].location, (0.9, -4.7, 21))
        self.atoms[1].rotate(math.pi, "y")
        self.assertEqual(self.atoms[1].location, (-1, 0, 0))
        self.atoms[1].rotate(90, "z", degrees=True)
        self.assertEqual(self.atoms[1].location, (0, -1, 0))
        self.atoms[1].rotate(25, "x", degrees=True, trim=2)
        self.assertEqual(self.atoms[1].location, (0, -0.91, -0.42))
        self.atoms[1].translate(1, 0.91, 0.42)
        self.assertEqual(self.atoms[1].location, (1, 0, 0))
        self.atoms[0].move_to(0, 0, 0)
        self.assertEqual(self.atoms[0].location, (0, 0, 0))

        # Other properties
        self.assertAlmostEqual(self.atoms[0].mass, 14, delta=0.05)
        self.assertAlmostEqual(self.atoms[-1].mass, 56, delta=0.5)
        self.assertEqual(self.atoms[0].distance_to(self.atoms[1]), 1)
        self.assertEqual(self.atoms[0].distance_to(self.atoms[2]), 2)
        self.assertEqual(self.atoms[0].distance_to(self.atoms[-1]), 12 ** 0.5)
        self.assertEqual(self.atoms[0].distance_to((10, 0, 0)), 10)

        # Awareness of environment
        for atom in self.atoms:
            self.assertIsNone(atom.molecule)
            self.assertIsNone(atom.residue)
            self.assertIsNone(atom.chain)
            self.assertIsNone(atom.model)
            self.assertEqual(atom.bonds, set())
            self.assertEqual(atom.bonded_atoms(), set())

        # Atom bonding
        self.atoms[2].bond_to(self.atoms[11])
        self.atoms[11].bond_to(self.atoms[20])
        self.assertEqual(self.atoms[2].bonded_atoms(), {self.atoms[11]})
        self.assertEqual(self.atoms[20].bonded_atoms(), {self.atoms[11]})
        self.assertEqual(
         self.atoms[11].bonded_atoms(), {self.atoms[2], self.atoms[20]}
        )
        self.assertEqual(self.atoms[11].bonded_atoms(name="C"), {self.atoms[2]})
        self.assertEqual(self.atoms[2].bond_with(self.atoms[11]).length, 1)

        self.atoms[20].bond_with(self.atoms[11]).destroy()
        self.atoms[2].unbond_from(self.atoms[11])
        self.assertEqual(self.atoms[2].bonded_atoms(), set())
        self.assertEqual(self.atoms[11].bonded_atoms(), set())
        self.assertEqual(self.atoms[20].bonded_atoms(), set())

        # Atom copying
        self.atoms[5].bond_to(self.atoms[1])
        new = self.atoms[5].copy()
        self.assertEqual(new.location, (2, 1, 0))
        self.assertEqual(new.id, 6)
        self.assertEqual(new.name, "OB")
        self.assertEqual(new.element, "O")
        self.assertEqual(new.charge, 0)
        self.assertEqual(new.bfactor, 1.5)
        self.assertEqual(new.bonds, set())
        self.assertEqual(new.bonded_atoms(), set())
        self.assertIsNot(self.atoms[5], new)


    def test_atoms_in_molecules(self):
        """Manipulation atoms in molecules and residues."""

        # Small molecule - atoms
        molecule = Molecule(*self.atoms[:3], id="A1000", name="ABC")
        molecule2 = Molecule(self.atoms[8], self.atoms[17], self.atoms[26])
        self.assertIn(self.atoms[0], molecule)
        self.assertIn(self.atoms[1], molecule)
        self.assertIn(self.atoms[2], molecule)
        self.assertNotIn(self.atoms[3], molecule)
        self.assertEqual(molecule.atoms(), set(self.atoms[:3]))
        self.assertEqual(molecule.atoms(element="C"), set(self.atoms[1:3]))
        self.assertIs(molecule.atom(1), self.atoms[0])
        self.assertIs(molecule.atom(id=2), self.atoms[1])
        self.assertIs(molecule.atom(3), self.atoms[2])
        self.assertIs(molecule.atom(name="CA"), self.atoms[1])
        molecule.add_atom(self.atoms[3])
        self.assertEqual(len(molecule.atoms()), 4)
        self.assertIn(self.atoms[3], molecule)
        self.assertIs(self.atoms[3].molecule, molecule)
        molecule.remove_atom(self.atoms[3])
        self.assertIs(self.atoms[3].molecule, None)
        self.assertEqual(len(molecule.atoms()), 3)
        self.assertNotIn(self.atoms[3], molecule)
        for pair in molecule.pairwise_atoms(element="C"):
            self.assertIn(self.atoms[1], pair)
            self.assertIn(self.atoms[2], pair)

        # Small molecule transformation
        molecule.translate(12, 13, 14)
        self.assertEqual(self.atoms[0].location, (12, 13, 14))
        self.assertEqual(self.atoms[1].location, (13, 13, 14))
        self.assertEqual(self.atoms[2].location, (14, 13, 14))
        molecule.translate(-1, 1.5, 9)
        self.assertEqual(self.atoms[0].location, (11, 14.5, 23))
        self.assertEqual(self.atoms[1].location, (12, 14.5, 23))
        self.assertEqual(self.atoms[2].location, (13, 14.5, 23))
        molecule.translate(-11, -14.5, -23)
        self.assertEqual(self.atoms[0].location, (0, 0, 0))
        self.assertEqual(self.atoms[1].location, (1, 0, 0))
        self.assertEqual(self.atoms[2].location, (2, 0, 0))
        molecule.rotate(10, "y", degrees=True, trim=2)
        self.assertEqual(self.atoms[0].location, (0, 0, 0))
        self.assertEqual(self.atoms[1].location, (0.98, 0, -0.17))
        self.assertEqual(self.atoms[2].location, (1.97, 0, -0.35))
        molecule2.rotate(-370, "z", degrees=True, trim=2)
        self.assertEqual(self.atoms[8].location, (2.32, 1.62, 0))
        self.assertEqual(self.atoms[17].location, (2.32, 1.62, 1))
        self.assertEqual(self.atoms[26].location, (2.32, 1.62, 2))
        molecule.rotate(-10, "y", degrees=True, trim=1)
        molecule2.rotate(370, "z", degrees=True, trim=1)
        self.assertEqual(self.atoms[0].location, (0, 0, 0))
        self.assertEqual(self.atoms[1].location, (1, 0, 0))
        self.assertEqual(self.atoms[2].location, (2, 0, 0))
        self.assertEqual(self.atoms[8].location, (2, 2, 0))
        self.assertEqual(self.atoms[17].location, (2, 2, 1))
        self.assertEqual(self.atoms[26].location, (2, 2, 2))

        # Molecule calculated properties
        self.assertAlmostEqual(molecule.mass, 38, delta=0.5)
        self.assertAlmostEqual(molecule.charge, 0.6, delta=0.5)
        self.assertEqual(molecule.formula, {"C": 2, "N": 1})
        self.assertEqual(molecule.center_of_mass, (0.9475124973374951, 0, 0))
        self.assertEqual(molecule.radius_of_gyration, 0.8181818896812696)

        # Molecule attributes
        self.assertEqual(molecule.id, "A1000")
        self.assertEqual(molecule.name, "ABC")
        self.assertIsNone(molecule2.id)
        self.assertIsNone(molecule2.name)
        molecule.name = "DEF"
        self.assertEqual(molecule.name, "DEF")
        self.assertIs(self.atoms[0].molecule, molecule)
        self.assertIs(self.atoms[1].molecule, molecule)
        self.assertIs(self.atoms[2].molecule, molecule)
        self.assertIsNone(molecule.model)
        self.assertIsNone(molecule2.model)

        # Molecule interactions
        self.assertEqual(molecule.pairing_with(molecule2), {
         self.atoms[0]: self.atoms[8],
         self.atoms[1]: self.atoms[17],
         self.atoms[2]: self.atoms[26],
        })
        self.atoms[26].element = "N"
        self.atoms[17].element = "C"
        self.atoms[8].element = "C"
        self.atoms[17].name = "CA"
        self.atoms[8].name = "C"
        self.assertEqual(molecule.pairing_with(molecule2), {
         self.atoms[0]: self.atoms[26],
         self.atoms[1]: self.atoms[17],
         self.atoms[2]: self.atoms[8],
        })
        molecule.superimpose_onto(molecule2)
        self.assertEqual(self.atoms[0].location, (2, 2, 2))
        self.assertEqual(self.atoms[1].location, (2, 2, 1))
        self.assertEqual(self.atoms[2].location, (2, 2, 0))
        self.assertEqual(self.atoms[26].location, (2, 2, 2))
        self.assertEqual(self.atoms[17].location, (2, 2, 1))
        self.assertEqual(self.atoms[8].location, (2, 2, 0))
        self.assertEqual(molecule.rmsd_with(molecule2), 0)
        self.atoms[0].move_to(0, 0, 0)
        self.atoms[1].move_to(0, 0, 1)
        self.atoms[2].move_to(0, 0, 2)
        self.assertAlmostEqual(molecule.rmsd_with(molecule2), 3.27, delta=0.05)
        self.assertEqual(molecule.rmsd_with(molecule2, superimpose=True), 0)
        self.assertEqual(self.atoms[0].location, (0, 0, 0))
        self.assertEqual(self.atoms[1].location, (0, 0, 1))
        self.assertEqual(self.atoms[2].location, (0, 0, 2))
        self.assertEqual(self.atoms[26].location, (2, 2, 2))
        self.assertEqual(self.atoms[17].location, (2, 2, 1))
        self.assertEqual(self.atoms[8].location, (2, 2, 0))

        # Residues
        res1 = Residue(*self.atoms[3:6], id="A1", name="VAL")
        res2 = Residue(*self.atoms[6:9], id="A2", name="ASP")
        res3 = Residue(*self.atoms[9:12], id="A2A", name="TRP")
        self.assertEqual(res1.atoms(), set(self.atoms[3:6]))
        self.assertIn(self.atoms[6], res2)
        self.assertEqual(res1.id, "A1")
        self.assertEqual(res2.id, "A2")
        self.assertEqual(res3.id, "A2A")
        self.assertEqual(res1.name, "VAL")
        res2.name = "GLU"
        self.assertEqual(res2.name, "GLU")
        res1.next = res2
        res3.previous = res2
        self.assertIs(res1.next, res2)
        self.assertIs(res2.next, res3)
        self.assertIs(res3.previous, res2)
        self.assertIs(res2.previous, res1)
        self.assertIsNone(res1.previous)
        self.assertIsNone(res3.next)
        res2.previous = None
        self.assertIsNone(res1.next)
        self.assertIsNone(res1.chain)
        self.assertIsNone(res2.chain)
        self.assertIsNone(res3.chain)
        self.assertEqual(res1.full_name, "valine")
        self.assertEqual(res2.full_name, "glutamic acid")
        self.assertEqual(res3.full_name, "tryptophan")
        self.assertEqual(res1.code, "V")
        self.assertEqual(res2.code, "E")
        self.assertEqual(res3.code, "W")
        res3.name = "XYZ"
        self.assertEqual(res3.full_name, "XYZ")
        self.assertEqual(res3.code, "X")

        # Copies
        copy = molecule.copy()
        self.assertEqual(len(copy.atoms()), 3)
        self.assertNotIn(self.atoms[0], copy)
        self.assertNotIn(self.atoms[1], copy)
        self.assertNotIn(self.atoms[2], copy)
        self.assertEqual(copy.id, "A1000")
        self.assertEqual(copy.name, "DEF")
        copy = res2.copy()
        self.assertEqual(len(copy.atoms()), 3)
        self.assertNotIn(self.atoms[6], copy)
        self.assertNotIn(self.atoms[7], copy)
        self.assertNotIn(self.atoms[8], copy)
        self.assertEqual(copy.id, "A2")
        self.assertEqual(copy.name, "GLU")


    def test_atoms_in_chains(self):
        """Manipilation of residues in chains"""

        # Make residues first
        res1 = Residue(*self.atoms[:3], id="A1", name="MET")
        res2 = Residue(*self.atoms[3:6], id="A2", name="TYR")
        res3 = Residue(*self.atoms[6:9], id="A3", name="CYS")
        res4 = Residue(*self.atoms[18:21], id="B1", name="MET")
        res5 = Residue(*self.atoms[21:24], id="B2", name="HIS")
        res6 = Residue(*self.atoms[24:27], id="B3", name="ILE")
        res1.next, res2.next, res4.next, res5.next = res2, res3, res5, res6

        # Chains
        chaina = Chain(res1, res2, res3, id="A")
        chainb = Chain(res4, res5, res6, id="B")
        self.assertEqual(chaina.atoms(), set(self.atoms[:9]))
        self.assertEqual(chainb.atoms(), set(self.atoms[18:]))
        self.assertEqual(chaina.residues(), (res1, res2, res3))
        self.assertEqual(chainb.residues(), (res4, res5, res6))
        self.assertEqual(chaina.residues(name="MET"), (res1,))
        res2.name = "MET"
        self.assertEqual(chaina.residues(name="MET"), (res1, res2))
        self.assertIs(chaina.residue(name="CYS"), res3)
        self.assertIs(chaina.residue("A1"), res1)
        self.assertEqual(chaina.length, 3)
        self.assertEqual(len(chainb), 3)
        self.assertIs(self.atoms[1].chain, chaina)
        self.assertIs(self.atoms[24].chain, chainb)
        self.assertIs(res6.chain, chainb)
        chainb.remove(res6)
        self.assertEqual(chainb.length, 2)
        self.assertIsNone(self.atoms[-1].chain)
        self.assertIsNone(res6.chain)
        chainb.add(res6)
        self.assertEqual(chainb.length, 3)
        self.assertIs(self.atoms[-1].chain, chainb)
        self.assertIs(res6.chain, chainb)
        self.assertEqual(chaina.sequence, "MMC")
        self.assertEqual(chainb.sequence, "MHI")

        # Complexes
        complex = Complex(chaina, chainb, id="1", name="HEAVY")
        self.assertEqual(complex.chains(), {chaina, chainb})
        for atom in self.atoms[:9] + self.atoms[18:]:
            self.assertIs(atom.complex, complex)


    def test_atoms_in_models(self):
        """Full model processing"""

        # Setup
        res1 = Residue(*self.atoms[:3], id="A1", name="MET")
        res2 = Residue(*self.atoms[3:6], id="A2", name="TYR")
        res3 = Residue(*self.atoms[6:9], id="A3", name="CYS")
        res4 = Residue(*self.atoms[9:12], id="B1", name="MET")
        res5 = Residue(*self.atoms[12:15], id="B2", name="HIS")
        res6 = Residue(*self.atoms[15:18], id="B3", name="ILE")
        res1.next, res2.next, res4.next, res5.next = res2, res3, res5, res6
        chaina = Chain(res1, res2, res3, id="A")
        chainb = Chain(res4, res5, res6, id="B")
        mol1 = Molecule(*self.atoms[18:21], id="A1000", name="XMP")
        mol2 = Molecule(*self.atoms[21:24], id="A1001", name="BIS")
        mol3 = Molecule(*self.atoms[24:27], id="A1002", name="BIS")
        complex = Complex(chaina, chainb, id="1", name="HEAVY")
        model = Model(complex, mol1, mol2, mol3)

        # Atoms in model
        for atom in self.atoms:
            self.assertIn(atom, model)
            self.assertIs(atom.model, model)
        self.assertEqual(model.atoms(), set(self.atoms))
        self.assertEqual(
         model.atoms(hydrogen=False), set(self.atoms) - set(self.atoms[16:18])
        )
        self.assertEqual(model.atoms(het=False), set(self.atoms[:18]))
        self.assertEqual(
         model.atoms(metal=False),
         set(self.atoms) - {self.atoms[-3], self.atoms[-1]}
        )
        self.assertEqual(self.atoms[0].nearby_atoms(0.5), set())
        self.assertEqual(self.atoms[0].nearby_atoms(1), {
         self.atoms[1], self.atoms[3], self.atoms[9]
        })
        self.assertEqual(model.atoms_in_sphere(0, 0, 0, 0.5), {self.atoms[0]})
        self.assertEqual(len(list(model.pairwise_atoms())), 351)
        self.assertAlmostEqual(model.mass, 479, delta=0.5)
        self.assertEqual(
         model.formula,
         {"C": 11, "N": 3, "O": 6, "Fe": 2, "P": 1, "S": 2, "H": 2}
        )

        # Residues in model
        self.assertEqual(model.residues(), {res1, res2, res3, res4, res5, res6})
        self.assertIs(model.residue("A2"), res2)
        self.assertIs(model.residue("B3"), res6)

        # Molecules in model
        self.assertEqual(model.molecules(), {chaina, chainb, mol1, mol2, mol3})
        self.assertEqual(model.molecules(generic=True), {mol1, mol2, mol3})
        mol1.name = "HOH"
        self.assertEqual(model.molecules(water=False, generic=True), {mol2, mol3})
        self.assertIs(model.molecule("A1002"), mol3)
        self.assertEqual(model.molecules(name="BIS"), {mol2, mol3})

        # Chains in model
        self.assertEqual(model.chains(), {chaina, chainb})
        self.assertIs(model.chain("B"), chainb)

        # Complexes in model
        self.assertEqual(model.complexes(name="HEAVY"), {complex})
        self.assertEqual(model.complex("1"), complex)

        # Model grid
        self.assertEqual(list(model.grid()), [(x, y, z)
         for x in range(3) for y in range(3) for z in range(3)])
        model.translate(-0.1, -0.1, -0.1)
        self.assertEqual(list(model.grid()), [(x, y, z)
         for x in range(-1, 3) for y in range(-1, 3) for z in range(-1, 3)])
        model.translate(-0.8, -0.8, -0.8)
        self.assertEqual(list(model.grid()), [(x, y, z)
         for x in range(-1, 3) for y in range(-1, 3) for z in range(-1, 3)])
        model.translate(1, 1, 1)
        self.assertEqual(list(model.grid()), [(x, y, z)
         for x in range(4) for y in range(4) for z in range(4)])
        model.translate(-0.1, -0.1, -0.1)
        self.assertEqual(list(model.grid(size=0.5)), [(x, y, z)
         for x in [0, 0.5, 1, 1.5, 2]
         for y in [0, 0.5, 1, 1.5, 2]
         for z in [0, 0.5, 1, 1.5, 2]])
        self.assertEqual(list(model.grid(size=0.5, margin=1.5)), [(x, y, z)
         for x in [-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5]
         for y in [-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5]
         for z in [-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5]])
        self.atoms[0].translate(-0.1, -0.5, -0.6)
        self.atoms[-1].translate(0.1, 0.5, 0.6)
        self.assertEqual(list(model.grid(size=0.25)), [(x, y, z)
         for x in [-0.25, 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25]
         for y in [-0.5, -0.25, 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5]
         for z in [-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75]])
        self.atoms[0].translate(0.1, 0.5, 0.6)
        self.atoms[-1].translate(-0.1, -0.5, -0.6)
        self.assertEqual(list(model.grid(size=5)), [(x, y, z)
         for x in [0, 5]
         for y in [0, 5]
         for z in [0, 5]])

        # Molecule sites
        site = mol2.site()
        self.assertEqual(site.residues(), {res3, res6, res5, res2})
        site = mol2.site(water=True)
        self.assertEqual(site.residues(), {res3, res6, res5, res2, mol1})
        site = mol2.site(cutoff=1)
        self.assertEqual(site.residues(), {res5})
        site = mol2.site(cutoff=1, water=True)
        self.assertEqual(site.residues(), {res5, mol1})
        for atom in mol1.atoms():
            atom.element = "C"
        site = mol2.site(cutoff=1, water=True, carbon=False)
        self.assertEqual(site.residues(), {res5})



class SavingTests(IntegratedTest):
    """Tests for saving structures."""

    def test_model_saving(self):

        model = Model()
        atom1 = Atom("N", 12.0, 11.5, 1.5)
        atom2 = Atom("C", 12.5, 10, 2, id=102, bfactor=22.1)
        model.add_atom(atom1)
        model.add_atom(atom2)

        model.save("tests/integration/files/model.xyz", description="Some atoms")
        new = atomium.xyz_from_file("tests/integration/files/model.xyz")
        self.assertEqual(new.title, "Some atoms")
        self.assertEqual(len(new.model.atoms()), 2)
        self.assertEqual(new.model.atom(element="N").location, (12, 11.5, 1.5))
        self.assertEqual(new.model.atom(element="C").location, (12.5, 10, 2))

        model.save("tests/integration/files/model.pdb", description="Some atoms")
        new = atomium.pdb_from_file("tests/integration/files/model.pdb")
        self.assertEqual(new.title, "Some atoms")
        self.assertEqual(len(new.model.atoms()), 2)
        self.assertEqual(new.model.atom(0).location, (12, 11.5, 1.5))
        self.assertEqual(new.model.atom(102).bfactor, 22.1)
