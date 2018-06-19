import math
from tests.integration.base import IntegratedTest
import atomium

class CreationTests(IntegratedTest):
    """Tests relating to creating and modifying structures."""

    def setUp(self):
        IntegratedTest.setUp(self)

        # Make 3x3 cube
        self.atoms = [
         atomium.Atom("N", id=1, name="N", bfactor=1.5),
         atomium.Atom("C", 1, 0, 0, id=2, name="CA", charge=0.9, bfactor=1.5),
         atomium.Atom("C", 2, 0, 0, id=3, name="C", charge=-0.3, bfactor=1.5),
         atomium.Atom("O", 0, 1, 0, id=4, name="OA", bfactor=1.5),
         atomium.Atom("C", 1, 1, 0, id=5, name="CB", charge=0.5, bfactor=1.5),
         atomium.Atom("O", 2, 1, 0, id=6, name="OB", bfactor=1.5),
         atomium.Atom("S", 0, 2, 0, id=7, name="S1", bfactor=1.5),
         atomium.Atom("C", 1, 2, 0, id=8, name="CC", bfactor=1.5),
         atomium.Atom("S", 2, 2, 0, id=9, name="S2", bfactor=1.5, anisotropy=[0.03, 0.4, 1.01, 0.2, 0.3, 1]),
         atomium.Atom("N", 0, 0, 1, id=10, name="N", bfactor=0.5),
         atomium.Atom("C", 1, 0, 1, id=11, name="CA", charge=0.5, bfactor=0.5),
         atomium.Atom("C", 2, 0, 1, id=12, name="C", bfactor=0.5),
         atomium.Atom("O", 0, 1, 1, id=13, name="OA", bfactor=0.5),
         atomium.Atom("C", 1, 1, 1, id=14, name="CB", charge=-2, bfactor=0.5),
         atomium.Atom("O", 2, 1, 1, id=15, name="OB", bfactor=0.5),
         atomium.Atom("P", 0, 2, 1, id=16, name="P1", bfactor=0.5),
         atomium.Atom("H", 1, 2, 1, id=17, name="CC", charge=0.5, bfactor=0.5),
         atomium.Atom("H", 2, 2, 1, id=18, name="P2", bfactor=0.5),
         atomium.Atom("N", 0, 0, 2, id=19, name="N", bfactor=1.5),
         atomium.Atom("C", 1, 0, 2, id=20, name="CA", bfactor=1.5),
         atomium.Atom("C", 2, 0, 2, id=21, name="CO", bfactor=1.5),
         atomium.Atom("O", 0, 1, 2, id=22, name="OA", bfactor=1.5),
         atomium.Atom("C", 1, 1, 2, id=23, name="CB", charge=0.5, bfactor=1.5),
         atomium.Atom("O", 2, 1, 2, id=24, name="OB", bfactor=1.5),
         atomium.Atom("FE", 0, 2, 2, name="F1", bfactor=1.5),
         atomium.Atom("C", 1, 2, 2,  name="CC", bfactor=1.5),
         atomium.Atom("FE", 2, 2, 2, name="F2", bfactor=1.5),
        ]


    def test_atom_creation(self):
        """Manipulation of pure atoms."""

        # Basic properties
        self.assertEqual(self.atoms[-1].element, "FE")
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

        # Awareness of environment
        for atom in self.atoms:
            self.assertIsNone(atom.ligand)
            self.assertIsNone(atom.residue)
            self.assertIsNone(atom.chain)
            self.assertIsNone(atom.model)
            self.assertEqual(atom.bonded_atoms, set())

        # Atom bonding
        self.atoms[2].bond_to(self.atoms[11])
        self.atoms[11].bond_to(self.atoms[20])
        self.assertEqual(self.atoms[2].bonded_atoms, {self.atoms[11]})
        self.assertEqual(self.atoms[20].bonded_atoms, {self.atoms[11]})
        self.assertEqual(
         self.atoms[11].bonded_atoms, {self.atoms[2], self.atoms[20]}
        )
        self.atoms[20].unbond_from(self.atoms[11])
        self.atoms[2].unbond_from(self.atoms[11])
        self.assertEqual(self.atoms[2].bonded_atoms, set())
        self.assertEqual(self.atoms[11].bonded_atoms, set())
        self.assertEqual(self.atoms[20].bonded_atoms, set())

        # Other properties
        self.assertAlmostEqual(self.atoms[0].mass, 14, delta=0.05)
        self.assertAlmostEqual(self.atoms[-1].mass, 56, delta=0.5)
        self.assertEqual(self.atoms[0].distance_to(self.atoms[1]), 1)
        self.assertEqual(self.atoms[0].distance_to(self.atoms[2]), 2)
        self.assertEqual(self.atoms[0].distance_to(self.atoms[-1]), 12 ** 0.5)
        self.assertEqual(self.atoms[0].distance_to((10, 0, 0)), 10)

        # Atom copying
        self.atoms[5].bond_to(self.atoms[1])
        new = self.atoms[5].copy()
        self.assertEqual(new.location, (2, 1, 0))
        self.assertEqual(new.id, 6)
        self.assertEqual(new.name, "OB")
        self.assertEqual(new.element, "O")
        self.assertEqual(new.charge, 0)
        self.assertEqual(new.bfactor, 1.5)
        self.assertEqual(new.bonded_atoms, set())
        self.assertIsNot(self.atoms[5], new)

        # Transformations
        self.atoms[0].translate(0.5, -1.9, 12)
        self.assertEqual(self.atoms[0].location, (0.5, -1.9, 12))
        self.atoms[0].translate(0.4, -2.8, 9)
        self.assertEqual(self.atoms[0].location, (0.9, -4.7, 21))
        self.atoms[1].rotate([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])
        self.assertEqual(self.atoms[1].location, (-1, 0, 0))
        self.atoms[1].rotate([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
        self.assertEqual(self.atoms[1].location, (0, -1, 0))
        self.atoms[1].rotate([[1, 0, 0], [0, 0.906, -0.423], [0, 0.423, 0.906]], trim=2)
        self.assertEqual(self.atoms[1].location, (0, -0.91, -0.42))
        self.atoms[1].translate(1, 0.91, 0.42)
        self.assertEqual(self.atoms[1].location, (1, 0, 0))
        self.atoms[0].move_to(0, 0, 0)
        self.assertEqual(self.atoms[0].location, (0, 0, 0))


    def test_simple_model_with_atoms(self):
        # Make model
        model = atomium.Model(*self.atoms)

        # Atoms are fine
        self.assertEqual(model.atoms(), set(self.atoms))
        for atom in self.atoms:
            self.assertIn(atom, model)
            self.assertIs(atom.model, model)
        model.remove(self.atoms[0])
        self.assertIsNone(self.atoms[0].model)
        self.assertNotIn(self.atoms[0], model)
        self.assertEqual(model.atoms(), set(self.atoms[1:]))
        model.add(self.atoms[0])
        self.assertIs(self.atoms[0].model, model)
        self.assertIn(self.atoms[0], model)
        self.assertEqual(model.atoms(), set(self.atoms))
        pairs = list(model.pairwise_atoms())
        for atom1 in self.atoms:
            for atom2 in self.atoms:
                if atom1 is not atom2:
                    self.assertIn({atom1, atom2}, pairs)
        self.assertEqual(model.atoms(element="FE"), {self.atoms[-3], self.atoms[-1]})
        self.assertEqual(model.atoms(name="CA"), {self.atoms[1], self.atoms[10], self.atoms[19]})
        self.assertEqual(model.atoms(is_metal=True), {self.atoms[-3], self.atoms[-1]})
        self.assertEqual(model.atoms(is_metal=False), set(self.atoms) - {self.atoms[-3], self.atoms[-1]})
        self.assertEqual(model.atoms(element__ne="H"), set(self.atoms) - {self.atoms[16], self.atoms[17]})

        # Transformation is fine
        model.translate(12, 13, 14)
        self.assertEqual(self.atoms[0].location, (12, 13, 14))
        self.assertEqual(self.atoms[2].location, (14, 13, 14))
        model.translate(-1, 1.5, 9)
        self.assertEqual(self.atoms[0].location, (11, 14.5, 23))
        self.assertEqual(self.atoms[2].location, (13, 14.5, 23))
        model.translate(-11, -14.5, -23)
        self.assertEqual(self.atoms[0].location, (0, 0, 0))
        self.assertEqual(self.atoms[2].location, (2, 0, 0))
        model.rotate([[0.985, 0, 0.174], [0, 1, 0], [-0.174, 0, 0.985]], trim=2)
        self.assertEqual(self.atoms[0].location, (0, 0, 0))
        self.assertEqual(self.atoms[2].location, (1.97, 0, -0.35))

        # Structure calculated properties
        self.assertAlmostEqual(model.mass, 479, delta=0.5)
        self.assertAlmostEqual(model.charge, 0.6, delta=0.5)
        self.assertEqual(model.formula, {"C": 11, "N": 3, "O": 6, "FE": 2, "H": 2, "S": 2, "P": 1})
        self.assertAlmostEqual(model.center_of_mass[0], 1.1017, delta=0.005)
        self.assertAlmostEqual(model.center_of_mass[1], 1.2479, delta=0.005)
        self.assertAlmostEqual(model.center_of_mass[2], 0.9213, delta=0.005)
        self.assertAlmostEqual(model.radius_of_gyration, 1.4413, delta=0.005)

        # Can get nearby atoms
        self.assertEqual(
         self.atoms[0].nearby_atoms(1), {self.atoms[1], self.atoms[3], self.atoms[9]}
        )
        self.assertEqual(
         self.atoms[0].nearby_atoms(1, element="C"), {self.atoms[1]}
        )

        # Model copying
        copy = model.copy()
        self.assertIsNot(model, copy)
        self.assertEqual(len(copy.atoms()), 27)
        self.assertFalse(model.atoms() & copy.atoms())


    def test_complex_model_with_structures(self):
        # Make ligands
        ligand1 = atomium.Ligand(*self.atoms[-9:-6], id="A100", name="BIF")
        ligand2 = atomium.Ligand(*self.atoms[-6:-3], id="B100", name="FOZ")
        ligand3 = atomium.Ligand(*self.atoms[-3:], id="C100", name="XMP")

        # Make residues
        res1 = atomium.Residue(*self.atoms[:2], id="A1", name="MET")
        res2 = atomium.Residue(*self.atoms[2:4], id="A2", name="TYR")
        res3 = atomium.Residue(*self.atoms[4:6], id="A3", name="VAL")
        res4 = atomium.Residue(*self.atoms[6:8], id="B1", name="HIS")
        res5 = atomium.Residue(*self.atoms[8:10], id="B2", name="HIS")
        res6 = atomium.Residue(*self.atoms[10:12], id="B3", name="PRO")
        res7 = atomium.Residue(*self.atoms[12:14], id="C1", name="LEU")
        res8 = atomium.Residue(*self.atoms[14:16], id="C2", name="ILE")
        res9 = atomium.Residue(*self.atoms[16:18], id="C3", name="SER")

        # Make chains and put in model
        res1.next, res2.next = res2, res3
        res4.next, res5.next = res5, res6
        res7.next, res8.next = res8, res9
        chaina = atomium.Chain(res1, res2, res3, ligand1, id="A")
        chainb = atomium.Chain(res4, res5, res6, ligand2, id="B")
        chainc = atomium.Chain(res7, res8, res9, ligand3, id="C", rep="AALIH")
        model = atomium.Model(chaina, chainb, chainc)

        # The model is fine
        self.assertEqual(model.atoms(), set(self.atoms))
        self.assertEqual(
         model.atoms(name_regex="(F|P2)"),
         {self.atoms[-3], self.atoms[-1], self.atoms[17]}
        )
        self.assertEqual(model.ligands(), {ligand1, ligand2, ligand3})
        self.assertEqual(model.ligand(name="BIF"), ligand1)
        self.assertEqual(
         model.residues(), {res1, res2, res3, res4, res5, res6, res7, res8, res9}
        )
        self.assertEqual(model.residues(name="HIS"), {res4, res5})
        self.assertEqual(model.residues(name_regex="HIS|MET"), {res1, res4, res5})
        self.assertEqual(model.residue(id="C2"), res8)
        self.assertEqual(model.chains(), {chaina, chainb, chainc})
        self.assertEqual(model.chain("A"), chaina)
        self.assertEqual(model.chain("B"), chainb)
        self.assertEqual(model.chains(id_regex="A|B"), {chaina, chainb})
        for atom in self.atoms: self.assertIs(atom.model, model)
        self.assertIs(ligand1.model, model)
        self.assertIs(res6.model, model)
        self.assertIs(chainc.model, model)

        # The chains are fine
        self.assertEqual(chaina.atoms(), set(self.atoms[:6] + self.atoms[-9:-6]))
        self.assertEqual(chainb.atoms(), set(self.atoms[6:12] + self.atoms[-6:-3]))
        self.assertEqual(chainc.atoms(), set(self.atoms[12:18] + self.atoms[-3:]))
        self.assertEqual(chaina.residues(), (res1, res2, res3))
        self.assertEqual(chainb.ligands(), {ligand2})
        for atom in self.atoms[:6] + self.atoms[-9:-6]:
            self.assertIs(atom.chain, chaina)
        for atom in self.atoms[12:18] + self.atoms[-3:]:
            self.assertIs(atom.chain, chainc)
        self.assertIs(res2.chain, chaina)
        self.assertIs(ligand1.chain, chaina)
        self.assertEqual(chaina.formula, {"C": 5, "O": 2, "N": 2})

        # The residues are fine
        self.assertEqual(res1.atoms(), set(self.atoms[:2]))
        self.assertEqual(res5.atoms(), set(self.atoms[8:10]))
        self.assertIs(self.atoms[10].residue, res6)

        # The ligands are fine
        self.assertEqual(ligand1.atoms(), set(self.atoms[-9:-6]))
        self.assertEqual(ligand2.atoms(), set(self.atoms[-6:-3]))
        self.assertEqual(ligand3.atoms(), set(self.atoms[-3:]))
        self.assertIs(self.atoms[-9].ligand, ligand1)
        self.assertIs(self.atoms[-5].ligand, ligand2)
        self.assertIs(self.atoms[-2].ligand, ligand3)
        self.assertIs(self.atoms[0].ligand, None)

        # The properties can be messed with
        res9.name = "HIS"
        self.assertEqual(model.residues(name="HIS"), {res4, res5, res9})
        self.assertEqual(chainb.residues(name="HIS"), (res4, res5))
        self.assertEqual(chaina.residues(name="HIS"), ())

        # Membership is fine
        self.assertIn(chaina, model)
        self.assertIn(ligand2, chainb)
        self.assertIn(res1, model)
        self.assertNotIn(res1, chainb)
        self.assertNotIn(res9, chaina)
        self.assertNotIn(res9, chainb)

        # Stuff can be added and removed
        chaina.remove(res2)
        self.assertIsNone(res2.chain)
        self.assertEqual(chaina.residues(), (res1, res3))
        chainb.add(ligand3)
        self.assertEqual(chainb.ligands(), {ligand3, ligand2})
        self.assertIs(ligand3.chain, chainb)

        # Specific properties - residues
        self.assertIs(res1.next, res2)
        self.assertIs(res2.next, res3)
        self.assertIsNone(res3.next)
        self.assertIs(res3.previous, res2)
        self.assertIs(res2.previous, res1)
        self.assertIsNone(res1.previous)
        self.assertEqual(res1.code, "M")
        self.assertEqual(res1.full_name, "methionine")
        res1.name = "SPLO"
        self.assertEqual(res1.code, "X")
        self.assertEqual(res1.full_name, "SPLO")

        # Specific properties - chains
        self.assertEqual(len(chainc), 3)
        self.assertEqual(chainc.length, 3)
        self.assertEqual(chainc[1], res8)
        self.assertEqual(list(chainc), [res7, res8, res9])
        for res in chainc: pass
        self.assertEqual(chainc.sequence, "LIH")
        self.assertEqual(chainc.rep_sequence, "AALIH")

        # Atom Nearby residues
        self.assertEqual(
         self.atoms[0].nearby_residues(1), {res1, res2, res5}
        )
        self.assertEqual(
         self.atoms[0].nearby_residues(2), {res1, res2, res3, res4, res5, res6, res7}
        )
        self.assertEqual(
         self.atoms[0].nearby_residues(2, ligands=True),
         {res1, res2, res3, res4, res5, res6, res7, ligand1}
        )

        # Structure Nearby residues
        self.assertEqual(
         res1.nearby_residues(1), {res2, res3, res5, res6}
        )
        self.assertEqual(
         res1.nearby_residues(1, element="C"), {res2, res3, res6}
        )

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

        # Molecule interactions
        self.assertEqual(ligand1.pairing_with(ligand2), {
         self.atoms[19]: self.atoms[22],
         self.atoms[20]: self.atoms[21],
         self.atoms[18]: self.atoms[23],
        })
        self.atoms[18].element, self.atoms[21].element = "NN"
        self.atoms[19].element, self.atoms[22].element = "OO"
        self.atoms[20].element, self.atoms[23].element = "PP"
        self.assertEqual(ligand1.pairing_with(ligand2), {
         self.atoms[18]: self.atoms[21],
         self.atoms[19]: self.atoms[22],
         self.atoms[20]: self.atoms[23],
        })
        ligand1.superimpose_onto(ligand2)
        self.assertEqual(self.atoms[18].location, self.atoms[21].location)
        self.assertEqual(self.atoms[19].location, self.atoms[22].location)
        self.assertEqual(self.atoms[20].location, self.atoms[23].location)
        self.assertEqual(ligand1.rmsd_with(ligand2), 0)
        self.atoms[18].move_to(0, 0, 0)
        self.atoms[19].move_to(0, 0, 1)
        self.atoms[20].move_to(0, 0, 2)
        self.assertAlmostEqual(ligand1.rmsd_with(ligand2), 2.08, delta=0.05)
        self.assertEqual(ligand1.rmsd_with(ligand2, superimpose=True), 0)
        self.assertEqual(self.atoms[18].location, (0, 0, 0))
        self.assertEqual(self.atoms[19].location, (0, 0, 1))
        self.assertEqual(self.atoms[20].location, (0, 0, 2))
        self.assertEqual(self.atoms[26].location, (2, 2, 2))
        self.assertEqual(self.atoms[17].location, (2, 2, 1))
        self.assertEqual(self.atoms[8].location, (2, 2, 0))



class SavingTests(IntegratedTest):
    """Tests for saving structures."""

    def test_model_saving(self):
        model = atomium.Model()
        atom1 = atomium.Atom("N", 12.0, 11.5, 1.5)
        atom2 = atomium.Atom("C", 12.5, 10, 2, id=102, bfactor=22.1)
        model.add(atom1)
        model.add(atom2)

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
