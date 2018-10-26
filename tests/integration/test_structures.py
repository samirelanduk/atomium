import math
from unittest import TestCase
import atomium

class CreationTests(TestCase):
    """Tests relating to creating and modifying structures."""

    def setUp(self):

        # Make 3x3 cube
        self.atoms = [
         atomium.Atom("N", 0, 0, 0, 1, "N", 0, 1.5, [0] * 6),
         atomium.Atom("C", 1, 0, 0, 2, "CA", 0.9, 1.5, [0] * 6),
         atomium.Atom("C", 2, 0, 0, 3, "C", -0.3, 1.5, [0] * 6),
         atomium.Atom("O", 0, 1, 0, 4, "OA", 0, 1.5, [0] * 6),
         atomium.Atom("C", 1, 1, 0, 5, "CB", 0.5, 1.5, [0] * 6),
         atomium.Atom("O", 2, 1, 0, 6, "OB", 0, 1.5, [0] * 6),
         atomium.Atom("S", 0, 2, 0, 7, "S1", 0, 1.5, [0] * 6),
         atomium.Atom("C", 1, 2, 0, 8, "CC", 0, 1.5, [0] * 6),
         atomium.Atom("S", 2, 2, 0, 9, "S2", 0, 1.5, [0.03, 0.4, 1.01, 0.2, 0.3, 1]),
         atomium.Atom("N", 0, 0, 1, 10, "N", 0, 0.5, [0] * 6),
         atomium.Atom("C", 1, 0, 1, 11, "CA", 0.5, 0.5, [0] * 6),
         atomium.Atom("C", 2, 0, 1, 12, "C", 0, 0.5, [0] * 6),
         atomium.Atom("O", 0, 1, 1, 13, "OA", 0, 0.5, [0] * 6),
         atomium.Atom("C", 1, 1, 1, 14, "CB", -2, 0.5, [0] * 6),
         atomium.Atom("O", 2, 1, 1, 15, "OB", 0, 0.5, [0] * 6),
         atomium.Atom("P", 0, 2, 1, 16, "P1", 0, 0.5, [0] * 6),
         atomium.Atom("H", 1, 2, 1, 17, "CC", 0.5, 0.5, [0] * 6),
         atomium.Atom("H", 2, 2, 1, 18, "P2", 0, 0.5, [0] * 6),
         atomium.Atom("N", 0, 0, 2, 19, "N", 0, 1.5, [0] * 6),
         atomium.Atom("C", 1, 0, 2, 20, "CA", 0, 1.5, [0] * 6),
         atomium.Atom("C", 2, 0, 2, 21, "CO", 0, 1.5, [0] * 6),
         atomium.Atom("O", 0, 1, 2, 22, "OA", 0, 1.5, [0] * 6),
         atomium.Atom("C", 1, 1, 2, 23, "CB", 0.5, 1.5, [0] * 6),
         atomium.Atom("O", 2, 1, 2, 24, "OB", 0, 1.5, [0] * 6),
         atomium.Atom("FE", 0, 2, 2, 25, "F1", 0, 1.5, [0] * 6),
         atomium.Atom("C", 1, 2, 2, 26, "CC", 0, 1.5, [0] * 6),
         atomium.Atom("FE", 2, 2, 2, 27, "F2", 0, 1.5, [0] * 6),
        ]


    def test_atom_creation(self):
        """Manipulation of pure atoms."""

        # Basic properties
        self.assertEqual(self.atoms[-1].element, "FE")
        self.assertEqual(self.atoms[-1].name, "F2")
        self.assertEqual(self.atoms[0].id, 1)
        self.assertEqual(self.atoms[13].charge, -2)
        self.assertEqual(self.atoms[13].location, (1, 1, 1))
        self.assertEqual(self.atoms[0].location, (0, 0, 0))
        self.assertEqual(self.atoms[13].bvalue, 0.5)
        self.assertEqual(self.atoms[13].anisotropy, [0] * 6)
        self.assertEqual(self.atoms[8].anisotropy, [0.03, 0.4, 1.01, 0.2, 0.3, 1])
        self.atoms[0].x = 4
        self.assertEqual(self.atoms[0].location, (4, 0, 0))
        self.atoms[0].x -= 4
        self.assertEqual(self.atoms[0].location, (0, 0, 0))

        # Awareness of environment
        for atom in self.atoms:
            self.assertIsNone(atom.structure)
            self.assertIsNone(atom.chain)
            self.assertIsNone(atom.model)

        # Other properties
        self.assertAlmostEqual(self.atoms[0].mass, 14, delta=0.05)
        self.assertAlmostEqual(self.atoms[-1].mass, 56, delta=0.5)
        self.assertEqual(self.atoms[0].distance_to(self.atoms[1]), 1)
        self.assertEqual(self.atoms[0].distance_to(self.atoms[2]), 2)
        self.assertEqual(self.atoms[0].distance_to(self.atoms[-1]), 12 ** 0.5)
        self.assertEqual(self.atoms[0].distance_to((10, 0, 0)), 10)

        # Atom copying
        new = self.atoms[5].copy()
        self.assertEqual(new.location, (2, 1, 0))
        self.assertEqual(new.id, 6)
        self.assertEqual(new.name, "OB")
        self.assertEqual(new.element, "O")
        self.assertEqual(new.charge, 0)
        self.assertEqual(new.bvalue, 1.5)
        self.assertIsNot(self.atoms[5], new)

        # Transformations
        self.atoms[0].translate(0.5, -1.9, 12)
        self.assertEqual(self.atoms[0].location, (0.5, -1.9, 12))
        self.atoms[0].translate(0.4, -2.8, 9)
        self.assertEqual(self.atoms[0].location, (0.9, -4.7, 21))
        self.atoms[1].transform([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])
        self.assertEqual(self.atoms[1].location, (-1, 0, 0))
        self.atoms[1].transform([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
        self.assertEqual(self.atoms[1].location, (0, -1, 0))
        self.atoms[1].transform([[1, 0, 0], [0, 0.906, -0.423], [0, 0.423, 0.906]], trim=2)
        self.assertEqual(self.atoms[1].location, (0, -0.91, -0.42))
        self.atoms[1].translate(1, 0.91, 0.42)
        self.assertEqual(self.atoms[1].location, (1, 0, 0))
        self.atoms[0].move_to(0, 0, 0)
        self.assertEqual(self.atoms[0].location, (0, 0, 0))


    def test_simple_structure_with_atoms(self):
        # Make mol
        mol = atomium.Ligand(*self.atoms)

        # Atoms are fine
        self.assertEqual(mol.atoms(), set(self.atoms))
        for atom in self.atoms:
            self.assertIn(atom, mol)
            self.assertIs(atom.structure, mol)
        mol.remove(self.atoms[0])
        self.assertIsNone(self.atoms[0].structure)
        self.assertNotIn(self.atoms[0], mol)
        self.assertEqual(mol.atoms(), set(self.atoms[1:]))
        mol.add(self.atoms[0])
        self.assertIs(self.atoms[0].structure, mol)
        self.assertIn(self.atoms[0], mol)
        self.assertEqual(mol.atoms(), set(self.atoms))
        pairs = list(mol.pairwise_atoms())
        for atom1 in self.atoms:
            for atom2 in self.atoms:
                if atom1 is not atom2:
                    self.assertIn({atom1, atom2}, pairs)
        self.assertEqual(mol.atoms(element="FE"), {self.atoms[-3], self.atoms[-1]})
        self.assertEqual(mol.atoms(name="CA"), {self.atoms[1], self.atoms[10], self.atoms[19]})
        self.assertEqual(mol.atoms(is_metal=True), {self.atoms[-3], self.atoms[-1]})
        self.assertEqual(mol.atoms(is_metal=False), set(self.atoms) - {self.atoms[-3], self.atoms[-1]})
        self.assertEqual(mol.atoms(element__ne="H"), set(self.atoms) - {self.atoms[16], self.atoms[17]})

        # Geometry is fine
        self.assertEqual(self.atoms[0].angle(self.atoms[1], self.atoms[3]), math.pi / 2)
        self.assertEqual(self.atoms[0].angle(self.atoms[1], self.atoms[2]), 0)
        self.assertAlmostEqual(
         self.atoms[0].angle(self.atoms[1], self.atoms[4]),
         math.pi / 4, delta=0.00001
        )
        self.assertAlmostEqual(
         self.atoms[8].angle(self.atoms[7], self.atoms[4]),
         math.pi / 4, delta=0.00001
        )

        # Transformation is fine
        mol.translate(12, 13, 14)
        self.assertEqual(self.atoms[0].location, (12, 13, 14))
        self.assertEqual(self.atoms[2].location, (14, 13, 14))
        mol.translate(-1, 1.5, 9)
        self.assertEqual(self.atoms[0].location, (11, 14.5, 23))
        self.assertEqual(self.atoms[2].location, (13, 14.5, 23))
        mol.translate(-11, -14.5, -23)
        self.assertEqual(self.atoms[0].location, (0, 0, 0))
        self.assertEqual(self.atoms[2].location, (2, 0, 0))
        mol.transform([[1, 0, 0], [0, 1, 0], [0, 0, 1]]) # identity
        self.assertEqual(self.atoms[0].location, (0, 0, 0))
        self.assertEqual(self.atoms[2].location, (2, 0, 0))
        mol.transform([[1, 0, 0], [0, 0, -1], [0, 1, 0]]) # 90 around x-axis
        self.assertEqual(self.atoms[0].location, (0, 0, 0))
        self.assertEqual(self.atoms[2].location, (2, 0, 0))
        self.assertEqual(self.atoms[3].location, (0, 0, 1))
        mol.transform([[1, 0, 0], [0, 0, 1], [0, -1, 0]]) # back again
        self.assertEqual(self.atoms[0].location, (0, 0, 0))
        self.assertEqual(self.atoms[2].location, (2, 0, 0))
        self.assertEqual(self.atoms[3].location, (0, 1, 0))
        mol.transform([[0.985, 0, 0.174], [0, 1, 0], [-0.174, 0, 0.985]], trim=2) # 10 around y-axis
        self.assertEqual(self.atoms[0].location, (0, 0, 0))
        self.assertEqual(self.atoms[2].location, (1.97, 0, -0.35))
        mol.rotate(math.pi * 2, "y", trim=2)
        self.assertEqual(self.atoms[0].location, (0, 0, 0))
        self.assertEqual(self.atoms[2].location, (1.97, 0, -0.35))

        # Structure calculated properties
        self.assertAlmostEqual(mol.mass, 479, delta=0.5)
        self.assertAlmostEqual(mol.charge, 0.6, delta=0.5)
        self.assertEqual(mol.formula, {"C": 11, "N": 3, "O": 6, "FE": 2, "H": 2, "S": 2, "P": 1})
        self.assertAlmostEqual(mol.center_of_mass[0], 1.1017, delta=0.005)
        self.assertAlmostEqual(mol.center_of_mass[1], 1.2479, delta=0.005)
        self.assertAlmostEqual(mol.center_of_mass[2], 0.9213, delta=0.005)
        self.assertAlmostEqual(mol.radius_of_gyration, 1.4413, delta=0.005)

        # Model copying
        copy = mol.copy()
        self.assertIsNot(mol, copy)
        self.assertEqual(len(copy.atoms()), 27)
        self.assertFalse(mol.atoms() & copy.atoms())


    def test_complex_model_with_structures(self):
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
        chaina = atomium.Chain(res1, res2, res3, id="A")
        chainb = atomium.Chain(res4, res5, res6, id="B")
        chainc = atomium.Chain(res7, res8, res9, id="C", sequence="AALIH")

        # Make ligands
        ligand1 = atomium.Ligand(*self.atoms[-9:-6], id="A100", name="BIF", chain=chaina)
        ligand2 = atomium.Ligand(*self.atoms[-6:-3], id="B100", name="FOZ", chain=chainb)
        ligand3 = atomium.Ligand(*self.atoms[-3:], id="C100", name="XMP", chain=chainc)

        model = atomium.Model(chaina, chainb, chainc, ligand1, ligand2, ligand3)

        # The model is fine
        self.assertEqual(model.atoms(), set(self.atoms))
        self.assertEqual(
         model.atoms(name__regex="(F|P2)"),
         {self.atoms[-3], self.atoms[-1], self.atoms[17]}
        )
        self.assertEqual(model.ligands(), {ligand1, ligand2, ligand3})
        self.assertEqual(model.ligand(name="BIF"), ligand1)
        self.assertEqual(
         model.residues(), {res1, res2, res3, res4, res5, res6, res7, res8, res9}
        )
        self.assertEqual(model.residues(name="HIS"), {res4, res5})
        self.assertEqual(model.residues(name__regex="HIS|MET"), {res1, res4, res5})
        self.assertEqual(model.residue(id="C2"), res8)
        self.assertEqual(model.chains(), {chaina, chainb, chainc})
        self.assertEqual(model.chain("A"), chaina)
        self.assertEqual(model.chain("B"), chainb)
        self.assertEqual(model.chains(id__regex="A|B"), {chaina, chainb})
        for atom in self.atoms: self.assertIs(atom.model, model)
        self.assertIs(ligand1.model, model)
        self.assertIs(res6.model, model)
        self.assertIs(chainc.model, model)

        # Can get nearby atoms
        self.assertEqual(
         model.atoms_in_sphere([0, 0, 0], 1), {self.atoms[0], self.atoms[1], self.atoms[3], self.atoms[9]}
        )
        self.assertEqual(
         self.atoms[0].nearby_atoms(1), {self.atoms[1], self.atoms[3], self.atoms[9]}
        )
        self.assertEqual(
         self.atoms[0].nearby_atoms(1, element="C"), {self.atoms[1]}
        )

        # The chains are fine
        self.assertEqual(chaina.id, "A")
        self.assertEqual(chaina.atoms(), set(self.atoms[:6]))
        self.assertEqual(chainb.atoms(), set(self.atoms[6:12]))
        self.assertEqual(chainc.atoms(), set(self.atoms[12:18]))
        self.assertEqual(chaina.residues(), (res1, res2, res3))
        self.assertEqual(chainb.ligands(), {ligand2})
        for atom in self.atoms[:6] + self.atoms[-9:-6]:
            self.assertIs(atom.chain, chaina)
        for atom in self.atoms[12:18] + self.atoms[-3:]:
            self.assertIs(atom.chain, chainc)
        self.assertIs(res2.chain, chaina)
        self.assertIs(ligand1.chain, chaina)
        self.assertEqual(chaina.formula, {"C": 3, "O": 2, "N": 1})

        # The residues are fine
        self.assertEqual(res1.atoms(), set(self.atoms[:2]))
        self.assertEqual(res5.atoms(), set(self.atoms[8:10]))
        self.assertIs(self.atoms[10].structure, res6)

        # The ligands are fine
        self.assertEqual(ligand1.atoms(), set(self.atoms[-9:-6]))
        self.assertEqual(ligand2.atoms(), set(self.atoms[-6:-3]))
        self.assertEqual(ligand3.atoms(), set(self.atoms[-3:]))
        self.assertEqual(ligand1.id, "A100")
        self.assertEqual(ligand1.name, "BIF")
        self.assertIs(self.atoms[-9].structure, ligand1)
        self.assertIs(self.atoms[-5].structure, ligand2)
        self.assertIs(self.atoms[-2].structure, ligand3)

        # The properties can be messed with
        res9.name = "HIS"
        self.assertEqual(model.residues(name="HIS"), {res4, res5, res9})
        self.assertEqual(chainb.residues(name="HIS"), (res4, res5))
        self.assertEqual(chaina.residues(name="HIS"), ())

        # Membership is fine
        self.assertIn(chaina, model)
        self.assertIn(ligand2, model)
        self.assertIn(res1, model)
        self.assertIn(self.atoms[10], chainb)
        self.assertNotIn(res1, chainb)
        self.assertNotIn(res9, chaina)
        self.assertNotIn(res9, chainb)

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
        self.assertEqual(chainc.sequence, "AALIH")

        # Atom Nearby residues
        self.assertEqual(
         self.atoms[0].nearby_structures(1, ligands=False), {res2, res5}
        )
        self.assertEqual(
         self.atoms[0].nearby_structures(2, ligands=False), {res2, res3, res4, res5, res6, res7}
        )
        self.assertEqual(
         self.atoms[0].nearby_structures(2,),
         {res2, res3, res4, res5, res6, res7, ligand1}
        )

        # Structure Nearby residues
        self.assertEqual(
         res1.nearby_atoms(1), set(self.atoms[2:5] + self.atoms[9:11])
        )
        self.assertEqual(
         res1.nearby_atoms(1, element="C"), {self.atoms[2], self.atoms[4], self.atoms[10]}
        )
        self.assertEqual(
         res1.nearby_structures(1, ligands=False), {res2, res3, res5, res6}
        )
        self.assertEqual(
         res1.nearby_structures(1, ligands=False, element="C"), {res2, res3, res6}
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
        self.assertEqual(ligand1.rmsd_with(ligand2), 0)

        # Chain copying
        chain = chaina.copy()
        self.assertIsNot(chain, chaina)
        self.assertEqual(len(chain.atoms()), 6)
        self.assertFalse(chain.atoms() & chaina.atoms())
        self.assertEqual(len(chain.residues()), 3)
        self.assertFalse(set(chain.residues()) & set(chaina.residues()))
        self.assertEqual(chain.sequence, chaina.sequence)
        self.assertEqual(chain._internal_id, chaina._internal_id)

        # Stuff can be added and removed
        model.remove(chaina)
        self.assertIsNone(chaina.model)
        self.assertIsNone(self.atoms[0].model)
        self.assertEqual(model.chains(), {chainb, chainc})
        model.add(chaina)
        self.assertIs(chaina.model, model)
        self.assertIs(self.atoms[0].model, model)
        self.assertEqual(model.chains(), {chaina, chainb, chainc})
