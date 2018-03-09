import math
from tests.integration.base import IntegratedTest
from atomium.structures import Model, Atom, Residue, Chain, Molecule
import atomium

class StructureTests(IntegratedTest):

    def test_structures(self):
        # Make atoms
        atoms = [
         Atom("N", 0, 0, 0, 1, "N", bfactor=1.5),
         Atom("C", 1, 0, 0, 2, "CA", bfactor=1.5),
         Atom("C", 2, 0, 0, 3, "C", bfactor=1.5),
         Atom("O", 0, 1, 0, 4, "OA", bfactor=1.5),
         Atom("C", 1, 1, 0, 5, "CB", charge=0.5, bfactor=1.5),
         Atom("O", 2, 1, 0, 6, "OB", bfactor=1.5),
         Atom("S", 0, 2, 0, 7, "S1", bfactor=1.5),
         Atom("C", 1, 2, 0, 8, "CC", bfactor=1.5),
         Atom("S", 2, 2, 0, 9, "S2", bfactor=1.5),
         Atom("N", 0, 0, 1, 10, "N", bfactor=0.5),
         Atom("C", 1, 0, 1, 11, "CA", charge=0.5, bfactor=0.5),
         Atom("C", 2, 0, 1, 12, "C", bfactor=0.5),
         Atom("O", 0, 1, 1, 13, "OA", bfactor=0.5),
         Atom("C", 1, 1, 1, 14, "CB", charge=-2, bfactor=0.5),
         Atom("O", 2, 1, 1, 15, "OB", bfactor=0.5),
         Atom("P", 0, 2, 1, 16, "P1", bfactor=0.5),
         Atom("H", 1, 2, 1, 17, "CC", charge=0.5, bfactor=0.5),
         Atom("H", 2, 2, 1, 18, "P2", bfactor=0.5),
         Atom("N", 0, 0, 2, 19, "N", bfactor=1.5),
         Atom("C", 1, 0, 2, 20, "CA", bfactor=1.5),
         Atom("C", 2, 0, 2, 21, "CO", bfactor=1.5),
         Atom("O", 0, 1, 2, 22, "OA", bfactor=1.5),
         Atom("C", 1, 1, 2, 23, "CB", charge=0.5, bfactor=1.5),
         Atom("O", 2, 1, 2, 24, "OB", bfactor=1.5),
         Atom("Fe", 0, 2, 2, 25, "F1", bfactor=1.5),
         Atom("C", 1, 2, 2, 26, "CC", bfactor=1.5),
         Atom("Fe", 2, 2, 2, 27, "F2", bfactor=1.5),
        ]
        self.assertEqual(atoms[-1].element(), "Fe")
        self.assertEqual(atoms[-1].name(), "F2")
        self.assertEqual(atoms[0].atom_id(), 1)
        self.assertEqual(atoms[13].charge(), -2)
        self.assertEqual(atoms[13].location(), (1, 1, 1))
        self.assertEqual(atoms[13].bfactor(), 0.5)
        for atom in atoms:
            self.assertIsNone(atom.molecule())
            self.assertIsNone(atom.residue())
            self.assertIsNone(atom.chain())
            self.assertIsNone(atom.model())
            self.assertEqual(atom.bonds(), set())
            self.assertEqual(atom.bonded_atoms(), set())
        atoms[0].x(0.1)
        self.assertEqual(atoms[0].x(), 0.1)
        atoms[0].round(0)
        self.assertEqual(atoms[0].x(), 0)
        atoms[13].translate(10, 20, 30)
        self.assertEqual(atoms[13].location(), (11, 21, 31))
        atoms[13].translate(-10, -20, -30)
        self.assertEqual(atoms[13].location(), (1, 1, 1))
        atoms[13].rotate(math.pi / 2, "x")
        atoms[13].round(8)
        self.assertEqual(atoms[13].location(), (1, -1, 1))
        atoms[13].rotate(math.pi / 2, "x")
        atoms[13].round(8)
        self.assertEqual(atoms[13].location(), (1, -1, -1))
        atoms[13].rotate(math.pi, "x")
        atoms[13].round(8)
        self.assertEqual(atoms[13].location(), (1, 1, 1))
        self.assertAlmostEqual(atoms[0].mass(), 14, delta=0.5)
        self.assertAlmostEqual(atoms[-1].mass(), 56, delta=0.5)
        self.assertEqual(atoms[0].distance_to(atoms[1]), 1)
        self.assertEqual(atoms[0].distance_to(atoms[2]), 2)
        self.assertEqual(atoms[0].distance_to(atoms[3]), 1)
        self.assertEqual(atoms[0].distance_to(atoms[4]), 2 ** 0.5)
        self.assertEqual(atoms[0].distance_to((10, 0, 0)), 10)

        # Bond atoms
        atoms[2].bond(atoms[11])
        atoms[11].bond(atoms[20])
        self.assertEqual(atoms[2].bonded_atoms(), set([atoms[11]]))
        self.assertEqual(atoms[20].bonded_atoms(), set([atoms[11]]))
        self.assertEqual(atoms[11].bonded_atoms(), set([atoms[2], atoms[20]]))
        self.assertEqual(atoms[11].bonded_atoms(name="C"), set([atoms[2]]))
        self.assertEqual(atoms[2].bond_with(atoms[11]).length(), 1)
        self.assertEqual(atoms[2].bond_with(atoms[11]).angle_with(
         atoms[20].bond_with(atoms[11])
        ), math.pi)
        atoms[20].bond_with(atoms[11]).destroy()
        atoms[2].unbond(atoms[11])
        self.assertEqual(atoms[2].bonded_atoms(), set())
        self.assertEqual(atoms[11].bonded_atoms(), set())
        self.assertEqual(atoms[20].bonded_atoms(), set())

        # Put atoms in model
        model = Model(*atoms)
        for atom in atoms:
            self.assertIs(atom.model(), model)
        self.assertEqual(model.atoms(), set(atoms))
        self.assertIn(atoms[-1], model)
        model.remove_atom(atoms[-1])
        self.assertNotIn(atoms[-1], model)
        self.assertIsNone(atoms[-1].model())
        self.assertEqual(model.atoms(), set(atoms[:-1]))
        model.add_atom(atoms[-1])
        self.assertIn(atoms[-1], model)
        self.assertIs(atoms[-1].model(), model)
        self.assertEqual(model.atoms(element="N"), set(atoms[::9]))
        self.assertEqual(len(model.atoms(hydrogen=False)), 25)
        self.assertEqual(model.atoms(element="fe"), set([atoms[-3]] + [atoms[-1]]))
        self.assertEqual(model.atoms(name="CA"), set(atoms[1::9]))
        self.assertEqual(model.atoms(metal=False), set(atoms) - set(atoms[-3::2]))
        self.assertIs(model.atom(1), atoms[0])
        self.assertEqual(model.atom(name="CA").name(), "CA")
        pairs = []
        for pair in model.pairwise_atoms():
            pairs.append(pair)
        self.assertEqual(len(pairs), 351)
        self.assertAlmostEqual(model.mass(), 479, delta=0.5)
        self.assertEqual(model.charge(), 0)
        atoms[13].charge(-3)
        self.assertEqual(model.charge(), -1)
        atoms[13].charge(-2)
        self.assertEqual(model.charge(), 0)
        self.assertEqual(
         model.formula(),
         {"C": 11, "N": 3, "O": 6, "Fe": 2, "P": 1, "S": 2, "H": 2}
        )
        self.assertAlmostEqual(model.center_of_mass()[0], 0.9249, delta=0.005)
        self.assertAlmostEqual(model.center_of_mass()[1], 1.2479, delta=0.005)
        self.assertAlmostEqual(model.center_of_mass()[2], 1.0993, delta=0.005)
        self.assertAlmostEqual(model.radius_of_gyration(), 1.4412, delta=0.005)

        # The model can be transformed correctly
        model.translate(10, 20, 30)
        self.assertEqual(atoms[0].location(), (10, 20, 30))
        self.assertEqual(atoms[13].location(), (11, 21, 31))
        model.translate(-10, -20, -30)
        self.assertEqual(atoms[0].location(), (0, 0, 0))
        self.assertEqual(atoms[13].location(), (1, 1, 1))
        model.rotate(math.pi / 2, "x"), model.round(8)
        self.assertEqual(atoms[0].location(), (0, 0, 0))
        self.assertEqual(atoms[13].location(), (1, -1, 1))
        model.rotate(math.pi / 4, "y"), model.round(6)
        self.assertEqual(atoms[0].location(), (0, 0, 0))
        self.assertEqual(atoms[13].location(), (1.414214, -1, 0))
        model.rotate(math.pi / 8, "z"), model.round(6)
        self.assertEqual(atoms[0].location(), (0, 0, 0))
        self.assertEqual(atoms[13].location(), (1.689247, -0.382683, 0))
        model.orient(atoms[-1])
        self.assertEqual(atoms[-1].location(), (0, 0, 0))
        self.assertNotEqual(atoms[0].location(), (0, 0, 0))
        model.orient(atoms[4], atom2=atoms[-1], axis="y"), model.round(8)
        self.assertEqual(atoms[4].location(), (0, 0, 0))
        self.assertEqual((atoms[-1].x(), atoms[-1].z()), (0, 0))
        model.orient(atoms[8], atom2=atoms[9], axis="z"), model.round(8)
        self.assertEqual(atoms[8].location(), (0, 0, 0))
        self.assertEqual((atoms[9].x(), atoms[9].y()), (0, 0))
        model.orient(atoms[19], atom2=atoms[1], axis="x"), model.round(8)
        self.assertEqual(atoms[19].location(), (0, 0, 0))
        self.assertEqual((atoms[1].y(), atoms[1].z()), (0, 0))
        model.orient(atoms[1], atom2=atoms[19], axis="x"), model.round(8)
        self.assertEqual(atoms[1].location(), (0, 0, 0))
        self.assertEqual((atoms[19].y(), atoms[19].z()), (0, 0))
        model.orient(atoms[25], atom2=atoms[19], axis="x", atom3=atoms[5], plane="xy")
        model.round(8)
        self.assertEqual(atoms[25].location(), (0, 0, 0))
        self.assertEqual((atoms[19].y(), atoms[19].z()), (0, 0))
        self.assertEqual(atoms[5].z(), 0)
        model.orient(atoms[5], atom2=atoms[25], axis="x", atom3=atoms[19], plane="yx")
        model.round(8)
        self.assertEqual(atoms[5].location(), (0, 0, 0))
        self.assertEqual((atoms[25].y(), atoms[25].z()), (0, 0))
        self.assertEqual(atoms[19].z(), 0)
        model.orient(atoms[5], atom2=atoms[25], axis="x", atom3=atoms[19], plane="xz")
        model.round(8)
        self.assertEqual(atoms[5].location(), (0, 0, 0))
        self.assertEqual((atoms[25].y(), atoms[25].z()), (0, 0))
        self.assertEqual(atoms[19].y(), 0)
        model.orient(atoms[21], atom2=atoms[2], axis="x", atom3=atoms[13], plane="zx")
        model.round(8)
        self.assertEqual(atoms[21].location(), (0, 0, 0))
        self.assertEqual((atoms[2].y(), atoms[2].z()), (0, 0))
        self.assertEqual(atoms[13].y(), 0)
        model.orient(atoms[21], atom2=atoms[2], axis="y", atom3=atoms[13], plane="xy")
        model.round(8)
        self.assertEqual(atoms[21].location(), (0, 0, 0))
        self.assertEqual((atoms[2].x(), atoms[2].z()), (0, 0))
        self.assertEqual(atoms[13].z(), 0)
        model.orient(atoms[11], atom2=atoms[4], axis="y", atom3=atoms[10], plane="yx")
        model.round(8)
        self.assertEqual(atoms[11].location(), (0, 0, 0))
        self.assertEqual((atoms[4].x(), atoms[4].z()), (0, 0))
        self.assertEqual(atoms[10].z(), 0)
        model.orient(atoms[11], atom2=atoms[4], axis="y", atom3=atoms[10], plane="yz")
        model.round(8)
        self.assertEqual(atoms[11].location(), (0, 0, 0))
        self.assertEqual((atoms[4].x(), atoms[4].z()), (0, 0))
        self.assertEqual(atoms[10].x(), 0)
        model.orient(atoms[1], atom2=atoms[2], axis="y", atom3=atoms[17], plane="zy")
        model.round(8)
        self.assertEqual(atoms[1].location(), (0, 0, 0))
        self.assertEqual((atoms[2].x(), atoms[2].z()), (0, 0))
        self.assertEqual(atoms[17].x(), 0)
        model.orient(atoms[1], atom2=atoms[2], axis="z", atom3=atoms[17], plane="xz")
        model.round(8)
        self.assertEqual(atoms[1].location(), (0, 0, 0))
        self.assertEqual((atoms[2].x(), atoms[2].y()), (0, 0))
        self.assertEqual(atoms[17].y(), 0)
        model.orient(atoms[11], atom2=atoms[0], axis="z", atom3=atoms[16], plane="zx")
        model.round(8)
        self.assertEqual(atoms[11].location(), (0, 0, 0))
        self.assertEqual((atoms[0].x(), atoms[0].y()), (0, 0))
        self.assertEqual(atoms[16].y(), 0)
        model.orient(atoms[11], atom2=atoms[0], axis="z", atom3=atoms[16], plane="yz")
        model.round(8)
        self.assertEqual(atoms[11].location(), (0, 0, 0))
        self.assertEqual((atoms[0].x(), atoms[0].y()), (0, 0))
        self.assertEqual(atoms[16].x(), 0)
        model.orient(atoms[10], atom2=atoms[9], axis="z", atom3=atoms[11], plane="zy")
        model.round(8)
        self.assertEqual(atoms[10].location(), (0, 0, 0))
        self.assertEqual((atoms[9].x(), atoms[9].y()), (0, 0))
        self.assertEqual(atoms[11].x(), 0)
        model.orient(atoms[0], atom2=atoms[1], axis="x", atom3=atoms[4], plane="xy")
        model.round(4)
        self.assertEqual(atoms[13].location(), (1, 1, 1))

        # The atoms work in the model
        self.assertEqual(atoms[0].nearby(0.5), set())
        self.assertEqual(atoms[0].nearby(1), set([atoms[1], atoms[3], atoms[9]]))
        self.assertEqual(atoms[0].nearby(1, name="CA"), set([atoms[1]]))
        self.assertEqual(atoms[0].nearby(1.5), set([
         atoms[1], atoms[3], atoms[9], atoms[4], atoms[10], atoms[12]
        ]))
        self.assertEqual(
         atoms[0].nearby(1.5), set([atoms[1], atoms[3], atoms[4], atoms[12], atoms[9], atoms[10]])
        )
        self.assertEqual(
         atoms[0].nearby(1.5, element="C"), set([atoms[1], atoms[4], atoms[10]])
        )
        self.assertEqual(
         atoms[0].nearby(1.5, name="CA"), set([atoms[1], atoms[10]])
        )
        self.assertEqual(model.atoms_in_sphere(0, 0, 0, 0.5), set([atoms[0]]))
        self.assertEqual(
         model.atoms_in_sphere(0, 0, 0, 1),
         set([atoms[0], atoms[1], atoms[3], atoms[9]])
        )
        self.assertEqual(model.atoms_in_sphere(2, 2, 2, 0.5), set([atoms[-1]]))
        self.assertEqual(model.atoms_in_sphere(2, 2, 2, 5), set(atoms))
        self.assertEqual(model.atoms_in_sphere(2, 2, 2, 5, hydrogen=False), set(atoms) - set(atoms[16:18]))

        # Grid
        self.assertEqual(model.grid(), tuple([(x, y, z)
         for x in range(3) for y in range(3) for z in range(3)]))
        model.translate(-0.1, -0.1, -0.1), model.round(6)
        self.assertEqual(model.grid(), tuple([(x, y, z)
         for x in range(-1, 3) for y in range(-1, 3) for z in range(-1, 3)]))
        model.translate(-0.8, -0.8, -0.8), model.round(6)
        self.assertEqual(model.grid(), tuple([(x, y, z)
         for x in range(-1, 3) for y in range(-1, 3) for z in range(-1, 3)]))
        model.translate(1, 1, 1), model.round(6)
        self.assertEqual(model.grid(), tuple([(x, y, z)
         for x in range(4) for y in range(4) for z in range(4)]))
        model.translate(-0.1, -0.1, -0.1), model.round(6)
        self.assertEqual(model.grid(size=0.5), tuple([(x, y, z)
         for x in [0, 0.5, 1, 1.5, 2]
         for y in [0, 0.5, 1, 1.5, 2]
         for z in [0, 0.5, 1, 1.5, 2]]))
        self.assertEqual(model.grid(size=0.5, margin=1.5), tuple([(x, y, z)
         for x in [-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5]
         for y in [-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5]
         for z in [-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5]]))
        atoms[0].translate(-0.1, -0.5, -0.6), model.round(6)
        atoms[-1].translate(0.1, 0.5, 0.6), model.round(6)
        self.assertEqual(model.grid(size=0.25), tuple([(x, y, z)
         for x in [-0.25, 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25]
         for y in [-0.5, -0.25, 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5]
         for z in [-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75]]))
        atoms[0].translate(0.1, 0.5, 0.6), model.round(6)
        atoms[-1].translate(-0.1, -0.5, -0.6), model.round(6)
        self.assertEqual(model.grid(size=5), tuple([(x, y, z)
         for x in [0, 5]
         for y in [0, 5]
         for z in [0, 5]]))


        # Residues can be made
        res1 = Residue(*atoms[:9], residue_id="A1", name="MET")
        res2 = Residue(*atoms[9:18], residue_id="A2", name="VAL")
        res3 = Residue(*atoms[18:], residue_id="A2A", name="MET")
        for atom in atoms[:9]:
            self.assertIs(atom.residue(), res1)
            self.assertIs(atom.molecule(), res1)
        for atom in atoms[9:18]:
            self.assertIs(atom.residue(), res2)
            self.assertIs(atom.molecule(), res2)
        for atom in atoms[18:27]:
            self.assertIs(atom.residue(), res3)
            self.assertIs(atom.molecule(), res3)
        self.assertIsNone(res1.chain())
        self.assertIsNone(res2.chain())
        self.assertIsNone(res3.chain())
        self.assertIs(res1.model(), model)
        self.assertIs(res2.model(), model)
        self.assertIs(res3.model(), model)
        res1.next(res2)
        res3.previous(res2)
        self.assertIs(res1.next(), res2)
        self.assertIs(res2.next(), res3)
        self.assertIs(res3.previous(), res2)
        self.assertIs(res2.previous(), res1)
        res1.next(None)
        self.assertIsNone(res2.previous())
        res1.next(res2)

        # Chains can be made
        chain = Chain(res1, res2, res3, chain_id="A")
        for atom in atoms:
            self.assertIs(atom.chain(), chain)
            self.assertIs(atom.molecule(), chain)
        self.assertIs(res1.chain(), chain)
        self.assertIs(res2.chain(), chain)
        self.assertIs(res3.chain(), chain)
        self.assertEqual(chain.residues(), (res1, res2, res3))
        self.assertEqual(chain.atoms(), set(atoms))
        self.assertIsNone(chain.residue("A4"))
        self.assertIs(chain.residue(name="VAL"), res2)
        self.assertIs(chain.residue("A2A"), res3)
        self.assertEqual(chain.residues(name="MET"), (res1, res3))
        self.assertEqual(len(chain), chain.length(), 3)
        for index, res in enumerate(chain):
            self.assertIs(res, [res1, res2, res3][index])
        chain.remove_residue(res3)
        self.assertEqual(len(chain), chain.length(), 2)
        self.assertIsNone(res3.chain())
        chain.add_residue(res3)
        self.assertEqual(len(chain), chain.length(), 3)
        self.assertIs(res3.chain(), chain)

        # Molecules can be made
        molatom1 = Atom("N", 1, -4, 0, 100, "NA")
        molatom2 = Atom("FE", 1, -4, 1, 101, "NB")
        molatom1.bond(molatom2)
        mol = Molecule(molatom1, molatom2, molecule_id="A1000", name="XOP")
        model.add_molecule(mol)
        self.assertEqual(mol.atoms(), set([molatom1, molatom2]))
        self.assertIs(molatom1.molecule(), mol)
        self.assertIs(molatom1.residue(), None)
        self.assertIs(molatom1.chain(), None)
        self.assertIs(molatom1.model(), model)
        self.assertIs(molatom2.molecule(), mol)
        self.assertIs(molatom2.residue(), None)
        self.assertIs(molatom2.chain(), None)
        self.assertIs(molatom2.model(), model)
        self.assertEqual(model.molecules(), set([chain, mol]))
        site = mol.site(main_chain=True)
        self.assertEqual(site.residues(), set([res1, res2]))
        self.assertEqual(site.atoms(), set(atoms[:18]))
        self.assertIs(site.ligand(), mol)
        site = mol.site()
        self.assertEqual(site.residues(), set())
        self.assertEqual(len(model.atoms()), 29)
        self.assertEqual(len(model.atoms(het=False)), 27)
        self.assertEqual(len(model.atoms(het=False, hydrogen=False)), 25)

        # Copies can be made
        copy = model.copy()
        self.assertEqual(len(copy.atoms()), 29)
        self.assertEqual(len(copy.atoms() | model.atoms()), 58)
        for atom in atoms:
            self.assertNotIn(atom, copy)
        copy.translate(100, 100, 100)
        self.assertEqual(copy.atom(name="F2").location(), (102, 102, 102))
        self.assertEqual(model.atom(name="F2").location(), (2, 2, 2))
        for atom in copy.atoms():
            self.assertIsNone(atom.atom_id())
