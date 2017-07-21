from base import IntegratedTest
from atomium.structures import Model, Atom, Residue, Chain
import atomium

class StructureTests(IntegratedTest):

    def test_can_manipualte_model(self):
        # Create a model
        model = Model()
        self.assertEqual(model.atoms(), set())
        self.assertEqual(model.mass(), 0)

        # Create some atoms
        atom1 = Atom("N", 12.0, 11.5, 1.5, atom_id=101)
        atom2 = Atom("C", 12.5, 10, 2, atom_id=102)

        # The atoms can work out their distance to each other
        self.assertAlmostEqual(atom1.distance_to(atom2), 1.658312, delta=0.0005)

        # Give the model some atoms
        model.add_atom(atom1)
        model.add_atom(atom2)

        # The model can give a full account of itself
        self.assertEqual(model.atoms(), set([atom1, atom2]))
        self.assertAlmostEqual(model.mass(), 26, delta=0.1)
        self.assertIn(atom1, model)
        self.assertIn(atom2, model)

        # Atom retrieval can pick by element and ID
        self.assertIs(model.atom(element="N"), atom1)
        self.assertIs(model.atom(element="C"), atom2)
        self.assertIs(model.atom(atom_id=101), atom1)
        self.assertIs(model.atom(atom_id=102), atom2)

        # Atoms can be bonded to each other
        self.assertEqual(atom1.bonded_atoms(), set())
        self.assertEqual(atom2.bonded_atoms(), set())
        atom1.bond(atom2)
        self.assertEqual(atom1.bonded_atoms(), set([atom2]))
        self.assertEqual(atom2.bonded_atoms(), set([atom1]))
        atom2.unbond(atom1)
        self.assertEqual(atom1.bonded_atoms(), set())
        self.assertEqual(atom2.bonded_atoms(), set())

        # The structure has a center of mass and radius of gyration
        self.assertAlmostEqual(model.center_of_mass()[0], 12.23, delta=0.005)
        self.assertAlmostEqual(model.center_of_mass()[1], 10.81, delta=0.005)
        self.assertAlmostEqual(model.center_of_mass()[2], 1.73, delta=0.005)
        self.assertAlmostEqual(model.radius_of_gyration(), 0.83, delta=0.005)

        # The model is fully transformable
        model.translate(-12, -11.5, -1.5)
        self.assertEqual((atom1.x(), atom1.y(), atom1.z()), (0, 0, 0))
        self.assertEqual((atom2.x(), atom2.y(), atom2.z()), (0.5, -1.5, 0.5))
        model.rotate("x", 180)
        self.assertEqual((atom1.x(), atom1.y(), atom1.z()), (0, 0, 0))

        # Atoms can be removed
        model.remove_atom(atom1)
        self.assertIn(atom2, model)
        self.assertNotIn(atom1, model)

        # Model can be saved and reloaded
        model.save("itests/files/model.xyz", description="Some atoms")
        new = atomium.xyz_from_file("itests/files/model.xyz")
        self.assertEqual(new.comment(), "Some atoms")
        self.assertEqual(len(new.model().atoms()), 1)
        self.assertEqual(new.model().atom().x(), 0.5)
        self.assertEqual(new.model().atom().y(), 1.5)
        self.assertEqual(new.model().atom().z(), -0.5)


    def test_can_manipualte_advanced_model(self):
        # Create 18 atoms
        atom1 = Atom("C", 0, 0, 0, atom_id=1, name="CA", charge=0.8)
        atom2 = Atom("C", 1, 0, 0.1, atom_id=2, name="CB")
        atom3 = Atom("N", 0.5, 1, 0, atom_id=3, name="N1")
        atom4 = Atom("C", 2, 0, 0, atom_id=4, name="CA")
        atom5 = Atom("C", 3, 0, -0.2, atom_id=5, name="CB")
        atom6 = Atom("N", 2.5, 1, 0, atom_id=6, name="N1")
        atom7 = Atom("C", 4, 0, 0, atom_id=7, name="CA")
        atom8 = Atom("C", 5, 0, 0, atom_id=8, name="CB")
        atom9 = Atom("N", 4.5, 1, 0, atom_id=9, name="N1")
        atom10 = Atom("C", 0, 0, 10, atom_id=10, name="CA")
        atom11 = Atom("C", 1, 0, 10.1, atom_id=11, name="CB")
        atom12 = Atom("N", 0.5, 1, 10, atom_id=12, name="N1")
        atom13 = Atom("C", 2, 0, 10, atom_id=13, name="CA")
        atom14 = Atom("C", 3, 0, -10.2, atom_id=14, name="CB")
        atom15 = Atom("N", 2.5, 1, 10, atom_id=15, name="N1")
        atom16 = Atom("C", 4, 0, 10, atom_id=16, name="CA")
        atom17 = Atom("C", 5, 0, 10, atom_id=17, name="CB")
        atom18 = Atom("N", 4.5, 1, 10, atom_id=18, name="N1", charge=-0.6)

        # Bond atoms together
        atom1.bond(atom2), atom1.bond(atom3)
        atom4.bond(atom3), atom4.bond(atom5), atom4.bond(atom6)
        atom7.bond(atom6), atom7.bond(atom8), atom7.bond(atom9)
        atom10.bond(atom11), atom10.bond(atom12)
        atom13.bond(atom12), atom13.bond(atom14), atom13.bond(atom15)
        atom16.bond(atom15), atom16.bond(atom17), atom16.bond(atom18)

        # Make residues
        residue1a = Residue(atom1, atom2, atom3, residue_id="A1", name="GLY")
        residue2a = Residue(atom4, atom5, atom6, residue_id="A2", name="HIS")
        residue3a = Residue(atom7, atom8, atom9, residue_id="A3", name="GLY")
        residue1b = Residue(atom10, atom11, atom12, residue_id="B1", name="ASP")
        residue2b = Residue(atom13, atom14, atom15, residue_id="B2", name="GLU")
        residue3b = Residue(atom16, atom17, atom18, residue_id="B3", name="CYS")

        # Check residues
        self.assertEqual(residue1a.atoms(), set([atom1, atom2, atom3]))
        self.assertEqual(residue1a.residue_id(), "A1")
        self.assertEqual(residue1a.molecule_id(), "A1")
        self.assertEqual(residue1a.name(), "GLY")
        self.assertEqual(residue1a.charge(), 0.8)
        residue1a.name("LYS")
        self.assertEqual(residue1a.name(), "LYS")
        residue1a.name("GLY")
        self.assertEqual(residue2a.charge(), 0)
        self.assertEqual(residue3b.charge(), -0.6)
        self.assertIs(atom1.residue(), residue1a)
        self.assertIs(atom14.residue(), residue2b)
        self.assertIs(atom1.molecule(), residue1a)
        self.assertIs(atom14.molecule(), residue2b)

        # Make chains
        chaina = Chain(residue1a, residue2a, residue3a, chain_id="A")
        chainb = Chain(residue1b, residue2b, residue3b, chain_id="B")

        # Check chains
        self.assertEqual(chaina.residues(), (residue1a, residue2a, residue3a))
        self.assertEqual(chainb.residues(), (residue1b, residue2b, residue3b))
        self.assertEqual(chaina.residues(name="GLY"), (residue1a, residue3a))
        self.assertEqual(chaina.residues(residue_id="A2"), (residue2a,))
        self.assertIs(residue1b.chain(), chainb)
        self.assertIs(atom15.chain(), chainb)
        self.assertIs(atom15.molecule(), chainb)

        # Make model
        model = Model()
        model.add_chain(chaina)
        model.add_chain(chainb)

        # Check model
        self.assertEqual(len(model.atoms()), 18)
        self.assertEqual(model.atom(element="C", name="CA").name(), "CA")
        self.assertEqual(model.residues(), set(
         [residue1a, residue1b, residue1c, residue2a, residue2b, residue2c]
        ))
        self.assertEqual(model.chains(), set([chaina, chainb]))
        self.assertIs(atom1.model(), model)
        self.assertIs(residue1a.model(), model)
        self.assertIs(chaina.model(), model)
        model.remove_chain(chainb)
        self.assertEqual(len(model.atoms()), 9)
        self.assertEqual(model.residues(), set(
         [residue1a, residue1b, residue1c]
        ))
        self.assertEqual(model.chains(), set([chaina]))
        chaina.remove_residue(chain1a)
        self.assertEqual(len(model.atoms()), 6)
        self.assertEqual(model.residues(), set([residue1b, residue1c]))

        # Make model from atoms
        model = Model(
         atom1, atom2, atom3, atom4, atom5, atom6, atom7, atom8, atom9,
         atom10, atom11, atom12, atom13, atom14, atom15, atom16, atom17, atom18
        )
        self.assertEqual(len(model.atoms()), 18)
        self.assertEqual(model.residues(), set(
         [residue1a, residue1b, residue1c, residue2a, residue2b, residue2c]
        ))
        self.assertEqual(model.chains(), set([chaina, chainb]))
