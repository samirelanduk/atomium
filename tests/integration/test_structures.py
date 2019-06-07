import math
from unittest import TestCase
import atomium

class Tests(TestCase):

    def test_structure_processing(self):
        # Create five atoms of a residue
        atom1 = atomium.Atom("N", 0, 0, 0, 1, "N", 0.5, 0.5, [0] * 6)
        atom2 = atomium.Atom("C", 1.5, 0, 0, 2, "CA", 0, 0.4, [0] * 6)
        atom3 = atomium.Atom("C", 1.5, 1.5, 0, 3, "CB", 0, 0.3, [1] * 6)
        atom4 = atomium.Atom("C", 3, 0, 0, 4, "C", 0, 0.2, [0] * 6)
        atom5 = atomium.Atom("O", 3, -1.5, 0, 5, "O", 0, 0.1, [0] * 6)

        # Check basic atom properties
        self.assertEqual(atom2.element, "C")
        self.assertEqual(atom1.location, (0, 0, 0))
        self.assertEqual(tuple(atom1), (0, 0, 0))
        self.assertEqual(atom2.name, "CA")
        self.assertEqual(atom1.charge, 0.5)
        self.assertEqual(atom2.charge, 0)
        self.assertEqual(atom2.bvalue, 0.4)
        self.assertEqual(atom3.anisotropy, [1] * 6)

        # Check can update some properties
        atom1.name = "HG"
        atom1.charge = 200
        atom1.bvalue = 20
        self.assertEqual(atom1.name, "HG")
        self.assertEqual(atom1.charge, 200)
        self.assertEqual(atom1.bvalue, 20)
        atom1.name = "N"
        atom1.charge = 0.5
        atom1.bvalue = 0.5

        # Check atoms are not part of any higher structures
        for atom in (atom1, atom2, atom3, atom4, atom5):
            self.assertIsNone(atom.het)
            self.assertIsNone(atom.chain)
            self.assertIsNone(atom.model)
            self.assertEqual(atom.bonded_atoms, set())

        # Check atoms' calculated properties
        self.assertAlmostEqual(atom5.mass, 16, delta=0.05)
        for atom in (atom1, atom2, atom3, atom4, atom5):
            self.assertFalse(atom.is_metal)
            self.assertFalse(atom.is_backbone) # Not yet
            self.assertFalse(atom.is_side_chain) # Not yet

        # Check atom magic methods
        self.assertEqual(list(atom5), [3, -1.5, 0])
        for a1 in (atom1, atom2, atom3, atom4, atom5):
            for a2 in (atom1, atom2, atom3, atom4, atom5):
                if a1 is a2:
                    self.assertEqual(a1, a2)
                else:
                    self.assertNotEqual(a1, a2)

        # Check atom safe methods
        self.assertEqual(atom1.distance_to(atom2), 1.5)
        self.assertEqual(atom1.distance_to(atom3), 4.5 ** 0.5)
        self.assertEqual(atom2.angle(atom3, atom4), math.pi / 2)
        for atom in (atom1, atom2, atom3, atom4, atom5):
            self.assertEqual(atom.nearby_atoms(5), set()) # Not without model
            self.assertEqual(atom.nearby_hets(5), set()) # Ditto

        # Check atom side effect methods
        atom2.translate(0, 0, 1)
        self.assertEqual(atom2.location, (1.5, 0, 1))
        atom2.translate(0, 0, -1)
        self.assertEqual(atom2.location, (1.5, 0, 0))
        atom2.transform([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])
        self.assertEqual(atom2.location, (-1.5, 0, 0))
        atom2.transform([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])
        self.assertEqual(atom2.location, (1.5, 0, 0))
        atom2.rotate(math.pi / 2, "y")
        self.assertEqual(atom2.location, (0, 0, -1.5))
        atom2.rotate(math.pi * 1.5, "y")
        self.assertEqual(atom2.location, (1.5, 0, 0))
        atom2.move_to(10, 10, 10)
        self.assertEqual(atom2.location, (10, 10, 10))
        atom2.move_to(1.5, 0, 0)

        # Bond atoms
        atom1.bond(atom2)
        atom2.bond(atom3)
        atom2.bond(atom4)
        atom4.bond(atom5)
        self.assertEqual(atom1.bonded_atoms, {atom2})
        self.assertEqual(atom2.bonded_atoms, {atom1, atom3, atom4})
        self.assertEqual(atom3.bonded_atoms, {atom2})
        self.assertEqual(atom4.bonded_atoms, {atom2, atom5})
        self.assertEqual(atom5.bonded_atoms, {atom4})

        # Check can copy atom
        copy = atom2.copy()
        self.assertEqual(copy.element, "C")
        self.assertEqual(copy.location, (1.5, 0, 0))
        self.assertEqual(copy.name, "CA")
        self.assertEqual(copy.id, 2)
        self.assertEqual(copy.charge, 0)
        self.assertEqual(copy.charge, 0)
        self.assertEqual(copy.bvalue, 0.4)
        self.assertEqual(copy.anisotropy, [0] * 6)
        self.assertIsNone(copy.het)
        self.assertIsNone(copy.chain)
        self.assertIsNone(copy.model)
        self.assertEqual(copy.bonded_atoms, set())
        self.assertEqual(atom2, copy)

        # Can copy atom with new ID
        copy = atom2.copy(id=10000)
        self.assertEqual(copy.id, 10000)

        # Create residue
        res1 = atomium.Residue(
         atom1, atom2, atom3, atom4, atom5, id="A5", name="AL"
        )

        # Check basic residue properties
        self.assertEqual(res1.id, "A5")
        self.assertEqual(res1.name, "AL")
        self.assertEqual(res1.code, "X")
        res1.name = "ALA"
        self.assertEqual(res1.name, "ALA")
        self.assertEqual(res1.code, "A")
        self.assertIsNone(res1.next)
        self.assertIsNone(res1.previous)
        self.assertIsNone(res1.chain)

        # Check residue and atoms
        for atom in (atom1, atom2, atom3, atom4, atom5):
            self.assertIs(atom.het, res1)
        self.assertTrue(atom1.is_backbone)
        self.assertTrue(atom3.is_side_chain)
        self.assertEqual(res1.atoms(), {atom1, atom2, atom3, atom4, atom5})
        self.assertEqual(res1.atoms(element="C"), {atom2, atom3, atom4})
        self.assertEqual(res1.atoms(name="O"), {atom5})
        self.assertEqual(res1.atoms(is_backbone=True), {atom1, atom2, atom4, atom5})
        self.assertEqual(res1.atoms(mass__gt=13), {atom1, atom5})
        self.assertEqual(res1.atoms(name__regex="N|O"), {atom1, atom5})

        # Check residue is container
        self.assertIn(atom1, res1)
        self.assertIn(copy, res1)

        # Check residue calculated properties
        self.assertEqual(res1.code, "A")
        self.assertEqual(res1.full_name, "alanine")
        self.assertIsNone(res1.model)
        self.assertAlmostEqual(res1.mass, 66, delta=0.05)
        self.assertEqual(res1.charge, 0.5)
        self.assertEqual(res1.formula, {"C": 3, "O": 1, "N": 1})
        self.assertAlmostEqual(res1.center_of_mass[0], 1.818, delta=0.001)
        self.assertAlmostEqual(res1.center_of_mass[1], -0.091, delta=0.001)
        self.assertEqual(res1.center_of_mass[2], 0)
        self.assertAlmostEqual(res1.radius_of_gyration, 1.473, delta=0.001)

        # Check residue safe methods
        self.assertEqual(len(tuple(res1.pairwise_atoms())), 10)
        self.assertEqual(res1.nearby_hets(10), set())
        self.assertEqual(res1.nearby_atoms(10), set())
        self.assertEqual(res1.nearby_chains(10), set())
        self.assertEqual(
         tuple(res1.create_grid(size=3)),
         ((0, -3, 0), (0, 0, 0), (0, 3, 0), (3, -3, 0), (3, 0, 0), (3, 3, 0))
        )
        self.assertEqual(res1.atoms_in_sphere((1.5, 0, 0), 1.5), {atom2, atom1, atom3, atom4})
        self.assertEqual(res1.atoms_in_sphere((1.5, 0, 0), 1.5, element="C"), {atom2, atom3, atom4})
        res1.check_ids()
        self.assertFalse(res1.helix)
        self.assertFalse(res1.strand)

        # Check residue side effect methods
        res1.translate(0, 0, 1)
        self.assertEqual(atom2.location, (1.5, 0, 1))
        res1.translate(0, 0, -1)
        self.assertEqual(atom2.location, (1.5, 0, 0))
        res1.transform([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])
        self.assertEqual(atom2.location, (-1.5, 0, 0))
        res1.transform([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])
        self.assertEqual(atom2.location, (1.5, 0, 0))
        res1.rotate(math.pi / 2, "y")
        self.assertEqual(atom2.location, (0, 0, -1.5))
        res1.rotate(math.pi * 1.5, "y")
        self.assertEqual(atom2.location, (1.5, 0, 0))

        # Can make copy of residue
        res_copy = res1.copy()
        self.assertEqual(res1, res_copy)
        self.assertEqual(res_copy.id, "A5")
        self.assertEqual(res_copy.name, "ALA")
        self.assertEqual(res1.pairing_with(res_copy), {
         atom1: res_copy.atom(1), atom2: res_copy.atom(2),
         atom3: res_copy.atom(3), atom4: res_copy.atom(4),
         atom5: res_copy.atom(5)
        })
        self.assertEqual(len(res1.atoms() | res_copy.atoms()), 10)
        self.assertEqual(res1.rmsd_with(res_copy), 0)
        res_copy.atom(1).translate(1)
        self.assertAlmostEqual(res1.rmsd_with(res_copy), 0.4, delta=0.001)

        # Can make copy of residue with new IDs
        res_copy = res1.copy(id="C5", atom_ids=lambda i: i * 100)
        self.assertEqual(res_copy.id, "C5")
        self.assertEqual({a.id for a in res_copy.atoms()}, {100, 200, 300, 400, 500})

        # Make more residues
        atom6 = atomium.Atom("N", 4.5, 0, 0, 6, "N", 0, 0.5, [0] * 6)
        atom7 = atomium.Atom("C", 6, 0, 0, 7, "CA", 0, 0.5, [0] * 6)
        atom8 = atomium.Atom("C", 6, -1.5, 0, 8, "CB", 0, 0.5, [0] * 6)
        atom9 = atomium.Atom("S", 6, -3, 0, 9, "S", 0, 0.5, [0] * 6)
        atom10 = atomium.Atom("C", 7.5, 0, 0, 10, "C", 0, 0.5, [0] * 6)
        atom11 = atomium.Atom("O", 7.5, 1.5, 0, 11, "O", 0, 0.5, [0] * 6)
        atom12 = atomium.Atom("N", 9, 0, 0, 12, "CA", 0, 0.5, [0] * 6)
        atom13 = atomium.Atom("C", 10.5, 0, 0, 13, "CB", 0, 0.5, [0] * 6)
        atom14 = atomium.Atom("C", 10.5, 1.5, 0, 14, "OG", 0, 0.5, [0] * 6)
        atom15 = atomium.Atom("O", 10.5, 3, 0, 15, "C", 0, 0.5, [0] * 6)
        atom16 = atomium.Atom("C", 12, 0, 0, 16, "C", 0, 0.5, [0] * 6)
        atom17 = atomium.Atom("O", 13.5, 1.25, 0, 17, "OX1", 0, 0.5, [0] * 6)
        atom18 = atomium.Atom("O", 13.5, -1.25, 0, 18, "OX2", 0, 0.5, [0] * 6)
        atom6.bond(atom7)
        atom6.bond(atom4)
        atom7.bond(atom8)
        atom7.bond(atom10)
        atom8.bond(atom9)
        atom10.bond(atom11)
        atom10.bond(atom12)
        atom12.bond(atom13)
        atom13.bond(atom14)
        atom13.bond(atom16)
        atom14.bond(atom15)
        atom16.bond(atom17)
        atom16.bond(atom18)
        res2 = atomium.Residue(
         atom6, atom7, atom8, atom9, atom10, atom11, id="A5A", name="CYS"
        )
        res3 = atomium.Residue(
         atom12, atom13, atom14, atom15, atom16, atom17, atom18, id="A6", name="SER"
        )

        # Connect residues
        res1.next = res2
        res3.previous = res2
        self.assertIs(res1.next, res2)
        self.assertIs(res2.next, res3)
        self.assertIs(res3.previous, res2)
        self.assertIs(res2.previous, res1)

        # Create chain
        chain1 = atomium.Chain(
         res1, res2, res3, id="A", sequence="MACSD", helices=((res1, res2),), strands=((res3,),)
        )
        self.assertEqual(chain1.id, "A")
        self.assertEqual(chain1.internal_id, "A")
        self.assertIsNone(chain1.name)
        self.assertIsNone(chain1.model)

        # Chain properties
        self.assertEqual(chain1.sequence, "MACSD")
        chain1.sequence = "MACSDA"
        self.assertEqual(chain1.sequence, "MACSDA")
        self.assertEqual(chain1.present_sequence, "ACS")
        self.assertEqual(chain1.helices[0], (res1, res2))
        self.assertEqual(chain1.strands[0], (res3,))

        # Check chain residues and atoms
        self.assertEqual(chain1.residues(), (res1, res2, res3))
        self.assertEqual(chain1.residues(mass__gt=80), (res2, res3))
        self.assertEqual(chain1.residues(mass__lt=98.1), (res1, res3))
        self.assertEqual(chain1.atoms(), {
         atom1, atom2, atom3, atom4, atom5, atom6, atom7, atom8, atom9, atom10,
         atom11, atom12, atom13, atom14, atom15, atom16, atom17, atom18
        })
        self.assertTrue(res1.helix)
        self.assertFalse(res1.strand)
        self.assertTrue(res2.helix)
        self.assertFalse(res2.strand)
        self.assertFalse(res3.helix)
        self.assertTrue(res3.strand)

        # Check chain magic methods
        self.assertEqual(chain1.length, 3)
        self.assertIs(chain1[0], res1)
        for res in chain1: self.assertIn(res, (res1, res2, res3))
        self.assertIn(res1, chain1)
        self.assertIn(atom10, chain1)

        # Check chain ligands and atoms
        self.assertEqual(chain1.ligands(), set())
        self.assertEqual(
         chain1.atoms(element__regex="O|S"),
         {atom5, atom9, atom11, atom15, atom17, atom18}
        )
        self.assertEqual(
         chain1.atoms(het__name="CYS"),
         {atom6, atom7, atom8, atom9, atom10, atom11}
        )

        # Make copy of chain
        chain2 = chain1.copy(
         id="B",
         residue_ids=lambda i: i.replace("A", "B"),
         atom_ids=lambda i: i * 100
        )
        self.assertEqual(chain1, chain2)
        self.assertEqual(chain2.id, "B")
        self.assertEqual(chain2.internal_id, "B")
        self.assertEqual([res.id for res in chain2], ["B5", "B5B", "B6"])
        self.assertEqual({a.id for a in chain2.atoms()}, set([x * 100 for x in range(1, 19)]))
        self.assertIs(chain2[0].next, chain2[1])
        self.assertIs(chain2[2].previous, chain2[1])
        self.assertEqual(chain2.helices, ((chain2[0], chain2[1]),))
        self.assertEqual(chain2.strands, ((chain2[2],),))

        # Move chain into place
        chain2.rotate(math.pi, "x")
        chain2.rotate(math.pi, "y")
        chain2.translate(12, -10.5)
        self.assertEqual(chain1.rmsd_with(chain2), 0)

        # Make ligand
        copper_atom = atomium.Atom("Cu", 6, -5.25, 2, 100, "Cu", 2, 0, [0] * 6)
        copper = atomium.Ligand(
         copper_atom, id="A100", internal_id="M", name="CU", chain=chain1, full_name="copper"
        )

        # Check ligand properties
        self.assertEqual(copper.id, "A100")
        self.assertEqual(copper.name, "CU")
        self.assertEqual(copper.full_name, "copper")
        copper.full_name = None
        self.assertEqual(copper.full_name, "CU")
        self.assertEqual(copper.internal_id, "M")
        self.assertIs(copper.chain, chain1)
        self.assertIsNone(copper.model)
        self.assertEqual(copper.atom(), copper_atom)
        self.assertIn(copper_atom, copper)
        self.assertFalse(copper.is_water)

        # Can make copy of ligand
        cu_copy = copper.copy()
        self.assertEqual(copper, cu_copy)
        self.assertEqual(cu_copy.id, "A100")
        self.assertEqual(cu_copy.name, "CU")
        self.assertEqual(len(cu_copy.atoms() | cu_copy.atoms()), 1)
        self.assertFalse(cu_copy.is_water)

        # Can make copy of ligand with new IDs
        cu_copy = copper.copy(id="C100", atom_ids=lambda i: i * 100)
        self.assertEqual(cu_copy.id, "C100")
        self.assertEqual(cu_copy.atom().id, 10000)

        # Create waters
        hoh1 = atomium.Ligand(
         atomium.Atom("O", 3, -3, 3, 500, "O", 0, 0, [0] * 6),
         id="A1000", name="HOH", water=True
        )
        self.assertTrue(hoh1.is_water)
        hoh2 = atomium.Ligand(
         atomium.Atom("O", 3, -9, -3, 500, "O", 0, 0, [0] * 6),
         id="B1000", name="HOH", water=True
        )
        self.assertTrue(hoh2.is_water)

        # Create model
        model = atomium.Model(chain1, chain2, copper, hoh1, hoh2)
        
        # Model properties
        self.assertIsNone(model.file)
        self.assertEqual(model.chains(), {chain1, chain2})
        self.assertEqual(model.ligands(), {copper})
        self.assertEqual(chain1.ligands(), {copper})
        self.assertEqual(model.waters(), {hoh1, hoh2})
        self.assertEqual(model.molecules(), {chain1, chain2, copper, hoh1, hoh2})
        self.assertEqual(model.residues(), set(chain1.residues() + chain2.residues()))
        self.assertEqual(model.residues(name="ALA"), {chain1[0], chain2[0]})
        self.assertEqual(model.ligand(), copper)
        self.assertEqual(model.atom(1), atom1)
        self.assertEqual(model.atom(name="N", het__name="ALA", chain__id="A"), atom1)

        # Everything points upwards correctly
        self.assertIs(atom1.model, model)
        self.assertIs(res1.model, model)
        self.assertIs(chain1.model, model)
        self.assertIs(copper.model, model)

        # Now that atoms are in a model, find nearby things
        self.assertEqual(atom2.nearby_atoms(1.5), {atom1, atom3, atom4})
        self.assertEqual(atom4.nearby_atoms(1.5), {atom2, atom5, atom6})
        self.assertEqual(atom4.nearby_atoms(1.5, het__name="CYS"), {atom6})
        self.assertEqual(atom4.nearby_hets(1.5), {res2})
        self.assertEqual(atom4.nearby_hets(9), {res2, res3, chain2[1], copper, hoh1})
        self.assertEqual(atom4.nearby_hets(9, ligands=False), {res2, res3, chain2[1]})
        self.assertEqual(atom4.nearby_hets(9, residues=False), {copper, hoh1})
        self.assertEqual(atom4.nearby_hets(9, residues=False, het__is_water=False), {copper})
        self.assertEqual(atom4.nearby_chains(9), {chain2})
        self.assertEqual(atom4.nearby_chains(9, chain__id="A"), set())
        self.assertEqual(res2.nearby_hets(3), {res1, res3})
        self.assertEqual(res2.nearby_hets(6), {res1, res3, hoh1, copper, chain2[1]})
        self.assertEqual(res2.nearby_hets(6, ligands=False), {res1, res3, chain2[1]})
        self.assertEqual(copper.nearby_chains(5), {chain2})
        self.assertEqual(chain2.nearby_chains(5), {chain1})

        # Dehydrate model
        model.dehydrate()
        self.assertEqual(model.waters(), set())
        self.assertEqual(model.ligands(), {copper})
        self.assertEqual(model.chains(), {chain1, chain2})
