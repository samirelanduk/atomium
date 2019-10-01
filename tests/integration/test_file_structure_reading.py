from datetime import date
import math
import atomium
from unittest import TestCase

class DeNovoStructureTests(TestCase):

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
        self.assertEqual(atom1.covalent_radius, 0.71)
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
        self.assertEqual(chain1.atoms(het__name="ALA"), {
         atom1, atom2, atom3, atom4, atom5
        })
        self.assertEqual(chain1.atoms(het__name="CYS"), {
         atom6, atom7, atom8, atom9, atom10, atom11
        })
        self.assertEqual(chain1.atoms(het__name__regex="CYS|ALA"), {
         atom1, atom2, atom3, atom4, atom5, atom6, atom7, atom8, atom9, atom10, atom11
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



class FileReadingTests(TestCase):

    def test_1lol(self):
        for e in ["cif", "mmtf", "pdb"]:
            f = atomium.open("tests/integration/files/1lol." + e)
            self.assertEqual(f.filetype, e)
            self.assertEqual(f.code, "1LOL")
            if e == "pdb":
                self.assertEqual(
                 f.title, "CRYSTAL STRUCTURE OF OROTIDINE MONOPHOSPHATE DECARBOXYLASE COMPLEX WITH XMP"
                )
            else:
                self.assertEqual(
                 f.title, "Crystal structure of orotidine monophosphate decarboxylase complex with XMP"
                )
            self.assertEqual(f.deposition_date, date(2002, 5, 6))
            self.assertEqual(f.classification, None if e == "mmtf" else "LYASE")
            self.assertEqual(f.keywords, [] if e == "mmtf" else ["TIM BARREL", "LYASE"] if e == "pdb" else ["TIM barrel", "LYASE"])
            self.assertEqual(f.authors, [] if e == "mmtf" else ["N.WU", "E.F.PAI"] if e == "pdb" else ["Wu, N.", "Pai, E.F."])
            self.assertEqual(f.technique, "X-RAY DIFFRACTION")
            missing_residues = [{"id": id, "name": name} for id, name in zip([
             "A.1", "A.2", "A.3", "A.4", "A.5", "A.6", "A.7", "A.8", "A.9", "A.10",
             "A.182", "A.183", "A.184", "A.185", "A.186", "A.187", "A.188", "A.189",
             "A.223", "A.224", "A.225", "A.226", "A.227", "A.228", "A.229", "B.1001",
             "B.1002", "B.1003", "B.1004", "B.1005", "B.1006", "B.1007", "B.1008",
             "B.1009", "B.1010", "B.1182", "B.1183", "B.1184", "B.1185", "B.1186"
            ], [
             "LEU", "ARG", "SER", "ARG", "ARG", "VAL", "ASP", "VAL", "MET", "ASP",
             "VAL", "GLY", "ALA", "GLN", "GLY", "GLY", "ASP", "PRO", "LYS", "ASP",
             "LEU", "LEU", "ILE", "PRO", "GLU", "LEU", "ARG", "SER", "ARG", "ARG",
             "VAL", "ASP", "VAL", "MET", "ASP", "VAL", "GLY", "ALA", "GLN", "GLY"
            ])]
            if e == "pdb":
                self.assertEqual(f.source_organism, "METHANOTHERMOBACTER THERMAUTOTROPHICUS STR. DELTA H")
                self.assertEqual(f.expression_system, "ESCHERICHIA COLI")
                self.assertEqual(f.missing_residues, missing_residues)
            else:
                self.assertEqual(
                 f.source_organism,
                 None if e == "mmtf" else "Methanothermobacter thermautotrophicus str. Delta H"
                )
                self.assertEqual(
                 f.expression_system, None if e == "mmtf" else "Escherichia coli"
                )
                self.assertEqual(f.missing_residues, [] if e == "mmtf" else missing_residues)
            self.assertEqual(f.resolution, 1.9)
            self.assertEqual(f.rvalue, 0.193)
            self.assertEqual(f.rfree, 0.229)
            self.assertEqual(f.assemblies, [{
             "id": 1,
             "software": None if e == "mmtf" else "PISA",
             "delta_energy": None if e == "mmtf" else -31.0,
             "buried_surface_area": None if e == "mmtf" else 5230,
             "surface_area": None if e == "mmtf" else 16550,
             "transformations": [{
              "chains": ["A", "B"] if e == "pdb" else ["A", "B", "C", "D", "E", "F", "G", "H"],
              "matrix": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
              "vector": [0.0, 0.0, 0.0]
             }]
            }])

            self.assertEqual(len(f.models), 1)
            model = f.model
            self.assertEqual(len(model.chains()), 2)
            self.assertIsInstance(model.chains(), set)
            self.assertEqual(len(model.ligands()), 4)
            self.assertIsInstance(model.ligands(), set)
            self.assertEqual(len(model.waters()), 180)
            self.assertIsInstance(model.waters(), set)
            self.assertEqual(len(model.molecules()), 186)
            self.assertIsInstance(model.molecules(), set)
            self.assertEqual(len(model.residues()), 418)
            self.assertIsInstance(model.residues(), set)
            self.assertEqual(len(model.atoms()), 3431)
            self.assertIsInstance(model.atoms(), set)
            self.assertEqual(len(model.chains(length__gt=200)), 2)
            self.assertEqual(len(model.chains(length__gt=210)), 1)
            self.assertEqual(len(model.ligands(name="XMP")), 2)
            self.assertEqual(len(model.residues(name="VAL")), 28)
            self.assertEqual(len(model.residues(name="CYS")), 6)
            self.assertEqual(len(model.residues(name__regex="CYS|VAL")), 34)
            self.assertAlmostEqual(
             model.mass, 46018.5, delta=0.005
            )

            chaina = model.chain("A")
            chainb = model.chain(id="B")
            self.assertIs(chaina.model, model)
            self.assertIs(chainb.model, model)
            self.assertEqual(chaina.id, "A")
            self.assertEqual(chaina.length, 204)
            self.assertEqual(chainb.length, 214)
            self.assertTrue(chaina.sequence.startswith("LRSRRVDVMDVMNRLILAMDL"))
            self.assertTrue(chaina.sequence.endswith("LADNPAAAAAGIIESIKDLLIPE"))
            self.assertTrue(chainb.sequence.startswith("LRSRRVDVMDVMNRLILAMDL"))
            self.assertTrue(chainb.sequence.endswith("LADNPAAAAAGIIESIKDLLIPE"))
            for res in chaina: self.assertIn(res, chaina)
            self.assertEqual(len(chaina.residues()), 204)
            self.assertIsInstance(chaina.residues(), tuple)
            self.assertEqual(len(chaina.ligands()), 2)
            self.assertIsInstance(chaina.ligands(), set)
            self.assertEqual(len(chaina.atoms()), 1557)
            self.assertIsInstance(chaina.atoms(), set)
            self.assertEqual(len(chainb.atoms()), 1634)
            self.assertIsInstance(chainb.atoms(), set)
            res = chaina.residue("A.13")
            self.assertTrue(res.helix)
            self.assertFalse(res.strand)
            res = chaina.residue("A.15")
            self.assertFalse(res.helix)
            self.assertTrue(res.strand)
            self.assertIs(res.chain, chaina)
            self.assertIs(res.model, model)
            self.assertEqual(res.name, "LEU")
            self.assertEqual(res.code, "L")
            self.assertEqual(res.full_name, "leucine")
            self.assertEqual(len(res.atoms()), 8)
            self.assertIsInstance(chaina.atoms(), set)
            self.assertEqual(len(res.atoms(element="C")), 6)
            self.assertEqual(len(res.atoms(element__regex="C|O")), 7)
            self.assertEqual(len(res.atoms(name__regex="^CD")), 2)
            self.assertIs(chaina[0], chaina.residue("A.11"))
            self.assertIs(res.next, chaina[5])
            self.assertIn(chaina.residue(name="GLN"), [chaina.residue("A.136"), chaina.residue("A.173")])

            lig = model.ligand(name="XMP")
            self.assertIs(lig.model, model)
            self.assertEqual(len(lig.atoms()), 24)
            self.assertEqual(lig.formula, {"C": 10, "O": 9, "N": 4, "P": 1})
            self.assertEqual(lig.full_name, "XANTHOSINE-5'-MONOPHOSPHATE")
            lig = model.ligand("A.5001")
            self.assertIs(lig.model, model)
            self.assertIs(lig.chain, chaina)
            self.assertEqual(len(lig.atoms()), 6)
            self.assertEqual(lig.mass, 80.0416)
            pairs = list(lig.pairwise_atoms())
            self.assertEqual(len(pairs), 15)
            for pair in pairs:
                pair = list(pair)
                self.assertTrue(0 < pair[0].distance_to(pair[1]), 5)
            hoh = model.water("A.3005")
            self.assertEqual(hoh.name, "HOH")
            self.assertIs(lig.model, model)
            self.assertIs(lig.chain, chaina)
            lig1, lig2 = model.ligands(name="XMP")
            self.assertAlmostEqual(lig1.rmsd_with(lig2), 0.133, delta=0.001)
            self.assertAlmostEqual(lig2.rmsd_with(lig1), 0.133, delta=0.001)

            atom = model.atom(934)
            self.assertEqual(atom.anisotropy, [0, 0, 0, 0, 0, 0])
            self.assertEqual(atom.element, "C")
            self.assertEqual(atom.name, "CA")
            self.assertEqual(atom.location, (4.534, 53.864, 43.326))
            self.assertEqual(atom.bvalue, 17.14)
            self.assertEqual(atom.charge, 0)
            self.assertAlmostEqual(atom.mass, 12, delta=0.1)
            self.assertIs(atom.chain, chaina)
            self.assertIs(atom.model, model)
            self.assertIs(atom.het, model.residue("A.131"))

            self.assertEqual(model.molecule("A"), chaina)
            self.assertEqual(model.molecule("A.5001"), lig)
            self.assertEqual(model.molecule("A.3005"), hoh)
            self.assertEqual(len(model.molecules(mass__gt=18)), 6)
            self.assertEqual(len(model.molecules(mass__gt=90)), 4)
            self.assertEqual(len(model.molecules(mass__gt=1000)), 2)
            self.assertEqual(len(model.molecules(mass__gt=90, mass__lt=1000)), 2)

            for optimise in [False, True]:
                if optimise: model.optimise_distances()
                atom = model.atom(1587 if e == "pdb" else 1586)
                four_angstrom = atom.nearby_atoms(cutoff=4)
                self.assertEqual(len(four_angstrom), 10)
                self.assertEqual(
                sorted([atom.id for atom in four_angstrom]),
                [n - (e != "pdb") for n in [1576, 1582, 1583, 1584, 1586, 1588, 1589, 1590, 1591, 2957]]
                )
                self.assertEqual(len(atom.nearby_atoms(cutoff=4, element="O")), 1)
                four_angstrom = model.atoms_in_sphere(atom.location, 4)
                self.assertEqual(len(four_angstrom), 11)
                self.assertEqual(
                sorted([atom.id for atom in four_angstrom]),
                [n - (e != "pdb") for n in [1576, 1582, 1583, 1584, 1586, 1587, 1588, 1589, 1590, 1591, 2957]]
                )
                self.assertEqual(len(model.atoms_in_sphere(atom.location, 4, element="O")), 1)
                self.assertEqual(len(model.atoms_in_sphere([10, 20, 30], 40)), 1281)
                self.assertEqual(len(model.atoms_in_sphere([10, 20, 30], 41)), 1360)
                self.assertEqual(len(model.atoms_in_sphere([10, 20, 30], 40, element="C")), 760)
                self.assertEqual(len(model.atoms_in_sphere([10, 20, 30], 39, element="C")), 711)

                atom = model.atom(905)
                self.assertEqual(len(atom.nearby_hets(5)), 9)
                self.assertEqual(len(atom.nearby_hets(5, ligands=False)), 7)
                self.assertEqual(len(atom.nearby_hets(5, het__is_water=False)), 8)
                self.assertEqual(len(atom.nearby_hets(5, residues=False)), 2)
                self.assertEqual(len(atom.nearby_hets(5, element="O")), 4)

            model.dehydrate()
            self.assertEqual(model.waters(), set())


    def test_5xme(self):
        for e in ["cif", "mmtf", "pdb"]:
            f = atomium.open("tests/integration/files/5xme." + e)
            self.assertEqual(f.resolution, None)
            models = f.models
            self.assertEqual(len(models), 10)
            self.assertIs(f.model, f.models[0])
            x_values = [
             33.969, 34.064, 37.369, 36.023, 35.245,
             35.835, 37.525, 35.062, 36.244, 37.677
            ]
            all_atoms = set()
            for x, model in zip(x_values, models):
                self.assertEqual(len(model.atoms()), 1827)
                all_atoms.update(model.atoms())
                atom = model.chain()[0].atom(name="N")
                self.assertEqual(atom.location[0], x)
            self.assertEqual(len(all_atoms), 18270)


    def test_1cbn(self):
        for e in ["cif", "mmtf", "pdb"]:
            f = atomium.open("tests/integration/files/1cbn." + e)
            chain = f.model.chain()
            residue1, residue2, residue3 = chain[:3]
            self.assertEqual(len(residue1.atoms()), 16)
            self.assertEqual(len(residue2.atoms()), 14)
            self.assertEqual(len(residue3.atoms()), 10)
            for residue in chain[:3]:
                for name in ["N", "C", "CA", "CB"]:
                    self.assertEqual(len(residue.atoms(name=name)), 1)


    def test_1xda(self):
        for e in ["cif", "mmtf", "pdb"]:
            f = atomium.open("tests/integration/files/1xda." + e)
            self.assertEqual(len(f.model.atoms()), 1842)
            self.assertEqual(len(f.model.atoms(is_metal=True)), 4)
            self.assertEqual(len(f.model.atoms(is_metal=False)), 1838)

            model = f.model
            self.assertEqual(len(model.atoms()), 1842)
            self.assertEqual(len(model.chains()), 8)
            self.assertEqual(len(model.ligands()), 16)

            model = f.generate_assembly(1)
            self.assertEqual(len(model.chains()), 2)
            self.assertEqual(set([c.id for c in model.chains()]), {"A", "B"})
            self.assertEqual(len(model.ligands()), 4)

            model = f.generate_assembly(2)
            self.assertEqual(len(model.chains()), 2)
            self.assertEqual(set([c.id for c in model.chains()]), {"C", "D"})
            self.assertEqual(len(model.ligands()), 4)

            model = f.generate_assembly(3)
            self.assertEqual(len(model.chains()), 2)
            self.assertEqual(set([c.id for c in model.chains()]), {"E", "F"})
            self.assertEqual(len(model.ligands()), 4)

            model = f.generate_assembly(4)
            self.assertEqual(len(model.chains()), 2)
            self.assertEqual(set([c.id for c in model.chains()]), {"G", "H"})
            self.assertEqual(len(model.ligands()), 4)

            model = f.generate_assembly(7)
            self.assertEqual(len(model.chains()), 6)
            self.assertEqual(set([c.id for c in model.chains()]), {"A", "B"})
            self.assertEqual(len(model.ligands()), 12)
            zn = model.atom(element="ZN")
            liganding_residues = zn.nearby_hets(3, is_metal=False, element__ne="CL")
            self.assertEqual(len(liganding_residues), 3)
            self.assertEqual(set([r.id for r in liganding_residues]), {"B.10"})
            self.assertEqual(set([r.name for r in liganding_residues]), {"HIS"})
            res1, res2, res3 = liganding_residues

            self.assertGreater(res1.atom(name="N").distance_to(res2.atom(name="N")), 10)
        

    def test_4opj(self):
        for e in ["cif", "mmtf", "pdb"]:
            f = atomium.open("tests/integration/files/4opj." + e)
            if e == "cif":
                self.assertEqual(
                 f.model.residue("B.6").full_name,
                 "(2R,3aS,4aR,5aR,5bS)-2-(6-amino-9H-purin-9-yl)-3a-hydroxyhexahydrocyclopropa[4,5]cyclopenta[1,2-b]furan-5a(4H)-yl dihydrogen phosphate"
                )
            elif e =="mmtf":
                self.assertEqual(f.model.residue("B.6").full_name, "TCY")
            else:
                self.assertEqual(
                 f.model.residue("B.6").full_name,
                 "(2R,3AS,4AR,5AR,5BS)-2-(6-AMINO-9H-PURIN-9-YL)-3A-HYDROXYHEXAHYDROCYCLOPROPA[4,5]CYCLOPENTA[1,2-B]FURAN-5A(4H)-YL DIHYDROGEN PHOSPHATE"
                )
