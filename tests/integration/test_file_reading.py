from datetime import date
import atomium
from unittest import TestCase

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
            res = chaina.residue("A.15")
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
            self.assertIs(atom.structure, model.residue("A.131"))

            self.assertEqual(model.molecule("A"), chaina)
            self.assertEqual(model.molecule("A.5001"), lig)
            self.assertEqual(model.molecule("A.3005"), hoh)
            self.assertEqual(len(model.molecules(mass__gt=18)), 6)
            self.assertEqual(len(model.molecules(mass__gt=90)), 4)
            self.assertEqual(len(model.molecules(mass__gt=1000)), 2)
            self.assertEqual(len(model.molecules(mass__gt=90, mass__lt=1000)), 2)

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

            atom = model.atom(905)
            self.assertEqual(len(atom.nearby_structures(5)), 8)
            self.assertEqual(len(atom.nearby_structures(5, ligands=False)), 7)
            self.assertEqual(len(atom.nearby_structures(5, waters=True)), 9)
            self.assertEqual(len(atom.nearby_structures(5, residues=False)), 1)
            self.assertEqual(len(atom.nearby_structures(5, element="O")), 3)

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
                self.assertEqual(atom.x, x)
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
        for e in ["cif"]:#, "mmtf", "pdb"]:
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
            liganding_residues = zn.nearby_structures(3, is_metal=False, element__ne="CL")
            self.assertEqual(len(liganding_residues), 3)
            self.assertEqual(set([r.id for r in liganding_residues]), {"B.10"})
            self.assertEqual(set([r.name for r in liganding_residues]), {"HIS"})
            res1, res2, res3 = liganding_residues

            self.assertGreater(res1.atom(name="N").distance_to(res2.atom(name="N")), 10)
