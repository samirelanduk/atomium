from datetime import date
from unittest import TestCase
import atomium

class DictParsingTests(TestCase):

    def test_1lol_mmcif(self):
        d = atomium.open("tests/integration/files/1lol.cif", dictionary=True)
        self.assertEqual(len(d.keys()), 66)



class ParsingTests(TestCase):

    def test_1lol_mmcif(self):
        # Open file
        pdb = atomium.open("tests/integration/files/1lol.cif")

        # File information
        self.assertEqual(pdb.name, "1LOL")
        self.assertEqual(len(pdb.source.keys()), 66)
        self.assertEqual(pdb.source["entry"][0]["id"], "1LOL")
        self.assertEqual(pdb.entry__id, "1LOL")
        self.assertEqual(pdb.title, "Crystal structure of orotidine monophosphate decarboxylase complex with XMP")
        self.assertEqual(pdb.deposition_date, date(2002, 5, 6))
        self.assertEqual(pdb.classification, "LYASE")
        self.assertEqual(pdb.keywords, ["TIM barrel", "LYASE"])
        self.assertEqual(pdb.authors, ["Wu, N.", "Pai, E.F."])
        self.assertEqual(pdb.technique, "X-RAY DIFFRACTION")
        self.assertEqual(pdb.source_organism, "Methanothermobacter thermautotrophicus str. Delta H")
        self.assertEqual(pdb.expression_system, "Escherichia coli")
        self.assertEqual(pdb.resolution, 1.9)
        self.assertEqual(pdb.rvalue, 0.193)
        self.assertEqual(pdb.rfree, 0.229)

        # More complex file information
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
        self.assertEqual(pdb.missing_residues, missing_residues)
        self.assertEqual(pdb.assemblies, [{
            "id": 1, "software": "PISA", "delta_energy": -31.0,
            "buried_surface_area":5230, "surface_area": 16550,
            "transformations": [{
                "chains": ["A", "B", "C", "D", "E", "F", "G", "H"],
                "matrix": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                "vector": [0.0, 0.0, 0.0]
            }]
        }])

        # Models
        self.assertEqual(len(pdb.models), 1)
        self.assertIs(pdb.model, pdb.models[0])
        self.assertEqual(str(pdb.model), "<Model (2 chains, 4 ligands)>")
        self.assertEqual(len(pdb.model.entities()), 4)
        self.assertEqual(len(pdb.model.chains()), 2)
        self.assertEqual(len(pdb.model.carbohydrates()), 0)
        self.assertEqual(len(pdb.model.ligands()), 4)
        self.assertEqual(len(pdb.model.waters()), 180)
        self.assertEqual(len(pdb.model.residues()), 418)
        self.assertEqual(len(pdb.model.atoms()), 3431)

        # Entities
        entity = pdb.model.entity("1")
        self.assertEqual(str(entity), "<Entity 1 (polymer)>")
        self.assertEqual(entity.id, "1")
        self.assertEqual(entity.type, "polymer")
        self.assertEqual(entity.name, "orotidine 5'-monophosphate decarboxylase")
        self.assertTrue(entity.sequence.startswith("LRSRRVDVMDVMNRLILAMDL"))
        self.assertTrue(entity.sequence.endswith("LADNPAAAAAGIIESIKDLLIPE"))
        self.assertEqual(len(entity.molecules()), 2)
        entity = pdb.model.entity("2")
        self.assertEqual(str(entity), "<Entity 2 (non-polymer)>")
        self.assertEqual(entity.id, "2")
        self.assertEqual(entity.type, "non-polymer")
        self.assertEqual(entity.name, "1,3-BUTANEDIOL")
        self.assertEqual(entity.sequence, "")
        self.assertEqual(len(entity.molecules()), 2)
        entity = pdb.model.entity("4")
        self.assertEqual(str(entity), "<Entity 4 (water)>")
        self.assertEqual(entity.id, "4")
        self.assertEqual(entity.type, "water")
        self.assertEqual(entity.name, "water")
        self.assertEqual(entity.sequence, "")
        self.assertEqual(len(entity.molecules()), 180)

        # Chains
        chain_a = pdb.model.chain("A")
        self.assertEqual(str(chain_a), "<Chain A (204 residues)>")
        self.assertEqual(chain_a.id, "A")
        self.assertEqual(chain_a.asym_id, "A")
        self.assertEqual(chain_a.auth_asym_id, "A")
        self.assertEqual(len(chain_a), 204)
        self.assertEqual(len(chain_a.residues()), 204)
        for res in chain_a: self.assertIn(res, chain_a)
        self.assertIs(chain_a[2], chain_a.residues()[2])
        self.assertTrue(chain_a.sequence.startswith("LRSRRVDVMDVMNRLILAMDL"))
        self.assertTrue(chain_a.sequence.endswith("LADNPAAAAAGIIESIKDLLIPE"))
        self.assertEqual(len(chain_a.ligands()), 2)
        self.assertEqual(len(chain_a.waters()), 96)
        self.assertIs(chain_a.model, pdb.model)
        self.assertIs(chain_a.entity, pdb.model.entity("1"))
        self.assertEqual(len(chain_a.atoms()), 1557)
        self.assertTrue(list(chain_a.atoms())[0] in chain_a)
        self.assertEqual(len(chain_a.helices), 11)
        self.assertEqual(chain_a.helices[0], chain_a.residues()[:3])
        self.assertEqual(len(chain_a.helices[-1]), 12)
        self.assertEqual(len(chain_a.strands), 9)
        self.assertEqual(chain_a.strands[0], chain_a.residues()[4:9])

        # Residues
        residue = pdb.model.residue("A.11")
        self.assertEqual(str(residue), "<Residue VAL (A.11)>")
        self.assertEqual(residue.id, "A.11")
        self.assertEqual(residue.name, "VAL")
        self.assertIs(residue.chain, chain_a)
        self.assertIsNone(residue.carbohydrate)
        self.assertIs(residue.model, pdb.model)
        self.assertEqual(len(residue.atoms()), 7)
        for atom in residue.atoms(): self.assertIn(atom, residue)

        # Ligands
        ligand = pdb.model.ligand("A.5001")
        self.assertEqual(str(ligand), "<Ligand BU2 (A.5001)>")
        self.assertFalse(ligand.is_water)
        self.assertEqual(ligand.id, "A.5001")
        self.assertEqual(ligand.name, "BU2")
        self.assertIs(ligand.entity, pdb.model.entity("2"))
        self.assertIs(ligand.chain, chain_a)
        self.assertIsNone(ligand.carbohydrate)
        self.assertIs(ligand.model, pdb.model)
        self.assertEqual(len(ligand.atoms()), 6)
        for atom in ligand.atoms(): self.assertIn(atom, ligand)

        # Waters
        for water in pdb.model.waters():
            self.assertEqual(len(water.atoms()), 1)
        water = pdb.model.water("A.3005")
        self.assertEqual(str(water), "<Water HOH (A.3005)>")
        self.assertTrue(water.is_water)
        self.assertIs(water.model, pdb.model)
        self.assertIs(water.entity, pdb.model.entity("4"))
        self.assertIsNone(water.carbohydrate)
        self.assertIs(water.chain, chain_a)
        self.assertIn(water.atom(), water)

        # Atoms 2649
        atom = pdb.model.atom(3231)
        self.assertEqual(str(atom), "<Atom 3231 (O5')>")
        self.assertEqual(atom.element, "O")
        self.assertEqual(atom.name, "O5'")
        self.assertEqual(atom.charge, 0)
        self.assertEqual(atom.bvalue, 34.58)
        self.assertEqual(list(atom.location), [-22.61, 62.264, 52.258])
        for i, value in enumerate(atom): self.assertEqual(value, atom.location[i])
        self.assertIsNone(atom.anisotropy)
        self.assertIs(atom.model, pdb.model)
        self.assertIs(atom.ligand, pdb.model.ligand("B.2002"))
        self.assertIsNone(atom.residue)
        self.assertIs(atom.chain, pdb.model.chain("B"))
        atom = pdb.model.atom(2649)
        self.assertEqual(atom.element, "N")
        self.assertEqual(atom.name, "NZ")
        self.assertIs(atom.model, pdb.model)
        self.assertIs(atom.residue, pdb.model.residue("B.1152"))
        self.assertIsNone(atom.ligand)
        self.assertIs(atom.chain, pdb.model.chain("B"))
        self.assertIsNone(atom.carbohydrate)


    def test_5xme_mmcif(self):
        # Parse
        pdb = atomium.open("tests/integration/files/5xme.cif")
        self.assertIsNone(pdb.resolution)

        # Models are correct
        self.assertEqual(len(pdb.models), 10)
        self.assertIs(pdb.model, pdb.models[0])
        x_values = [
            33.969, 34.064, 37.369, 36.023, 35.245,
            35.835, 37.525, 35.062, 36.244, 37.677
        ]
        all_atoms = set()
        for x, model in zip(x_values, pdb.models):
            self.assertEqual(len(model.atoms()), 1827)
            all_atoms.update(model.atoms())
            res = model.chain()[0]
            atom = res.atom(name="N")
            self.assertEqual(atom.location[0], x)
        self.assertEqual(len(all_atoms), 18270)
    

    def test_1cbn(self):
        # Multiple occupancy is handled
        pdb = atomium.open("tests/integration/files/1cbn.cif")
        chain = pdb.model.chain()
        residue1, residue2, residue3 = chain[:3]
        self.assertEqual(len(residue1.atoms()), 16)
        self.assertEqual(len(residue2.atoms()), 14)
        self.assertEqual(len(residue3.atoms()), 10)
        for residue in chain[:3]:
            for name in ["N", "C", "CA", "CB"]:
                self.assertEqual(len(residue.atoms(name=name)), 1)
    

    def test_3jbp(self):
        # Multi character secondary structure
        pdb = atomium.open("tests/integration/files/3jbp.cif")
        self.assertEqual(len(pdb.model.chain("TA").helices), 4)
    

    def test_6xlu(self):
        # Carbs parsed
        pdb = atomium.open("tests/integration/files/6xlu.cif")

        # Model
        self.assertEqual(len(pdb.model.chains()), 3)
        self.assertEqual(len(pdb.model.carbohydrates()), 15)
        self.assertEqual(len(pdb.model.residues()), 3178)
        self.assertEqual(len(pdb.model.ligands()), 32)
        self.assertEqual(len(pdb.model.atoms()), 25948)
        self.assertEqual(len(pdb.model.entities()), 4)

        # Entities
        entity = pdb.model.entity("2")
        self.assertEqual(str(entity), "<Entity 2 (branched)>")
        self.assertEqual(entity.id, "2")
        self.assertEqual(entity.type, "branched")
        self.assertEqual(entity.name, "2-acetamido-2-deoxy-beta-D-glucopyranose-(1-4)-2-acetamido-2-deoxy-beta-D-glucopyranose")
        self.assertEqual(len(entity.molecules()), 15)
        entity = pdb.model.entity("2")

        # Carb
        carb = pdb.model.carbohydrate("D")
        self.assertEqual(str(carb), "<Carbohydrate D (2 residues)>")
        self.assertEqual(carb.id, "D")
        self.assertEqual(carb.asym_id, "D")
        self.assertEqual(carb.auth_asym_id, "D")
        self.assertEqual(len(carb.residues()), 2)
        self.assertEqual(len(carb.ligands()), 0)
        self.assertEqual(len(carb.waters()), 0)
        self.assertIs(carb.model, pdb.model)
        self.assertEqual(len(carb.atoms()), 28)
        self.assertTrue(list(carb.atoms())[0] in carb)
    

    def test_1xda_assembly_5(self):
        # Parse PDBe assembly files
        pdb = atomium.open("tests/integration/files/1xda-assembly-5.cif")
        self.assertEqual(len(pdb.model.entities(type="polymer")), 2)
        self.assertEqual(len(pdb.model.entities(type="non-polymer")), 4)
        self.assertEqual(pdb.model.entity("1").name, "FATTY ACID ACYLATED INSULIN")
        self.assertEqual(pdb.model.entity("3").name, "PHENOL")
        self.assertEqual(len(pdb.model.chains()), 12)