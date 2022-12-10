from datetime import date
from unittest import TestCase
import atomium

class DictParsingTests(TestCase):

    def test_1lol_mmcif(self):
        d = atomium.open("tests/integration/files/1lol.cif", dictionary=True)
        self.assertEqual(len(d.keys()), 66)
        self.assertEqual(d["entry"], [{"id": "1LOL"}])
        self.assertEqual(d["database_2"], [
            {"database_id": "PDB", "database_code": "1LOL"},
            {"database_id": "RCSB", "database_code": "RCSB016137"},
            {"database_id": "WWPDB", "database_code": "D_1000016137"},
        ])
        self.assertEqual(
            d["pdbx_database_related"][0]["details"],
            "orotidine monophosphate decarboxylase complexed with UMP"
        )
        self.assertEqual(d["entity"][2]["pdbx_description"], "XANTHOSINE-5'-MONOPHOSPHATE")
        self.assertEqual(
            d["entity_poly"][0]["pdbx_seq_one_letter_code"],
            "LRSRRVDVMDVMNRLILAMDLMNRDDALRVTGEVREYIDTVKIGYPLVLSEGMDIIAEFRKRFG" +
            "CRIIADFKVADIPETNEKICRATFKAGADAIIVHGFPGADSVRACLNVAEEMGREVFLLTEMSH" +
            "PGAEMFIQGAADEIARMGVDLGVKNYVGPSTRPERLSRLREIIGQDSFLISPGVGAQGGDPGET" +
            "LRFADAIIVGRSIYLADNPAAAAAGIIESIKDLLIPE"
        )
        self.assertEqual(
            d["struct_biol"][0]["details"],
            "The biologically functional unit is a dimer composed of the two monomers in the asymmetric unit."
        )
        self.assertEqual(
            d["pdbx_database_remark"][0]["text"], "\n".join([
                "sequence",
                "Author states that although residues 1 and 1001 are MET",
                "and residues 101 and 1101 are Arg according to the",
                "SwissProt entry, residues 1 and 1001 were LEU and residues",
                "101 and 1101 were Pro in the original construct cloned",
                "of MT genomic dna.",
            ])
        )
    

    def test_1lol_mmcif_compressed(self):
        d = atomium.open("tests/integration/files/1lol.cif.gz", dictionary=True)
        self.assertEqual(len(d.keys()), 66)
        self.assertEqual(d["entry"], [{"id": "1LOL"}])
    

    def test_2bfb(self):
        d = atomium.open("tests/integration/files/2bfb.cif", dictionary=True)
        self.assertEqual(d["pdbx_database_related"][0], {
            "db_name": "PDB", "db_id": "1DTW", "content_type": "unspecified",
            "details": "HUMAN BRANCHED-CHAIN ALPHA-KETO ACID DEHYDROGENASE"
        })
        self.assertEqual(d["pdbx_database_related"][1], {
            "db_name": "PDB", "db_id": "1OLS", "content_type": "unspecified",
            "details": "ROLES OF HIS291-ALPHA AND HIS146-BETA' IN THE REDUCTIVE ACYLATION REACTION CATALYZED BY HUMAN BRANCHED-CHAIN ALPHA-KETOACID DEHYDROGENASE"
        })
    

    def test_6uaj(self):
        d = atomium.open("tests/integration/files/6uaj.cif", dictionary=True)
        self.assertEqual(d["pdbx_struct_assembly_gen"], [{
            "assembly_id": "1", "oper_expression": "1",
            "asym_id_list": "A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,AA,BA,CA,DA,EA,FA,GA,HA,IA,JA,KA,LA,MA,NA,OA,PA,QA,RA,SA,TA,UA,VA"
        }])
        
    



class ParsingTests(TestCase):

    def test_1lol(self):
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
            "A.223", "A.224", "A.225", "A.226", "A.227", "A.228", "A.229", "B.1",
            "B.2", "B.3", "B.4", "B.5", "B.6", "B.7", "B.8",
            "B.9", "B.10", "B.182", "B.183", "B.184", "B.185", "B.186"
        ], [
            "LEU", "ARG", "SER", "ARG", "ARG", "VAL", "ASP", "VAL", "MET", "ASP",
            "VAL", "GLY", "ALA", "GLN", "GLY", "GLY", "ASP", "PRO", "LYS", "ASP",
            "LEU", "LEU", "ILE", "PRO", "GLU", "LEU", "ARG", "SER", "ARG", "ARG",
            "VAL", "ASP", "VAL", "MET", "ASP", "VAL", "GLY", "ALA", "GLN", "GLY"
        ])]
        self.assertEqual(pdb.missing_residues, missing_residues)

        # Model
        self.assertEqual(len(pdb.models), 1)
        self.assertIs(pdb.model, pdb.models[0])
        self.assertEqual(str(pdb.model), "<Model (2 polymers, 4 non-polymers)>")
        self.assertEqual(len(pdb.model.molecules()), 8)
        self.assertEqual(len(pdb.model.polymers()), 2)
        self.assertEqual(len(pdb.model.branched_polymers()), 0)
        self.assertEqual(len(pdb.model.non_polymers()), 4)
        self.assertEqual(len(pdb.model.waters()), 2)
        self.assertEqual(len(pdb.model.residues()), 418)
        self.assertEqual(len(pdb.model.atoms()), 3431)

        # Molecule A
        mol_a = pdb.model.polymer("A")
        self.assertEqual(str(mol_a), "<Polymer A (204 residues)>")
        self.assertEqual(mol_a.id, "A")
        self.assertEqual(mol_a.auth_id, "A")
        self.assertEqual(len(mol_a), 204)
        self.assertEqual(len(mol_a.residues()), 204)
        for res in mol_a: self.assertIn(res, mol_a)
        self.assertIs(mol_a[2], mol_a.residues()[2])
        self.assertIs(mol_a.model, pdb.model)
        self.assertEqual(len(mol_a.atoms()), 1557)
        self.assertTrue(list(mol_a.atoms())[0] in mol_a)
        self.assertEqual(len(mol_a.helices), 11)
        self.assertEqual(mol_a.helices[0], mol_a.residues()[:3])
        self.assertEqual(len(mol_a.helices[-1]), 12)
        self.assertEqual(len(mol_a.strands), 9)
        self.assertEqual(mol_a.strands[0], mol_a.residues()[4:9])
        self.assertEqual(mol_a.entity_id, "1")
        self.assertEqual(mol_a.entity_name, "orotidine 5'-monophosphate decarboxylase")
        self.assertEqual(len(mol_a.sequence), 229)
        self.assertTrue(mol_a.sequence.startswith("LRSRRVDVMDVMNRLILAMDL"))
        self.assertTrue(mol_a.sequence.endswith("LADNPAAAAAGIIESIKDLLIPE"))
        self.assertIn(mol_a, pdb.model)

        # Molecule B
        mol_b = pdb.model.polymer("B")
        self.assertEqual(str(mol_b), "<Polymer B (214 residues)>")
        self.assertEqual(mol_b.id, "B")
        self.assertEqual(mol_b.auth_id, "B")
        self.assertEqual(len(mol_b), 214)
        self.assertEqual(len(mol_b.residues()), 214)
        for res in mol_b: self.assertIn(res, mol_b)
        self.assertIs(mol_b[2], mol_b.residues()[2])
        self.assertIs(mol_b.model, pdb.model)
        self.assertEqual(len(mol_b.atoms()), 1634)
        self.assertTrue(list(mol_b.atoms())[0] in mol_b)
        self.assertEqual(len(mol_b.helices), 11)
        self.assertEqual(mol_b.helices[0], mol_b.residues()[:3])
        self.assertEqual(len(mol_b.helices[-1]), 16)
        self.assertEqual(len(mol_b.strands), 9)
        self.assertEqual(mol_b.strands[0], mol_b.residues()[4:8])
        self.assertEqual(mol_b.entity_id, "1")
        self.assertEqual(mol_b.entity_name, "orotidine 5'-monophosphate decarboxylase")
        self.assertEqual(len(mol_b.sequence), 229)
        self.assertTrue(mol_b.sequence.startswith("LRSRRVDVMDVMNRLILAMDL"))
        self.assertTrue(mol_b.sequence.endswith("LADNPAAAAAGIIESIKDLLIPE"))
        self.assertIn(mol_b, pdb.model)

        # Molecule C
        mol_c = pdb.model.non_polymer("C")
        self.assertEqual(str(mol_c), "<NonPolymer C (BU2)>")
        self.assertEqual(mol_c.id, "C")
        self.assertEqual(mol_c.auth_id, "A")
        self.assertEqual(mol_c.name, "BU2")
        self.assertIs(mol_c.model, pdb.model)
        self.assertEqual(len(mol_c.atoms()), 6)
        for atom in mol_c.atoms(): self.assertIn(atom, mol_c)
        self.assertEqual(mol_c.entity_id, "2")
        self.assertEqual(mol_c.entity_name, "1,3-BUTANEDIOL")
        self.assertIn(mol_c, pdb.model)

        # Molecule D
        mol_d = pdb.model.non_polymer("D")
        self.assertEqual(str(mol_d), "<NonPolymer D (XMP)>")
        self.assertEqual(mol_d.id, "D")
        self.assertEqual(mol_d.auth_id, "A")
        self.assertEqual(mol_d.name, "XMP")
        self.assertIs(mol_d.model, pdb.model)
        self.assertEqual(len(mol_d.atoms()), 24)
        for atom in mol_d.atoms(): self.assertIn(atom, mol_d)
        self.assertEqual(mol_d.entity_id, "3")
        self.assertEqual(mol_d.entity_name, "XANTHOSINE-5'-MONOPHOSPHATE")
        self.assertIn(mol_d, pdb.model)

        # Molecule E
        mol_e = pdb.model.non_polymer("E")
        self.assertEqual(str(mol_e), "<NonPolymer E (BU2)>")
        self.assertEqual(mol_e.id, "E")
        self.assertEqual(mol_e.auth_id, "B")
        self.assertEqual(mol_e.name, "BU2")
        self.assertIs(mol_e.model, pdb.model)
        self.assertEqual(len(mol_e.atoms()), 6)
        for atom in mol_e.atoms(): self.assertIn(atom, mol_e)
        self.assertEqual(mol_e.entity_id, "2")
        self.assertEqual(mol_e.entity_name, "1,3-BUTANEDIOL")
        self.assertIn(mol_e, pdb.model)

        # Molecule F
        mol_f = pdb.model.non_polymer("F")
        self.assertEqual(str(mol_f), "<NonPolymer F (XMP)>")
        self.assertEqual(mol_f.id, "F")
        self.assertEqual(mol_f.auth_id, "B")
        self.assertEqual(mol_f.name, "XMP")
        self.assertIs(mol_f.model, pdb.model)
        self.assertEqual(len(mol_f.atoms()), 24)
        for atom in mol_f.atoms(): self.assertIn(atom, mol_f)
        self.assertEqual(mol_f.entity_id, "3")
        self.assertEqual(mol_f.entity_name, "XANTHOSINE-5'-MONOPHOSPHATE")
        self.assertIn(mol_f, pdb.model)

        # Molecule G
        mol_g = pdb.model.water("G")
        self.assertEqual(str(mol_g), "<Water G (HOH)>")
        self.assertEqual(mol_g.id, "G")
        self.assertEqual(mol_g.auth_id, "A")
        self.assertEqual(mol_g.name, "HOH")
        self.assertIs(mol_g.model, pdb.model)
        self.assertEqual(len(mol_g.atoms()), 96)
        for atom in mol_g.atoms(): self.assertIn(atom, mol_g)
        self.assertEqual(mol_g.entity_id, "4")
        self.assertEqual(mol_g.entity_name, "water")
        self.assertIn(mol_g, pdb.model)

        # Molecule H
        mol_h = pdb.model.water("H")
        self.assertEqual(str(mol_h), "<Water H (HOH)>")
        self.assertEqual(mol_h.id, "H")
        self.assertEqual(mol_h.auth_id, "B")
        self.assertEqual(mol_h.name, "HOH")
        self.assertIs(mol_h.model, pdb.model)
        self.assertEqual(len(mol_h.atoms()), 84)
        for atom in mol_h.atoms(): self.assertIn(atom, mol_h)
        self.assertEqual(mol_h.entity_id, "4")
        self.assertEqual(mol_h.entity_name, "water")
        self.assertIn(mol_h, pdb.model)

        # Residue A.11
        residue = pdb.model.residue("A.11")
        self.assertEqual(str(residue), "<Residue VAL (A.11)>")
        self.assertEqual(residue.id, "A.11")
        self.assertEqual(residue.number, 11)
        self.assertEqual(residue.name, "VAL")
        self.assertIs(residue.polymer, mol_a)
        self.assertIsNone(residue.branched_polymer)
        self.assertIs(residue.model, pdb.model)
        self.assertEqual(len(residue.atoms()), 7)
        for atom in residue.atoms(): self.assertIn(atom, residue)
        self.assertIn(residue, pdb.model)

        # Atom 3231
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
        self.assertIsNone(atom.polymer)
        self.assertIsNone(atom.branched_polymer)
        self.assertIs(atom.non_polymer, mol_f)
        self.assertIsNone(atom.water)
        self.assertIsNone(atom.residue)
        self.assertIn(atom, pdb.model)

        # Atom 2649
        atom = pdb.model.atom(2649)
        self.assertEqual(atom.element, "N")
        self.assertEqual(atom.name, "NZ")
        self.assertIs(atom.model, pdb.model)
        self.assertIs(atom.polymer, mol_b)
        self.assertIsNone(atom.branched_polymer)
        self.assertIsNone(atom.non_polymer)
        self.assertIsNone(atom.water)
        self.assertIs(atom.residue, pdb.model.residue("B.152"))
        self.assertIn(atom, pdb.model)

        # Atom 3323
        atom = pdb.model.atom(3323)
        self.assertEqual(atom.element, "O")
        self.assertEqual(atom.name, "O")
        self.assertIs(atom.model, pdb.model)
        self.assertIsNone(atom.polymer,)
        self.assertIsNone(atom.branched_polymer)
        self.assertIsNone(atom.non_polymer)
        self.assertIs(atom.water, mol_g)
        self.assertIsNone(atom.residue)
        self.assertIn(atom, pdb.model)
    

    def test_2igd(self):
        # Anisotropy
        pdb = atomium.open("tests/integration/files/2igd.cif")
        atom = pdb.model.atom(1)
        self.assertEqual(atom.anisotropy, [
            0.1085, 0.125, 0.1612, -0.0074, -0.0061, 0.0108
        ])


    def test_5xme(self):
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
            res = model.polymer()[0]
            atom = res.atom(name="N")
            self.assertEqual(atom.location[0], x)
        self.assertEqual(len(all_atoms), 18270)
    

    def test_1cbn(self):
        # Multiple occupancy is handled
        pdb = atomium.open("tests/integration/files/1cbn.cif")
        polymer = pdb.model.polymer()
        residue1, residue2, residue3 = polymer[:3]
        self.assertEqual(len(residue1.atoms()), 16)
        self.assertEqual(len(residue2.atoms()), 14)
        self.assertEqual(len(residue3.atoms()), 10)
        for residue in polymer[:3]:
            for name in ["N", "C", "CA", "CB"]:
                self.assertEqual(len(residue.atoms(name=name)), 1)
    

    def test_3jbp(self):
        # Multi character secondary structure
        pdb = atomium.open("tests/integration/files/3jbp.cif")
        self.assertEqual(len(pdb.model.polymer("TA").helices), 4)
    

    def test_6xlu(self):
        # Carbs parsed
        pdb = atomium.open("tests/integration/files/6xlu.cif")

        # Model
        self.assertEqual(len(pdb.model.polymers()), 3)
        self.assertEqual(len(pdb.model.branched_polymers()), 15)
        self.assertEqual(len(pdb.model.residues()), 3208)
        self.assertEqual(len(pdb.model.non_polymers()), 32)
        self.assertEqual(len(pdb.model.atoms()), 25948)

        # Carb
        carb = pdb.model.branched_polymer("D")
        self.assertIn(carb, pdb.model)
        self.assertEqual(str(carb), "<BranchedPolymer D (2 residues)>")
        self.assertEqual(carb.id, "D")
        self.assertEqual(carb.auth_id, "D")
        self.assertEqual(len(carb.residues()), 2)
        for residue in carb.residues(): self.assertIn(residue, carb)
        self.assertIs(carb.model, pdb.model)
        self.assertEqual(len(carb.atoms()), 28)
        self.assertTrue(list(carb.atoms())[0] in carb)
        self.assertEqual(carb.entity_id, "2")
        self.assertEqual(carb.entity_name, "2-acetamido-2-deoxy-beta-D-glucopyranose-(1-4)-2-acetamido-2-deoxy-beta-D-glucopyranose")

        # Residue
        residue = pdb.model.residue("D.2")
        self.assertEqual(residue.id, "D.2")
        self.assertEqual(residue.number, 2)
        self.assertIs(residue.branched_polymer, carb)

    

    def test_1xda_assembly_5(self):
        # Parse PDBe assembly files
        pdb = atomium.open("tests/integration/files/1xda-assembly-5.cif")
        self.assertEqual(pdb.resolution, 1.8)
        self.assertIsNone(pdb.rvalue)

        self.assertEqual(len(pdb.model.polymers()), 12)
        self.assertEqual(len(pdb.model.non_polymers()), 24)
        self.assertEqual(pdb.model.polymer("E").entity_name, "FATTY ACID ACYLATED INSULIN")
        self.assertEqual(pdb.model.non_polymer("Q").entity_name, "PHENOL")
        self.assertEqual(len(set(m.entity_id for m in pdb.model.molecules())), 7)
    

    def test_12ca(self):
        # Values not ?
        pdb = atomium.open("tests/integration/files/12ca.cif")
        self.assertIsNone(pdb.expression_system)



class NetworkTests(TestCase):

    def test_fetching_3nir_mmcif_rcsb(self):
        pdb = atomium.fetch("3nir")
        self.assertEqual(pdb.title, "Crystal structure of small protein crambin at 0.48 A resolution")
    

    def test_fetching_3nir_mmcif_pdbe(self):
        pdb = atomium.fetch("https://www.ebi.ac.uk/pdbe/static/entry/3nir_updated.cif")
        self.assertEqual(pdb.title, "Crystal structure of small protein crambin at 0.48 A resolution")