from datetime import datetime
import atomium
from tests.integration.base import IntegratedTest

class PdbReadingTests(IntegratedTest):

    def test_can_read_pdb(self):
        pdb = atomium.pdb_from_file("tests/integration/files/1lol.pdb")
        self.assertEqual(pdb.code, "1LOL")
        self.assertEqual(
         pdb.title,
         "CRYSTAL STRUCTURE OF OROTIDINE MONOPHOSPHATE DECARBOXYLASE COMPLEX WITH XMP"
        )
        self.assertEqual(pdb.deposition_date, datetime(2002, 5, 6).date())
        self.assertEqual(pdb.organism, "METHANOTHERMOBACTER THERMAUTOTROPHICUS STR. DELTA H")
        self.assertEqual(pdb.expression_system, "ESCHERICHIA COLI")
        self.assertEqual(pdb.technique, "X-RAY DIFFRACTION")
        self.assertEqual(pdb.classification, "LYASE")
        self.assertEqual(pdb.resolution, 1.9)
        self.assertEqual(pdb.rfactor, 0.193)
        self.assertEqual(pdb.rfree, 0.229)
        self.assertEqual(pdb.rcount, 1583)
        self.assertEqual(pdb.keywords, ["TIM BARREL", "LYASE"])
        self.assertEqual(pdb.biomolecules, [{
         "id": 1,
         "software": "PISA",
         "delta_energy": -31.0,
         "buried_surface_area": 5230,
         "surface_area": 16550,
         "transformations": [{
          "chains": ["A", "B"],
          "matrix": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
          "vector": [0.0, 0.0, 0.0]
         }]
        }])

        # Atoms are correct
        model = pdb.model
        self.assertEqual(len(model.atoms()), 3431)
        atom = model.atom(2934)
        self.assertEqual(atom.anisotropy, [0, 0, 0, 0, 0, 0])
        self.assertEqual(
         model.atom(1).anisotropy,
         [0.2406, 0.1892, 0.1614, 0.0198, 0.0519, -0.0328]
        )
        self.assertEqual(atom.element, "N")
        self.assertEqual(atom.name, "NE")
        self.assertEqual(atom.location, (-20.082, 79.647, 41.645))
        self.assertEqual(atom.bfactor, 35.46)
        self.assertEqual(atom.charge, 0)
        self.assertAlmostEqual(
         model.mass, 46018.5, delta=0.005
        )
        self.assertEqual(atom.model, model)

        # Chains are correct
        self.assertEqual(len(model.chains()), 2)
        for chain in model.chains():
            self.assertIs(chain.model, model)
            self.assertIsNone(chain.name)
            self.assertIn(chain, model)
        chaina, chainb = model.chain(id="A"), model.chain(id="B")
        self.assertIs(atom.chain, chainb)

        # Residues are correct
        self.assertEqual(chaina[0].name, "VAL")
        self.assertEqual(chaina[0].full_name, "valine")
        self.assertEqual(chaina[0].next.name, "MET")
        self.assertEqual(chaina[0].next.full_name, "methionine")
        self.assertEqual(chaina[-1].name, "ILE")
        self.assertEqual(chaina[-1].full_name, "isoleucine")
        self.assertEqual(chaina[-1].previous.name, "SER")
        self.assertEqual(len(chaina.residues(name="ASN")), 6)
        for residue in chaina:
            self.assertIs(residue.model, model)
            self.assertIs(residue.chain, chaina)
            self.assertIn(residue, chaina)
            self.assertIn(residue, model)
        for residue in chainb:
            self.assertIs(residue.model, model)
            self.assertIs(residue.chain, chainb)
            self.assertIn(residue, chainb)
            self.assertIn(residue, model)
        self.assertTrue(chaina.sequence.startswith("VMNRLILAMDLMNRDDALRVTGEVR"))
        self.assertTrue(chaina.sequence.endswith("ADAIIVGRSIYLADNPAAAAAGIIESI"))
        self.assertTrue(chainb.sequence.startswith("VMNRLILAMDLMNRDDALRVTGEVR"))
        self.assertTrue(chainb.sequence.endswith("RSIYLADNPAAAAAGIIESIKDLLIPE"))
        self.assertTrue(chaina.rep_sequence.startswith("LRSRRVDVMDVMNRLILAMDL"))
        self.assertTrue(chaina.rep_sequence.endswith("LADNPAAAAAGIIESIKDLLIPE"))
        self.assertTrue(chainb.rep_sequence.startswith("LRSRRVDVMDVMNRLILAMDL"))
        self.assertTrue(chainb.rep_sequence.endswith("LADNPAAAAAGIIESIKDLLIPE"))
        residue = chaina.residue("A13")
        self.assertEqual(residue.name, "ASN")
        self.assertEqual(len(residue.atoms()), 8)
        self.assertEqual(len(residue.atoms(element="O")), 2)
        for atom in residue.atoms():
            self.assertIs(atom.residue, residue)
            self.assertIs(atom.chain, chaina)
            self.assertIs(atom.model, model)
        gly = chaina.residue(id="A32")
        pairs = list(gly.pairwise_atoms())
        self.assertEqual(len(pairs), 6)
        for pair in pairs:
            pair = list(pair)
            self.assertTrue(0 < pair[0].distance_to(pair[1]), 5)

        # Ligands are correct
        self.assertEqual(len(model.ligands()), 184)
        self.assertEqual(len(model.ligands(water=False)), 4)
        self.assertEqual(len(model.ligands(name="XMP")), 2)
        self.assertEqual(len(model.ligands(name="BU2")), 2)
        self.assertEqual(len(model.ligands(name="HOH")), 180)
        mol = model.ligand(id="A2001")
        self.assertIs(mol.model, model)
        self.assertEqual(mol.name, "XMP")
        self.assertEqual(len(mol.atoms()), 24)
        self.assertEqual(mol.formula, {"C": 10, "O": 9, "N": 4, "P": 1})
        mol1, mol2 = model.ligand("A5001"), model.ligand("B5002")
        self.assertEqual(mol1.pairing_with(mol2), {
         model.atom(3194): model.atom(3224),
         model.atom(3195): model.atom(3225),
         model.atom(3196): model.atom(3226),
         model.atom(3197): model.atom(3227),
         model.atom(3198): model.atom(3228),
         model.atom(3199): model.atom(3229),
        })
        self.assertIn(mol1, chaina)
        self.assertIn(mol2, chainb)
        self.assertEqual(mol1.rmsd_with(mol1), 0)
        self.assertAlmostEqual(mol1.rmsd_with(mol2), 23.5, delta=0.5)

        # Bonding is correct
        residue = chaina[0]
        n = residue.atom(name="N")
        ca = residue.atom(name="CA")
        c = residue.atom(name="C")
        o = residue.atom(name="O")
        cb = residue.atom(name="CB")
        cg1 = residue.atom(name="CG1")
        cg2 = residue.atom(name="CG2")
        next_atom = chaina[1].atom(name="N")
        self.assertEqual(n.bonded_atoms, set([ca]))
        self.assertEqual(ca.bonded_atoms, set([n, c, cb]))
        self.assertEqual(c.bonded_atoms, set([ca, o, next_atom]))
        self.assertEqual(o.bonded_atoms, set([c]))
        self.assertEqual(cb.bonded_atoms, set([ca, cg1, cg2]))
        self.assertEqual(cg1.bonded_atoms, set([cb]))
        self.assertEqual(cg2.bonded_atoms, set([cb]))
        res181 = chaina.residue("A181")
        res190 = chaina.residue("A190")
        c2, n2 = res181.atom(name="C"), res190.atom(name="N")
        self.assertNotIn(c2, n2.bonded_atoms)
        self.assertNotIn(n2, c2.bonded_atoms)
        self.assertIn(
         model.atom(3194), model.atom(3195).bonded_atoms
        )

        # Can get atoms in cutoff distance
        atom = model.atom(1587)
        four_angstrom = atom.nearby_atoms(cutoff=4)
        self.assertEqual(len(four_angstrom), 10)
        self.assertEqual(
         sorted([atom.id for atom in four_angstrom]),
         [1576, 1582, 1583, 1584, 1586, 1588, 1589, 1590, 1591, 2957]
        )
        self.assertEqual(len(atom.nearby_atoms(cutoff=4, element="O")), 1)
        four_angstrom = model.atoms_in_sphere(*atom.location, radius=4)
        self.assertEqual(len(four_angstrom), 11)
        self.assertEqual(
         sorted([atom.id for atom in four_angstrom]),
         [1576, 1582, 1583, 1584, 1586, 1587, 1588, 1589, 1590, 1591, 2957]
        )
        self.assertEqual(len(model.atoms_in_sphere(*atom.location, radius=4, element="O")), 1)

        # 'Bindsites' are correct
        nearby = model.ligand("A5001").nearby_residues(4)
        self.assertEqual(nearby, set([
         model.residue("A42"), model.residue("A70"), model.residue("A72"),
         model.residue("A96"), model.residue("A123"), model.residue("A155")
        ]))
        nearby = model.ligand("A5001").nearby_residues(4, ligands=True)
        self.assertEqual(nearby, set([
         model.residue("A42"), model.residue("A70"), model.residue("A72"),
         model.residue("A96"), model.residue("A123"), model.residue("A155"),
         model.ligand("A3015"), model.ligand("A2001")
        ]))


    def test_can_read_multi_model_pdbs(self):
        pdb = atomium.pdb_from_file("tests/integration/files/5xme.pdb")
        self.assertEqual(pdb.resolution, None)
        self.assertEqual(pdb.biomolecules, [{
         "id": 1,
         "software": None,
         "delta_energy": None,
         "buried_surface_area": None,
         "surface_area": None,
         "transformations": [{
          "chains": ["A"],
          "matrix": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
          "vector": [0.0, 0.0, 0.0]
         }]
        }])

        models = pdb.models
        self.assertEqual(len(models), 10)
        self.assertIs(pdb.model, pdb.models[0])

        x_values = [
         33.969, 34.064, 37.369, 36.023, 35.245,
         35.835, 37.525, 35.062, 36.244, 37.677
        ]
        all_atoms = set()
        for x, model in zip(x_values, models):
            self.assertEqual(len(model.atoms()), 1827)
            all_atoms.update(model.atoms())
            atom = model.atom(1)
            self.assertEqual(atom.x, x)
            self.assertEqual(len(atom.bonded_atoms), 1)
        self.assertEqual(len(all_atoms), 18270)


    def test_can_read_alt_loc_pdbs(self):
        pdb = atomium.pdb_from_file("tests/integration/files/1cbn.pdb")
        chain = pdb.model.chain()
        residue1, residue2, residue3 = chain[:3]

        # Residues have the correct number of atoms
        self.assertEqual(len(residue1.atoms()), 16)
        self.assertEqual(len(residue2.atoms()), 14)
        self.assertEqual(len(residue3.atoms()), 10)

        for residue in chain[:3]:
            for name in ["N", "C", "CA", "CB"]:
                self.assertEqual(len(residue.atoms(name=name)), 1)


    def test_can_create_biological_assemblies(self):
        pdb = atomium.pdb_from_file("tests/integration/files/1xda.pdb")

        # The PDB has the correct instructions for creating assemblies
        data = [
         [1, 1720, 3980, -7, [["A,B", 0, 0]]],
         [2, 1870, 4400, -2, [["C,D", 0, 0]]],
         [3, 1160, 4110, -11, [["E,F", 0, 0]]],
         [4, 1650, 4240, -7, [["G,H", 0, 0]]],
         [5, 21680, 12240, -332, [
          ["E,F,G,H", 0, 0],
          ["E,F,G,H", [[-0.5, -0.866025, 0.0], [0.866025, -0.5, 0.0], [0, 0, 1]], 0],
          ["E,F,G,H", [[-0.5, 0.866025, 0.0], [-0.866025, -0.5, 0.0], [0, 0, 1]], 0]]
         ],
         [6, 23900, 12440, -287, [
          ["A,B,C,D", 0, 0],
          ["A,B,C,D", [[-0.5, -0.866025, 0.0], [0.866025, -0.5, 0.0], [0, 0, 1]], 0],
          ["A,B,C,D", [[-0.5, 0.866025, 0.0], [-0.866025, -0.5, 0.0], [0, 0, 1]], 0]]
         ],
         [7, 7540, 9770, -137, [
          ["A,B", 0, 0],
          ["A,B", [[-0.5, -0.866025, 0.0], [0.866025, -0.5, 0.0], [0, 0, 1]], 0],
          ["A,B", [[-0.5, 0.866025, 0.0], [-0.866025, -0.5, 0.0], [0, 0, 1]], 0]]
         ],
         [8, 5500, 10530, -156, [
          ["E,F", 0, 0],
          ["E,F", [[-0.5, -0.866025, 0.0], [0.866025, -0.5, 0.0], [0, 0, 1]], 0],
          ["E,F", [[-0.5, 0.866025, 0.0], [-0.866025, -0.5, 0.0], [0, 0, 1]], 0]]
         ],
         [9, 7450, 10440, -149, [
          ["G,H", 0, 0],
          ["G,H", [[-0.5, -0.866025, 0.0], [0.866025, -0.5, 0.0], [0, 0, 1]], 0],
          ["G,H", [[-0.5, 0.866025, 0.0], [-0.866025, -0.5, 0.0], [0, 0, 1]], 0]]
         ],
         [10, 8150, 10880, -129, [
          ["C,D", 0, 0],
          ["C,D", [[-0.5, -0.866025, 0.0], [0.866025, -0.5, 0.0], [0, 0, 1]], 0],
          ["C,D", [[-0.5, 0.866025, 0.0], [-0.866025, -0.5, 0.0], [0, 0, 1]], 0]]
         ],
         [11, 4530, 6630, -25, [["E,F,G,H", 0, 0]]],
         [12, 5190, 6770, -14, [["A,B,C,D", 0, 0]]],
        ]
        self.assertEqual(pdb.biomolecules, [{
         "id": d[0],
         "software": "PISA",
         "buried_surface_area": d[1],
         "surface_area": d[2],
         "delta_energy": d[3],
         "transformations": [{
          "chains": t[0].split(","),
          "matrix": t[1] if t[1] else [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
          "vector": [0.0, 0.0, 0.0]
         } for t in d[4]]
        } for d in data])

        # Different assemblies can be generated
        self.assertEqual(len(pdb.model.chains()), 8)

        model = pdb.generate_assembly(1)
        self.assertEqual(len(model.chains()), 2)
        self.assertEqual(len(model.residues()), 50)
        self.assertEqual(len(model.ligands()), 44)
        for atom in model.atoms():
            self.assertIn(atom.chain, model.chains())
        self.assertFalse(model.chains() & pdb.model.chains())
        self.assertFalse(model.residues() & pdb.model.residues())
        self.assertFalse(model.ligands() & pdb.model.ligands())

        model = pdb.generate_best_assembly()
        self.assertEqual(len(model.chains()), 12)
        self.assertEqual(len(model.residues()), 300)
        self.assertEqual(len(model.ligands()), 267)
        self.assertFalse(model.chains() & pdb.model.chains())
        self.assertFalse(model.residues() & pdb.model.residues())
        self.assertFalse(model.ligands() & pdb.model.ligands())
        self.assertEqual(pdb.best_assembly["id"], 5)


    def test_can_read_pdb_data(self):
        data_file = atomium.pdb_data_from_file(
         "tests/integration/files/1lol.pdb"
        )
        self.assertEqual(data_file["code"], "1LOL")
        self.assertEqual(len(data_file["models"][0]), 2)
        self.assertEqual(
         data_file["models"][0][0]["residues"][0]["atoms"][0],
         {
          "atom_id": 1, "atom_name": "N", "alt_loc": None,
          "residue_name": "VAL", "full_id": "A11",
          "chain_id": "A", "residue_id": 11, "insert_code": "",
          "x": 3.696, "y": 33.898, "z": 63.219,
          "occupancy": 1.0, "temp_factor": 21.50,
          "element": "N", "charge": 0.0,
          "anisotropy": [0.2406, 0.1892, 0.1614, 0.0198, 0.0519, -0.0328]
         }
        )
        self.assertEqual(len(data_file["connections"]), 60)
        self.assertEqual(data_file["connections"][0], {
         "atom": 3194, "bond_to": [3195, 3196]
        })


    def test_can_fetch_pdb(self):
        pdb = atomium.fetch("1h4W", pdbe=True)
        model = pdb.model

        self.assertEqual(len(model.chains()), 1)
        residue = model.residue("A221A")
        self.assertEqual(residue.next.id, "A222")
        self.assertEqual(residue.previous.id, "A221")


    def test_can_fetch_pdb_data(self):
        data_file = atomium.fetch_data("1lol", pdbe=True)
        self.assertEqual(data_file["code"], "1LOL")
        self.assertEqual(len(data_file["models"][0]), 2)
        self.assertEqual(len(data_file["connections"]), 60)
        self.assertEqual(data_file["connections"][0], {
         "atom": 3194, "bond_to": [3195, 3196]
        })



class PdbSavingTests(IntegratedTest):

    def test_can_save_single_model_pdb(self):
        pdb = atomium.pdb_from_file("tests/integration/files/1lol.pdb")

        # Can be saved
        pdb.code = "9SAM"
        pdb.deposition_date = datetime(1990, 9, 1).date()
        pdb.title = (
         "FOR IN THAT SLEEP OF DEATH, WHAT DREAMS MAY COME, WHEN WE HAVE " +
         "SHUFFLED OFF THIS MORTAL COIL MUST GIVE US PAUSE"
        )
        pdb.organism = "HOMO SAPIENS"
        pdb.expression_system = "COOL NEW ORGANISM"
        pdb.technique = "MEME DIFFRACTION"
        pdb.model.chain("B").rep_sequence = ""
        while pdb.keywords: pdb.keywords.pop()
        for k in ["AMAZING", "SUPERB", "WOW"]: pdb.keywords.append(k)
        pdb.save("tests/integration/files/1LOL2.pdb")
        self.check_files_the_same("1LOL2.pdb", "1lol_output.pdb")

        new = atomium.pdb_from_file("tests/integration/files/1LOL2.pdb")
        model = new.model
        self.assertAlmostEqual(
         model.mass, 46018.5, delta=0.005
        )


    def test_can_save_structures(self):
        pdb = atomium.pdb_from_file("tests/integration/files/1lol.pdb")

        # Save chains
        for chain in pdb.model.chains():
            chain.save("tests/integration/files/chain{}.pdb".format(chain.id))

        self.check_files_the_same("chainA.pdb", "chaina_output.pdb")
        self.check_files_the_same("chainB.pdb", "chainb_output.pdb")
        new = atomium.pdb_from_file("tests/integration/files/chainA.pdb")
        model = new.model
        self.assertEqual(len(model.chains()), 1)

        # Save molecules
        pdb.model.ligand("A5001").save("tests/integration/files/5001.pdb")
        self.check_files_the_same("5001.pdb", "5001_output.pdb")
        new = atomium.pdb_from_file("tests/integration/files/5001.pdb")
        model = new.model
        self.assertEqual(len(model.atoms()), 6)


    def test_can_save_multi_model_pdb(self):
        pdb = atomium.pdb_from_file("tests/integration/files/5xme.pdb")
        while pdb.keywords: pdb.keywords.pop()
        for k in ["INTEGRAL"] * 10: pdb.keywords.append(k)

        pdb.save("tests/integration/files/5XME2.pdb")
        self.check_files_the_same("5XME2.pdb", "5xme_output.pdb")
        new = atomium.pdb_from_file("tests/integration/files/5XME2.pdb")
        models = new.models
        self.assertEqual(len(models), 10)
        x_values = [
         33.969, 34.064, 37.369, 36.023, 35.245,
         35.835, 37.525, 35.062, 36.244, 37.677
        ]
        for x, model in zip(x_values, models):
            self.assertEqual(len(model.atoms()), 1827)


    def test_can_save_alt_loc_pdbs(self):
        pdb = atomium.pdb_from_file("tests/integration/files/1cbn.pdb")

        pdb.save("tests/integration/files/1CBN2.pdb")
        self.check_files_the_same("1CBN2.pdb", "1cbn_output.pdb")

        new = atomium.pdb_from_file("tests/integration/files/1CBN2.pdb")
        chain = pdb.model.chain()
        residue1, residue2, residue3 = chain[:3]
        self.assertEqual(len(residue1.atoms()), 16)
        self.assertEqual(len(residue2.atoms()), 14)
        self.assertEqual(len(residue3.atoms()), 10)


    def test_can_save_multi_assembly_pdbs(self):
        pdb = atomium.pdb_from_file("tests/integration/files/1xda.pdb")

        pdb.save("tests/integration/files/1XDA2.pdb")
        self.check_files_the_same("1XDA2.pdb", "1xda_output.pdb")

        new = atomium.pdb_from_file("tests/integration/files/1XDA2.pdb")
