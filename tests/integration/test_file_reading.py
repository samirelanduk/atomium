from datetime import date
from tests.integration.base import IntegratedTest
import atomium

class FileDictReadingTests(IntegratedTest):

    def test_1lol_pdb(self):
        d = atomium.open("tests/integration/files/1lol.pdb", file_dict=True)
        self.assertEqual(d["HEADER"], [
         "    LYASE                                   06-MAY-02   1LOL"
        ])
        self.assertEqual(len(d["ATOM"]), 3191)


    def test_glucose_xyz(self):
        d = atomium.open("tests/integration/files/glucose.xyz", file_dict=True)
        self.assertEqual(d["header_lines"], ["12", "glucose from 2gbp"])
        self.assertEqual(len(d["atom_lines"]), 12)


    def test_1lol_mmcif(self):
        d = atomium.open("tests/integration/files/1lol.cif", file_dict=True)
        self.assertEqual(d["entry"], [{"id": "1LOL"}])
        self.assertEqual(d["audit_author"], [
         {"name": "Wu, N.", "pdbx_ordinal": "1"},
         {"name": "Pai, E.F.", "pdbx_ordinal": "2"}
        ])
        entity = d["entity"]
        self.assertEqual(len(entity), 4)
        self.assertEqual(
         entity[0]["pdbx_description"], "orotidine 5'-monophosphate decarboxylase"
        )
        self.assertEqual(entity[1]["type"], "non-polymer")
        self.assertTrue(d["citation"][0]["title"].startswith("Crystal"))
        self.assertTrue(d["citation"][0]["title"].endswith("decarboxylase."))



class DataDictReadingTests(IntegratedTest):

    def test_1lol_pdb(self):
        d = atomium.open("tests/integration/files/1lol.pdb", data_dict=True)
        self.assertEqual(set(d.keys()), {
         "description", "experiment", "quality", "geometry", "models"
        })
        self.assertEqual(d["description"], {
         "code": "1LOL",
         "title": "CRYSTAL STRUCTURE OF OROTIDINE MONOPHOSPHATE DECARBOXYLASE COMPLEX WITH XMP",
         "deposition_date": date(2002, 5, 6),
         "classification": "LYASE",
         "keywords": ["TIM BARREL", "LYASE"],
         "authors": ["N.WU", "E.F.PAI"]
        })
        self.assertEqual(d["experiment"], {
         "technique": "X-RAY DIFFRACTION",
         "source_organism": "METHANOTHERMOBACTER THERMAUTOTROPHICUS STR. DELTA H",
         "expression_system": "ESCHERICHIA COLI"
        })
        self.assertEqual(d["quality"], {
         "resolution": 1.9, "rvalue": 0.193, "rfree": 0.229
        })
        self.assertEqual(d["geometry"], {"assemblies": [{
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
        }]})

        self.assertEqual(len(d["models"]), 1)
        model = d["models"][0]
        self.assertEqual(len(model["chains"]), 2)
        self.assertEqual(model["chains"][0]["id"], "A")
        self.assertEqual(model["chains"][1]["id"], "B")
        for chain in model["chains"]:
            self.assertEqual(len(chain["full_sequence"]), 229)
        self.assertEqual(len(model["residues"]), 418)
        self.assertEqual(model["residues"][0], {
         "id": "A11", "name": "VAL", "chain_id": "A"
        })
        self.assertEqual(model["residues"][-1], {
         "id": "B1229", "name": "GLU", "chain_id": "B"
        })
        self.assertEqual(len(model["ligands"]), 184)
        self.assertEqual(model["ligands"][0], {
         "id": "A5001", "name": "BU2", "chain_id": "A"
        })
        self.assertEqual(model["ligands"][-1], {
         "id": "B3180", "name": "HOH", "chain_id": "B"
        })
        self.assertEqual(len(model["atoms"]), 3431)
        self.assertEqual(model["atoms"][0], {
         "id": 1,
         "element": "N",
         "name": "N",
         "x": 3.696,
         "y": 33.898,
         "z": 63.219,
         "bfactor": 21.5,
         "charge": -1,
         "residue_id": 11,
         "residue_name": "VAL",
         "residue_insert": "",
         "chain_id": "A",
         "occupancy": 1.0,
         "alt_loc": None,
         "anisotropy": [0.2406, 0.1892, 0.1614, 0.0198, 0.0519, -0.0328],
         "polymer": True,
         "full_res_id": "A11"
        })
        self.assertEqual(len(model["connections"]), 60)
        self.assertEqual(
         model["connections"][0], {"atom": 3194, "bond_to": [3195, 3196]}
        )
        self.assertEqual(
         model["connections"][-1], {"atom": 3253, "bond_to": [3252]}
        )


    def test_5xme_pdb(self):
        d = atomium.open("tests/integration/files/5xme.pdb", data_dict=True)
        self.assertEqual(d["experiment"]["technique"], "SOLUTION NMR")
        self.assertEqual(d["quality"], {
         "resolution": None, "rvalue": None, "rfree": None
        })
        self.assertEqual(len(d["models"]), 10)
        for model in d["models"][1:]:
            self.assertEqual(model["chains"], d["models"][0]["chains"])
            self.assertEqual(model["residues"], d["models"][0]["residues"])
            self.assertEqual(model["ligands"], d["models"][0]["ligands"])
            self.assertEqual(len(model["atoms"]), len(d["models"][0]["atoms"]))
            self.assertNotEqual(model["atoms"], d["models"][0]["atoms"])
            self.assertEqual(model["connections"], d["models"][0]["connections"])


    def test_1cbn_pdb(self):
        d = atomium.open("tests/integration/files/1cbn.pdb", data_dict=True)
        atoms = d["models"][0]["atoms"]
        self.assertEqual(len(atoms), 777)
        self.assertEqual(atoms[0], {
         "id": 1,
         "element": "N",
         "name": "N",
         "x": 16.864,
         "y": 14.059,
         "z": 3.442,
         "bfactor": 6.22,
         "charge": 0,
         "residue_id": 1,
         "residue_name": "THR",
         "residue_insert": "",
         "chain_id": "A",
         "occupancy": 0.8,
         "alt_loc": "A",
         "anisotropy": [0, 0, 0, 0, 0, 0],
         "polymer": True,
         "full_res_id": "A1"
        })


    def test_1xda_pdb(self):
        d = atomium.open("tests/integration/files/1xda.pdb", data_dict=True)
        self.assertEqual(len(d["geometry"]["assemblies"]), 12)
        self.assertEqual(d["geometry"]["assemblies"][0], {
         "id": 1,
         "software": "PISA",
         "delta_energy": -7.0,
         "buried_surface_area": 1720.0,
         "surface_area": 3980.0,
         "transformations": [{
          "chains": ["A", "B"],
          "matrix": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
          "vector": [0.0, 0.0, 0.0]
         }]
        })
        self.assertEqual(d["geometry"]["assemblies"][4], {
         "id": 5,
         "software": "PISA",
         "delta_energy": -332.0,
         "buried_surface_area": 21680.0,
         "surface_area": 12240.0,
         "transformations": [{
          "chains": ["E", "F", "G", "H"],
          "matrix": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
          "vector": [0.0, 0.0, 0.0]
         }, {
          "chains": ["E", "F", "G", "H"],
          "matrix": [[-0.5, -0.866025, 0.0], [0.866025, -0.5, 0.0], [0.0, 0.0, 1.0]],
          "vector": [0.0, 0.0, 0.0]
         }, {
          "chains": ["E", "F", "G", "H"],
          "matrix": [[-0.5, 0.866025, 0.0], [-0.866025, -0.5, 0.0], [0.0, 0.0, 1.0]],
          "vector": [0.0, 0.0, 0.0]
         }]
        })


    def test_glucose_xyz(self):
        d = atomium.open("tests/integration/files/glucose.xyz", data_dict=True)
        self.assertEqual(d["description"]["title"], "glucose from 2gbp")
        self.assertEqual(len(d["models"]), 1)
        model = d["models"][0]
        for key in ["residues", "ligands", "chains", "connections"]:
            self.assertEqual(model[key], [])
        self.assertEqual(len(model["atoms"]), 12)
        self.assertEqual(model["atoms"][0], {
         "id": 0,
         "element": "C",
         "name": None,
         "x": 38.553,
         "y": 30.4,
         "z": 50.259,
         "bfactor": None,
         "charge": 0,
         "residue_id": None,
         "residue_name": None,
         "residue_insert": "",
         "chain_id": None,
         "occupancy": 1,
         "alt_loc": None,
         "anisotropy": [],
         "polymer": False,
         "full_res_id": None
        })


    def test_1lol_mmcif(self):
        d = atomium.open("tests/integration/files/1lol.cif", data_dict=True)
        self.assertEqual(set(d.keys()), {
         "description", "experiment", "quality", "geometry", "models"
        })
        self.assertEqual(d["description"], {
         "code": "1LOL",
         "title": "Crystal structure of orotidine monophosphate decarboxylase complex with XMP",
         "deposition_date": date(2002, 5, 6),
         "classification": "LYASE",
         "keywords": ["TIM barrel", "LYASE"],
         "authors": ["Wu, N.", "Pai, E.F."]
        })
        self.assertEqual(d["experiment"], {
         "technique": "X-RAY DIFFRACTION",
         "source_organism": "Methanothermobacter thermautotrophicus str. Delta H",
         "expression_system": "Escherichia coli"
        })
        self.assertEqual(d["quality"], {
         "resolution": 1.9, "rvalue": 0.193, "rfree": 0.229
        })
        self.assertEqual(d["geometry"], {"assemblies": [{
         "id": 1,
         "software": "PISA",
         "delta_energy": -31.0,
         "buried_surface_area": 5230,
         "surface_area": 16550,
         "transformations": [{
          "chains": ["A", "B", "C", "D", "E", "F", "G", "H"],
          "matrix": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
          "vector": [0.0, 0.0, 0.0]
         }]
        }]})

        self.assertEqual(len(d["models"]), 1)
        model = d["models"][0]
        self.assertEqual(len(model["chains"]), 2)
        self.assertEqual(model["chains"][0]["id"], "A")
        self.assertEqual(model["chains"][1]["id"], "B")
        for chain in model["chains"]:
            self.assertEqual(len(chain["full_sequence"]), 229)
        self.assertEqual(len(model["residues"]), 418)
        self.assertEqual(model["residues"][0], {
         "id": "A11", "name": "VAL", "chain_id": "A"
        })
        self.assertEqual(model["residues"][-1], {
         "id": "B1229", "name": "GLU", "chain_id": "B"
        })
        self.assertEqual(len(model["ligands"]), 184)
        self.assertEqual(model["ligands"][0], {
         "id": "A5001", "name": "BU2", "chain_id": "A"
        })
        self.assertEqual(model["ligands"][-1], {
         "id": "B3180", "name": "HOH", "chain_id": "B"
        })
        self.assertEqual(len(model["atoms"]), 3431)
        self.assertEqual(model["atoms"][0], {
         "id": 1,
         "element": "N",
         "name": "N",
         "x": 3.696,
         "y": 33.898,
         "z": 63.219,
         "bfactor": 21.5,
         "charge": -1.0,
         "residue_id": 11,
         "residue_name": "VAL",
         "residue_insert": "",
         "chain_id": "A",
         "occupancy": 1.0,
         "alt_loc": None,
         "anisotropy": [0.2406, 0.1892, 0.1614, 0.0198, 0.0519, -0.0328],
         "polymer": True,
         "full_res_id": "A11"
        })
        self.assertEqual(len(model["connections"]), 0)


    def test_5xme_mmcif(self):
        d = atomium.open("tests/integration/files/5xme.cif", data_dict=True)
        self.assertEqual(d["experiment"]["technique"], "SOLUTION NMR")
        self.assertEqual(d["quality"], {
         "resolution": None, "rvalue": None, "rfree": None
        })
        self.assertEqual(len(d["models"]), 10)
        for model in d["models"][1:]:
            self.assertEqual(model["chains"], d["models"][0]["chains"])
            self.assertEqual(model["residues"], d["models"][0]["residues"])
            self.assertEqual(model["ligands"], d["models"][0]["ligands"])
            self.assertEqual(len(model["atoms"]), len(d["models"][0]["atoms"]))
            self.assertNotEqual(model["atoms"], d["models"][0]["atoms"])
            self.assertEqual(model["connections"], d["models"][0]["connections"])


    def test_1cbn_mmcif(self):
        d = atomium.open("tests/integration/files/1cbn.cif", data_dict=True)
        atoms = d["models"][0]["atoms"]
        self.assertEqual(len(atoms), 777)
        self.assertEqual(atoms[0], {
         "id": 1,
         "element": "N",
         "name": "N",
         "x": 16.864,
         "y": 14.059,
         "z": 3.442,
         "bfactor": 6.22,
         "charge": 0,
         "residue_id": 1,
         "residue_name": "THR",
         "residue_insert": "",
         "chain_id": "A",
         "occupancy": 0.8,
         "alt_loc": "A",
         "anisotropy": [0, 0, 0, 0, 0, 0],
         "polymer": True,
         "full_res_id": "A1"
        })


    def test_1xda_mmcif(self):
        d = atomium.open("tests/integration/files/1xda.cif", data_dict=True)
        self.assertEqual(len(d["geometry"]["assemblies"]), 12)
        self.assertEqual(d["geometry"]["assemblies"][0], {
         "id": 1,
         "software": "PISA",
         "delta_energy": -7.0,
         "buried_surface_area": 1720.0,
         "surface_area": 3980.0,
         "transformations": [{
          "chains": ["A", "B", "I", "J", "K", "L", "Y", "Z"],
          "matrix": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
          "vector": [0.0, 0.0, 0.0]
         }]
        })
        self.assertEqual(d["geometry"]["assemblies"][4], {
         "id": 5,
         "software": "PISA",
         "delta_energy": -332.0,
         "buried_surface_area": 21680.0,
         "surface_area": 12240.0,
         "transformations": [{
          "chains": ["E", "F", "G", "H", "Q", "R", "S", "T", "U", "V", "W", "X", "CA", "DA", "EA", "FA"],
          "matrix": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
          "vector": [0.0, 0.0, 0.0]
         }, {
          "chains": ["E", "F", "G", "H", "Q", "R", "S", "T", "U", "V", "W", "X", "CA", "DA", "EA", "FA"],
          "matrix": [[-0.5, -0.8660254038, 0.0], [0.8660254038, -0.5, 0.0], [0.0, 0.0, 1.0]],
          "vector": [0.0, 0.0, 0.0]
         }, {
          "chains": ["E", "F", "G", "H", "Q", "R", "S", "T", "U", "V", "W", "X", "CA", "DA", "EA", "FA"],
          "matrix": [[-0.5, 0.8660254038, 0.0], [-0.8660254038, -0.5, 0.0], [0.0, 0.0, 1.0]],
          "vector": [0.0, 0.0, 0.0]
         }]
        })


    def test_1ej6_mmcif(self):
        d = atomium.open("tests/integration/files/1ej6.cif", data_dict=True)
        self.assertEqual(len(d["geometry"]["assemblies"]), 6)
        self.assertEqual(d["geometry"]["assemblies"][1], {
         "id": 2,
         "software": None,
         "delta_energy": None,
         "buried_surface_area": None,
         "surface_area": None,
         "transformations": [{
          "chains": ["A", "B", "C", "D", "E", "F"],
          "matrix": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
          "vector": [-0.0, 0.0, -0.0]
         }]
        })
        self.assertEqual(len(d["geometry"]["assemblies"][0]["transformations"]), 60)
        self.assertEqual(len(d["geometry"]["assemblies"][2]["transformations"]), 5)
        self.assertEqual(len(d["geometry"]["assemblies"][3]["transformations"]), 6)
        self.assertEqual(len(d["geometry"]["assemblies"][4]["transformations"]), 1)
        self.assertEqual(len(d["geometry"]["assemblies"][5]["transformations"]), 6)



class FileReadingTests(IntegratedTest):

    def test_1lol_pdb(self):
        pdb = atomium.open("tests/integration/files/1lol.pdb")
        self.assertEqual(pdb.filetype, "pdb")

        self.assertEqual(pdb.code, "1LOL")
        self.assertEqual(pdb.title, "CRYSTAL STRUCTURE OF OROTIDINE MONOPHOSPHATE DECARBOXYLASE COMPLEX WITH XMP")
        self.assertEqual(pdb.deposition_date, date(2002, 5, 6))
        self.assertEqual(pdb.classification, "LYASE")
        self.assertEqual(pdb.keywords, ["TIM BARREL", "LYASE"])
        self.assertEqual(pdb.authors, ["N.WU", "E.F.PAI"])

        self.assertEqual(pdb.technique, "X-RAY DIFFRACTION")
        self.assertEqual(pdb.source_organism, "METHANOTHERMOBACTER THERMAUTOTROPHICUS STR. DELTA H")
        self.assertEqual(pdb.expression_system, "ESCHERICHIA COLI")

        self.assertEqual(pdb.resolution, 1.9)
        self.assertEqual(pdb.rvalue, 0.193)
        self.assertEqual(pdb.rfree, 0.229)

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


    def test_5xme_pdb(self):
        pdb = atomium.open("tests/integration/files/5xme.pdb")
        self.assertEqual(pdb.resolution, None)
        self.assertEqual(pdb.assemblies, [{
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


    def test_1cbn_pdb(self):
        pdb = atomium.open("tests/integration/files/1cbn.pdb")
        chain = pdb.model.chain()
        residue1, residue2, residue3 = chain[:3]
        self.assertEqual(len(residue1.atoms()), 16)
        self.assertEqual(len(residue2.atoms()), 14)
        self.assertEqual(len(residue3.atoms()), 10)
        for residue in chain[:3]:
            for name in ["N", "C", "CA", "CB"]:
                self.assertEqual(len(residue.atoms(name=name)), 1)


    def test_1xda_pdb(self):
        pdb = atomium.open("tests/integration/files/1xda.pdb")
        self.assertEqual(len(pdb.model.chains()), 8)
        self.assertEqual(len(pdb.model.residues()), 200)
        self.assertEqual(len(pdb.model.ligands()), 170)
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


    def test_glucose_xyz(self):
        xyz = atomium.open("tests/integration/files/glucose.xyz")
        self.assertEqual(xyz.title, "glucose from 2gbp")
        model = xyz.model
        self.assertEqual(len(model.atoms()), 12)
        self.assertEqual(len(model.atoms(element="C")), 6)
        self.assertEqual(len(model.atoms(element="O")), 6)
        self.assertEqual(len(model.atoms(element_regex="[CO]")), 12)
        self.assertAlmostEqual(model.mass, 168, delta=0.5)
        for atom in model.atoms():
            self.assertAlmostEqual(
             atom.mass, {"C": 12, "O": 16}[atom.element], delta=0.2
            )
            self.assertIs(atom.model, model)
        self.assertFalse(model.chains())
        self.assertFalse(model.ligands())
        self.assertFalse(model.residues())


    def test_1lol_mmcif(self):
        pdb = atomium.open("tests/integration/files/1lol.cif")
        self.assertEqual(pdb.filetype, "cif")
        self.assertEqual(pdb.title, "Crystal structure of orotidine monophosphate decarboxylase complex with XMP")
        self.assertEqual(pdb.deposition_date, date(2002, 5, 6))
        self.assertEqual(pdb.technique, "X-RAY DIFFRACTION")
        self.assertEqual(pdb.resolution, 1.9)
        self.assertEqual(pdb.rvalue, 0.193)
        self.assertEqual(pdb.rfree, 0.229)
        model = pdb.model
        self.assertEqual(len(model.atoms()), 3431)
        atom = model.atom(2934)
        self.assertEqual(atom.anisotropy, [0, 0, 0, 0, 0, 0])
        self.assertEqual(
         model.atom(1).anisotropy,
         [0.2406, 0.1892, 0.1614, 0.0198, 0.0519, -0.0328]
        )
        self.assertEqual(atom.element, "C")
        self.assertEqual(atom.name, "CZ")
        self.assertAlmostEqual(
         model.mass, 46018.5, delta=0.005
        )
        self.assertEqual(atom.model, model)


    def test_5xme_mmcif(self):
        pdb = atomium.open("tests/integration/files/5xme.cif")
        models = pdb.models
        self.assertEqual(len(models), 10)
        self.assertIs(pdb.model, pdb.models[0])
        all_atoms = set()
        for model in models:
            self.assertEqual(len(model.atoms()), 1827)
            all_atoms.update(model.atoms())
        self.assertEqual(len(all_atoms), 18270)


    def test_1cbn_mmcif(self):
        pdb = atomium.open("tests/integration/files/1cbn.cif")
        chain = pdb.model.chain()
        residue1, residue2, residue3 = chain[:3]
        self.assertEqual(len(residue1.atoms()), 16)
        self.assertEqual(len(residue2.atoms()), 14)
        self.assertEqual(len(residue3.atoms()), 10)
        for residue in chain[:3]:
            for name in ["N", "C", "CA", "CB"]:
                self.assertEqual(len(residue.atoms(name=name)), 1)


    def test_1xda_mmcif(self):
        pdb = atomium.open("tests/integration/files/1xda.cif")
        self.assertEqual(len(pdb.model.chains()), 8)
        self.assertEqual(len(pdb.model.residues()), 200)
        self.assertEqual(len(pdb.model.ligands()), 170)
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
