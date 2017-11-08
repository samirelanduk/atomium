from datetime import datetime
import atomium
from tests.integration.base import IntegratedTest

class PdbReadingTests(IntegratedTest):

    def test_can_read_pdb(self):
        pdb = atomium.pdb_from_file("tests/integration/files/1lol.pdb")
        self.assertEqual(pdb.code(), "1LOL")
        self.assertEqual(
         pdb.title(),
         "CRYSTAL STRUCTURE OF OROTIDINE MONOPHOSPHATE DECARBOXYLASE COMPLEX WITH XMP"
        )
        self.assertEqual(pdb.deposition_date(), datetime(2002, 5, 6).date())
        self.assertEqual(pdb.resolution(), 1.9)
        self.assertEqual(pdb.organism(), "METHANOTHERMOBACTER THERMAUTOTROPHICUS STR. DELTA H")
        self.assertEqual(pdb.expression_system(), "ESCHERICHIA COLI")

        # Atoms are correct
        model = pdb.model()
        self.assertEqual(len(model.atoms()), 3431)
        atom = model.atom(atom_id=2934)
        self.assertEqual(atom.element(), "N")
        self.assertEqual(atom.name(), "NE")
        self.assertEqual(atom.location(), (-20.082, 79.647, 41.645))
        self.assertEqual(atom.bfactor(), 35.46)
        self.assertAlmostEqual(
         model.mass(), 46018.5, delta=0.005
        )

        # Chains are correct
        self.assertEqual(len(model.chains()), 2)
        for chain in model.chains():
            self.assertIs(chain.model(), model)
            self.assertIsNone(chain.name())
        chaina, chainb = model.chain(chain_id="A"), model.chain(chain_id="B")

        # Residues are correct
        self.assertEqual(chaina[0].name(), "VAL")
        self.assertEqual(chaina[0].next().name(), "MET")
        self.assertEqual(chaina[-1].name(), "ILE")
        self.assertEqual(chaina[-1].previous().name(), "SER")
        self.assertEqual(len(chaina.residues(name="ASN")), 6)
        for residue in chaina:
            self.assertIs(residue.model(), model)
            self.assertIs(residue.chain(), chaina)
        for residue in chainb:
            self.assertIs(residue.model(), model)
            self.assertIs(residue.chain(), chainb)
        residue = chaina.residue(name="ASN")
        self.assertEqual(residue.residue_id(), "A13")
        self.assertEqual(len(residue.atoms()), 8)
        self.assertEqual(len(residue.atoms(exclude="O")), 6)
        self.assertEqual(len(residue.atoms(element="O")), 2)
        for atom in residue.atoms():
            self.assertIs(atom.residue(), residue)
            self.assertIs(atom.molecule(), chaina)
            self.assertIs(atom.chain(), chaina)
            self.assertIs(atom.model(), model)
        gly = chaina.residue(residue_id="A32")
        pairs = list(gly.pairwise_atoms())
        self.assertEqual(len(pairs), 6)
        for pair in pairs:
            self.assertTrue(0 < pair[0].distance_to(pair[1]), 5)

        # Molecules are correct
        self.assertEqual(len(model.molecules()), 186)
        self.assertEqual(len(model.molecules(generic=True)), 184)
        self.assertEqual(len(model.molecules(water=False)), 6)
        self.assertEqual(len(model.molecules(water=False, generic=True)), 4)
        self.assertEqual(len(model.molecules(name="XMP")), 2)
        self.assertEqual(len(model.molecules(name="BU2")), 2)
        self.assertEqual(len(model.molecules(name="HOH")), 180)
        mol = model.molecule(molecule_id="A2001")
        self.assertIs(mol.model(), model)
        self.assertEqual(mol.name(), "XMP")
        self.assertEqual(len(mol.atoms()), 24)

        # Bindsites are correct
        site = model.molecule("A5001").site()
        self.assertIs(site.ligand(), model.molecule("A5001"))
        self.assertEqual(site.residues(), set([
         model.residue("A42"), model.residue("A70"), model.residue("A72"),
         model.residue("A96"), model.residue("A123"), model.residue("A155")
        ]))

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
        self.assertEqual(n.bonded_atoms(), set([ca]))
        self.assertEqual(len(n.bonds()), 1)
        self.assertEqual(ca.bonded_atoms(), set([n, c, cb]))
        self.assertEqual(c.bonded_atoms(), set([ca, o, next_atom]))
        self.assertEqual(o.bonded_atoms(), set([c]))
        self.assertEqual(cb.bonded_atoms(), set([ca, cg1, cg2]))
        self.assertEqual(cg1.bonded_atoms(), set([cb]))
        self.assertEqual(cg2.bonded_atoms(), set([cb]))
        res181 = chaina.residue(residue_id="A181")
        res190 = chaina.residue(residue_id="A190")
        c2, n2 = res181.atom(name="C"), res190.atom(name="N")
        self.assertNotIn(c2, n2.bonded_atoms())
        self.assertNotIn(n2, c2.bonded_atoms())
        self.assertIn(
         model.atom(atom_id=3194), model.atom(atom_id=3195).bonded_atoms()
        )
        bond1 = n.bond_with(ca)
        bond2 = ca.bond_with(c)
        self.assertAlmostEqual(
         bond1.angle_with(bond2, degrees=True), 109.474, delta=0.005
        )

        # Can get atoms in cutoff distance
        atom = model.atom(1587)
        four_angstrom = atom.nearby(cutoff=4)
        self.assertEqual(len(four_angstrom), 10)
        self.assertEqual(
         sorted([atom.atom_id() for atom in four_angstrom]),
         [1576, 1582, 1583, 1584, 1586, 1588, 1589, 1590, 1591, 2957]
        )
        self.assertEqual(len(atom.nearby(cutoff=4, element="O")), 1)
        self.assertEqual(len(atom.nearby(cutoff=4, exclude="C")), 4)


    def test_can_read_multi_model_pdbs(self):
        pdb = atomium.pdb_from_file("tests/integration/files/5xme.pdb")
        self.assertEqual(pdb.resolution(), 0)

        models = pdb.models()
        self.assertEqual(len(models), 10)
        self.assertIs(pdb.model(), pdb.models()[0])

        x_values = [
         33.969, 34.064, 37.369, 36.023, 35.245,
         35.835, 37.525, 35.062, 36.244, 37.677
        ]
        all_atoms = set()
        for x, model in zip(x_values, models):
            self.assertEqual(len(model.atoms()), 1827)
            all_atoms.update(model.atoms())
            atom = model.atom(1)
            self.assertEqual(atom.x(), x)
            self.assertEqual(len(atom.bonded_atoms()), 1)
        self.assertEqual(len(all_atoms), 18270)


    def test_can_read_alt_loc_pdbs(self):
        pdb = atomium.pdb_from_file("tests/integration/files/1cbn.pdb")
        chain = pdb.model().chain()
        residue1, residue2, residue3 = chain[:3]

        # Residues have the correct number of atoms
        self.assertEqual(len(residue1.atoms()), 16)
        self.assertEqual(len(residue2.atoms()), 14)
        self.assertEqual(len(residue3.atoms()), 10)

        for residue in chain[:3]:
            for name in ["N", "C", "CA", "CB"]:
                self.assertEqual(len(residue.atoms(name=name)), 1)


    def test_can_read_pdb_data(self):
        data_file = atomium.pdb_data_from_file(
         "tests/integration/files/1lol.pdb"
        )
        self.assertEqual(data_file["code"], "1LOL")
        self.assertEqual(len(data_file["models"][0]["chains"]), 2)
        self.assertEqual(
         data_file["models"][0]["chains"][0]["residues"][0]["atoms"][0],
         {
          "atom_id": 1, "atom_name": "N", "alt_loc": None,
          "residue_name": "VAL", "full_id": "A11",
          "chain_id": "A", "residue_id": 11, "insert_code": "",
          "x": 3.696, "y": 33.898, "z": 63.219,
          "occupancy": 1.0, "temp_factor": 21.50,
          "element": "N", "charge": 0.0,
         }
        )
        self.assertEqual(len(data_file["connections"]), 60)
        self.assertEqual(data_file["connections"][0], {
         "atom": 3194, "bond_to": [3195, 3196]
        })


    def test_can_fetch_pdb(self):
        pdb = atomium.fetch("1h4W", pdbe=True)
        model = pdb.model()

        self.assertEqual(len(model.chains()), 1)
        residue = model.residue(residue_id="A221A")
        self.assertEqual(residue.next().residue_id(), "A222")
        self.assertEqual(residue.previous().residue_id(), "A221")


    def test_can_fetch_pdb_data(self):
        data_file = atomium.fetch_data("1lol", pdbe=True)
        self.assertEqual(data_file["code"], "1LOL")
        self.assertEqual(len(data_file["models"][0]["chains"]), 2)
        self.assertEqual(
         data_file["models"][0]["chains"][0]["residues"][0]["atoms"][0],
         {
          "atom_id": 1, "atom_name": "N", "alt_loc": None,
          "residue_name": "VAL", "full_id": "A11",
          "chain_id": "A", "residue_id": 11, "insert_code": "",
          "x": 3.696, "y": 33.898, "z": 63.219,
          "occupancy": 1.0, "temp_factor": 21.50,
          "element": "N", "charge": 0.0,
         }
        )
        self.assertEqual(len(data_file["connections"]), 60)
        self.assertEqual(data_file["connections"][0], {
         "atom": 3194, "bond_to": [3195, 3196]
        })



class PdbSavingTests(IntegratedTest):

    def test_can_save_single_model_pdb(self):
        pdb = atomium.pdb_from_file("tests/integration/files/1lol.pdb")

        # Can be saved
        pdb.code("9SAM")
        pdb.deposition_date(datetime(1990, 9, 1).date())
        pdb.title(
         "FOR IN THAT SLEEP OF DEATH WHAT DREAMS MAY COME WHEN WE HAVE " +
         "SHUFFLED OFF THIS MORTAL COIL MUST GIVE US PAUSE"
        )
        pdb.save("tests/integration/files/1LOL2.pdb")
        with open("tests/integration/files/1LOL2.pdb") as f:
            new = [l.strip() for l in f.readlines() if l.strip()]
        with open("tests/integration/files/1lol_output.pdb") as f:
            ref = [l.strip() for l in f.readlines() if l.strip()]
        for new_line, ref_line in zip(new[3:], ref[3:]):
            self.assertEqual(new_line, ref_line)
        self.assertIn("9SAM", new[0])
        self.assertIn("90", new[0])
        self.assertIn("DREAMS", new[1])
        self.assertIn("PAUSE", new[2])
        new = atomium.pdb_from_file("tests/integration/files/1LOL2.pdb")
        model = new.model()
        self.assertAlmostEqual(
         model.mass(), 46018.5, delta=0.005
        )


    def test_can_save_structures(self):
        pdb = atomium.pdb_from_file("tests/integration/files/1lol.pdb")

        # Save chains
        for chain in pdb.model().chains():
            chain.save("tests/integration/files/chain{}.pdb".format(chain.chain_id()))

        with open("tests/integration/files/chainA.pdb") as f:
            new = [l.strip() for l in f.readlines() if l.strip()]
        with open("tests/integration/files/chaina_output.pdb") as f:
            ref = [l.strip() for l in f.readlines() if l.strip()]
        for new_line, ref_line in zip(new, ref):
            self.assertEqual(new_line, ref_line)
        new = atomium.pdb_from_file("tests/integration/files/chainA.pdb")
        model = new.model()
        self.assertEqual(len(model.chains()), 1)

        with open("tests/integration/files/chainB.pdb") as f:
            new = [l.strip() for l in f.readlines() if l.strip()]
        with open("tests/integration/files/chainb_output.pdb") as f:
            ref = [l.strip() for l in f.readlines() if l.strip()]
        for new_line, ref_line in zip(new, ref):
            self.assertEqual(new_line, ref_line)
        new = atomium.pdb_from_file("tests/integration/files/chainB.pdb")
        model = new.model()
        self.assertEqual(len(model.chains()), 1)

        # Save molecules
        pdb.model().molecule("A5001").save("tests/integration/files/5001.pdb")
        with open("tests/integration/files/5001.pdb") as f:
            new = [l.strip() for l in f.readlines() if l.strip()]
        with open("tests/integration/files/5001_output.pdb") as f:
            ref = [l.strip() for l in f.readlines() if l.strip()]
        for new_line, ref_line in zip(new, ref):
            self.assertEqual(new_line, ref_line)
        new = atomium.pdb_from_file("tests/integration/files/5001.pdb")
        model = new.model()
        self.assertEqual(len(model.atoms()), 6)


    def test_can_save_multi_model_pdb(self):
        pdb = atomium.pdb_from_file("tests/integration/files/5xme.pdb")

        pdb.save("tests/integration/files/5XME2.pdb")
        with open("tests/integration/files/5XME2.pdb") as f:
            new = [l.strip() for l in f.readlines() if l.strip()]
        with open("tests/integration/files/5xme_output.pdb") as f:
            ref = [l.strip() for l in f.readlines() if l.strip()]
        for new_line, ref_line in zip(new, ref):
            self.assertEqual(new_line, ref_line)
        new = atomium.pdb_from_file("tests/integration/files/5XME2.pdb")
        models = new.models()
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
        with open("tests/integration/files/1CBN2.pdb") as f:
            new = [l.strip() for l in f.readlines() if l.strip()]
        with open("tests/integration/files/1cbn_output.pdb") as f:
            ref = [l.strip() for l in f.readlines() if l.strip()]
        for new_line, ref_line in zip(new, ref):
            self.assertEqual(new_line, ref_line)

        new = atomium.pdb_from_file("tests/integration/files/1CBN2.pdb")
        chain = pdb.model().chain()
        residue1, residue2, residue3 = chain[:3]
        self.assertEqual(len(residue1.atoms()), 16)
        self.assertEqual(len(residue2.atoms()), 14)
        self.assertEqual(len(residue3.atoms()), 10)
