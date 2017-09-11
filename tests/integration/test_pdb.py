import atomium
from tests.integration.base import IntegratedTest

class PdbReadingTests(IntegratedTest):

    def test_can_read_pdb_file(self):
        pdb = atomium.pdb_from_file("tests/integration/files/1lol.pdb")

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
        self.assertEqual(len(residue.atoms(element="O")), 2)
        for atom in residue.atoms():
            self.assertIs(atom.residue(), residue)
            self.assertIs(atom.molecule(), chaina)
            self.assertIs(atom.chain(), chaina)
            self.assertIs(atom.model(), model)

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

        # Can be saved
        pdb.save("tests/integration/files/1LOL2.pdb")
        with open("tests/integration/files/1LOL2.pdb") as f:
            new = [l.strip() for l in f.readlines() if l.strip()]
        with open("tests/integration/files/1lol_output.pdb") as f:
            ref = [l.strip() for l in f.readlines() if l.strip()]
        self.assertEqual(new, ref)
        new = atomium.pdb_from_file("tests/integration/files/1LOL2.pdb")
        model = new.model()
        self.assertAlmostEqual(
         model.mass(), 46018.5, delta=0.005
        )


    def test_can_read_multi_model_pdbs(self):
        pdb = atomium.pdb_from_file("tests/integration/files/5xme.pdb")

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

        pdb.save("tests/integration/files/5XME2.pdb")
        with open("tests/integration/files/5XME2.pdb") as f:
            new = [l.strip() for l in f.readlines() if l.strip()]
        with open("tests/integration/files/5xme_output.pdb") as f:
            ref = [l.strip() for l in f.readlines() if l.strip()]
        self.assertEqual(new, ref)
        new = atomium.pdb_from_file("tests/integration/files/5XME2.pdb")
        models = new.models()
        self.assertEqual(len(models), 10)
        for x, model in zip(x_values, models):
            self.assertEqual(len(model.atoms()), 1827)





class PdbFetchingTests(IntegratedTest):

    def test_can_fetch_file(self):
        pdb = atomium.fetch("1h4W", pdbe=True)
        model = pdb.model()

        self.assertEqual(len(model.chains()), 1)
        residue = model.residue(residue_id="A221A")
        self.assertEqual(residue.next().residue_id(), "A222")
        self.assertEqual(residue.previous().residue_id(), "A221")
