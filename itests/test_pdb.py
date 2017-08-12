import atomium
from base import IntegratedTest

class PdbReadingTests(IntegratedTest):

    def test_can_read_pdb_file(self):
        pdb = atomium.pdb_from_file("itests/files/1lol.pdb")

        # Atoms are correct
        model = pdb.model()
        self.assertEqual(len(model.atoms()), 3431)
        atom = model.atom(atom_id=2934)
        self.assertEqual(atom.element(), "N")
        self.assertEqual(atom.name(), "NE")
        self.assertEqual(atom.location(), (-20.082, 79.647, 41.645))
        self.assertAlmostEqual(
         model.mass(), 2 * 24994.8 + 2 * 90.1 + 2 * 365.2 + 180 * 18, delta=0.005
        )

        # Chains are correct
        self.assertEqual(len(model.chains()), 2)
        for chain in model.chains():
            self.assertIs(chain.model(), model)
            self.assertIsNone(chain.name())
            self.assertEqual(chain.length(), 229)
            self.assertAlmostEqual(chain.mass(), 24994.8, delta=0.005)
        chaina, chainb = model.chain(chain_id="A", chain_id="B")

        # Residues are correct
        self.assertEqual(chaina[0].name(), "LEU")
        self.assertEqual(chaina[0].next().name(), "ARG")
        self.assertEqual(chaina[-1].name(), "GLU")
        self.assertEqual(chaina[-1].previous().name(), "PRO")
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
        for atom in residue:
            self.assertIs(atom.residue(), residue)
            self.assertIs(atom.molecule(), residue)
            self.assertIs(atom.chain(), chain)
            self.assertIs(atom.model(), model)

        # Molecules are correct
        self.assertEqual(len(model.molecules()), 186)
        self.assertEqual(len(model.molecules(filter=True)), 180)
        self.assertEqual(len(model.molecules(water=False)), 6)
        self.assertEqual(len(model.molecules(water=False, filter=True)), 4)
        self.assertEqual(len(model.molecules(molecule_name="XMP")), 2)
        self.assertEqual(len(model.molecules(molecule_name="BU2")), 2)
        self.assertEqual(len(model.molecules(molecule_name="HOH")), 180)
        mol = model.molecules(molecule_id="A2001")
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
        next_atom = chain[1].atom(name="N")
        self.assertEqual(n.bonded_atoms(), set([ca]))
        self.assertEqual(ca.bonded_atoms(), set([n, c, cb]))
        self.assertEqual(c.bonded_atoms(), set([ca, o, next_atom]))
        self.assertEqual(o.bonded_atoms(), set([c]))
        self.assertEqual(cb.bonded_atoms(), set([ca, cg1, cg2]))
        self.assertEqual(cg1.bonded_atoms(), set([cb]))
        self.assertEqual(cg2.bonded_atoms(), set([cb]))
        bond1 = n.bond_with(ca)
        bond2 = ca.bond_with(c)
        self.assertAlmostEqual(bond1.angle_with(bond2), 109.474, delta=0.005)
        self.assertIn(
         model.atom(atom_id=48), model.atom(atom_id=1082).bonded_atoms()
        )



class PdbFetchingTests(IntegratedTest):

    def test_can_fetch_file(self):
        pdb = atomium.fetch("1h4W")
        model = pdb.model()

        self.assertEqual(len(model.chains()), 1)
        residue = model.residue(residue_id="A221A")
        self.assertEqual(residue.next().residue_id(), "A222")
        self.assertEqual(residue.previous().residue_id(), "A221")
