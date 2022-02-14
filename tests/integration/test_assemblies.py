from unittest import TestCase
import atomium

class SimpleAssemblyGenerationTests(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.pdb = atomium.open("tests/integration/files/1xda.cif")


    def test_assembly_1(self):
        model = self.pdb.generate_assembly(1)
        self.assertEqual(len(model.polymers()), 2)
        self.assertEqual(set([p.id for p in model.polymers()]), {"A", "B"})
        self.assertEqual(len(model.non_polymers()), 4)


    def test_assembly_2(self):
        model = self.pdb.generate_assembly(2)
        self.assertEqual(len(model.polymers()), 2)
        self.assertEqual(set([p.id for p in model.polymers()]), {"C", "D"})
        self.assertEqual(len(model.non_polymers()), 4)


    def test_assembly_3(self):
        model = self.pdb.generate_assembly(3)
        self.assertEqual(len(model.polymers()), 2)
        self.assertEqual(set([p.id for p in model.polymers()]), {"E", "F"})
        self.assertEqual(len(model.non_polymers()), 4)
    

    def test_assembly_4(self):
        model = self.pdb.generate_assembly(4)
        self.assertEqual(len(model.polymers()), 2)
        self.assertEqual(set([p.id for p in model.polymers()]), {"G", "H"})
        self.assertEqual(len(model.non_polymers()), 4)
    

    def test_assembly_5(self):
        model = self.pdb.generate_assembly(1)
        self.assertEqual(len(model.polymers()), 2)
        self.assertEqual(set([p.id for p in model.polymers()]), {"A", "B"})
        self.assertEqual(len(model.non_polymers()), 4)
    

    def test_assembly_6(self):
        model = self.pdb.generate_assembly(1)
        self.assertEqual(len(model.polymers()), 2)
        self.assertEqual(set([p.id for p in model.polymers()]), {"A", "B"})
        self.assertEqual(len(model.non_polymers()), 4)
    

    def test_assembly_7(self):
        model = self.pdb.generate_assembly(7)
        self.assertEqual(len(model.polymers()), 6)
        self.assertEqual(set([p.id for p in model.polymers()]), {"A", "B"})
        self.assertEqual(len(model.non_polymers()), 12)
        zn = model.atom(element="ZN")
        liganding_residues = zn.nearby_residues(3, is_metal=False, element__ne="CL")
        self.assertEqual(len(liganding_residues), 3)
        self.assertEqual(set([r.id for r in liganding_residues]), {"B.10"})
        self.assertEqual(set([r.name for r in liganding_residues]), {"HIS"})
        res1, res2, res3 = liganding_residues
        self.assertGreater(res1.atom(name="N").distance_to(res2.atom(name="N")), 10)