import atomium
import os
from unittest import TestCase

class SavingTest(TestCase):

    def setUp(self):
        self.files_at_start = os.listdir("tests/integration/files")


    def tearDown(self):
        files_at_end = os.listdir("tests/integration/files")
        to_remove = [f for f in files_at_end if f not in self.files_at_start]
        for f in to_remove:
            os.remove("tests/integration/files/" + f)


    def check_file_saving(self, filename):
        f = atomium.open("tests/integration/files/" + filename)
        f.model.save("tests/integration/files/saved_" + filename)
        f2 = atomium.open("tests/integration/files/saved_" + filename)
        self.assertEqual(f.model, f2.model)
        self.assertEqual(len(f.model.chains()), len(f2.model.chains()))
        for chain1, chain2 in zip(sorted(f.model.chains(), key=lambda c: c.id),
         sorted(f2.model.chains(), key=lambda c: c.id)):
            self.assertEqual(chain1.sequence, chain2.sequence)
            self.assertEqual(chain1.id, chain2.id)
            self.assertEqual(chain1, chain2)
            for res1, res2 in zip(sorted(chain1.residues(), key=lambda r: r.id),
             sorted(chain2.residues(), key=lambda r: r.id)):
                self.assertEqual(res1.id, res2.id)
                self.assertEqual(res1.name, res2.name)
                self.assertEqual(res1, res2)
        for lig1, lig2 in zip(sorted(f.model.ligands(), key=lambda c: c.id),
         sorted(f2.model.ligands(), key=lambda c: c.id)):
            self.assertEqual(lig1.name, lig2.name)
            self.assertEqual(lig1.id, lig2.id)
            self.assertEqual(lig1, lig2)



class MmcifFileSavingTests(SavingTest):

    def test_can_save_1lol(self):
        self.check_file_saving("1lol.cif")


    def test_can_save_1cbn(self):
        self.check_file_saving("1cbn.cif")


    def test_can_save_1m4x(self):
        self.check_file_saving("1m4x.cif")


    def test_can_save_1xda(self):
        self.check_file_saving("1xda.cif")


    def test_can_save_5xme(self):
        self.check_file_saving("5xme.cif")


    def test_can_save_4y60(self):
        self.check_file_saving("4y60.cif")


    def test_chain(self):
        f = atomium.open("tests/integration/files/1lol.cif")
        f.model.chain("A").save("tests/integration/files/chaina.cif")
        chain = atomium.open("tests/integration/files/chaina.cif").model
        self.assertEqual(f.model.chain("A"), chain)


    def test_biological_assembly_warns_on_saving(self):
        f = atomium.open("tests/integration/files/1xda.cif")
        model = f.generate_assembly(5)
        with self.assertWarns(Warning):
            model.save("tests/integration/files/assembly.cif")



class MmtfFileSavingTests(SavingTest):

    def test_can_save_1lol(self):
        self.check_file_saving("1lol.mmtf")


    def test_can_save_1cbn(self):
        self.check_file_saving("1cbn.mmtf")


    def test_can_save_1m4x(self):
        self.check_file_saving("1m4x.mmtf")


    def test_can_save_1xda(self):
        self.check_file_saving("1xda.mmtf")


    def test_can_save_5xme(self):
        self.check_file_saving("5xme.mmtf")


    def test_can_save_4y60(self):
        self.check_file_saving("4y60.mmtf")


    def test_chain(self):
        f = atomium.open("tests/integration/files/1lol.mmtf")
        f.model.chain("A").save("tests/integration/files/chaina.mmtf")
        chain = atomium.open("tests/integration/files/chaina.mmtf").model
        self.assertEqual(f.model.chain("A"), chain)


    def test_biological_assembly_warns_on_saving(self):
        f = atomium.open("tests/integration/files/1xda.cif")
        model = f.generate_assembly(5)
        with self.assertWarns(Warning):
            model.save("tests/integration/files/assembly.cif")



class PdbFileSavingTests(SavingTest):

    def test_can_save_1lol(self):
        self.check_file_saving("1lol.pdb")


    def test_can_save_1cbn(self):
        self.check_file_saving("1cbn.pdb")


    def test_can_save_1m4x(self):
        self.check_file_saving("1m4x.pdb")


    def test_can_save_1xda(self):
        self.check_file_saving("1xda.pdb")


    def test_can_save_5xme(self):
        self.check_file_saving("5xme.pdb")


    def test_can_save_4y60(self):
        self.check_file_saving("4y60.pdb")
    

    def test_can_save_1d5t(self):
        self.check_file_saving("1d5t.pdb")
    

    def test_can_save_1grm(self):
        self.check_file_saving("1grm.pdb")
        f = atomium.open("tests/integration/files/1grm.pdb")
        f.model.save("tests/integration/files/saved_1grm.pdb")
        with open("tests/integration/files/1grm.pdb") as f:
            old_text = f.read()
            old_text = old_text[:old_text.find("ENDMDL")]
            old_remark_count = old_text.count("HETATM")
        with open("tests/integration/files/saved_1grm.pdb") as f:
            new_remark_count = f.read().count("HETATM")
        self.assertEqual(old_remark_count, new_remark_count)


    def test_chain(self):
        f = atomium.open("tests/integration/files/1lol.pdb")
        f.model.chain("A").save("tests/integration/files/chaina.pdb")
        chain = atomium.open("tests/integration/files/chaina.pdb").model
        self.assertEqual(f.model.chain("A"), chain)


    def test_biological_assembly_warns_on_saving(self):
        f = atomium.open("tests/integration/files/1xda.pdb")
        model = f.generate_assembly(5)
        with self.assertWarns(Warning):
            model.save("tests/integration/files/assembly.pdb")
