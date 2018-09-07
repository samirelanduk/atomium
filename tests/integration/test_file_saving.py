from datetime import datetime
import atomium
from tests.integration.base import IntegratedTest

class PdbSavingTests(IntegratedTest):

    def test_1lol_pdb(self):
        pdb = atomium.open("tests/integration/files/1lol.pdb")

        # Can be saved
        pdb.code = "9SAM"
        pdb.deposition_date = datetime(1990, 9, 1).date()
        pdb.title = (
         "FOR IN THAT SLEEP OF DEATH, WHAT DREAMS MAY COME, WHEN WE HAVE " +
         "SHUFFLED OFF THIS MORTAL COIL MUST GIVE US PAUSE"
        )
        pdb.source_organism = "HOMO SAPIENS"
        pdb.expression_system = "COOL NEW ORGANISM"
        pdb.technique = "MEME DIFFRACTION"
        pdb.model.chain("B").rep_sequence = ""
        while pdb.keywords: pdb.keywords.pop()
        for k in ["AMAZING", "SUPERB", "WOW"]: pdb.keywords.append(k)
        pdb.save("tests/integration/files/1LOL2.pdb")
        self.check_files_the_same("1LOL2.pdb", "1lol_output.pdb")

        new = atomium.open("tests/integration/files/1LOL2.pdb")
        model = new.model
        self.assertAlmostEqual(
         model.mass, 46018.5, delta=0.005
        )


    def test_can_save_multi_model_pdb(self):
        pdb = atomium.open("tests/integration/files/5xme.pdb")
        while pdb.keywords: pdb.keywords.pop()
        for k in ["INTEGRAL"] * 10: pdb.keywords.append(k)

        pdb.save("tests/integration/files/5XME2.pdb")
        self.check_files_the_same("5XME2.pdb", "5xme_output.pdb")
        new = atomium.open("tests/integration/files/5XME2.pdb")
        models = new.models
        self.assertEqual(len(models), 10)
        x_values = [
         33.969, 34.064, 37.369, 36.023, 35.245,
         35.835, 37.525, 35.062, 36.244, 37.677
        ]
        for x, model in zip(x_values, models):
            self.assertEqual(len(model.atoms()), 1827)


    def test_can_save_alt_loc_pdbs(self):
        pdb = atomium.open("tests/integration/files/1cbn.pdb")

        pdb.save("tests/integration/files/1CBN2.pdb")
        self.check_files_the_same("1CBN2.pdb", "1cbn_output.pdb")

        new = atomium.open("tests/integration/files/1CBN2.pdb")
        chain = pdb.model.chain()
        residue1, residue2, residue3 = chain[:3]
        self.assertEqual(len(residue1.atoms()), 16)
        self.assertEqual(len(residue2.atoms()), 14)
        self.assertEqual(len(residue3.atoms()), 10)


    def test_can_save_multi_assembly_pdbs(self):
        pdb = atomium.open("tests/integration/files/1xda.pdb")

        pdb.save("tests/integration/files/1XDA2.pdb")
        self.check_files_the_same("1XDA2.pdb", "1xda_output.pdb")

        new = atomium.open("tests/integration/files/1XDA2.pdb")



class XyzSavingTests(IntegratedTest):

    def test_can_save_xyz_file(self):
        # Open and save file
        xyz = atomium.open("tests/integration/files/glucose.xyz")
        xyz.save("tests/integration/files/glucose2.xyz")

        # The saved xyz is correct
        with open("tests/integration/files/glucose2.xyz") as f:
            new = [l.strip() for l in f.readlines()]
        with open("tests/integration/files/glucose.xyz") as f:
            old = [l.strip() for l in f.readlines()]
        self.assertEqual(old[:-12], new[:-12])
        self.assertEqual(set(old[-12:]), set(new[-12:]))
        new = atomium.open("tests/integration/files/glucose2.xyz")
        model = xyz.model
        self.assertAlmostEqual(model.mass, 168, delta=0.5)



class StructureSavingTests(IntegratedTest):

    def test_can_save_structures_to_pdb(self):
        pdb = atomium.open("tests/integration/files/1lol.pdb")

        # Save chains
        for chain in pdb.model.chains():
            chain.save("tests/integration/files/chain{}.pdb".format(chain.id))

        self.check_files_the_same("chainA.pdb", "chaina_output.pdb")
        self.check_files_the_same("chainB.pdb", "chainb_output.pdb")
        new = atomium.open("tests/integration/files/chainA.pdb")
        model = new.model
        self.assertEqual(len(model.chains()), 1)

        # Save molecules
        pdb.model.ligand("A:5001").save("tests/integration/files/5001.pdb")
        self.check_files_the_same("5001.pdb", "5001_output.pdb")
        new = atomium.open("tests/integration/files/5001.pdb")
        model = new.model
        self.assertEqual(len(model.atoms()), 6)
