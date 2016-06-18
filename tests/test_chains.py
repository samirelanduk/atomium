import math
from unittest import TestCase
from omnicanvas.canvas import Canvas
import omnicanvas.graphics
from molecupy import exceptions
from molecupy.structures import PdbChain, ResiduicSequence, PdbAtom, PdbResidue
from molecupy.structures import PdbAlphaHelix, PdbBetaStrand

class ChainTest(TestCase):

    def setUp(self):
        self.atom1 = PdbAtom(1.0, 1.0, 1.0, "H", 1, "H1")
        self.atom2 = PdbAtom(1.0, 1.0, 2.0, "C", 2, "CA")
        self.atom3 = PdbAtom(1.0, 1.0, 3.0, "O", 3, "OX1")
        self.residue1 = PdbResidue("A1", "ARG", self.atom1, self.atom2, self.atom3)
        self.atom4 = PdbAtom(1.0, 1.0, 4.0, "H", 4, "H1")
        self.atom5 = PdbAtom(1.0, 1.0, 5.0, "C", 5, "CA")
        self.atom6 = PdbAtom(1.0, 1.0, 6.0, "O", 6, "OX1")
        self.residue2 = PdbResidue("A2", "HST", self.atom4, self.atom5, self.atom6)
        self.atom7 = PdbAtom(1.0, 1.0, 7.0, "H", 7, "H1")
        self.atom8 = PdbAtom(1.0, 1.0, 8.0, "C", 8, "CA")
        self.atom9 = PdbAtom(1.0, 1.0, 9.0, "O", 9, "OX1")
        self.residue3 = PdbResidue("A3", "TRP", self.atom7, self.atom8, self.atom9)


    def check_valid_chain(self, chain):
        self.assertIsInstance(chain, PdbChain)
        self.assertIsInstance(chain, ResiduicSequence)
        self.assertIsInstance(chain.chain_id, str)
        self.assertIsInstance(chain.residues, list)
        self.assertIsInstance(chain.alpha_helices, set)
        self.assertIsInstance(chain.beta_strands, set)
        for residue in chain.residues:
            self.assertIs(residue.chain, chain)
        chain.model
        self.assertIsInstance(chain.missing_residues, list)
        self.assertRegex(str(chain), r"<Chain [A-Z] \((\d+) residues\)>")



class ChainCreationTests(ChainTest):

    def test_can_create_chain(self):
        chain = PdbChain("A", self.residue1, self.residue2, self.residue3)
        self.check_valid_chain(chain)


    def test_chain_id_must_be_str(self):
        with self.assertRaises(TypeError):
            chain = PdbChain(1, self.residue1, self.residue2, self.residue3)
        with self.assertRaises(TypeError):
            chain = PdbChain(None, self.residue1, self.residue2, self.residue3)



class ChainSequenceTests(ChainTest):

    def test_can_get_sequence_with_missing_residues(self):
        chain = PdbChain("A", self.residue1, self.residue2, self.residue3)
        chain.missing_residues = ["A2A", "A2B", "A3A", "A4", "A5", "A5A"]
        self.assertEqual(
         chain.get_all_residue_ids(),
         ["A1", "A2", "A2A", "A2B", "A3", "A3A", "A4", "A5", "A5A"]
        )


    def test_can_get_missing_residues_when_residue_ids_are_messed_up(self):
        self.residue1.residue_id = "A2"
        self.residue2.residue_id = "A100"
        self.residue3.residue_id = "A4"
        chain = PdbChain("A", self.residue1, self.residue2, self.residue3)
        chain.missing_residues = ["A1", "A5"]
        self.assertEqual(
         chain.get_all_residue_ids(),
         ["A1", "A2", "A100", "A4", "A5"]
        )


    def test_can_order_negative_residue_ids(self):
        chain = PdbChain("A", self.residue1, self.residue2, self.residue3)
        chain.missing_residues = ["A-2", "A-1", "A0", "A4"]
        self.assertEqual(
         chain.get_all_residue_ids(),
         ["A-2", "A-1", "A0", "A1", "A2", "A3", "A4"]
        )


class ChainSecondaryStructureTests(ChainTest):

    def setUp(self):
        ChainTest.setUp(self)
        self.chain = PdbChain("A", self.residue1, self.residue2, self.residue3)
        self.helix1 = PdbAlphaHelix("AH", self.residue1, self.residue2)
        self.helix2 = PdbAlphaHelix("AG", self.residue3)
        self.strand1 = PdbBetaStrand(1, 0, self.residue1, self.residue2)
        self.strand2 = PdbBetaStrand(2, 1, self.residue3)


    def test_can_get_helices_by_id(self):
        self.assertEqual(
         self.chain.get_alpha_helix_by_id("AH"),
         self.helix1
        )
        self.assertEqual(
         self.chain.get_alpha_helix_by_id("AG"),
         self.helix2
        )
        self.assertEqual(
         self.chain.get_alpha_helix_by_id("AF"),
         None
        )


    def test_can_only_search_helices_by_str_id(self):
        with self.assertRaises(TypeError):
            self.chain.get_alpha_helix_by_id(None)


    def test_can_get_strands_by_id(self):
        self.assertEqual(
         self.chain.get_beta_strand_by_id(1),
         self.strand1
        )
        self.assertEqual(
         self.chain.get_beta_strand_by_id(2),
         self.strand2
        )
        self.assertEqual(
         self.chain.get_beta_strand_by_id(3),
         None
        )


    def test_can_only_search_strands_by_int_id(self):
        with self.assertRaises(TypeError):
            self.chain.get_beta_strand_by_id(None)



class ChainMatrixTests(ChainTest):

    def setUp(self):
        self.atom1 = PdbAtom(3.696, 33.898, 63.219, "N", 1, "N")
        self.atom2 = PdbAtom(3.198, 33.218, 61.983, "C", 2, "CA")
        self.atom3 = PdbAtom(3.914, 31.863, 61.818, "O", 3, "C")
        self.residue1 = PdbResidue("A11", "VAL", self.atom1, self.atom2, self.atom3)
        self.atom4 = PdbAtom(3.155, 30.797, 61.557, "N", 4, "N")
        self.atom5 = PdbAtom(3.728, 29.464, 61.400, "C", 5, "CA")
        self.atom6 = PdbAtom(4.757, 29.459, 60.275, "O", 6, "C")
        self.residue2 = PdbResidue("A12", "MET", self.atom4, self.atom5, self.atom6)
        self.atom7 = PdbAtom(5.980, 29.039, 60.600, "N", 7, "N")
        self.atom8 = PdbAtom(7.092, 28.983, 59.649, "C", 8, "CA")
        self.atom9 = PdbAtom(7.092, 30.310, 58.987, "O", 9, "C")
        self.residue3 = PdbResidue("A13", "ASN", self.atom7, self.atom8, self.atom9)
        self.chain = PdbChain("A", self.residue1, self.residue2, self.residue3)
        self.chain.missing_residues = ["A12A", "A14"]


    def get_cells(self, matrix):
        return [g for g in matrix.graphics if
         isinstance(g, omnicanvas.graphics.Rectangle) and g.width == g.height]


    def check_valid_matrix(self, matrix):
        self.assertIsInstance(matrix, Canvas)
        cells = self.get_cells(matrix)
        dimension = ((math.sqrt((8 * len(cells)) + 1)) - 1) / 2
        self.assertEqual(dimension, int(dimension))
        top_row = cells[:int(dimension)]
        self.assertEqual(len(set([cell.y for cell in top_row])), 1)
        furthest_left = min([cell.x for cell in cells])
        self.assertEqual(
         len([cell.x for cell in cells if cell.x == furthest_left]),
         len(top_row)
        )
        for cell in cells:
            self.assertEqual(cell.line_width, 0)


    def test_can_generate_basic_matrix(self):
        matrix = self.chain.generate_residue_distance_matrix()
        self.check_valid_matrix(matrix)


    def test_matrix_cells_are_correct(self):
        matrix = self.chain.generate_residue_distance_matrix()
        self.check_valid_matrix(matrix)
        cells = self.get_cells(matrix)
        self.assertEqual(len(cells), 10)


    def test_missing_residues_produce_white_cells(self):
        matrix = self.chain.generate_residue_distance_matrix()
        self.check_valid_matrix(matrix)
        cells = self.get_cells(matrix)
        for cell in [cells[0], cells[1], cells[2], cells[3], cells[6], cells[7], cells[8]]:
            self.assertEqual(cell.fill_color, "#FFFFFF")
        for cell in [cells[4], cells[5], cells[9]]:
            self.assertNotEqual(cell.fill_color, "#FFFFFF")


    def test_present_residues_produce_sort_of_greenish_cells(self):
        matrix = self.chain.generate_residue_distance_matrix()
        self.check_valid_matrix(matrix)
        cells = self.get_cells(matrix)
        for cell in [cells[0], cells[2], cells[5]]:
            self.assertEqual(cell.fill_color[3], "F")


    def test_can_vary_color_range(self):
        matrix = self.chain.generate_residue_distance_matrix(close_color=0, far_color=120)
        self.check_valid_matrix(matrix)
        cells = self.get_cells(matrix)
        for cell in [cells[0], cells[2], cells[5]]:
            self.assertEqual(cell.fill_color[1], "F")


    def test_can_vary_cutoff(self):
        matrix = self.chain.generate_residue_distance_matrix(cutoff=0.1)
        self.check_valid_matrix(matrix)
        cells = self.get_cells(matrix)
        for cell in [cells[0], cells[2], cells[5]]:
            self.assertEqual(cell.fill_color[1], "F")


    def test_can_vary_dimension(self):
        matrix = self.chain.generate_residue_distance_matrix(dimension=70)
        self.check_valid_matrix(matrix)
        self.assertEqual(matrix.width, 70)


    def test_javascript_correct(self):
        matrix = self.chain.generate_residue_distance_matrix()
        for cell in self.get_cells(matrix):
            self.assertIn("onmouseover", cell.data.keys())
            self.assertIn("onmouseleave", cell.data.keys())
            self.assertIn("data", cell.data.keys())
            self.assertEqual(cell.data["onmouseover"], "cellHovered(this)")
            self.assertEqual(cell.data["onmouseleave"], "cellLeft(this)")


    def test_can_supply_subsequence(self):
        PdbAlphaHelix("1", self.residue1)
        matrix = self.chain.generate_residue_distance_matrix(subsequence=(self.residue1, self.residue2))
