import math
from omnicanvas.canvas import Canvas
import omnicanvas.graphics
from unittest import TestCase
import unittest.mock
from molecupy.structures import Chain, ResiduicSequence, Residue, PdbAtom, Atom

class ChainTest(TestCase):

    def setUp(self):
        self.residues = [unittest.mock.Mock(spec=Residue) for _ in range(10)]



class ChainCreationTests(ChainTest):

    def test_can_create_chain(self):
        chain = Chain("A", *self.residues)
        self.assertIsInstance(chain, ResiduicSequence)
        self.assertEqual(chain._chain_id, "A")
        self.assertEqual(chain._residues, self.residues)


    def test_chain_updates_residues(self):
        for residue in self.residues:
            residue._chain = None
        chain = Chain("A", *self.residues)
        for residue in self.residues:
            self.assertIs(residue._chain, chain)


    def test_chain_id_must_be_str(self):
        with self.assertRaises(TypeError):
            Chain(200, *self.residues)


    def test_chain_repr(self):
        chain = Chain("A", *self.residues)
        self.assertEqual(str(chain), "<Chain A (10 residues)>")



class ChainPropertyTests(ChainTest):

    def test_chain_properties(self):
        chain = Chain("A", *self.residues)
        self.assertEqual(chain.chain_id(), "A")


    def test_can_add_residues_and_update_them(self):
        chain = Chain("A", *self.residues)
        residue = unittest.mock.Mock(spec=Residue)
        residue._chain = None
        chain.add_residue(residue)
        self.assertIs(chain.residues()[-1], residue)
        self.assertIs(residue._chain, chain)


    def test_can_remove_residues_and_update_them(self):
        chain = Chain("A", *self.residues)
        chain.remove_residue(self.residues[5])
        self.assertNotIn(self.residues[5], chain.residues())
        self.assertIs(self.residues[5]._chain, None)



class ChainMatrixTests(ChainTest):

    def setUp(self):
        self.atom1 = PdbAtom(3.696, 33.898, 63.219, "N", 1, "N")
        self.atom2 = PdbAtom(3.198, 33.218, 61.983, "C", 2, "CA")
        self.atom3 = PdbAtom(3.914, 31.863, 61.818, "O", 3, "C")
        self.residue1 = Residue("A11", "VAL", self.atom1, self.atom2, self.atom3)
        self.atom4 = PdbAtom(3.155, 30.797, 61.557, "N", 4, "N")
        self.atom5 = PdbAtom(3.728, 29.464, 61.400, "C", 5, "CA")
        self.atom6 = PdbAtom(4.757, 29.459, 60.275, "O", 6, "C")
        self.residue2 = Residue("A12", "MET", self.atom4, self.atom5, self.atom6)
        self.atomA = Atom("C", 1001, "CA")
        self.residue3 = Residue("A12A", "TRP", self.atomA)
        self.atom7 = PdbAtom(5.980, 29.039, 60.600, "N", 7, "N")
        self.atom8 = PdbAtom(7.092, 28.983, 59.649, "C", 8, "CA")
        self.atom9 = PdbAtom(7.092, 30.310, 58.987, "O", 9, "C")
        self.residue4 = Residue("A13", "ASN", self.atom7, self.atom8, self.atom9)
        self.atomB = Atom("C", 1002, "CA")
        self.residue5 = Residue("14", "TRP", self.atomA)
        self.chain = Chain(
         "A", self.residue1, self.residue2, self.residue3, self.residue4, self.residue5
        )


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


    def test_colors_must_be_numeric(self):
        with self.assertRaises(TypeError):
            self.chain.generate_residue_distance_matrix(close_color="red")
        with self.assertRaises(TypeError):
            self.chain.generate_residue_distance_matrix(far_color="red")


    def test_colors_must_be_between_0_and_360(self):
        with self.assertRaises(ValueError):
            self.chain.generate_residue_distance_matrix(close_color=-1)
        with self.assertRaises(ValueError):
            self.chain.generate_residue_distance_matrix(close_color=360)
        with self.assertRaises(ValueError):
            self.chain.generate_residue_distance_matrix(far_color=-1)
        with self.assertRaises(ValueError):
            self.chain.generate_residue_distance_matrix(far_color=360)


    def test_can_vary_cutoff(self):
        matrix = self.chain.generate_residue_distance_matrix(cutoff=0.1)
        self.check_valid_matrix(matrix)
        cells = self.get_cells(matrix)
        for cell in [cells[0], cells[2], cells[5]]:
            self.assertEqual(cell.fill_color[1], "F")


    def test_cutoff_must_be_numeric(self):
        with self.assertRaises(TypeError):
            self.chain.generate_residue_distance_matrix(cutoff=None)


    def test_can_vary_dimension(self):
        matrix = self.chain.generate_residue_distance_matrix(dimension=70)
        self.check_valid_matrix(matrix)
        self.assertEqual(matrix.width, 70)


    def test_dimension_must_be_int(self):
        with self.assertRaises(TypeError):
            self.chain.generate_residue_distance_matrix(dimension=0.5)


    def test_javascript_correct(self):
        matrix = self.chain.generate_residue_distance_matrix()
        for cell in self.get_cells(matrix):
            self.assertIn("onmouseover", cell.data.keys())
            self.assertIn("onmouseleave", cell.data.keys())
            self.assertIn("data", cell.data.keys())
            self.assertEqual(cell.data["onmouseover"], "cellHovered(this)")
            self.assertEqual(cell.data["onmouseleave"], "cellLeft(this)")


    def test_can_supply_subsequence(self):
        matrix = self.chain.generate_residue_distance_matrix(
         subsequence=(self.residue1, self.residue2)
        )


    def test_subsequence_must_be_sequence_of_residues(self):
        with self.assertRaises(TypeError):
            matrix = self.chain.generate_residue_distance_matrix(
             subsequence=(self.residue1, "residue")
            )


    def test_subsequence_must_be_sequence_of_two_residues(self):
        with self.assertRaises(ValueError):
            matrix = self.chain.generate_residue_distance_matrix(
             subsequence=(self.residue1,)
            )


    def test_subsequence_residues_must_be_in_chain(self):
        with self.assertRaises(ValueError):
            matrix = self.chain.generate_residue_distance_matrix(
             subsequence=(self.residue1, Residue(
              "A4", "RES", PdbAtom(1.0, 1.0, 1.0, "H", 1, "H")
             ))
            )
