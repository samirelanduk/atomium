import math
import omnicanvas
from .molecules import Residue

def generate_residue_distance_matrix(self, dimension=700, close_color=120,
     far_color=0, cutoff=40, subsequence=None):
        # Validation
        if not isinstance(close_color, int):
            raise TypeError("close_color must be int, not '%s'" % str(close_color))
        if not isinstance(far_color, int):
            raise TypeError("far_color must be int, not '%s'" % str(far_color))
        if not 0 <= close_color < 360:
            raise ValueError("close_color must be between 0 and 360, not %i" % close_color)
        if not 0 <= far_color < 360:
            raise ValueError("far_color must be between 0 and 360, not %i" % far_color)
        if not isinstance(cutoff, int) and not isinstance(cutoff, float):
            raise TypeError("cutoff must be numeric, not '%s'" % str(dimension))
        if not isinstance(dimension, int):
            raise TypeError("dimension must be int, not '%s'" % str(dimension))
        if subsequence:
            if not isinstance(subsequence, list) and not isinstance(subsequence, tuple):
                raise TypeError("subsequence must be a list or tuple")
            for residue in subsequence:
                if not isinstance(residue, Residue):
                    raise TypeError("Only Residues can be given for subsequence")
                if residue not in self.residues():
                    raise ValueError(
                     "%s does not have subsequence residue %s" % (str(self), str(residue))
                    )
            if len(subsequence) != 2:
                raise ValueError(
                 "Only two residues can be used to define a subsequence, not %s" % str(subsequence)
                )
        # Set up canvas
        matrix = omnicanvas.Canvas(dimension, dimension)
        residues = self.residues()
        residue_ids = [res.residue_id() for res in residues]

        # Set up parameters
        padding_proportion = 0.09
        padding = padding_proportion * dimension
        plot_dimension = dimension - (2 * padding)
        chain_length = len(residues)
        cell_dimension = plot_dimension / chain_length
        plot_width = dimension - (2 * padding)
        bar_width = 4
        bar_left = (dimension / 2) - (bar_width / 2) + 5
        hypoteneuse = math.sqrt((plot_width ** 2) + (plot_width ** 2))
        bar_top = (dimension / 2) - (hypoteneuse / 2) + 5
        diagonal_chunk = hypoteneuse / len(residues)
        chain_color = 80
        helix_color = 325
        strand_color = 182
        tick = 0
        if len(residues) >= 10000:
            tick = 5000
        elif len(residues) >= 1000:
            tick = 500
        elif len(residues) >= 100:
            tick = 50
        elif len(residues) >= 10:
            tick = 5
        else:
            tick = 1

        # Calculate distances
        distances = []
        for index, residue2 in enumerate(residues[1:][::-1]):
            row = []
            for residue1 in residues[:0 - (index + 1)]:
                row.append({
                 "residue1": residue1,
                 "residue2": residue2,
                 "distance": residue1.alpha_carbon().distance_to(residue2.alpha_carbon()) if not residue1.is_missing() and not residue2.is_missing() else None
                })
            distances.append(row)

        # Add cells
        for row_index, row in enumerate(distances):
            for cell_index, cell in enumerate(row):
                color = "#FFFFFF"
                if cell["distance"] is not None:
                    fraction = cell["distance"] / cutoff if cell["distance"] <= cutoff else 1
                    if far_color >= close_color:
                        distance_from_start = fraction * (far_color - close_color)
                        color = close_color + distance_from_start
                    else:
                        distance_from_start = fraction * (close_color - far_color)
                        color = close_color - distance_from_start
                    color = omnicanvas.hsl_to_rgb(color, 100, 50)

                matrix.add_rectangle(
                 padding + (cell_index * cell_dimension),
                 padding + (row_index * cell_dimension),
                 cell_dimension + 0.5,
                 cell_dimension + 0.5,
                 fill_color=color,
                 line_width=0,
                 data={
                  "onmouseover": "cellHovered(this)",
                  "onmouseleave": "cellLeft(this)",
                  "data": "%s,%i,%s,%i,%.2f" % (
                   cell["residue1"].residue_name(),
                   cell_index + 1,
                   cell["residue2"].residue_name(),
                   len(residues) - row_index,
                   cell["distance"] if cell["distance"] is not None else 0.0
                  )
                 }
                )

        # Subsequence
        if subsequence:
            x = (residues.index(subsequence[0]) * cell_dimension) + padding
            y = dimension - (((residues.index(subsequence[1]) + 1) * cell_dimension) + padding)
            matrix.add_line(
             x, dimension - padding, x, y,
             line_width=1.5,
             line_style=".."
            )
            matrix.add_line(
             dimension - padding, y, x, y,
             line_width=1.5,
             line_style=".."
            )

        # Add gridlines, border and labels
        residue_number = 0
        while residue_number <= len(residues) - 1:
            x = padding + (residue_number * cell_dimension)
            y = dimension - x
            matrix.add_line(
             x, padding, x, dimension - padding
            )
            matrix.add_line(
             padding, y, dimension - padding, y
            )
            if residue_number != len(residues) - 1:
                matrix.add_text(
                 x + (0.5 * cell_dimension), padding * 0.75, str(residue_number + 1)
                )
            if residue_number != 0:
                matrix.add_text(
                 padding - 2, y - (0.5 * cell_dimension), str(residue_number + 1),
                 horizontal_align="left"
                )
            residue_number += tick
        matrix.add_text(dimension / 2, padding * 0.35, "Residue 1")
        matrix.add_text(
         padding * 0.35, dimension / 2, "Residue 2", rotation=(
          padding * 0.35, dimension / 2, 270
         )
        )
        matrix.add_polygon(
         dimension - padding, padding,
         dimension - padding, dimension - padding,
         padding, dimension - padding,
         line_width=2,
         line_color="#FFFFFF"
        )
        matrix.add_line(
         padding, padding, dimension - padding, padding, line_width=2
        )
        matrix.add_line(
         padding, padding, padding, dimension - padding, line_width=2
        )
        matrix.add_line(
         dimension - padding, padding,  padding, dimension - padding, line_width=2
        )

        # Add secondary structure
        matrix.add_rectangle(
         bar_left, bar_top, bar_width, hypoteneuse,
         rotation=((dimension / 2) + 5, (dimension / 2) + 5, 45),
         line_width=0,
         fill_color=omnicanvas.hsl_to_rgb(chain_color, 100, 50)
        )
        for helix in self.alpha_helices():
            start = (len(residues) - residue_ids.index(
             helix.residues()[-1].residue_id()
            )) - 1
            end = (len(residues) - residue_ids.index(
             helix.residues()[0].residue_id()
            )) - 1
            matrix.add_rectangle(
             bar_left - 1, bar_top + (diagonal_chunk * start),
             bar_width + 2, diagonal_chunk * ((end - start) + 1),
             fill_color=omnicanvas.hsl_to_rgb(helix_color, 100, 50),
             line_width=0,
             rotation=((dimension / 2) + 5, (dimension / 2) + 5, 45)
            )
        for strand in self.beta_strands():
            start = (len(residues) - residue_ids.index(
             strand.residues()[-1].residue_id()
            )) - 1
            end = (len(residues) - residue_ids.index(
             strand.residues()[0].residue_id()
            )) - 1
            matrix.add_rectangle(
             bar_left - 1, bar_top + (diagonal_chunk * start),
             bar_width + 2, diagonal_chunk * (end - start),
             fill_color=omnicanvas.hsl_to_rgb(strand_color, 100, 50),
             line_width=0,
             rotation=((dimension / 2) + 5, (dimension / 2) + 5, 45)
            )

        # Add legend
        legend_dimension = plot_width * 0.4
        legend_left = padding + (dimension / 2)
        legend_top = dimension - (padding + legend_dimension)
        scale_width = legend_dimension * 0.8
        scale_height = legend_dimension * 0.1
        scale_left = legend_left + (0.1 * legend_dimension)
        scale_right = legend_left + (0.9 * legend_dimension)
        scale_top = legend_top + (0.1 * legend_dimension)
        scale_bottom = legend_top + (0.2 * legend_dimension)
        scale_label_y = legend_top + (0.05 * legend_dimension)
        number_label_y = legend_top + (0.28 * legend_dimension)
        white_top = legend_top + (0.36 * legend_dimension)
        helix_top = legend_top + (0.5 * legend_dimension)
        helix_bottom = legend_top + (0.6 * legend_dimension)
        helix_width = legend_dimension * 0.4
        helix_left = scale_left
        helix_right = helix_left + helix_width
        strand_top = legend_top + (0.7 * legend_dimension)
        strand_bottom = legend_top + (0.8 * legend_dimension)
        bindseq_top = legend_top + (0.9 * legend_dimension)
        bindseq_bottom = legend_top + (1.0 * legend_dimension)
        x_pixels = range(math.floor(scale_left), math.ceil(scale_right))

        matrix.add_text(
         scale_left + (scale_width / 2),
         scale_label_y,
         "Distance (&#8491;ngstroms)",
         font_size=int(scale_width / 9)
        )
        for x_pixel in x_pixels:
            color = 0
            fraction = (x_pixel - scale_left) / (scale_right - scale_left)
            if far_color >= close_color:
                distance_from_start = fraction * (far_color - close_color)
                color = int(close_color + distance_from_start)
            else:
                distance_from_start = fraction * (close_color - far_color)
                color = int(close_color - distance_from_start)
            matrix.add_rectangle(
             x_pixel - 1, scale_top, 2, scale_bottom - scale_top,
             fill_color=omnicanvas.hsl_to_rgb(color, 100, 50),
             line_width=0
            )
        matrix.add_text(
         scale_left, number_label_y, "0",
         font_size=int(scale_width / 10)
        )
        matrix.add_text(
         scale_right, number_label_y, "%i+" % cutoff,
         font_size=int(scale_width / 10)
        )
        matrix.add_text(
         scale_right - (scale_width * 0.5), white_top,
         "White areas: PDB missing residues",
         font_size=int(scale_width / 15),
         horizontal_align="center"
        )
        matrix.add_rectangle(
         helix_left, ((helix_bottom + helix_top) / 2) - ((bar_width / 2) + 0),
         helix_width, bar_width,
         fill_color=omnicanvas.hsl_to_rgb(chain_color, 100, 50),
         line_width=0
        )
        matrix.add_rectangle(
         helix_left + (0.1 * helix_width),
         ((helix_bottom + helix_top) / 2) - ((bar_width / 2) + 1),
         helix_width * 0.8,
         bar_width + 2,
         fill_color=omnicanvas.hsl_to_rgb(helix_color, 100, 50),
         line_width=0
        )
        matrix.add_rectangle(
         helix_left, ((strand_bottom + strand_top) / 2) - ((bar_width / 2) + 0),
         helix_width, bar_width,
         fill_color=omnicanvas.hsl_to_rgb(chain_color, 100, 50),
         line_width=0
        )
        matrix.add_rectangle(
         helix_left + (0.1 * helix_width),
         ((strand_bottom + strand_top) / 2) - ((bar_width / 2) + 1),
         helix_width * 0.8,
         bar_width + 2,
         fill_color=omnicanvas.hsl_to_rgb(strand_color, 100, 50),
         line_width=0
        )
        matrix.add_text(
         helix_right + (legend_dimension * 0.1),
         ((helix_bottom + helix_top) / 2),
         "&#945;-helix",
         font_size=int(scale_width / 10),
         horizontal_align="right"
        )
        matrix.add_text(
         helix_right + (legend_dimension * 0.1),
         ((strand_bottom + strand_top) / 2),
         "&#946;-strand",
         font_size=int(scale_width / 10),
         horizontal_align="right"
        )
        if subsequence:
            matrix.add_line(
             helix_left, ((bindseq_bottom + bindseq_top) / 2) - ((bar_width / 2) + 0),
             helix_right, ((bindseq_bottom + bindseq_top) / 2) - ((bar_width / 2) + 0),
             line_width=1.5,
             line_style=".."
            )
            matrix.add_text(
             helix_right + (legend_dimension * 0.1),
             ((bindseq_bottom + bindseq_top) / 2),
             "Bind Sequence",
             font_size=int(scale_width / 10),
             horizontal_align="right"
            )

        return matrix
