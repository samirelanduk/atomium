def align(long, short):
    matrix = create_nw_matrix(long, short)
    back_directions = calculate_matrix_scores(matrix, long, short)
    return calculate_aligment_indices(long, short, back_directions)
    

def create_nw_matrix(long, short):
    matrix = [[
        None for _ in range(len(long) + 1)
    ] for _ in range(len(short) + 1)]
    matrix[0] = [0 - n for n in range(len(matrix[0]))]
    for n, row in enumerate(matrix):
        row[0] = 0 - n
    return matrix


def calculate_matrix_scores(matrix, long, short, match=1, mismatch=-1, gap=-5):
    back_directions = {r: {
        c: [1] if r == 0 else [3] if c == 0 else [] for c in range(len(matrix[0]))
    } for r in range(len(matrix))}
    for row_num, row in enumerate(matrix):
        for col_num, cell in enumerate(row):
            if cell is None:
                short_res = short[row_num - 1]
                long_res = long[col_num - 1]
                score = match if short_res == long_res else mismatch
                corner = matrix[row_num - 1][col_num - 1] + score
                left = matrix[row_num][col_num - 1] + gap
                top = matrix[row_num - 1][col_num] + gap
                biggest = max(corner, left, top)
                matrix[row_num][col_num] = biggest
                if left == biggest: back_directions[row_num][col_num].append(1)
                if corner == biggest: back_directions[row_num][col_num].append(2)
                if top == biggest: back_directions[row_num][col_num].append(3)
    return back_directions


def calculate_aligment_indices(long, short, back_directions):
    row_num, col_num = len(short), len(long)
    path = [[row_num, col_num]]
    line1, line2 = [], []
    while path[-1] != [0, 0]:
        direction = back_directions[path[-1][0]][path[-1][1]][0]
        if direction == 1:
            line1.append(long[col_num - 1])
            line2.append("-")
        if direction == 2:
            line1.append(long[col_num - 1])
            line2.append(short[row_num - 1])
        if direction == 3:
            line1.append("-")
            line2.append(short[row_num - 1])
        if direction in (2, 3): row_num -= 1
        if direction in (1, 2): col_num -= 1
        path.append([row_num, col_num])
    indexes = []
    for n, char in enumerate(line2[::-1]):
        if char != "-": indexes.append(n)
    return indexes