from collections import deque
import re

def mmcif_string_to_mmcif_dict(filestring):
    """Takes a .cif filestring and turns into a ``dict`` which represents its
    table structure. Only lines which aren't empty and which don't begin with
    ``#`` are used.

    Multi-line strings are consolidated onto one line, and the whole thing is
    then split into the blocks that will become table lists. At the end, quote
    marks are removed from any string which retains them.

    :param str filestring: the .cif filestring to process.
    :rtype: ``dict``"""

    lines = deque(filter(lambda l: l and l[0] != "#", filestring.split("\n")))
    lines = consolidate_strings(lines)
    blocks = mmcif_lines_to_mmcif_blocks(lines)
    mmcif_dict = {}
    for block in blocks:
        if block["lines"][0] == "loop_":
            mmcif_dict[block["category"]] = loop_block_to_list(block)
        else:
            mmcif_dict[block["category"]] = non_loop_block_to_list(block)
    strip_quotes(mmcif_dict)
    return mmcif_dict


def consolidate_strings(lines):
    """Generally, .cif files have a one file line to one table row
    correspondence. Sometimes however, a string cell is given a line of its own,
    breaking the row over several lines. This function takes the lines of a .cif
    file and puts all table rows on a single line.

    :param deque lines: the .cif file lines.
    :rtype: ``deque``"""

    new_lines = deque()
    while lines:
        line = lines.popleft()
        if line.startswith(";"):
            string = [line[1:].strip()]
            while not lines[0].startswith(";"):
                string.append(lines.popleft())
            lines.popleft()
            join = "\x00" if len(string) > 1 and " " in "".join(string) else ""
            string = join.join(string).replace('"', "\x1a")
            new_lines.append("\"{}\"".format(string))
        else:
            new_lines.append(line)
    return new_lines


def mmcif_lines_to_mmcif_blocks(lines):
    """A .cif file is ultimately a list of tables. This function takes a list of
    .cif file lines and splits them into these table blocks. Each block will be
    a ``dict`` containing a category name and a list of lines.

    :param deque lines: the .cif file lines.
    :rtype: ``list``"""

    category = None
    block, blocks = [], []
    while lines:
        line = lines.popleft()
        if line.startswith("data_"): continue
        if line.startswith("_"):
            line_category = line.split(".")[0]
            if line_category != category:
                if category:
                    blocks.append({"category": category[1:], "lines": block})
                category = line_category
                block = []
        if line.startswith("loop_"):
            if category:
                blocks.append({"category": category[1:], "lines": block})
            category = lines[0].split(".")[0]
            block = []
        block.append(line)
    if block: blocks.append({"category": category[1:], "lines": block})
    return blocks


def non_loop_block_to_list(block):
    """Takes a simple block ``dict`` with no loop and turns it into a table
    ``list``.

    :param dict block: the .cif block to process.
    :rtype: ``list``"""

    d = {}
    for index in range(len(block["lines"]) - 1):
        if block["lines"][index + 1][0] != "_":
            block["lines"][index] += " " + block["lines"][index + 1]
    block["lines"] = [l for l in block["lines"] if l[0] == "_"]
    for line in block["lines"]:
        name = line.split(".")[1].split()[0]
        value = line
        if line.startswith("_"):
            value = " ".join(line.split()[1:])
        d[name] = value
    return [d]


def loop_block_to_list(block):
    """Takes a loop block ``dict`` where the initial lines are table headers and
    turns it into a table ``list``. Sometimes a row is broken over several lines
    so this function deals with that too.

    :param dict block: the .cif block to process.
    :rtype: ``list``"""

    names, lines = [], []
    body_start = 0
    for index, line in enumerate(block["lines"][1:], start=1):
        if not line.startswith("_" + block["category"]):
            body_start = index
            break
    names = [l.split(".")[1].rstrip() for l in block["lines"][1:body_start]]
    lines = [split_values(l) for l in block["lines"][body_start:]]
    l = []
    for n in range(len(lines) - 1):
        while n < len(lines) - 1 and\
         len(lines[n]) + len(lines[n + 1]) <= len(names):
            lines[n] += lines[n + 1]
            lines.pop(n + 1)
    for line in lines:
        l.append({
            name: value for name, value in zip(names, line)
        })
    return l


def split_values(line):
    """The body of a .cif table is a series of lines, with each cell divided by
    whitespace. This function takes a string line and breaks it into cells.

    There are a few peculiarities to handle. Sometimes a cell is a string
    enclosed in quote marks, and spaces within this string obviously shouldn't
    be used to break the line. This function handles all of that.

    :param str line: the .cif line to split.
    :rtype: ``list``"""

    if not re.search("[\'\"]", line): return line.split()
    chars = deque(line.strip())
    values, value, in_string = [], [], False
    while chars:
        char = chars.popleft()
        if char == " " and not in_string:
            values.append(value)
            value = []
        elif char in "'\"":
            if in_string and chars and chars[0] != " ":
                value.append(char)
            else:
                in_string = not in_string
        else:
            value.append(char)
    values.append(value)
    return ["".join(v) for v in values if v]


def strip_quotes(mmcif_dict):
    """Goes through each table in the mmcif ``dict`` and removes any unneeded
    quote marks from the cells.

    :param dict mmcif_dict: the almost finished .mmcif dictionary to clean."""

    for table in mmcif_dict.values():
        for row in table:
            for k, value in row.items():
                for char in "'\"":
                    if value[0] == char and value[-1] == char:
                        row[k] = value[1:-1]
                    row[k] = row[k].replace("\x1a", '"').replace("\x1b", "'").replace("\x00", "\n")


def save_mmcif_dict(mmcif_dict, path):
    code = mmcif_dict.get("entry", [{"id": ""}])[0]["id"]
    lines = [f"data_{code}" if code else "data_XXXX"]
    for category_name, rows in mmcif_dict.items():
        lines.append("#")
        if len(rows) == 1:
            lines += get_non_list_lines(category_name, rows)
        else:
            lines += get_list_lines(category_name, rows)
    lines.append("#")
    with open(path, "w") as f:
        f.write("\n".join(lines))


def get_non_list_lines(name, rows):
    lines = []
    keys = list(rows[0].keys())
    max_length = max(len(k) for k in keys)
    for key in keys:
        value = format_value(rows[0][key])
        if value[0] == ";":
            lines.append(f"_{name}.{key.ljust(max_length)}")
            lines.append(value)
        else:
            lines.append(f"_{name}.{key.ljust(max_length)} {value}")
    return lines


def get_list_lines(name, rows):
    lines = []
    lines.append("loop_")
    keys = list(rows[0].keys())
    for key in keys:
        lines.append(f"_{name}.{key}")
    grid = [[ format_value(row[key]) for key in keys] for row in rows]
    max_lengths = [max([len(row[col_num]) for row in grid])
        for col_num in range(len(keys))]
    for row in grid:
        line = []
        for i, cell in enumerate(row):
            line.append(cell.ljust(max_lengths[i]))
        lines.append(" ".join(line))
    return lines


def format_value(s):
    if "\n" in s:
        lines = s.splitlines()
        lines[0] = f";{lines[0]}"
        lines.append(";")
        return "\n".join(lines)
    if " " in s:
        if "'" in s:
            return f'"{s}"'
        else:
            return f"'{s}'"
    elif "'" in s:
        return f'"{s}"'
    elif '"' in s:
        return f"'{s}'"
    return s