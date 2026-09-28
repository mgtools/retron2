#!/usr/bin/env python3
"""
ct2dot.py — Convert an RNA secondary structure CT file to dot-bracket format.

CT file column layout (1-indexed):
  1: nucleotide index
  2: base (letter)
  3: index - 1  (previous)
  4: index + 1  (next)
  5: paired partner index (0 = unpaired)
  6: natural numbering index

Usage:
  python ct2dot.py input.ct [output.dot]

If output path is omitted the result is printed to stdout.
"""

import sys
import os


def parse_ct(path: str):
    """
    Parse a CT file and return (header_line, sequence, pairs_dict).

    pairs_dict maps each 1-based nucleotide index to its partner index
    (0 means unpaired).
    """
    with open(path) as fh:
        lines = fh.readlines()

    if not lines:
        raise ValueError("CT file is empty.")

    # First line: "  N  ENERGY = ...  date"
    header = lines[0].rstrip()
    parts = header.split()
    try:
        length = int(parts[0])
    except (IndexError, ValueError):
        raise ValueError(f"Cannot read sequence length from header: {header!r}")

    sequence = []
    pairs = {}          # index (1-based) → partner index (0 = unpaired)

    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue
        cols = line.split()
        if len(cols) < 5:
            continue
        try:
            idx     = int(cols[0])
            base    = cols[1]
            partner = int(cols[4])
        except ValueError:
            continue
        sequence.append(base)
        pairs[idx] = partner

    if len(sequence) != length:
        # Warn but continue; use however many nucleotides we actually found.
        sys.stderr.write(
            f"Warning: header says {length} nucleotides but {len(sequence)} "
            "data rows were found. Proceeding with actual count.\n"
        )

    return header, "".join(sequence), pairs


def pairs_to_dotbracket(pairs: dict) -> str:
    """
    Convert a pairs dict {1-based index → partner} to a dot-bracket string.

    Nested base pairs  → standard ( )
    Pseudoknots        → [ ] then { } then < > (successive nesting levels)
    """
    n = max(pairs.keys())
    db = ["."] * (n + 1)   # index 0 unused; positions 1..n

    # Collect canonical pairs (i < j only, skip unpaired)
    canonical = sorted(
        (min(i, p), max(i, p))
        for i, p in pairs.items()
        if p != 0 and i < p
    )

    # Assign bracket levels using an interval-nesting check
    OPEN  = ["(", "[", "{", "<"]
    CLOSE = [")", "]", "}", ">"]
    MAX_LEVELS = len(OPEN)

    # For each pair find the deepest nesting level whose stack doesn't conflict
    level_stacks: list[list[int]] = [[] for _ in range(MAX_LEVELS)]

    def conflicts(j_open: int, level: int) -> bool:
        """Does adding a pair whose opening bracket is at j_open conflict with level's stack?"""
        # A conflict exists if any open bracket on the stack is between
        # j_open and the partner of the bracket on the stack (i.e. crossing).
        # Actually we just check: is j_open inside any currently open pair at this level?
        for open_pos in level_stacks[level]:
            # open_pos is an already-opened-but-not-yet-closed index at this level
            # We need the closing position, which we store separately.
            pass
        return False  # placeholder; real check below

    # Simpler greedy: assign each pair a level based on whether it crosses
    # any pair already assigned to that level.
    assigned: list[tuple[int, int, int]] = []  # (i, j, level)

    def crosses(i1: int, j1: int, i2: int, j2: int) -> bool:
        return (i1 < i2 < j1 < j2) or (i2 < i1 < j2 < j1)

    for (i, j) in canonical:
        placed = False
        for lvl in range(MAX_LEVELS):
            # Check if (i, j) crosses any pair already at this level
            if not any(crosses(i, j, ai, aj) for ai, aj, al in assigned if al == lvl):
                assigned.append((i, j, lvl))
                placed = True
                break
        if not placed:
            sys.stderr.write(
                f"Warning: pair ({i},{j}) exceeds {MAX_LEVELS} nesting levels "
                "and will be shown as '.'\n"
            )

    for (i, j, lvl) in assigned:
        db[i] = OPEN[lvl]
        db[j] = CLOSE[lvl]

    return "".join(db[1:])   # return 1..n as a 0-based string


def ct_to_dot(ct_path: str, dot_path: str | None = None) -> None:
    header, sequence, pairs = parse_ct(ct_path)
    dotbracket = pairs_to_dotbracket(pairs)

    # The CT header starts with the sequence length followed by the rest
    # (e.g. "  196  ENERGY = 16.4  date").  The dot-bracket convention
    # omits that leading count, so strip it here.
    header_parts = header.split(None, 1)          # ["196", "ENERGY = ..."]
    header_label = header_parts[1] if len(header_parts) > 1 else header

    lines_out = [
        f">{header_label}",
        sequence,
        dotbracket,
    ]
    output = "\n".join(lines_out) + "\n"

    if dot_path:
        with open(dot_path, "w") as fh:
            fh.write(output)
        print(f"Written to {dot_path}")
    else:
        sys.stdout.write(output)


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    ct_path  = sys.argv[1]
    dot_path = sys.argv[2] if len(sys.argv) >= 3 else None

    if not os.path.isfile(ct_path):
        sys.stderr.write(f"Error: file not found: {ct_path}\n")
        sys.exit(1)

    ct_to_dot(ct_path, dot_path)


if __name__ == "__main__":
    main()
