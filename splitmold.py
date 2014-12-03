#!/usr/bin/env python
# Copyright (c) 2014 Mathias Laurin, 3-clause BSD License
"""Select a range of atoms and their normal modes in a Molden file.

Example:
    Keep the first 6 atoms:
        # splitmold.py 1 7 data.molden > smaller.molden

"""
from __future__ import print_function, absolute_import
import os
import sys
from coroutines import parse_section, periodic_split, append, printer


def molden_splitter(beg, end, filename):
    """Write Molden file containing `filename`'s atom indexes [beg end[.

    Attributes:
        beg (int): Index of first atom (starts at 1).
        end (int): Index of last-plus-one atom.
        filename: Path to the Molden file to split.

    """
    beg -= 1
    end -= 1
    # null printer
    null = printer(open(os.devnull, "w"))
    # First pass: count atoms
    atom_list = []
    parser = parse_section("[Atoms]", append(atom_list), null)
    with open(filename, "rb") as input:
        for line in input:
            parser.send(line)
    atom_count = len(atom_list)
    # Check indexes
    if not (0 <= beg < end <= atom_count):
        raise IndexError("%s expects Molden indexes starting at 1."
                         % sys.argv[0])
    # Second pass: strip file
    stripper = periodic_split(range(beg, end),
                              atom_count,
                              printer(),
                              null)
    parser = parse_section(
        "[Atoms]",
        stripper,
        parse_section("vibration", stripper, printer()))
    with open(filename, "rb") as input:
        for line in input:
            parser.send(line)
    null.close()


if __name__ == "__main__":
    beg = int(sys.argv[1])
    end = int(sys.argv[2])
    filename = sys.argv[3]
    molden_splitter(beg, end, filename)
