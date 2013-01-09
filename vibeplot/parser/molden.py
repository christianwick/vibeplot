# Copyright (c) 2011, Mathias Laurin
#
# Permission to use, copy, modify, and/or distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.


import numpy as np
import vibeplot.chemistry as chem
import vibeplot.parser as cor


def pm_to_au(x):
    return float(x) * 0.01890

def au_to_pm(x):
    return float(x) * 52.918

@cor.coroutine
def to_atom(target):
    while True:
        line = (yield)
        symbol, index, mass, x, y, z = line.split()
        symbol = symbol.capitalize()
        atom = chem.Atom(**dict(symbol=symbol, index=int(index) - 1,
                           mass=int(mass)))
        atom.x = au_to_pm(x)
        atom.y = au_to_pm(y)
        atom.z = au_to_pm(z)
        target.send(atom)

@cor.coroutine
def to_peak(target):
    while True:
        freq = (yield)
        target.send((freq, 1.0))

@cor.coroutine
def parse_normal_coords(target):
    while True:
        line = (yield)
        line = line.strip()
        if line.startswith("vibration"):
            header = line
        else:
            xyz = np.fromstring(line, sep=" ")
            xyz = np.vectorize(au_to_pm)(xyz)
            target.send((header, xyz))

@cor.coroutine
def parse_geometries(target):
    while True:
        line = (yield)
        line = line.split()
        if len(line) == 1:
            nb = int(line[0])
            counter = 0
            xyz = ""
        elif line:
            counter += 1
            line = " ".join(line[1:])
            xyz = " ".join((xyz, line))
            if counter == nb:
                xyz = np.fromstring(xyz, sep=" ")
                xyz = np.vectorize(au_to_pm)(xyz)
                xyz = np.reshape(xyz, (-1, 3))
                target.send(xyz)

def molden_parser(atom_list, spk_dct, vib_dct):
    """
    :arg atom_list: list of atoms for graph
    :arg pl_list: list of peaks for spectrum
    :arg vib_dct: dict of normal coords {"vibration `n`": "[coords]"}
    :return: coroutine
    """
    mapping = {"[Atoms]": to_atom(cor.append(atom_list)),
               "[FREQ]": cor.to_type(float, to_peak(cor.dct_set(spk_dct))),
               "[FR-NORM-COORD": parse_normal_coords(cor.extend(vib_dct))}
    return cor.broadcast([cor.parse_section(header, coroutine)
                          for header, coroutine in mapping.iteritems()])

def load_molden(filename):
    a_list = []
    spk = {}
    v_dct = {}

    parser = molden_parser(a_list, spk, v_dct)
    with open(filename, "r") as f:
        for line in f:
            parser.send(line)

    delta = len(v_dct) - len(spk)  # number of 0-freq vibrations
    nc_list = [v_dct[key] for key in sorted(v_dct)][delta:]

    for freq, normal_coords in zip(sorted(spk), nc_list):
        for atom, delta_coords in zip(a_list, normal_coords):
            atom.normal_coords[freq] = atom.coords + delta_coords

    return chem.Molecule(a_list), spk

if __name__ == "__main__":
    print(repr(load_molden("tmp/molden.input")))

