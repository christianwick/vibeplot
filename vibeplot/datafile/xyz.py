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

from base import AbstractDataFile
from vibeplot.chemistry import *


class XYZAtom(object):

    def __init__(self):
        self.symbol = None
        self.index = None
        self.xmin = self.xmax = None
        self.ymin = self.ymax = None
        self.zmin = self.zmax = None

    @property
    def x(self):
        return 0.5 * (self.xmin + self.xmax)

    @x.setter
    def x(self, x):
        if self.xmin is None and self.xmax is None:
            self.xmin = self.xmax = x
        elif x < self.xmin:
            self.xmin = x
        elif x > self.xmax:
            self.xmax = x

    @property
    def y(self):
        return 0.5 * (self.ymin + self.ymax)

    @y.setter
    def y(self, y):
        if self.ymin is None and self.ymax is None:
            self.ymin = self.ymax = y
        elif y < self.ymin:
            self.ymin = y
        elif y > self.ymax:
            self.ymax = y

    @property
    def z(self):
        return 0.5 * (self.zmin + self.zmax)

    @z.setter
    def z(self, z):
        if self.zmin is None and self.zmax is None:
            self.zmin = self.zmax = z
        elif z < self.zmin:
            self.zmin = z
        elif z > self.zmax:
            self.zmax = z

    @property
    def xyz(self):
        return np.array([self.x, self.y, self.z])

    @property
    def dx(self):
        return self.xmax - self.xmin

    @property
    def dy(self):
        return self.ymax - self.ymin

    @property
    def dz(self):
        return self.zmax - self.zmin

    @property
    def dxdydz(self):
        return np.array([self.dx, self.dy, self.dz])


def angstrom_to_pm(x):
    return 100.0 * x


class XYZFile(AbstractDataFile):
    """ Parse animated (frame by frame) XYZ format """

    def __init__(self, filename='file.xyz'):
        self._atoms = {}
        super(XYZFile, self).__init__(filename)

    def parse(self):
        with open(self.filename, 'r') as f:
            while True:
                # header
                try:
                    n_atoms = int(f.readline().strip())
                except ValueError:
                    break
                f.readline()
                # body
                for n in range(n_atoms):
                    atom = self._atoms.get(n, XYZAtom())
                    symbol, x, y, z = f.readline().split()
                    atom.symbol = symbol
                    atom.index = n
                    atom.x = angstrom_to_pm(float(x))
                    atom.y = angstrom_to_pm(float(y))
                    atom.z = angstrom_to_pm(float(z))
                    self._atoms[n] = atom

            
        freq = 1000
        for atom in self._atoms.itervalues():
            a = Atom(symbol=atom.symbol,
                     index=atom.index,
                     xyz=atom.xyz
                    )
            a.normal_coords[freq] = atom.xyz + atom.dxdydz
            self.atoms.append(a)
        self.vibrations.append(Vibration(freq))
