# Copyright (c) 2011, 2012 Mathias Laurin
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


class Atom(object):

    __slots__ = """symbol index data_index mass
                   coords projected_coordinates
                   color label radius normal_coords""".split()

    X, Y, Z = 0, 1, 2

    vdw_radius = dict(H=120, He=140, Li=182, C=170, N=155, O=152, F=147,
                      Ne=154, Na=227, Si=210, P=180, S=180, Cl=175, Ar=188,
                      K=275, Ni=163, Cu=140, Zn=139, Ga=187, As=185, Se=190,
                      Br=185, Kr=202,)
    '''*class variable* van der Waal radius
    from webelements.com_
    
    .. _webelements.com: http://www.webelements.com/
    '''

    cov_radius = dict(H=32, Ne=71, F=72, O=73, N=75, C=77, B=82, Be=90, He=93,
                      Ar=98, Cl=99, S=102, P=106, Si=111, Kr=112, Br=114,
                      Ni=115, Se=116, Co=116, Cu=117, Fe=117, Mn=117, Al=118,
                      Cr=118, As=120, Ge=122, V=122, Li=123, Rh=125, Ru=125,
                      Zn=125, Ga=126, Os=126, Ir=127, Tc=127, Re=128, Pd=128,
                      W=130, Pt=130, Mo=130, Xe=131, Ti=132, I=133, Ta=134,
                      Nb=134, Ag=134, Au=134, Te=136, Mg=136, Sn=141, Sb=141,
                      U=142, In=144, Sc=144, Hf=144, Zr=145, At=145, Bi=146,
                      Po=147, Pb=147, Cd=148, Tl=148, Hg=149, Na=154, Tm=156,
                      Lu=156, Er=157, Ho=158, Dy=159, Tb=159, Gd=161, Y=162,
                      Sm=162, Pm=163, Nd=164, Th=165, Ce=165, Pr=165, La=169,
                      Yb=174, Ca=174, Eu=185, Sr=191, Ba=198, K=203, Rb=216,
                      Cs=235,)

    '''*class variable* covalent radius in pm from environmentalchemistry.com_
    
    .. _environmentalchemistry.com:
       http://environmentalchemistry.com/yogi/periodic/covalentradius.html
    '''

    colors = dict(H='0.5',  # gray50
                  C='k',  # black
                  N='b',  # blue
                  O='r',  # red
                  S='y',  # yellow
                  P='m',  # purple
                  F='g',  # halogens, green
                  Cl='g',
                  Br='g',
                  I='g',
                  At='g',)
    '''*class variable*'''


    def __init__(self,
                 symbol='',
                 index=0,
                 mass=1,
                 radius=None,
                 xyz=None,
                 color=None):
        self.symbol = symbol
        # keep two indexes: Python starts indexing from 0
        # but other softwares (ex Molden) index from 1
        self.data_index = index
        self.index = index
        self.mass = mass
        self.radius = (Atom.cov_radius.get(symbol, 150) if radius is None
                       else radius)
        # coordinates in pm
        self.coords = xyz if xyz is not None else np.array([0.0, 0.0, 0.0])
        # normal_coords are coordinates after displacement
        self.normal_coords = {}

        self.projected_coordinates = np.array([0.0, 0.0])
        self.color = (Atom.colors.get(symbol, "black") if color is None
                      else color)

    def __repr__(self):
        return "%s(symbol=%s, index=%r, mass=%r, radius=%r, xyz=%r, color=%r)"\
                % (self.__class__.__name__, self.symbol, self.index, self.mass,
                   self.radius, self.coords, self.color)

    @property
    def label(self):
        return "%s (%i)" % (self.symbol, self.index)

    @property
    def x(self):
        return self.coords[Atom.X]

    @x.setter
    def x(self, val):
        self.coords[Atom.X] = val

    @property
    def y(self):
        return self.coords[Atom.Y]

    @y.setter
    def y(self, val):
        self.coords[Atom.Y] = val

    @property
    def z(self):
        return self.coords[Atom.Z]

    @z.setter
    def z(self, val):
        self.coords[Atom.Z] = val

    @property
    def x_paper(self):
        return self.projected_coordinates[Atom.X]

    @x_paper.setter
    def x_paper(self, val):
        self.projected_coordinates[Atom.X] = val

    @property
    def y_paper(self):
        return self.projected_coordinates[Atom.Y]

    @y_paper.setter
    def y_paper(self, val):
        self.projected_coordinates[Atom.Y] = val



class Molecule(object):

    THRESHOLD = 1e-3  # neglect value if below

    def __init__(self, graph, vibrations):
        self.vibrations = vibrations
        self.graph = graph  # hold molecule as an undirected graph
        self.projection = None

    def __repr__(self):
        return "%s()" % self.__class__.__name__

    def dump_graph(self):
        return [(k.index, [a.index for a in v])
                for k, v in self.graph.iteritems()]

    @property
    def atoms(self):
        return sorted(self.graph, key=lambda a: a.index)

    def count_bonds(self):
        count = 0
        for children in self.graph.itervalues():
            count += len(children)
        return count / 2



class Vibration(object):

    __slots__ = ['frequency', 'intensity']

    bond_colors = arc_colors = ('b', 'r')
    oop_colors = ('g', 'y')
    threshold = 0.01     # 0.01 == 1%
    angle_threshold = 1.0  # degree

    def __init__(self, freq=0.0, intensity=1.0):
        self.frequency = freq
        self.intensity = intensity

    def __repr__(self):
        return "%s(freq=%r, intensity=%r)" % (self.__class__.__name__,
                                              self.frequency, self.intensity)

    @property
    def key(self):
        return self.mk_key(self.frequency)

    @classmethod
    def mk_key(self, val):
        # helper for sorting freq as directory keys without rounding problem
        return int(round(val * 100.0))

