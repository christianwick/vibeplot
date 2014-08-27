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


import matplotlib as mpl
from matplotlib.patches import Circle, Arc
from matplotlib.collections import PatchCollection, PathCollection
from matplotlib.path import Path

import numpy as np

import utils.vecalgebra as va
from utils import Plane
from utils.graph import make_paths

from chemistry import *


class MoleculePlotter(object):
    """
    do the maths and return matplotlib artists to be presented in a figure
    """

    def __init__(self):
        self.molecule = None
        self.fig = mpl.figure.Figure()
        self.fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
        self.axes = self.fig.add_subplot(111)
        for loc, spine in self.axes.spines.items():
            spine.set_visible(False)
        self.__initAxes()
        self.show_atom_index = False
        self.black_and_white = False
        self.fontsize = 12
        self.linewidth = 1.0
        self.padding = 0.3
        self._animated = []

    def __initAxes(self):
        self.axes.clear()
        self.axes.set_xticks(())
        self.axes.set_yticks(())
        self.axes.hold(False)
        self._animated = []  # cleared

    def plot_molecule(self):
        if self.molecule is None:
            return
        self.__initAxes()
        self.axes.add_collection(self.get_bond_collection(zorder=10,
                                                          lw=self.linewidth))
        for label in self.atom_label_text_iter(
                self.show_atom_index,
                zorder=100,
                size=self.fontsize):
            self.axes.text(*label)
        # padding
        self.axes.axis("image")
        self.axes.axis("equal")
        xmin, xmax, ymin, ymax = self.axes.axis()
        self.axes.axis((xmin - self.padding, xmax + self.padding,
                        ymin - self.padding, ymax + self.padding))

    def save_molecule(self, filename):
        # savefig does not save animated artists
        for artist in self._animated:
            artist.set_animated(False)
        self.fig.savefig(filename, dpi=300)
        for artist in self._animated:
            artist.set_animated(True)

    def atom_label_text_iter(self, show_index=False, **kwargs):
        """return list usable as argument to mpl.axes.text"""
        box_props = dict(boxstyle='round', facecolor='white', edgecolor='none')
        for atom in self.molecule.graph:
            x = atom.x_paper
            y = atom.y_paper
            label = atom.label if show_index else atom.symbol
            kw = dict(horizontalalignment="center",
                      verticalalignment="center",
                      bbox=box_props)
            if not self.black_and_white:
                kw["color"] = atom.color
            kw.update(kwargs)
            yield [x, y, label, kw]

    def get_atom_collection(self, radius=None, **kwargs):
        """return mpl.collection of atoms"""
        col = []
        colors = []
        for atom in self.molecule.graph:
            color = Atom.colors.get(atom.color, 'k')
            colors.append(color)
            if radius is None:
                radius = atom.radius
            circle = Circle(atom.projected_coordinates, radius)
            col.append(circle)
        kw = {'facecolors': colors, 'edgecolors': colors}
        kw.update(kwargs)
        return PatchCollection(col, **kw)

    def get_bond_collection(self, **kwargs):
        """return mpl.collection of bonds"""
        col = []
        codes = [Path.MOVETO,
                 Path.LINETO,  # segment
                ]
        for path in make_paths(self.molecule.graph, order=2):
            atom1, atom2 = path
            verts = atom1.projected_coordinates, atom2.projected_coordinates
            segment = Path(verts, codes)
            col.append(segment)
        kw = {'edgecolors': 'k'}
        kw.update(kwargs)
        return PathCollection(col, **kw)


class VibrationPlotter(MoleculePlotter):

    def __init__(self):
        super(VibrationPlotter, self).__init__()
        self.scaling_factor = 10
        self.oop_curve_type = 4

    def plot_vibration(self, vibration):
        while self._animated:
            artist = self._animated.pop()
            try:
                artist.remove()
            except ValueError:
                pass
        self._animated = [
            self.get_bondlength_change_collection(
                vibration, factor=self.scaling_factor, zorder=20,
                lw=self.linewidth),
            self.get_angle_change_collection(
                vibration, factor=self.scaling_factor, zorder=25,
                lw=self.linewidth),
            self.get_oop_angle_change_collection(
                vibration, factor=self.scaling_factor, zorder=50,
                lw=self.linewidth,
                CURVE_TYPE=self.oop_curve_type)]
        for collection in self._animated:
            collection.set_animated(True)
            self.axes.add_collection(collection)
            self.axes.draw_artist(collection)

    def get_bondlength_change_collection(self, vib, factor=10.0, **kwargs):
        """return mpl.collection of bondlength change markers (lines)"""

        def bond_length_change(vib, path):
            atom1, atom2 = path
            # distance at equilibrium
            equil  = va.distance(atom1.coords, atom2.coords)
            # distance after displacement
            length = va.distance(atom1.normal_coords[vib.frequency],
                                atom2.normal_coords[vib.frequency])
            # normalized difference
            return (length - equil) / equil

        codes = [Path.MOVETO,
                 Path.LINETO,
                ]
        col = []
        amp = []
        colors = []
        for path in make_paths(self.molecule.graph, order=2):
            atom1, atom2 = path
            amplitude = bond_length_change(vib, path)
            if abs(amplitude) < Vibration.threshold: continue
            amp.append(abs(amplitude) * factor)
            colors.append(Vibration.bond_colors[0 if amplitude < 0.0 else 1])
            verts = atom1.projected_coordinates, atom2.projected_coordinates
            segment = Path(verts, codes)
            col.append(segment)
        lw = 0.0 if not col else np.array(amp) * kwargs.pop("lw", 1.0)
        kw = {'edgecolors': colors, 'linewidths': lw}
        kw.update(kwargs)
        return PathCollection(col, **kw)

    def get_angle_change_collection(self, vib, factor=10.0, **kwargs):
        """return mpl.collection of angle change markers (arcs)"""

        def angle_change(vib, path):
            atom1, vertex, atom2 = path
            # angle at equilibrium
            vector1 = atom1.coords - vertex.coords
            vector2 = atom2.coords - vertex.coords
            equil = va.dangle(vector1, vector2)
            # angle after displacement
            vector1 = (atom1.normal_coords[vib.frequency] -
                    vertex.normal_coords[vib.frequency])
            vector2 = (atom2.normal_coords[vib.frequency] -
                    vertex.normal_coords[vib.frequency])
            angle = va.dangle(vector1, vector2)
            # difference
            return angle - equil

        col = []
        colors = []
        for path in make_paths(self.molecule.graph, order=3):
            atom1, vertex, atom2 = path
            xy = vertex.projected_coordinates
            amplitude = angle_change(vib, path)
            if abs(amplitude) < Vibration.angle_threshold: continue
            width = height = abs(amplitude) / 180.0 * factor
            angle = 0.0
            d1 = atom1.projected_coordinates - vertex.projected_coordinates
            d2 = atom2.projected_coordinates - vertex.projected_coordinates
            theta1 = va.dangle2d(np.array([1.0, 0.0]), d1)
            theta2 = va.dangle2d(np.array([1.0, 0.0]), d2)
            # always plot smaller arc [ 0.0, 180.0 [
            if (theta2 - theta1 + 360.0) % 360.0 > 180.0:
                theta2, theta1 = theta1, theta2
            color = Vibration.arc_colors[0 if amplitude < 0.0 else 1]
            colors.append(color)
            arc = Arc(xy, width, height, angle, theta1, theta2)
            col.append(arc)
        kw = {'edgecolors': colors, 'facecolors': 'none'}
        kw.update(kwargs)
        return PatchCollection(col, **kw)

    def get_oop_angle_change_collection(self, vib, factor=10.0,
                                        CURVE_TYPE=4, **kwargs):
        """return mpl.collection of bond torsion markers (beziers)

        possible values to CURVE_TYPE are 3 and 4"""

        def oop_angle_change(vib, path):
            atom_ext1, atom_mid1, atom_mid2, atom_ext2 = path
            # dihedral angle at equilibrium
            plane1 = Plane(
                atom_mid1.coords,
                atom_mid2.coords,
                atom_ext1.coords)
            plane2 = Plane(
                atom_mid1.coords,
                atom_mid2.coords,
                atom_ext2.coords)
            equil = plane1.dihedral_angle(plane2)
            # dihedral angle after displacement
            plane1 = Plane(
                atom_mid1.normal_coords[vib.frequency],
                atom_mid2.normal_coords[vib.frequency],
                atom_ext1.normal_coords[vib.frequency])
            plane2 = Plane(
                atom_mid1.normal_coords[vib.frequency],
                atom_mid2.normal_coords[vib.frequency],
                atom_ext2.normal_coords[vib.frequency])
            oop = plane1.dihedral_angle(plane2)
            # difference
            return oop - equil

        CURVE_TYPE_3, CURVE_TYPE_4 = 3, 4
        col = []
        edgecolors = []
        if CURVE_TYPE is CURVE_TYPE_3:
            codes = [Path.MOVETO,
                     Path.CURVE3,
                     Path.CURVE3,
                    ]
        elif CURVE_TYPE is CURVE_TYPE_4:
            codes = [Path.MOVETO,
                     Path.CURVE4,
                     Path.CURVE4,
                     Path.CURVE4,
                    ]
        for path in make_paths(self.molecule.graph, order=4):
            atom_pos = [atom.projected_coordinates for atom in path]
            amplitude = oop_angle_change(vib, path)
            # amplitude in degree
            damplitude = va.rad2deg(amplitude)
            if damplitude > 180.0:
                damplitude -= 360.0
            if abs(damplitude) < Vibration.angle_threshold: continue
            intensity = abs(amplitude) / np.pi * 0.25 * factor
            p2 = 0.5 * (atom_pos[2] + atom_pos[1])  # middle
            p1 = intensity * (atom_pos[0] - p2) + p2
            p3 = intensity * (atom_pos[3] - p2) + p2
            color = Vibration.oop_colors[0 if amplitude < 0.0 else 1]
            if CURVE_TYPE is CURVE_TYPE_3:
                verts = [p1, p2, p3]
            elif CURVE_TYPE is CURVE_TYPE_4:
                verts = [p1, atom_pos[1], atom_pos[2], p3]
            curve = Path(verts, codes)
            col.append(curve)
            edgecolors.append(color)
        kw = {'edgecolors': edgecolors, 'facecolors': 'none'}
        kw.update(kwargs)
        return PathCollection(col, **kw)


class SpectrumPlotter(object):

    def __init__(self):
        self.spectrum = None
        self.fig = mpl.figure.Figure()
        self.axes = self.fig.add_subplot(111)
        self.axes.hold(False)
        self.__initAxes()

    def __initAxes(self):
        self.axes.clear()
        self.axes.plot((0.0, 0.0), (0.0, 1.0), color="r", lw=2.0)
        # settings do not stick after axes.clear()
        self.axes.set_xlabel("Wavenumber [cm$^{-1}$]")
        self.axes.axis([0, 4000, 0, 1.2])
        self.axes.set_yticks(())
        self.fig.tight_layout()

    def plot_spectrum(self):
        """Plot spectrum on self.axes."""
        self.__initAxes()
        self.axes.add_collection(self.get_spectrum_collection())

    def mark_line(self, marked):
        marker = self.axes.lines[0]
        marker.set_xdata(marked)
        self.axes.draw_artist(marker)

    def get_spectrum_collection(self, **kwargs):
        """return mpl.collection of lines representing the calculated
        spectrum"""
        codes = [Path.MOVETO,
                 Path.LINETO,
                ]
        col = []
        for vib in self.spectrum.itervalues():
            verts = [(vib.frequency, 0.0), (vib.frequency, vib.intensity)]
            col.append(Path(verts, codes))
        kw = {}
        kw.update(kwargs)
        return PathCollection(col, **kw)


