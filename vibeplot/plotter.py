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


class MoleculePlotter(object):
    """Use `matplotlib` to draw a molecule.
    
    Attributes
    ----------
    molecule : `vibeplot.chemistry.Molecule`
    fig : `matplotlib.Figure`
    axes : `matplotlib.Axes`
    show_atom_index : bool
        If True, the index of the atom is written next to its symbol.
    black_and_white : bool
        It True, no colors are used.
    fontsize : int
        The font size used to write the atomic labels.
    linewidth : float
        The linewidth used to draw the molecule.
    padding : float

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
        """Convenience function used to plot the molecule."""
        self.__initAxes()
        if (self.molecule is None) or (self.linewidth == 0.0): return
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
        """Convenience function used to save the drawing."""
        for artist in self._animated:
            artist.set_animated(False)
        self.fig.savefig(filename, dpi=300)
        for artist in self._animated:
            artist.set_animated(True)

    def atom_label_text_iter(self, show_index=False, **kwargs):
        """Generate atomic labels.

        Returns
        -------
        x, y : float
            Position of the label.
        label : string

        kw : dict
            Other arguments to pass to `self.axes.text`.

        """
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
        """Return `PathCollection` representing atoms as circles."""
        col = []
        colors = []
        for atom in self.molecule.graph:
            color = atom.colors.get(atom.color, 'k')
            colors.append(color)
            if radius is None:
                radius = atom.radius
            circle = Circle(atom.projected_coordinates, radius)
            col.append(circle)
        kw = {'facecolors': colors, 'edgecolors': colors}
        kw.update(kwargs)
        return PatchCollection(col, **kw)

    def get_bond_collection(self, **kwargs):
        """Return `PathCollection` representing atomic bonds as segments."""
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
    """Use `matplotlib` to draw the vibration markers.
    
    Attributes
    ----------
    scaling_factor : float
        Scale the amplitude of the markers.
    oop_curve_type : {3, 4}
        Use either 3-points or 4-points bezier to represent bond
        torsions.
    bond_colors, arc_colors, oop_colors : tuple of two matplotlib colors
    threshold : float
        If the amplitude of the change is below the threshold, the
        marker is not drawn.
        
    """
    def __init__(self):
        super(VibrationPlotter, self).__init__()
        self.scaling_factor = 10
        self.oop_curve_type = 4
        self.bond_colors = self.arc_colors = ("b", "r")
        self.oop_colors = ("g", "y")
        self.threshold = 0.0

    def plot_vibration(self, freq, **kwargs):
        """Convenience function used to plot the markers.

        Parameters
        ----------
        freq : float
            The frequency of the vibration to show.
        kwargs : dict
            Keyword arguments forwarded to matplotlib.

        """
        while self._animated:
            artist = self._animated.pop()
            try:
                artist.remove()
            except ValueError:
                pass
        linewidth = self.linewidth if self.linewidth != 0.0 else 1.0
        self._animated = [
            self.get_bondlength_change_collection(
                freq, factor=self.scaling_factor, zorder=20,
                lw=linewidth, **kwargs),
            self.get_angle_change_collection(
                freq, factor=self.scaling_factor, zorder=25,
                lw=linewidth, **kwargs),
            self.get_oop_angle_change_collection(
                freq, factor=self.scaling_factor, zorder=50,
                lw=linewidth,
                CURVE_TYPE=self.oop_curve_type, **kwargs)]
        for collection in self._animated:
            self.axes.add_collection(collection)
            self.axes.draw_artist(collection)

    def get_bondlength_change_collection(self, freq, factor=10.0, **kwargs):
        """Return PathCollection of bondlength change markers (lines)"""

        def bond_length_change(freq, path):
            atom1, atom2 = path
            # distance at equilibrium
            equil  = va.distance(atom1.coords, atom2.coords)
            # distance after displacement
            length = va.distance(atom1.normal_coords[freq],
                                atom2.normal_coords[freq])
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
            amplitude = bond_length_change(freq, path)
            if abs(amplitude) <= self.threshold: continue
            amp.append(abs(amplitude) * 5.0 * factor)
            colors.append(self.bond_colors[0 if amplitude < 0.0 else 1])
            verts = atom1.projected_coordinates, atom2.projected_coordinates
            segment = Path(verts, codes)
            col.append(segment)
        lw = 0.0 if not col else np.array(amp) * kwargs.pop("lw", 1.0)
        kw = {'edgecolors': colors, 'linewidths': lw}
        kw.update(kwargs)
        return PathCollection(col, **kw)

    def get_angle_change_collection(self, freq, factor=10.0, **kwargs):
        """Return `PathCollection` of angle change markers (arcs)"""

        def angle_change(freq, path):
            atom1, vertex, atom2 = path
            # angle at equilibrium
            vector1 = atom1.coords - vertex.coords
            vector2 = atom2.coords - vertex.coords
            equil = va.dangle(vector1, vector2)
            # angle after displacement
            vector1 = (atom1.normal_coords[freq] -
                    vertex.normal_coords[freq])
            vector2 = (atom2.normal_coords[freq] -
                    vertex.normal_coords[freq])
            angle = va.dangle(vector1, vector2)
            # difference
            return angle - equil

        col = []
        colors = []
        for path in make_paths(self.molecule.graph, order=3):
            atom1, vertex, atom2 = path
            xy = vertex.projected_coordinates
            amplitude = angle_change(freq, path)
            if abs(amplitude) <= self.threshold: continue
            width = height = abs(amplitude) / 180.0 * factor
            angle = 0.0
            d1 = atom1.projected_coordinates - vertex.projected_coordinates
            d2 = atom2.projected_coordinates - vertex.projected_coordinates
            theta1 = va.dangle2d(np.array([1.0, 0.0]), d1)
            theta2 = va.dangle2d(np.array([1.0, 0.0]), d2)
            # always plot smaller arc [ 0.0, 180.0 [
            if (theta2 - theta1 + 360.0) % 360.0 > 180.0:
                theta2, theta1 = theta1, theta2
            color = self.arc_colors[0 if amplitude < 0.0 else 1]
            colors.append(color)
            arc = Arc(xy, width, height, angle, theta1, theta2)
            col.append(arc)
        kw = {'edgecolors': colors, 'facecolors': 'none'}
        kw.update(kwargs)
        return PatchCollection(col, **kw)

    def get_oop_angle_change_collection(self, freq, factor=10.0,
                                        CURVE_TYPE=4, **kwargs):
        """Return `PathCollection` of bond torsion markers (beziers)."""

        def oop_angle_change(freq, path):
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
                atom_mid1.normal_coords[freq],
                atom_mid2.normal_coords[freq],
                atom_ext1.normal_coords[freq])
            plane2 = Plane(
                atom_mid1.normal_coords[freq],
                atom_mid2.normal_coords[freq],
                atom_ext2.normal_coords[freq])
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
            amplitude = oop_angle_change(freq, path)
            # amplitude in degree
            damplitude = va.rad2deg(amplitude)
            if damplitude > 180.0:
                damplitude -= 360.0
            if abs(damplitude) <= self.threshold: continue
            intensity = abs(amplitude) / np.pi * 0.25 * factor
            p2 = 0.5 * (atom_pos[2] + atom_pos[1])  # middle
            p1 = intensity * (atom_pos[0] - p2) + p2
            p3 = intensity * (atom_pos[3] - p2) + p2
            color = self.oop_colors[0 if amplitude < 0.0 else 1]
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
    """Use `matplotlib` to draw a spectrum as dirac vectors.

    Attributes
    ----------
    spectrum : dict
        Mapping intensities to frequencies.
    fig : `matplotlib.Figure`
    axes : `matplotlib.Axes`

    """

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
        """Convenience function to plot the spectrum."""
        self.__initAxes()
        self.axes.add_collection(self.get_spectrum_collection())

    def mark_line(self, marked):
        """Place marker at the frequency `marked`."""
        marker = self.axes.lines[0]
        marker.set_xdata(marked)
        self.axes.draw_artist(marker)

    def get_spectrum_collection(self, **kwargs):
        """Return `PathCollection` representing the spectrum."""
        codes = [Path.MOVETO,
                 Path.LINETO,
                ]
        col = []
        for frequency, intensity in self.spectrum.iteritems():
            verts = [(frequency, 0.0), (frequency, intensity)]
            col.append(Path(verts, codes))
        kw = {}
        kw.update(kwargs)
        return PathCollection(col, **kw)


