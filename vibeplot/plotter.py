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

"""
vibeplot.plotter module
=======================

This module contains the classes used to generate the plots.

"""

import logging
logging.basicConfig()

import matplotlib as mpl
from matplotlib.patches import Circle, Arc
from matplotlib.collections import PatchCollection, PathCollection
from matplotlib.path import Path
from matplotlib.lines import Line2D

import openbabel as ob
import numpy as np

import vibeplot.utils.vecalgebra as va
import vibeplot.utils.broaden as broaden


logger = logging.getLogger("vibeplot")


def _coords(atom):
    """Returns a numpy array with the coordinates of the atom."""
    return np.array((atom.GetX(), atom.GetY(), atom.GetZ()))


class MoleculePlotter(object):
    """Use :mod:`matplotlib` to draw a molecule.

    Attributes
    ----------
    fig : :class:`~matplotlib.figure.Figure`
    axes : :class:`~matplotlib.axes.Axes`
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
        self._molecule = ob.OBMol()
        self._molecule2D = ob.OBMol()
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

    def _2Dcoords(self, atom):
        atom2D = self._molecule2D.GetAtom(atom.GetIdx())
        assert(atom2D.GetZ() == 0.0)
        return np.array([atom2D.GetX(), atom2D.GetY()])

    @property
    def molecule(self):
        """`openbabel.OBMol` molecule with original 3D coordinates."""
        return self._molecule

    @molecule.setter
    def molecule(self, molecule):
        """`openbabel.OBMol` molecule with original 3D coordinates."""
        self._molecule = molecule
        self._molecule2D = ob.OBMol(molecule)
        gen2D = ob.OBOp.FindType("gen2D")
        if not gen2D:
            raise NameError("name 'gen2D' is not defined")
        gen2D.Do(self._molecule2D)
        assert(not self._molecule2D.Has3D())
        assert(self._molecule.NumAtoms() == self._molecule2D.NumAtoms())

    def plot_molecule(self):
        """Convenience function used to plot the molecule."""
        self.__initAxes()
        if (not self.molecule.NumAtoms() and not self.linewidth): return
        for collection in (
                self.get_bond_collection(zorder=10, lw=self.linewidth),):
            self.axes.add_collection(collection)
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
        """
        Convenience function used to save the drawing of the molecule.
        """
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
        etab = ob.OBElementTable()
        for atom in ob.OBMolAtomIter(self.molecule):
            x, y = self._2Dcoords(atom)
            label = (etab.GetSymbol(atom.GetAtomicNum()) if not show_index else
                     "%s(%i)" % (etab.GetSymbol(atom.GetAtomicNum()),
                                 atom.GetIdx()))
            kw = dict(horizontalalignment="center",
                      verticalalignment="center",
                      bbox=box_props)
            if not self.black_and_white:
                kw["color"] = etab.GetRGB(atom.GetAtomicNum())
            kw.update(kwargs)
            yield [x, y, label, kw]

    def get_atom_collection(self, **kwargs):
        """
        Return :class:`~matplotlib.collections.PathCollection`
        representing atoms as circles.

        """
        col = []
        colors = []
        etab = ob.OBElementTable()
        for atom in ob.OBMolAtomIter(self.molecule):
            colors.append(etab.GetRGB(atom.GetAtomicNum()))
            radius = etab.GetCovalentRad(atom.GetAtomicNum())
            circle = Circle(self._2Dcoords(atom), radius)
            col.append(circle)
        kw = {'facecolors': colors, 'edgecolors': colors}
        kw.update(kwargs)
        return PatchCollection(col, **kw)

    def get_bond_collection(self, **kwargs):
        """
        Return :class:`~matplotlib.collections.PathCollection`
        representing atomic bonds as segments.

        """
        col = []
        codes = [Path.MOVETO, Path.LINETO,]  # segment
        for obbond in ob.OBMolBondIter(self.molecule):
            atom1, atom2 = (self.molecule.GetAtom(obbond.GetBeginAtomIdx()),
                            self.molecule.GetAtom(obbond.GetEndAtomIdx()))
            verts = self._2Dcoords(atom1), self._2Dcoords(atom2)
            segment = Path(verts, codes)
            col.append(segment)
        kw = {'edgecolors': 'k'}
        kw.update(kwargs)
        return PathCollection(col, **kw)


class VibrationPlotter(MoleculePlotter):
    """Use :mod:`matplotlib` to draw the vibration markers.

    Attributes
    ----------
    normal_coordinates : [[`openbabel.vector3`]]
        Indexed with vibration_number and atom_index.
    scaling_factor : float
        Scale the amplitude of the markers.
    oop_curve_type : {3, 4}
        Use either 3- or 4-points bezier to represent bond torsions.
    bond_colors, arc_colors, oop_colors : tuple
        Two matplotlib colors.
    threshold : float
        Do not show marker if amplitude below `threshold`.

    """
    def __init__(self):
        super(VibrationPlotter, self).__init__()
        self.normal_coordinates = [[None]]
        self.scaling_factor = 100
        self.oop_curve_type = 4
        self.bond_colors = self.arc_colors = ("b", "r")
        self.oop_colors = ("g", "y")
        self.threshold = 0.0

    def plot_vibration(self, index, **kwargs):
        """Convenience function used to plot the markers.

        Parameters
        ----------
        index : int
            Index in `vibrations` list.
        kwargs : dict
            Keyword arguments forwarded to `matplotlib`.

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
                index, factor=self.scaling_factor, zorder=20,
                lw=linewidth, **kwargs),
            self.get_angle_change_collection(
                index, factor=self.scaling_factor, zorder=25,
                lw=linewidth, **kwargs),
            self.get_oop_angle_change_collection(
                index, factor=self.scaling_factor, zorder=50,
                lw=linewidth,
                CURVE_TYPE=self.oop_curve_type, **kwargs)]
        for collection in self._animated:
            self.axes.add_collection(collection)
            self.axes.draw_artist(collection)

    def _to_normal_coordinates(self, atom, index):
        def au2angstrom(x):
            """Convert `x` from atomic units to Angstrom."""
            return 0.529177249 * x

        def au2ar(vec):
            return au2angstrom(_coords(vec))

        nc = (_coords(atom.GetVector()) +
              au2ar(self.normal_coordinates[index][atom.GetIdx() - 1]))
        atomnc = ob.OBAtom()
        atomnc.Duplicate(atom)
        atomnc.SetVector(ob.vector3(*nc))
        return atomnc

    def get_bondlength_change_collection(self, index, factor=10.0, **kwargs):
        """
        Return :class:`~matplotlib.collections.PathCollection` of
        bondlength change markers (lines).

        """
        codes = [Path.MOVETO, Path.LINETO]
        col = []
        amp = []
        colors = []
        for obbond in ob.OBMolBondIter(self.molecule):
            atom1, atom2 = (self.molecule.GetAtom(obbond.GetBeginAtomIdx()),
                            self.molecule.GetAtom(obbond.GetEndAtomIdx()))
            atom1nc, atom2nc = [self._to_normal_coordinates(atom, index)
                                for atom in atom1, atom2]
            if obbond.GetLength() == 0.0:
                logger.error(
                    "Bond between %i and %i with length %.1f ignored."
                    % (atom1.GetIdx(), atom2.GetIdx(), obbond.GetLength()))
                continue
            amplitude = ((atom1.GetDistance(atom2) -
                          atom1nc.GetDistance(atom2nc)) /
                         obbond.GetLength())
            if abs(amplitude) <= self.threshold: continue
            amp.append(abs(amplitude) * 5.0 * factor)
            colors.append(self.bond_colors[0 if amplitude < 0.0 else 1])

            verts = (self._2Dcoords(atom1), self._2Dcoords(atom2))
            segment = Path(verts, codes)
            col.append(segment)
        lw = 0.0 if not col else np.array(amp) * kwargs.pop("lw", 1.0)
        kw = {'edgecolors': colors, 'linewidths': lw}
        kw.update(kwargs)
        return PathCollection(col, **kw)

    def get_angle_change_collection(self, index, factor=10.0, **kwargs):
        """
        Return :class:`~matplotlib.collections.PathCollection` of angle
        change markers (arcs).

        """
        col = []
        colors = []
        for angle in ob.OBMolAngleIter(self.molecule):
            vertex, atom1, atom2 = [self.molecule.GetAtom(idx + 1)
                                    for idx in angle]
            vertexnc, atom1nc, atom2nc = [
                self._to_normal_coordinates(atom, index)
                for atom in vertex, atom1, atom2]
            amplitude = (atom1nc.GetAngle(vertexnc, atom2nc) - 
                         atom1.GetAngle(vertex, atom2))
            if abs(amplitude) <= self.threshold: continue
            width = height = abs(amplitude) / 180.0 * factor

            d1, d2 = (self._2Dcoords(atom1) - self._2Dcoords(vertex),
                      self._2Dcoords(atom2) - self._2Dcoords(vertex))
            theta1 = va.dangle2d(np.array([1.0, 0.0]), d1)
            theta2 = va.dangle2d(np.array([1.0, 0.0]), d2)
            # always plot smaller arc [ 0.0, 180.0 [
            if (theta2 - theta1 + 360.0) % 360.0 > 180.0:
                theta2, theta1 = theta1, theta2
            color = self.arc_colors[0 if amplitude < 0.0 else 1]
            colors.append(color)
            arc = Arc(self._2Dcoords(vertex),
                      width, height, 0.0, theta1, theta2)
            col.append(arc)
        kw = {'edgecolors': colors, 'facecolors': 'none'}
        kw.update(kwargs)
        return PatchCollection(col, **kw)

    def get_oop_angle_change_collection(self, index, factor=10.0,
                                        CURVE_TYPE=4, **kwargs):
        """
        Return :class:`~matplotlib.collections.PathCollection` of bond
        torsion markers (beziers).

        """
        CURVE_TYPE_3, CURVE_TYPE_4 = 3, 4
        col = []
        edgecolors = []
        if CURVE_TYPE is CURVE_TYPE_3:
            codes = [Path.MOVETO, Path.CURVE3, Path.CURVE3, ]
        elif CURVE_TYPE is CURVE_TYPE_4:
            codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, ]
        for torsion in ob.OBMolTorsionIter(self.molecule):
            atoms = [self.molecule.GetAtom(idx + 1)
                     for idx in torsion]
            atomsnc = [self._to_normal_coordinates(atom, index)
                       for atom in atoms]
            teq = self.molecule.GetTorsion(*atoms)
            tnc = self.molecule.GetTorsion(*atomsnc)
            amplitude = (tnc - teq + 360.0) % 360.0
            if amplitude > 180.0:
                    amplitude -= 360.0
            if abs(amplitude) <= self.threshold: continue
            intensity = abs(amplitude / 720.0) * factor

            a, b, c, d = [self._2Dcoords(atom) for atom in atoms]
            p2 = 0.5 * (b + c)  # middle
            p1 = intensity * (a - p2) + p2
            p3 = intensity * (d - p2) + p2
            color = self.oop_colors[0 if amplitude < 0.0 else 1]
            if CURVE_TYPE is CURVE_TYPE_3:
                verts = [p1, p2, p3]
            elif CURVE_TYPE is CURVE_TYPE_4:
                verts = [p1, b, c, p3]
            curve = Path(verts, codes)
            col.append(curve)
            edgecolors.append(color)
        kw = {'edgecolors': edgecolors, 'facecolors': 'none'}
        kw.update(kwargs)
        return PathCollection(col, **kw)


class SpectrumPlotter(object):
    """Use :mod:`matplotlib` to draw a spectrum as dirac vectors.

    Attributes
    ----------
    vibrations : `openbabel.OBVibrationData`
        Mapping intensities to frequencies.
    fig : :class:`~matplotlib.figure.Figure`
    axes : :class:`~matplotlib.axes.Axes`
    broadening : {None, "lorentzian", "gaussian"}
        Choose function to broaden the spectrum.

    """
    def __init__(self):
        self.vibrations = ob.OBVibrationData()
        self.fig = mpl.figure.Figure()
        self.axes = self.fig.add_subplot(111)
        self.axes.hold(False)
        self.broadening = None
        self.width = 8.0
        self.__initAxes()

    def __initAxes(self):
        self.axes.clear()
        self.axes.plot((0.0, 0.0), (0.0, 1.0), color="r", lw=2.0)
        # settings do not stick after axes.clear()
        self.axes.set_xlabel("Wavenumber [cm$^{-1}$]")
        self.axes.axis([0, 4000, 0, 1.0])
        self.axes.set_yticks(())
        self.fig.tight_layout()

    def plot_spectrum(self):
        """Convenience function to plot the spectrum."""
        self.__initAxes()
        self.axes.add_collection(self.get_spectrum_collection(color='0.30'))
        self.axes.add_line(self.get_broaden(linewidth=1.0, color='k'))

    def mark_line(self, marked):
        """Place marker at the frequency `marked`."""
        marker = self.axes.lines[0]
        marker.set_xdata(marked)
        self.axes.draw_artist(marker)

    def get_broaden(self, **kwargs):
        """
        Return :class:`~matplotlib.lines.Line2D` of the broadened
        spectrum.

        """
        if self.broadening is None or self.broadening == "none":
            return Line2D([0.0], [0.0])
        xmin, xmax = self.axes.get_xlim()
        spkx, spky = broaden.broaden(
            self.vibrations.GetFrequencies(), self.vibrations.GetIntensities(),
            width=self.width, xmin=xmin, xmax=xmax,
            fun=dict(lorentzian=broaden.lorentzian,
                     gaussian=broaden.gaussian)[self.broadening])
        if spky.any():
            spky /= spky.max()
        return Line2D(spkx, spky, **kwargs)

    def save_spectrum(self, filename):
        """Save broadened spectrum to file."""
        np.savetxt(filename, self.get_broaden().get_xydata())

    def get_spectrum_collection(self, **kwargs):
        """
        Return :class:`~matplotlib.collections.PathCollection`
        representing the spectrum.
 
        """
        codes = [Path.MOVETO,
                 Path.LINETO,
                ]
        col = []
        frequencies = self.vibrations.GetFrequencies()
        intensities = np.array(self.vibrations.GetIntensities())
        if intensities.any():
            intensities /= intensities.max()
        if len(intensities) != len(frequencies):
            intensities = [1.0] * len(frequencies)
        for frequency, intensity in zip(frequencies, intensities):
            verts = [(frequency, 0.0), (frequency, intensity)]
            col.append(Path(verts, codes))
        kw = {}
        kw.update(kwargs)
        return PathCollection(col, **kw)


