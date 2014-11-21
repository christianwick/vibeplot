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

from itertools import chain

from matplotlib.patches import Circle, Arc
from matplotlib.collections import PatchCollection, PathCollection
from matplotlib.path import Path
from matplotlib.lines import Line2D
from matplotlib.text import Text

import openbabel as ob
import numpy as np

import vibeplot.utils.vecalgebra as va
import vibeplot.utils.broaden as broaden


logger = logging.getLogger("vibeplot")


class AtomText(Text):

    def __init__(self, x, y, text, index, color, **kwargs):
        super(AtomText, self).__init__(x, y, text, color, **kwargs)
        self.__text = text
        self.__index = index
        self.__color = color

    def set_black_labels(self, black=True):
        self.set_color("black" if black else self.__color)

    def show_index(self, show=True):
        self.set_text("%s(%i)" % (self.__text, self.__index)
                      if show else self.__text)


class MoleculePlotter(object):
    """Use :mod:`matplotlib` to draw a molecule.

    Attributes
    ----------
    axes : mpl.Axes
    oop_curve_type : {3, 4}
        Use either 3- or 4-points bezier to represent bond torsions.
    bond_colors, arc_colors, oop_colors : tuple
        Two matplotlib colors.

    """
    def __init__(self, axes):
        super(MoleculePlotter, self).__init__()
        self.axes = axes
        self.axes.set_xticks(())
        self.axes.set_yticks(())
        for __, spine in self.axes.spines.items():
            spine.set_visible(False)
        self.axes.figure.tight_layout()
        self.draw = self.axes.figure.canvas.draw_idle
        self.oop_curve_type = 4
        self.bond_colors = self.arc_colors = ("b", "r")
        self.oop_colors = ("g", "y")
        self.molecule = ob.OBMol()
        self._molecule2D = ob.OBMol()
        self._vib_data = ob.OBVibrationData()
        self._vib_data_lx = []
        self._mol_bonds = None
        self._mol_atoms = None
        self._mol_labels = []
        self._vib_bonds = None
        self._vib_angles = None
        self._vib_oop = None

    def _2Dcoords(self, atom):
        atom2D = self._molecule2D.GetAtom(atom.GetIdx())
        assert(atom2D.GetZ() == 0.0)
        return np.array([atom2D.GetX(), atom2D.GetY()])

    def _to_normal_coordinates(self, atom, index):

        def ar(vec):
            """Returns a numpy array with the coordinates of the vector."""
            return np.array((vec.GetX(), vec.GetY(), vec.GetZ()))

        atomnc = ob.OBAtom()
        atomnc.Duplicate(atom)
        nc = ar(atom) + self._vib_data_lx[index][atom.GetIdx() - 1]
        atomnc.SetVector(*nc)
        return atomnc

    @staticmethod
    def _sdg(molecule):
        molecule2D = ob.OBMol(molecule)
        gen2D = ob.OBOp.FindType("gen2D")
        if not gen2D:
            raise NameError("name 'gen2D' is not defined")
        gen2D.Do(molecule2D)
        assert(not molecule2D.Has3D())
        assert(molecule.NumAtoms() == molecule2D.NumAtoms())
        return molecule2D

    def _add_atom_labels(self, zorder=100, **kwargs):
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
            kw = dict(horizontalalignment="center",
                      verticalalignment="center",
                      bbox=box_props)
            kwargs.update(kw)
            label = AtomText(x, y,
                             etab.GetSymbol(atom.GetAtomicNum()), atom.GetIdx(),
                             etab.GetRGB(atom.GetAtomicNum()),
                             zorder=zorder,
                             **kwargs)
            self._mol_labels.append(label)
            self.axes.add_artist(label)

    def _add_atom_collection(self, zorder=100, **kwargs):
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
        kwargs.update(kw)
        self._mol_atoms = PatchCollection(col, zorder=zorder, **kwargs)
        self.axes.add_collection(self._mol_atoms)

    def _add_bond_collection(self, zorder=10, **kwargs):
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
        kwargs.update(kw)
        self._mol_bonds = PathCollection(col, zorder=zorder, **kwargs)
        self.axes.add_collection(self._mol_bonds)

    def _add_bondlength_change_collection(
            self, index, threshold=0, zorder=20, **kwargs):
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
                                for atom in (atom1, atom2)]
            if obbond.GetLength() == 0.0:
                logger.error(
                    "Bond between %i and %i with length %.1f ignored."
                    % (atom1.GetIdx(), atom2.GetIdx(), obbond.GetLength()))
                continue
            amplitude = atom1.GetDistance(atom2) - atom1nc.GetDistance(atom2nc)
            if abs(amplitude * 100) <= threshold: continue
            amp.append(abs(amplitude * 50))
            colors.append(self.bond_colors[0 if amplitude < 0.0 else 1])

            verts = (self._2Dcoords(atom1), self._2Dcoords(atom2))
            segment = Path(verts, codes)
            col.append(segment)
        lw = 0.0 if not col else np.array(amp) * kwargs.pop("lw", 1.0)
        kw = {'edgecolors': colors, 'linewidths': lw}
        kwargs.update(kw)
        self._vib_bonds = PathCollection(col, zorder=zorder, **kwargs)
        self.axes.add_collection(self._vib_bonds)

    def _add_angle_change_collection(
            self, index, threshold=0, zorder=25, **kwargs):
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
                for atom in (vertex, atom1, atom2)]
            amplitude = (atom1nc.GetAngle(vertexnc, atom2nc) - 
                         atom1.GetAngle(vertex, atom2))
            if abs(amplitude) <= threshold: continue
            width = height = abs(amplitude) / 20

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
        kwargs.update(kw)
        self._vib_angles = PatchCollection(col, zorder=zorder, **kwargs)
        self.axes.add_collection(self._vib_angles)

    def _add_oop_angle_change_collection(
            self, index, threshold=0, CURVE_TYPE=4, zorder=50,
            **kwargs):
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
            if abs(amplitude) <= threshold: continue
            intensity = abs(amplitude / 40)

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
        kwargs.update(kw)
        self._vib_oop = PathCollection(col, zorder=zorder, **kwargs)
        self.axes.add_collection(self._vib_oop)

    def clear(self):
        for artist in self.axes.get_children():
            try:
                artist.remove()
            except NotImplementedError:
                pass
        self._mol_labels = []
        self._mol_atoms = self._mol_bonds = None
        self._vib_bonds = self._vib_angles = self._vib_oop = None

    def set_molecule(self, molecule):
        """`openbabel.OBMol` molecule with original 3D coordinates."""
        self.molecule = molecule
        self._molecule2D = self._sdg(molecule)
        self.clear()

    def set_vibration_data(self, vib_data):
        self._vib_data = vib_data
        # profiling shows performance increase:
        self._vib_data_lx = np.array(
            [[(vec.GetX(), vec.GetY(), vec.GetZ()) for vec in row]
             for row in vib_data.GetLx()], dtype=float)
        self._vib_data_lx *= 0.529177249  # to angstroem

    def draw_molecule(self, padding=0.3, lw=1.0, fontsize=12.0):
        for artist in chain((self._mol_atoms, self._mol_bonds),
                            self._mol_labels):
            if artist:
                artist.remove()
        self._add_bond_collection(lw=lw)
        self._add_atom_labels(fontsize=fontsize)
        self.axes.ignore_existing_data_limits = True
        xmin, xmax, ymin, ymax = self.axes.axis("image")
        self.axes.axis((xmin - padding, xmax + padding,
                        ymin - padding, ymax + padding))
        self.draw()

    def draw_vibration(self, row, bl_filter, angle_filter, torsion_filter):
        for artist in (self._vib_bonds, self._vib_angles, self._vib_oop):
            if artist:
                artist.remove()
        if row is -1: return
        self._add_bondlength_change_collection(row, bl_filter)
        self._add_angle_change_collection(row, angle_filter)
        self._add_oop_angle_change_collection(row, torsion_filter)
        self.draw()

    def show_atom_index(self, show=True):
        for artist in self._mol_labels:
            artist.show_index(show)
        self.draw()

    def set_black_labels(self, black=True):
        for artist in self._mol_labels:
            artist.set_black_labels(black)
        self.draw()

    def set_fontsize(self, fontsize):
        for artist in self._mol_labels:
            artist.set_fontsize(fontsize)
        self.draw()

    def set_linewidth(self, linewidth):
        self._mol_bonds.set_linewidth(float(linewidth))
        self.draw()


class SpectrumPlotter(object):
    """Use :mod:`matplotlib` to draw a spectrum as dirac vectors.

    Attributes
    ----------
    axes : mpl.Axes
    function : {None, "lorentzian", "gaussian"}
        Choose function to broaden the spectrum.

    """
    def __init__(self, axes):
        super(SpectrumPlotter, self).__init__()
        self.axes = axes
        self.axes.set_xlabel("Wavenumber [cm$^{-1}$]")
        self.axes.axis([0, 4000, 0, 1])
        self.axes.set_yticks(())
        self.axes.figure.tight_layout()
        self.draw = self.axes.figure.canvas.draw_idle
        self.needle, = self.axes.plot((0.0, 0.0), (0.0, 1.0),
                                      color="r", lw=2.0)
        self.broadening = Line2D([], [], linewidth=1.0, color="k")
        self.axes.add_line(self.needle)
        self.axes.add_line(self.broadening)
        self._broadening_function = None
        self._fwhm = 8.0
        self._vib_data = ob.OBVibrationData()
        self._spectrum = None

    def _add_spectrum_collection(self, **kwargs):
        """
        Return :class:`~matplotlib.collections.PathCollection`
        representing the spectrum.

        """
        codes = [Path.MOVETO,
                 Path.LINETO,
                ]
        col = []
        frequencies = self._vib_data.GetFrequencies()
        intensities = np.array(self._vib_data.GetIntensities())
        if intensities.any():
            intensities /= intensities.max()
        if len(intensities) != len(frequencies):
            intensities = [1.0] * len(frequencies)
        for frequency, intensity in zip(frequencies, intensities):
            verts = [(frequency, 0.0), (frequency, intensity)]
            col.append(Path(verts, codes))
        kw = {}
        kwargs.update(kw)
        self._spectrum = PathCollection(col, **kwargs)
        self.axes.add_collection(self._spectrum)

    def clear(self):
        for collection in self.axes.collections:
            collection.remove()
        self._spectrum = None

    def set_vibration_data(self, vib_data):
        self._vib_data = vib_data
        self.clear()
        self.set_vibration("")

    def draw_spectrum(self):
        self._add_spectrum_collection(color="0.30")
        self.update_broaden()
        self.draw()

    def set_vibration(self, freq):
        try:
            self.needle.set_xdata(float(freq))
        except ValueError:
            # freq not convertible to float
            self.needle.set_visible(False)
        else:
            self.needle.set_visible(True)
            self.draw()

    def set_fwhm(self, fwhm):
        self._fwhm = fwhm
        self.update_broaden()
        self.draw()

    def set_broadening_function(self, function_name):
        self._broadening_function = dict(
            lorentzian=broaden.lorentzian,
            gaussian=broaden.gaussian).get(function_name)
        self.update_broaden()
        self.draw()

    def update_broaden(self, **kwargs):
        """
        Return :class:`~matplotlib.lines.Line2D` of the broadened
        spectrum.

        """
        if self._broadening_function:
            spkx, spky = broaden.broaden(
                self._vib_data.GetFrequencies(),
                self._vib_data.GetIntensities(),
                width=self._fwhm, xmin=0.0, xmax=4000.0,
                fun=self._broadening_function)
            self.broadening.set_data(spkx, spky/spky.max()
                                     if spky.any() else spky)
            self.broadening.set_visible(True)
        else:
            self.broadening.set_visible(False)

    def save_spectrum(self, filename):
        """Save broadened spectrum to file."""
        np.savetxt(filename, self.broadening.get_xydata())

