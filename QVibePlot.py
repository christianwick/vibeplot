#!/usr/bin/env python
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
Qt4 graphical user interface (GUI) to vibeplot.
"""

import sip
for qtype in "QString QTextStream QVariant".split():
    sip.setapi(qtype, 2)

import sys
import os.path
from functools import partial
# Import Qt modules
from PyQt4 import QtGui, QtCore
Qt = QtCore.Qt
# Import Matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas

import vibeplot.plotter as plotter
from vibeplot.datafile import *
import vibeplot.chemistry as chem
import vibeplot.sdg as sdg


class SpectrumMpl(FigureCanvas):
    """Widget holding the spectrum.
    
    It presents a Qt-styled interface to `vibeplot.SpectrumPlotter`."""

    def __init__(self, parent=None):
        self._SpectrumPlotter = plotter.SpectrumPlotter()
        self.axes = self._SpectrumPlotter.axes
        self._background = None
        FigureCanvas.__init__(self, self._SpectrumPlotter.fig)
        self.setParent(parent)

    def setParent(self, parent):
        self._parent = parent
        if parent is None: return
        color = parent.palette().brush(QtGui.QPalette.Window).color()
        self._SpectrumPlotter.fig.set_facecolor("#%X%X%X" % (color.red(),
                                                             color.green(),
                                                             color.blue()))

    def parent(self):
        return self._parent

    def drawSpectrum(self):
        self._SpectrumPlotter.plot_spectrum()
        self.draw()
        self._background = self.copy_from_bbox(self.axes.bbox)
        self._oldSize = self.axes.bbox.width, self.axes.bbox.height

    def setSpectrum(self, spectrum):
        self._SpectrumPlotter.spectrum = spectrum

    def setMarker(self, frequency):
        if not frequency: return
        frequency = float(frequency)
        self._handleResize()
        self.restore_region(self._background)
        self._SpectrumPlotter.mark_line(frequency)
        self.blit(self.axes.bbox)

    def _handleResize(self):
        if self._oldSize != (self.axes.bbox.width, self.axes.bbox.height):
            self.drawSpectrum()


class MoleculeMpl(FigureCanvas):
    """Widget holding the molecule.
    
    It presents a Qt-styled interface to `vibeplot.VibrationPlotter`."""

    def __init__(self, parent=None):
        self._moleculePlotter = plotter.VibrationPlotter()
        self.axes = self._moleculePlotter.axes
        FigureCanvas.__init__(self, self._moleculePlotter.fig)
        self.setParent(parent)
        self._frequency = None
        self._background = None
        self._image_filename = None
        self.oop_curve_type = 4

    def setParent(self, parent):
        self._parent = parent
        if parent is None: return
        color = parent.palette().brush(QtGui.QPalette.Window).color()
        self._moleculePlotter.fig.set_facecolor("#%X%X%X" % (color.red(),
                                                             color.green(),
                                                             color.blue()))

    def parent(self):
        return self._parent

    def _handleResize(self):
        if self._oldSize != (self.axes.bbox.width, self.axes.bbox.height):
            self.redraw()

    def redraw(self):
        self.drawMolecule()
        if self._frequency is not None:
            self.drawVibration(self._frequency)

    def drawMolecule(self):
        self._moleculePlotter.plot_molecule()
        self.draw()
        self._background = self.copy_from_bbox(self.axes.bbox)
        self._oldSize = self.axes.bbox.width, self.axes.bbox.height

    def performSdg(self, algo, tolerance=1.1):
        sdg.sdg(self._moleculePlotter.molecule, algo, tolerance)
        try:
            self.drawMolecule()
        except ValueError:
            if tolerance < 2.0:
                tolerance += 0.1
                self.performSdg(algo, tolerance)
            else:
                raise

    def drawVibration(self, frequency):
        if not frequency: return
        self._image_filename = None  # reset
        self._frequency = float(frequency)
        self._handleResize()
        self.restore_region(self._background)
        vibration = self._moleculePlotter.molecule.vibrations[
            chem.Vibration.mk_key(self._frequency)]
        self._moleculePlotter.plot_vibration(vibration)
        self.blit(self.axes.bbox)

    def setMolecule(self, molecule):
        self._moleculePlotter.molecule = molecule
        self.performSdg(algo=sdg.oasa.denovo)
        self._frequency = None

    def setScalingFactor(self, scalingFactor):
        self._moleculePlotter.scaling_factor = float(scalingFactor)
        self.drawVibration()

    def scalingFactor(self):
        return self._moleculePlotter.scaling_factor

    def showAtomIndex(self, show=True):
        self._moleculePlotter.show_atom_index = show
        self.redraw()

    def isShowAtomIndex(self):
        return self._moleculePlotter.show_atom_index

    def setAllBlackAtomLabels(self, black=True):
        self._moleculePlotter.black_and_white = black
        self.redraw()

    def isAllBlackAtomLabels(self):
        return self._moleculePlotter.black_and_white

    def setFontSize(self, fontSize):
        self._moleculePlotter.fontsize = int(fontSize)
        self.redraw()

    def fontSize(self):
        return self._moleculePlotter.fontsize

    def setLineWidth(self, lw):
        self._moleculePlotter.linewidth = float(lw)
        self.redraw()

    def lineWidth(self):
        return self._moleculePlotter.linewidth

    def save_image(self, saveas=False):
        if (self._image_filename is None) or (saveas is True):
            filename = QtGui.QFileDialog.getSaveFileName(
                self,
                u"Save image",
                QtCore.QDir.homePath(),
                ";;".join(("pdf files (*.pdf)",
                           "raster images (*.png *.jpeg *.tiff)",
                           "vector images (*.pdf *.eps *.ps)",
                           "all files (*)",)))
            if not filename: return
            if "." not in filename:
                filename += ".pdf"
            self._image_filename = filename
        self._moleculePlotter.save_molecule(self._image_filename)
        self.redraw()


class MainWindow(QtGui.QMainWindow):

    def __init__(self):
        QtGui.QMainWindow.__init__(self)

        self.__initUI()
        self.__initMenuBar()

    def __initUI(self):
        self.setCentralWidget(QtGui.QWidget(self))
        mainLayout = QtGui.QGridLayout(self.centralWidget())
        self.centralWidget().setLayout(mainLayout)

        self.molecule_window = MoleculeMpl(self)
        self.molecule_window.setMinimumHeight(400)
        self.molecule_window.setMinimumWidth(400)
        mainLayout.addWidget(self.molecule_window, 0, 1)

        self.spectrum_window = SpectrumMpl(self)
        self.spectrum_window.setMinimumHeight(200)
        mainLayout.addWidget(self.spectrum_window, 1, 1)

        self._controlBox = QtGui.QGroupBox("Image", self)
        layout = QtGui.QFormLayout(self._controlBox)
        self._controlBox.setLayout(layout)

        fontSizeCombo = QtGui.QComboBox(self._controlBox)
        fontSizeCombo.addItems("6 8 10 12 14 18 20 24 32".split())
        currentIndex = \
                fontSizeCombo.findText(str(self.molecule_window.fontSize()))
        if currentIndex == -1:
            currentIndex = fontSizeCombo.findText("1.0")
        fontSizeCombo.setCurrentIndex(currentIndex)
        fontSizeCombo.currentIndexChanged[str].connect(
            partial(self.molecule_window.setFontSize))

        lineWidthCombo = QtGui.QComboBox(self._controlBox)
        lineWidthCombo.addItems("0.2 0.5 1.0 2.0 4.0".split())
        currentIndex = \
                lineWidthCombo.findText(str(self.molecule_window.lineWidth()))
        if currentIndex == -1:
            currentIndex = lineWidthCombo.findText("12")
        lineWidthCombo.setCurrentIndex(currentIndex)
        lineWidthCombo.currentIndexChanged[str].connect(
            partial(self.molecule_window.setLineWidth))

        colorLabelCheckBox = QtGui.QCheckBox(self._controlBox)
        colorLabelCheckBox.setCheckState(
            Qt.Checked if self.molecule_window.isAllBlackAtomLabels() else
            Qt.Unchecked)
        colorLabelCheckBox.stateChanged.connect(
            partial(self.molecule_window.setAllBlackAtomLabels))

        for label, editor in (("font size", fontSizeCombo),
                              ("linewidth", lineWidthCombo),
                              ("black labels", colorLabelCheckBox),
                             ):
            layout.addRow(QtGui.QLabel(label), editor)
        mainLayout.addWidget(self._controlBox, 0, 0,)

        self.frequency_list = QtGui.QListWidget(self)
        self.frequency_list.setSortingEnabled(True)
        self.frequency_list.currentTextChanged.connect(partial(
            self.spectrum_window.setMarker))
        self.frequency_list.currentTextChanged.connect(partial(
            self.molecule_window.drawVibration))

        mainLayout.addWidget(self.frequency_list, 1, 0)
        mainLayout.setColumnStretch(1, 1)

    def __initMenuBar(self):
        menuBar = QtGui.QMenuBar(self)
        self._fileMenu = QtGui.QMenu("File")
        for text, shortcut, callback in (
                ("Open", QtGui.QKeySequence.Open,
                 self._loadFile),
                ("Save image", QtGui.QKeySequence.Save,
                 self.molecule_window.save_image),
                ("Save image as...", "",
                 partial(self.molecule_window.save_image, dict(saveas=True))),
                ("Quit", QtGui.QKeySequence.Quit, self.close)):
            action = QtGui.QAction(self._fileMenu)
            action.setText(text)
            action.setShortcut(shortcut)
            action.triggered.connect(callback)
            self._fileMenu.addAction(action)
        menuBar.addMenu(self._fileMenu)

        self._viewMenu = QtGui.QMenu("View")
        for text, callback in (("Atom index",
                                partial(self.molecule_window.showAtomIndex)),):
            action = QtGui.QAction(self._viewMenu)
            action.setText(text)
            action.triggered.connect(callback)
            action.setCheckable(True)
            self._viewMenu.addAction(action)
        menuBar.addMenu(self._viewMenu)

        self._sdgMenu = QtGui.QMenu("SDG")
        for text, algo in (
                ("RDKit - 2D->3D conversion", sdg.rdkit.conversion),
                ("RDKit - de novo", sdg.rdkit.denovo),
                ("oasa - de novo", sdg.oasa.denovo)):
            callback = partial(self.molecule_window.performSdg, algo, 1.1)
            action = QtGui.QAction(self._sdgMenu)
            action.setText(text)
            action.triggered.connect(callback)
            self._sdgMenu.addAction(action)
        menuBar.addMenu(self._sdgMenu)

        self._helpMenu = QtGui.QMenu("help")
        for text, callback in (
                ("About", partial(
                    QtGui.QMessageBox.about,
                    self, "About QVibePlot", " ".join(("""
                    QVibePlot visualizes vibrational analysis performed by
                    density functional theory calculations (DFT) in terms of
                    changes of internal coordinates.
                    <p>QVibePlot is written in Python and depends on matplotlib
                    for the graphics and numpy for the maths. The GUI is
                    written using PyQt4.
                    <p>Copyright (c) 2011, 2012
                    <a href="mailto:Mathias.Laurin+vibeplot@gmail.com"> Mathias
                    Laurin</a>, FA Universitaet Erlangen Nuernberg.
                    <p>QVibePlot 0.13 is available under the Python Software
                    Foundation License.

                    """).splitlines()))),
                ("About Qt", partial(QtGui.QMessageBox.aboutQt, self))):
            action = QtGui.QAction(self._helpMenu)
            action.setText(text)
            action.triggered.connect(callback)
            self._helpMenu.addAction(action)
        menuBar.addMenu(self._helpMenu)
        self.setMenuBar(menuBar)

    def setWindowTitle(self, filename=None):
        super(MainWindow, self).setWindowTitle(
            'QVibeplot' if filename is None else
            '%s - QVibeplot' % os.path.basename(filename))
 
    def _loadFile(self, checked=None):
        if checked is None: return
        filename = QtGui.QFileDialog.getOpenFileName(
            self,
            u"Open file",
            QtCore.QDir.homePath(),
            ";;".join((
                "molden files (*.molden *.mold *.molf *input *.moinput)",
                "all files (*)")))
        if not filename: return
        self.setWindowTitle(filename)
        extension = os.path.splitext(filename)[1]
        filetype = {'.molden': MoldenFile}
        datafile = filetype.get(extension, MoldenFile)(filename)
        self.molecule = datafile.parse()
        self.spectrum_window.setSpectrum(self.molecule.vibrations)
        self.molecule_window.setMolecule(self.molecule)
        # reset
        self.frequency_list.clear()
        # show data
        self.spectrum_window.drawSpectrum()
        # populate frequency_list
        for vib in self.molecule.vibrations.itervalues():
            item = QtGui.QListWidgetItem()
            item.setData(Qt.DisplayRole, vib.frequency)
            self.frequency_list.addItem(item)


def main():
    app = QtGui.QApplication(sys.argv)
    window = MainWindow()
    window.show()

    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
