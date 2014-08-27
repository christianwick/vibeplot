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
import openbabel as ob

import vibeplot.plotter as plotter


class SpectrumMpl(FigureCanvas):
    """Widget holding the spectrum."""

    def __init__(self, parent=None):
        self._SpectrumPlotter = plotter.SpectrumPlotter()
        self.axes = self._SpectrumPlotter.axes
        self._oldSize = None
        self._background = None
        FigureCanvas.__init__(self, self._SpectrumPlotter.fig)
        self.setParent(parent)
        self.setContextMenuPolicy(Qt.ActionsContextMenu)

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

    def setVibrations(self, vibrations):
        self._SpectrumPlotter.vibrations = vibrations

    def setMarker(self, frequency):
        if not frequency: return
        frequency = float(frequency)
        self._handleResize()
        self.restore_region(self._background)
        self._SpectrumPlotter.mark_line(frequency)
        self.blit(self.axes.bbox)

    def setBroadeningFunction(self, function_name):
        self._SpectrumPlotter.broadening = function_name
        self.drawSpectrum()

    def setFwhm(self, fwhm):
        self._SpectrumPlotter.width = fwhm
        self.drawSpectrum()

    def saveSpectrum(self, filename):
        self._SpectrumPlotter.save_spectrum(filename)

    def _handleResize(self):
        if self._oldSize != (self.axes.bbox.width, self.axes.bbox.height):
            self.drawSpectrum()


class MoleculeMpl(FigureCanvas):
    """Widget holding the molecule."""

    def __init__(self, parent=None):
        self._moleculePlotter = plotter.VibrationPlotter()
        self.axes = self._moleculePlotter.axes
        FigureCanvas.__init__(self, self._moleculePlotter.fig)
        self.setParent(parent)
        self._row = -1
        self._oldSize = None
        self._background = None
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
        self.drawVibration(self._row)

    def drawMolecule(self):
        self._moleculePlotter.plot_molecule()
        self.draw()
        self._background = self.copy_from_bbox(self.axes.bbox)
        self._oldSize = self.axes.bbox.width, self.axes.bbox.height

    def drawVibration(self, row):
        if row == -1: return
        self._row = row
        self._handleResize()
        self.restore_region(self._background)
        self._moleculePlotter.plot_vibration(row, animated=True)
        self.blit(self.axes.bbox)

    def setMolecule(self, molecule):
        self._moleculePlotter.molecule = molecule
        self._row = -1

    def setNormalCoords(self, normCoords):
        self._moleculePlotter.normal_coordinates = normCoords

    def setScalingFactor(self, scalingFactor):
        self._moleculePlotter.scaling_factor = float(scalingFactor)
        self.redraw()

    def scalingFactor(self):
        return self._moleculePlotter.scaling_factor

    def setThreshold(self, threshold):
        self._moleculePlotter.threshold = float(threshold)
        self.redraw()

    def threshold(self):
        return self._moleculePlotter.threshold

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

    def saveImage(self, filename):
        self._moleculePlotter.save_molecule(filename)
        self.redraw()


class MainWindow(QtGui.QMainWindow):

    def __init__(self):
        QtGui.QMainWindow.__init__(self)

        self.__initUI()
        self.__initMenuBar()
        self.__initStatusBar()

        self._settings = QtCore.QSettings("Mathias Laurin", "QVibePlot")
        if not self._settings.contains("imageFile") or \
           not self._settings.contains("dataFile"):
            self._settings.setValue("imageFile", QtCore.QDir.homePath())
            self._settings.setValue("dataFile", QtCore.QDir.homePath())

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

        controlBoxLayout = QtGui.QVBoxLayout()
        mainLayout.addLayout(controlBoxLayout, 0, 0)

        self._moleculeControlBox = QtGui.QGroupBox("Molecule display")
        moleculeControlLayout = QtGui.QFormLayout(self._moleculeControlBox)
        self._moleculeControlBox.setLayout(moleculeControlLayout)

        self._vibrationControlBox = QtGui.QGroupBox("Vibration display")
        vibrationControlLayout = QtGui.QFormLayout(self._vibrationControlBox)
        self._vibrationControlBox.setLayout(vibrationControlLayout)

        self._spectrumControlBox = QtGui.QGroupBox("Spectrum display")
        spectrumControlLayout = QtGui.QFormLayout(self._spectrumControlBox)
        self._spectrumControlBox.setLayout(spectrumControlLayout)

        fontSizeCombo = QtGui.QComboBox(self._moleculeControlBox)
        fontSizeCombo.addItems("6 8 10 12 14 18 20 24 32".split())
        currentIndex = \
                fontSizeCombo.findText(str(self.molecule_window.fontSize()))
        if currentIndex == -1:
            currentIndex = fontSizeCombo.findText("12")
        fontSizeCombo.setCurrentIndex(currentIndex)
        fontSizeCombo.currentIndexChanged[str].connect(
            partial(self.molecule_window.setFontSize))

        lineWidthCombo = QtGui.QComboBox(self._moleculeControlBox)
        lineWidthCombo.addItems("0.0 0.2 0.5 1.0 2.0 4.0".split())
        currentIndex = \
                lineWidthCombo.findText(str(self.molecule_window.lineWidth()))
        if currentIndex == -1:
            currentIndex = lineWidthCombo.findText("1.0")
        lineWidthCombo.setCurrentIndex(currentIndex)
        lineWidthCombo.currentIndexChanged[str].connect(
            partial(self.molecule_window.setLineWidth))

        colorLabelCheckBox = QtGui.QCheckBox(self._moleculeControlBox)
        colorLabelCheckBox.setCheckState(
            Qt.Checked if self.molecule_window.isAllBlackAtomLabels() else
            Qt.Unchecked)
        colorLabelCheckBox.stateChanged.connect(
            partial(self.molecule_window.setAllBlackAtomLabels))

        scalingFactorSpinBox = QtGui.QSpinBox(self._vibrationControlBox)
        scalingFactorSpinBox.setRange(0, 500)
        scalingFactorSpinBox.setSingleStep(5)
        scalingFactorSpinBox.setValue(self.molecule_window.scalingFactor())
        scalingFactorSpinBox.valueChanged.connect(
            self.molecule_window.setScalingFactor)

        thresholdComboBox = QtGui.QComboBox(self._vibrationControlBox)
        thresholdComboBox.addItems("0.0 0.1 0.2 1.0 2.5 5.0 10.0".split())
        thresholdComboBox.setToolTip("percent bond length/degree angles")
        currentIndex = thresholdComboBox.findText(
            str(self.molecule_window.threshold()))
        if currentIndex == -1:
            currentIndex = thresholdComboBox.findText("1.0")
        thresholdComboBox.setCurrentIndex(currentIndex)
        thresholdComboBox.currentIndexChanged[str].connect(
            partial(self.molecule_window.setThreshold))

        broadeningComboBox = QtGui.QComboBox(self._spectrumControlBox)
        broadeningComboBox.addItems("none lorentzian gaussian".split())
        broadeningComboBox.setCurrentIndex(0)
        broadeningComboBox.currentIndexChanged[str].connect(
            self.spectrum_window.setBroadeningFunction)

        fwhmDoubleSpinBox = QtGui.QDoubleSpinBox(self._spectrumControlBox)
        fwhmDoubleSpinBox.setValue(8.0)
        fwhmDoubleSpinBox.valueChanged.connect(
            self.spectrum_window.setFwhm)

        saveSpectrumAction = QtGui.QAction("save spectrum",
                                           self.spectrum_window)
        saveSpectrumAction.triggered.connect(self._saveSpectrum)
        self.spectrum_window.addAction(saveSpectrumAction)

        for label, editor in (("font size", fontSizeCombo),
                              ("linewidth", lineWidthCombo),
                              ("black labels", colorLabelCheckBox),
                             ):
            moleculeControlLayout.addRow(QtGui.QLabel(label), editor)

        for label, editor in (("zooming factor", scalingFactorSpinBox),
                              ("threshold", thresholdComboBox),):
            vibrationControlLayout.addRow(QtGui.QLabel(label), editor)

        for label, editor in (("broadening", broadeningComboBox),
                              ("FWHM", fwhmDoubleSpinBox)):
            spectrumControlLayout.addRow(QtGui.QLabel(label), editor)

        controlBoxLayout.addWidget(self._moleculeControlBox)
        controlBoxLayout.addWidget(self._vibrationControlBox)
        controlBoxLayout.addWidget(self._spectrumControlBox)
        controlBoxLayout.addStretch()

        self.frequency_list = QtGui.QListWidget(self)
        self.frequency_list.setSortingEnabled(True)
        self.frequency_list.currentTextChanged.connect(partial(
            self.spectrum_window.setMarker))
        self.frequency_list.currentRowChanged.connect(partial(
            self.molecule_window.drawVibration))

        mainLayout.addWidget(self.frequency_list, 1, 0)
        mainLayout.setColumnStretch(1, 1)

    def __initMenuBar(self):
        menuBar = self.menuBar()
        self._fileMenu = QtGui.QMenu("File")
        for text, shortcut, callback in (
                ("Open", QtGui.QKeySequence.Open,
                 self._loadFile),
                ("Save image", QtGui.QKeySequence.Save,
                 self._saveImage),
                ("Save image as...", "",
                 partial(self._saveImage, dict(saveas=True))),
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
                    <p>Copyright (c) 2011-2013
                    <a href="mailto:Mathias.Laurin+vibeplot@gmail.com"> Mathias
                    Laurin</a>
                    <p>QVibePlot 0.14 is available under the Python Software
                    Foundation License.

                    """).splitlines()))),
                ("About Qt", partial(QtGui.QMessageBox.aboutQt, self))):
            action = QtGui.QAction(self._helpMenu)
            action.setText(text)
            action.triggered.connect(callback)
            self._helpMenu.addAction(action)
        menuBar.addMenu(self._helpMenu)

    def __initStatusBar(self):
        self.statusBar()

    def setWindowTitle(self, text=""):
        super(MainWindow, self).setWindowTitle(
            'QVibeplot' if not text else
            '%s - QVibeplot' % os.path.basename(text))
 
    def _loadFile(self, checked=None):
        if checked is None: return
        dataFile = QtCore.QFileInfo(self._settings.value("dataFile"))
        filename = QtGui.QFileDialog.getOpenFileName(
            self,
            u"Open file",
            dataFile.path() if dataFile.isFile() else dataFile.filePath(),
            ";;".join((
                "molden (*.molden *.mold *.molf)",
                "Gaussian (*.gal *.g92 *.g94 *.g98 *.g03 *.g09)",
                "ACES output (*.acesout)",
                "GAMESS-UK (*.gukout)",
                "NWChem output (*.nwo)",
                "VASP (*.contcar *.popscar *.vasp)",
                "all files (*)")))
        if filename:
            dataFile.setFile(filename)
            self._settings.setValue("dataFile", dataFile.filePath())
        else:
            return
        obconv = ob.OBConversion()
        obconv.SetInFormat(str(dataFile.suffix()))
        obconv.SetOutFormat("smi")
        mol = ob.OBMol()
        obconv.ReadFile(mol, str(filename))
        if not mol.NumAtoms():
            self.statusBar().showMessage(
                "Extension or file format unknown, see http://openbabel.org.")
        self.setWindowTitle(obconv.WriteString(mol))
        vibData = (ob.toVibrationData(mol.GetData(ob.VibrationData))
                   if mol.HasData(ob.VibrationData) else
                   ob.OBVibrationData())
        self.spectrum_window.setVibrations(vibData)
        self.molecule_window.setMolecule(mol)
        self.molecule_window.setNormalCoords(vibData.GetLx())
        # reset
        self.frequency_list.clear()
        imageFile = QtCore.QFileInfo(self._settings.value("imageFile"))
        if imageFile.isFile():
            self._settings.setValue("imageFile", imageFile.path())
        # show data
        self.spectrum_window.drawSpectrum()
        self.molecule_window.drawMolecule()
        # populate frequency_list
        for freq in vibData.GetFrequencies():
            item = QtGui.QListWidgetItem()
            item.setData(Qt.DisplayRole, freq)
            self.frequency_list.addItem(item)

    def _saveImage(self, saveas=False):
        imageFile = QtCore.QFileInfo(self._settings.value("imageFile"))
        if saveas or not imageFile.isFile():
            filename = QtGui.QFileDialog.getSaveFileName(
                self,
                u"Save image",
                imageFile.path() \
                        if imageFile.isFile() else imageFile.filePath(),
                ";;".join(("pdf files (*.pdf)",
                           "raster images (*.png *.jpeg *.tiff)",
                           "vector images (*.pdf *.eps *.ps)",
                           "all files (*)",)))
            if not filename: return
            if "." not in filename:
                filename += ".pdf"
            imageFile = QtCore.QFileInfo(filename)
            self._settings.setValue("imageFile", imageFile.filePath())
        self.molecule_window.saveImage(imageFile.filePath())

    def _saveSpectrum(self):
        imageFile = QtCore.QFileInfo(self._settings.value("imageFile"))
        filename = QtGui.QFileDialog.getSaveFileName(
            self,
            u"Save spectrum values",
            imageFile.path() \
                    if imageFile.isFile() else imageFile.filePath(),
            "plain text (*.txt)")
        if not filename: return
        if "." not in filename:
            filename += ".txt"
        self.spectrum_window.saveSpectrum(filename)


def main():
    app = QtGui.QApplication(sys.argv)
    window = MainWindow()
    window.show()

    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
