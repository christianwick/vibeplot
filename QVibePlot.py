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

"""Qt4 graphical user interface (GUI) to vibeplot."""

import sip
for qtype in "QString QTextStream QVariant".split():
    sip.setapi(qtype, 2)

import sys
import os.path
from glob import glob
from functools import partial

# Import Qt modules
from PyQt4 import QtGui, QtCore, QtSvg
Qt = QtCore.Qt
# Import Matplotlib
import matplotlib as mpl
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
import openbabel as ob

from qvibeplot_ui import Ui_MainWindow
import vibeplot.plotter as plotter


class QVibeplot(QtGui.QMainWindow, Ui_MainWindow):

    def __init__(self):
        super(QVibeplot, self).__init__()
        self.setupUi(self)

        self.toolbar = NavigationToolbar(self.moleculeCanvas, self)
        self.rightHLayout.insertWidget(
            self.rightHLayout.indexOf(self.moleculeCanvas) + 1, self.toolbar)

        self.svgWidget = QtSvg.QSvgWidget()
        palette = QtGui.QPalette(self.svgWidget.palette())
        palette.setColor(palette.Window, QtGui.QColor("white"))
        self.svgWidget.setPalette(palette)

        self.moleculePlotter = plotter.VibrationPlotter()
        self.spectrumPlotter = plotter.SpectrumPlotter()
        self.moleculeAxes = self.moleculeCanvas.figure.add_subplot(111)
        self.spectrumAxes = self.spectrumCanvas.figure.add_subplot(111)
        self.spectrumIndex, = self.spectrumAxes.plot((0.0, 0.0), (0.0, 1.0),
                                                     color="r", lw=2.0)

        self.clearMoleculeAxes()
        self.clearSpectrumAxes()

        self._settings = QtCore.QSettings("Mathias Laurin", "QVibePlot")
        if not self._settings.contains("imageFile") or \
           not self._settings.contains("dataFile"):
            self._settings.setValue("imageFile", QtCore.QDir.homePath())
            self._settings.setValue("dataFile", QtCore.QDir.homePath())

        # Set defaults
        self.fontSizeComboBox.setCurrentIndex(
            self.fontSizeComboBox.findText(str(self.moleculePlotter.fontsize)))
        self.lineWidthComboBox.setCurrentIndex(
            self.lineWidthComboBox.findText(str(self.moleculePlotter.linewidth)))
        self.colorLabelCheckBox.setCheckState(
            Qt.Checked if self.moleculePlotter.black_and_white else Qt.Unchecked)
        self.scalingFactorSpinBox.setValue(self.moleculePlotter.scaling_factor)
        self.thresholdComboBox.setCurrentIndex(
            self.thresholdComboBox.findText(str(self.moleculePlotter.threshold)))
        # Connect widgets
        self.fontSizeComboBox.currentIndexChanged[str].connect(self.setFontSize)
        self.lineWidthComboBox.currentIndexChanged[str].connect(self.setLineWidth)
        self.colorLabelCheckBox.stateChanged.connect(self.setAllBlackAtomLabels)
        self.scalingFactorSpinBox.valueChanged.connect(self.setScalingFactor)
        self.thresholdComboBox.currentIndexChanged[str].connect(self.setThreshold)
        self.broadeningComboBox.currentIndexChanged[str].connect(
            self.setBroadeningFunction)
        self.fwhmDoubleSpinBox.valueChanged.connect(self.setFwhm)
        self.frequencyList.currentTextChanged.connect(self.setMarker)
        self.frequencyList.currentRowChanged.connect(self.setVibration)
        # Create actions
        self.spectrumCanvas.setContextMenuPolicy(Qt.ActionsContextMenu)
        self.spectrumCanvas.addAction(
            QtGui.QAction(u"save spectrum", self.spectrumCanvas,
                          triggered=self._saveSpectrum))
        # File menu
        self.fileMenu.addActions((
            QtGui.QAction(u"Open", self.fileMenu,
                          shortcut=QtGui.QKeySequence.Open,
                          triggered=self._loadFile),
            QtGui.QAction(u"Save image", self.fileMenu,
                          shortcut=QtGui.QKeySequence.Save,
                          triggered=self._saveImage),
            QtGui.QAction(u"Save image as...", self.fileMenu,
                          shortcut=QtGui.QKeySequence.SaveAs,
                          triggered=self._saveImageAs),
            QtGui.QAction(u"Quit", self.fileMenu,
                          shortcut=QtGui.QKeySequence.Quit,
                          triggered=self.close),
        ))
        # File > OrbiMol menu
        self.orbiMolDbMenu = QtGui.QMenu(u"OrbiMol DB molecules")
        self.orbiMolDbMenu.addActions([
            QtGui.QAction(
                os.path.splitext(os.path.basename(filename))[0],  # text
                self.orbiMolDbMenu,
                triggered=partial(self._loadFile, filename, "g03"))
            for filename in glob("data/orbimol/*.freq")])
        self.orbiMolDbMenu.addSeparator()
        self.orbiMolDbMenu.addAction(QtGui.QAction(
            u"About OrbiMol", self.orbiMolDbMenu,
            triggered=partial(QtGui.QMessageBox.about, self, u"About OrbiMol",
                              " ".join((
            u"""
            <p>OrbiMol is a free molecular orbital database by Patrick
            Chaquin and Franck Fuster. Laboratoire de Chimie
            Th&eacute;orique, UPMC Univ Paris 06, UMR CNRS 7616, Paris.</p>

            <p>For more information, see
            <a href="http://www.lct.jussieu.fr/pagesperso/orbimol/">OrbiMol
            </a> or Chaquin, P.; Fuster, F. Enseigner la chimie organique
            avec les orbitales Pr&eacute;sentation d'une base de
            donn&eacute;es d'orbitales mol&eacute;culaires. <i>L'Act.
            Chim.</i> <b>2012</b>, <i>369</i>, 37-44.

            """).splitlines()))))
        self.fileMenu.addSeparator()
        self.fileMenu.addMenu(self.orbiMolDbMenu)
        # View menu
        self.viewMenu.addActions((
            QtGui.QAction(u"Atom index", self.viewMenu, checkable=True,
                          triggered=self.showAtomIndex),
            QtGui.QAction(u"Show skeleton", self.viewMenu,
                          triggered=self.svgWidget.show),
        ))
        # Help menu
        self.helpMenu.addActions((
            QtGui.QAction(u"About", self.helpMenu, triggered=partial(
                    QtGui.QMessageBox.about,
                    self, "About QVibePlot", " ".join((
            u"""
            QVibePlot visualizes vibrational analysis performed by
            density functional theory calculations (DFT) in terms of
            changes of internal coordinates.
            <p>QVibePlot is written in Python and depends on matplotlib
            for the graphics and numpy for the maths. The GUI is
            written using PyQt4.</p>
            <p>Copyright (c) 2011-2013
            <a href="mailto:Mathias.Laurin+vibeplot@gmail.com"> Mathias
            Laurin</a></p>
            <p>QVibePlot 0.14 is available under the modified BSD
            License.</p>
            <p>Support the program by citing: Laurin, M.  QVibeplot: A
            Program To Visualize Molecular Vibrations in Two
            Dimensions. <i>J. Chem. Educ.</i> <b>2013</b>
            <a href="http://dx.doi.org/10.1021/ed300554z">DOI:
            10.1021/ed300554z</a>.</p>

            """).splitlines()))),
            QtGui.QAction(u"About Qt", self.helpMenu,
                          triggered=partial(QtGui.QMessageBox.aboutQt, self)),
            QtGui.QAction(u"About Open Babal", self.helpMenu,
                          triggered=partial(QtGui.QMessageBox.about, self,
                                            u"About Open Babel", " ".join((
            u"""
            <P>This program uses Open Babel.</P>
            <P>Open Babel is a chemical toolbox designed to speak the
            many languages of chemical data. It's an open,
            collaborative project allowing anyone to search, convert,
            analyze, or store data from molecular modeling, chemistry,
            solid-state materials, biochemistry, or related areas.</P>
            <P>Open Babel is released under the GNU GPL.</P>
            <P>See <a href="http://openbabel.org">openbabel.org</a> for
            more information.</P>

            """).splitlines()))),
            QtGui.QAction(u"About Matplotlib", self.helpMenu,
                          triggered=partial(QtGui.QMessageBox.about, self,
                                            u"About Matplotlib", " ".join((
            u"""
            <P>This program uses Matplotlib {0}.</P>
            <P>Matplotlib is a python 2D plotting library which
            produces publication quality figures in a variety of
            hardcopy formats and interactive environments across
            platforms.</P>
            <P>Matplotlib is released under the
            <a href="http://matplotlib.org/users/license.html">
            Matplotlib License</a>.</P>
            <P>See <a href="http://matplotlib.org">matplotlib.org</a>
            for more information.</P>

            """.format(mpl.__version__)).splitlines())))
        ))

    def _loadFile(self, filename=None, inFormat=None):
        if not filename:
            dataFile = self._settings.value("dataFile")
            filename = QtGui.QFileDialog.getOpenFileName(
                self,
                u"Open file",
                os.path.dirname(dataFile)
                if os.path.isfile(dataFile) else dataFile,
                ";;".join((
                    " ".join(("Common formats (",
                              "*.moldem *.mold *.molf",
                              "*.gal *.g92 *.g94 *.g98 *.g03 *.g09",
                              "*.acesout *.gukout *.nwo",
                              "CONTCAR POSCAR *.vasp", ")")),
                    "molden (*.molden *.mold *.molf)",
                    "Gaussian (*.gal *.g92 *.g94 *.g98 *.g03 *.g09)",
                    "ACES output (*.acesout)",
                    "GAMESS-UK (*.gukout)",
                    "NWChem output (*.nwo)",
                    "VASP (CONTCAR POSCAR *.vasp)",
                    "all files (*)")))
        if filename:
            self._settings.setValue("dataFile", os.path.dirname(filename))
        else:
            return

        if inFormat is None:
            inFormat = str(os.path.splitext(filename)[1][1:])
        if not inFormat and os.path.basename(filename).lower() in (
                "poscar", "contcar"):
            inFormat = "vasp"

        # load data
        mol = ob.OBMol()
        obconv = ob.OBConversion()
        obconv.SetInFormat(inFormat)
        obconv.ReadFile(mol, str(filename))
        if not mol.NumAtoms():
            self.statusBar().showMessage(
                "".join((
                    "Extension or file format '%s' unknown, ",
                    "see http://openbabel.org for the list of ",
                    "supported files.")) % inFormat)
        vibData = (ob.toVibrationData(mol.GetData(ob.VibrationData))
                   if mol.HasData(ob.VibrationData) else
                   ob.OBVibrationData())
        self.clearMoleculeAxes()
        self.clearSpectrumAxes()
        self.spectrumPlotter.vibrations = vibData
        self.spectrumPlotter.plot_spectrum(self.spectrumAxes)
        self.moleculePlotter.molecule = mol
        self.moleculePlotter.plot_molecule(self.moleculeAxes)
        self.moleculePlotter.normal_coordinates = vibData.GetLx()

        # reset
        self.frequencyList.clear()
        imageFile = QtCore.QFileInfo(self._settings.value("imageFile"))
        if imageFile.isFile():
            self._settings.setValue("imageFile", imageFile.path())

        # show data
        self.spectrumCanvas.draw()
        self.moleculeCanvas.draw()

        # populate frequencyList
        for freq in vibData.GetFrequencies():
            item = QtGui.QListWidgetItem()
            item.setData(Qt.DisplayRole, freq)
            self.frequencyList.addItem(item)

        # window title
        obconv.SetOutFormat("smi")
        self.setWindowTitle(obconv.WriteString(mol))

        # SVG representation
        obconv.SetOutFormat("svg")
        obconv.AddOption("C", obconv.OUTOPTIONS)  # implicit carbons
        obconv.AddOption("d", obconv.OUTOPTIONS)  # no molecule name
        obconv.AddOption("d", obconv.GENOPTIONS)  # implicit hydrogens
        self.svgWidget.load(QtCore.QByteArray(obconv.WriteString(mol)))

    def clearMoleculeAxes(self):
        for __, spine in self.moleculeAxes.spines.items():
            spine.set_visible(False)
        self.moleculeAxes.clear()
        self.moleculeAxes.set_xticks(())
        self.moleculeAxes.set_yticks(())
        self.moleculeAxes.hold(False)

    def clearSpectrumAxes(self):
        self.spectrumAxes.clear()
        self.spectrumAxes.add_line(self.spectrumIndex)
        # settings do not stick after axes.clear()
        self.spectrumAxes.set_xlabel("Wavenumber [cm$^{-1}$]")
        self.spectrumAxes.axis([0, 4000, 0, 1])
        self.spectrumAxes.set_yticks(())
        self.spectrumCanvas.figure.tight_layout()

    def _saveImage(self):
        imageFile = QtCore.QFileInfo(self._settings.value("imageFile"))
        if not imageFile.isFile():
            self._saveImageAs()
            return
        self.moleculePlotter.save_molecule(imageFile.filePath())
        self.moleculeCanvas.draw()

    def _saveImageAs(self):
        imageFile = QtCore.QFileInfo(self._settings.value("imageFile"))
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
        self._saveImage()

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
        self.spectrumPlotter.save_spectrum(filename)

    def setWindowTitle(self, text=""):
        super(QVibeplot, self).setWindowTitle(
            'QVibeplot' if not text else
            '%s - QVibeplot' % os.path.basename(text))

    def setBroadeningFunction(self, function_name):
        self.spectrumPlotter.function = function_name
        self.spectrumPlotter.update_broaden()
        self.spectrumCanvas.draw()

    def setFwhm(self, fwhm):
        self.spectrumPlotter.width = fwhm
        self.spectrumPlotter.update_broaden()
        self.spectrumCanvas.draw()

    def setMarker(self, frequency):
        try:
            frequency = float(frequency)
        except ValueError:
            # frequency == u""
            return
        self.spectrumIndex.set_xdata(frequency)
        self.spectrumAxes.draw_artist(self.spectrumIndex)
        self.spectrumCanvas.draw()

    def setVibration(self, row):
        self.clearMoleculeAxes()
        self.moleculePlotter.plot_molecule(self.moleculeAxes)
        self.moleculePlotter.plot_vibration(self.moleculeAxes, row)
        self.moleculeCanvas.draw()

    def setScalingFactor(self, scalingFactor):
        self.moleculePlotter.scaling_factor = float(scalingFactor)
        self.moleculeCanvas.draw()

    def setThreshold(self, threshold):
        self.moleculePlotter.threshold = float(threshold)
        self.moleculeCanvas.draw()

    def showAtomIndex(self, show=True):
        self.moleculePlotter.show_atom_index = show
        self.moleculeCanvas.draw()

    def setAllBlackAtomLabels(self, black=True):
        self.moleculePlotter.black_and_white = black
        self.moleculeCanvas.draw()

    def setFontSize(self, fontSize):
        self.moleculePlotter.fontsize = int(fontSize)
        self.moleculeCanvas.draw()

    def setLineWidth(self, lw):
        self.moleculePlotter.linewidth = float(lw)
        self.moleculeCanvas.draw()


def main():
    app = QtGui.QApplication(sys.argv)
    app.lastWindowClosed.connect(app.quit)
    window = QVibeplot()
    window.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
