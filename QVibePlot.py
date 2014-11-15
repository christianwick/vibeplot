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
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
import openbabel as ob

import vibeplot.plotter as plotter


class MplCanvas(FigureCanvas):

    def __init__(self, figure, parent=None):
        super(MplCanvas, self).__init__(figure)
        self.setParent(parent)
        if not hasattr(figure, "ax"):
            figure.ax = figure.add_subplot(111)
        self._oldSize = None
        self._background = None

    def setParent(self, parent):
        super(MplCanvas, self).setParent(parent)
        if parent:
            color = parent.palette().brush(QtGui.QPalette.Window).color()
            self.figure.set_facecolor("#%X%X%X" % (color.red(),
                                                   color.green(),
                                                   color.blue()))

    def restore_background(self):
        self.restore_region(self._background)

    def draw(self):
        super(MplCanvas, self).draw()
        self._oldSize = self.figure.ax.bbox.width, self.figure.ax.bbox.height
        self._background = self.copy_from_bbox(self.figure.ax.bbox)

    def _handleResize(self):
        if self._oldSize != (self.figure.ax.bbox.width,
                             self.figure.ax.bbox.height):
            self.draw()


class MoleculeCanvas(MplCanvas):
    """Widget holding the molecule."""

    def __init__(self, figure, parent=None):
        super(MoleculeCanvas, self).__init__(figure, parent)
        self._row = -1
        self.oop_curve_type = 4

    def draw(self):
        self.figure.plot_molecule()
        super(MoleculeCanvas, self).draw()
        self.drawVibration(self._row)

    def drawVibration(self, row):
        if row == -1: return
        self._row = row
        self._handleResize()
        self.restore_background()
        self.figure.plot_vibration(row, animated=True)
        self.blit(self.figure.ax.bbox)


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
        mainLayout = QtGui.QHBoxLayout(self.centralWidget())
        self.centralWidget().setLayout(mainLayout)
        leftLayout = QtGui.QVBoxLayout()
        rightLayout = QtGui.QVBoxLayout()
        mainLayout.addLayout(leftLayout)
        mainLayout.addLayout(rightLayout)

        self.molecule_figure = plotter.VibrationFigure()
        self.molecule_window = MoleculeCanvas(self.molecule_figure, self)
        self.molecule_window.setMinimumHeight(400)
        self.molecule_window.setMinimumWidth(400)
        rightLayout.addWidget(self.molecule_window)
        self.toolbar = NavigationToolbar(self.molecule_window, self)
        rightLayout.addWidget(self.toolbar)

        self.spectrum_figure = plotter.SpectrumFigure()
        self.spectrum_window = MplCanvas(self.spectrum_figure, self)
        self.spectrum_window.setMinimumHeight(200)
        rightLayout.addWidget(self.spectrum_window)

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
                fontSizeCombo.findText(str(self.molecule_figure.fontsize))
        if currentIndex == -1:
            currentIndex = fontSizeCombo.findText("12")
        fontSizeCombo.setCurrentIndex(currentIndex)
        fontSizeCombo.currentIndexChanged[str].connect(self.setFontSize)

        lineWidthCombo = QtGui.QComboBox(self._moleculeControlBox)
        lineWidthCombo.addItems("0.0 0.2 0.5 1.0 2.0 4.0".split())
        currentIndex = \
                lineWidthCombo.findText(str(self.molecule_figure.linewidth))
        if currentIndex == -1:
            currentIndex = lineWidthCombo.findText("1.0")
        lineWidthCombo.setCurrentIndex(currentIndex)
        lineWidthCombo.currentIndexChanged[str].connect(self.setLineWidth)

        colorLabelCheckBox = QtGui.QCheckBox(self._moleculeControlBox)
        colorLabelCheckBox.setCheckState(
            Qt.Checked if self.molecule_figure.black_and_white else
            Qt.Unchecked)
        colorLabelCheckBox.stateChanged.connect(self.setAllBlackAtomLabels)

        scalingFactorSpinBox = QtGui.QSpinBox(self._vibrationControlBox)
        scalingFactorSpinBox.setRange(0, 500)
        scalingFactorSpinBox.setSingleStep(5)
        scalingFactorSpinBox.setValue(self.molecule_figure.scaling_factor)
        scalingFactorSpinBox.valueChanged.connect(self.setScalingFactor)

        thresholdComboBox = QtGui.QComboBox(self._vibrationControlBox)
        thresholdComboBox.addItems("0.0 0.1 0.2 1.0 2.5 5.0 10.0".split())
        thresholdComboBox.setToolTip("percent bond length/degree angles")
        currentIndex = thresholdComboBox.findText(
            str(self.molecule_figure.threshold))
        if currentIndex == -1:
            currentIndex = thresholdComboBox.findText("1.0")
        thresholdComboBox.setCurrentIndex(currentIndex)
        thresholdComboBox.currentIndexChanged[str].connect(self.setThreshold)

        broadeningComboBox = QtGui.QComboBox(self._spectrumControlBox)
        broadeningComboBox.addItems("none lorentzian gaussian".split())
        broadeningComboBox.setCurrentIndex(0)
        broadeningComboBox.currentIndexChanged[str].connect(
            self.setBroadeningFunction)

        fwhmDoubleSpinBox = QtGui.QDoubleSpinBox(self._spectrumControlBox)
        fwhmDoubleSpinBox.setValue(8.0)
        fwhmDoubleSpinBox.valueChanged.connect(self.setFwhm)

        saveSpectrumAction = QtGui.QAction("save spectrum",
                                           self.spectrum_window)
        saveSpectrumAction.triggered.connect(self._saveSpectrum)
        self.spectrum_window.setContextMenuPolicy(Qt.ActionsContextMenu)
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

        leftLayout.addWidget(self._moleculeControlBox)
        leftLayout.addWidget(self._vibrationControlBox)
        leftLayout.addWidget(self._spectrumControlBox)

        self.frequency_list = QtGui.QListWidget(self)
        self.frequency_list.setSortingEnabled(True)
        self.frequency_list.currentTextChanged.connect(self.setMarker)
        self.frequency_list.currentRowChanged.connect(partial(
            self.molecule_window.drawVibration))

        self.svgWidget = QtSvg.QSvgWidget()
        palette = QtGui.QPalette(self.svgWidget.palette())
        palette.setColor(palette.Window, QtGui.QColor("white"))
        self.svgWidget.setPalette(palette)

        leftLayout.addWidget(self.frequency_list)

    def __initMenuBar(self):
        menuBar = self.menuBar()

        # File menu
        self._fileMenu = QtGui.QMenu("File")
        menuBar.addMenu(self._fileMenu)

        self._openAction = QtGui.QAction("Open", self._fileMenu)
        self._openAction.setShortcut(QtGui.QKeySequence.Open)
        self._openAction.triggered.connect(self._loadFile)
        self._fileMenu.addAction(self._openAction)

        # File > OrbiMol menu
        self._orbiMolDbMenu = QtGui.QMenu("OrbiMol DB molecules")
        self._fileMenu.addMenu(self._orbiMolDbMenu)

        for filename in glob("data/orbimol/*.freq"):
            action = QtGui.QAction(
                os.path.splitext(os.path.basename(filename))[0],
                self._orbiMolDbMenu)
            action.triggered.connect(partial(
                self._loadFile, filename, "g03"))
            self._orbiMolDbMenu.addAction(action)

        self._orbiMolDbMenu.addSeparator()
        self._aboutOrbiMol = QtGui.QAction(
            "About OrbiMol...", self._orbiMolDbMenu)
        self._aboutOrbiMol.triggered.connect(partial(
            QtGui.QMessageBox.about,
            self, "About OrbiMol", " ".join((
                """
                <p>OrbiMol is a free molecular orbital database by Patrick
                Chaquin and Franck Fuster. Laboratoire de Chimie
                Th&eacute;orique, UPMC Univ Paris 06, UMR CNRS 7616, Paris.</p>

                <p>For more information, see
                <a href="http://www.lct.jussieu.fr/pagesperso/orbimol/">OrbiMol
                </a> or Chaquin, P.; Fuster, F. Enseigner la chimie organique
                avec les orbitales Pr&eacute;sentation d'une base de
                donn&eacute;es d'orbitales mol&eacute;culaires. <i>L'Act.
                Chim.</i> <b>2012</b>, <i>369</i>, 37-44.

                """).splitlines())))
        self._orbiMolDbMenu.addAction(self._aboutOrbiMol)

        # File menu
        self._saveImageAction = QtGui.QAction("Save image", self._fileMenu)
        self._saveImageAction.setShortcut(QtGui.QKeySequence.Save)
        self._saveImageAction.triggered.connect(self._saveImage)
        self._fileMenu.addAction(self._saveImageAction)

        self._saveImageAsAction = QtGui.QAction(
            "Save image as...", self._fileMenu)
        self._saveImageAsAction.setShortcut(QtGui.QKeySequence.SaveAs)
        self._saveImageAsAction.triggered.connect(partial(
            self._saveImage, dict(saveas=True)))
        self._fileMenu.addAction(self._saveImageAsAction)

        self._quitAction = QtGui.QAction("Quit", self._fileMenu)
        self._quitAction.setShortcut(QtGui.QKeySequence.Quit)
        self._quitAction.triggered.connect(self.close)
        self._fileMenu.addAction(self._quitAction)

        # View menu
        self._viewMenu = QtGui.QMenu("View")
        menuBar.addMenu(self._viewMenu)

        self._atomIndexAction = QtGui.QAction("Atom index", self._viewMenu)
        self._atomIndexAction.setCheckable(True)
        self._atomIndexAction.triggered.connect(self.showAtomIndex)
        self._viewMenu.addAction(self._atomIndexAction)

        self._showSkeletonAction = QtGui.QAction(
            "Show skeleton", self._viewMenu)
        self._showSkeletonAction.triggered.connect(partial(
            self.svgWidget.show), False)
        self._viewMenu.addAction(self._showSkeletonAction)

        # Help menu
        self._helpMenu = QtGui.QMenu("help")
        menuBar.addMenu(self._helpMenu)

        self._aboutAction = QtGui.QAction("About", self._helpMenu)
        self._aboutAction.triggered.connect(partial(
                    QtGui.QMessageBox.about,
                    self, "About QVibePlot", " ".join(("""
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

                    """).splitlines())))
        self._helpMenu.addAction(self._aboutAction)

        self._aboutQtAction = QtGui.QAction("About Qt", self._helpMenu)
        self._aboutQtAction.triggered.connect(partial(
            QtGui.QMessageBox.aboutQt, self))
        self._helpMenu.addAction(self._aboutQtAction)

        self._aboutOpenBabelAction = QtGui.QAction("About Open Babel", self._helpMenu)
        self._aboutOpenBabelAction.triggered.connect(partial(
                    QtGui.QMessageBox.about,
                    self, "About Open Babel", " ".join(("""
                    <P>This program uses Open Babel.</P>
                    <P>Open Babel is a chemical toolbox designed to speak the
                    many languages of chemical data. It's an open,
                    collaborative project allowing anyone to search, convert,
                    analyze, or store data from molecular modeling, chemistry,
                    solid-state materials, biochemistry, or related areas.</P>
                    <P>Open Babel is released under the GNU GPL.</P>
                    <P>See <a href="http://openbabel.org">openbabel.org</a> for
                    more information.</P>

                    """).splitlines())))
        self._helpMenu.addAction(self._aboutOpenBabelAction)

        self._aboutMplAction = QtGui.QAction("About Matplotlib", self._helpMenu)
        self._aboutMplAction.triggered.connect(partial(
                    QtGui.QMessageBox.about,
                    self, "About Matplotlib", " ".join(("""
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
        self._helpMenu.addAction(self._aboutMplAction)

    def __initStatusBar(self):
        self.statusBar()

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
        self.spectrum_figure.vibrations = vibData
        self.spectrum_figure.plot_spectrum()
        self.molecule_figure.molecule = mol
        self.molecule_figure.normal_coordinates = vibData.GetLx()

        # reset
        self.frequency_list.clear()
        imageFile = QtCore.QFileInfo(self._settings.value("imageFile"))
        if imageFile.isFile():
            self._settings.setValue("imageFile", imageFile.path())

        # show data
        self.spectrum_window.draw()
        self.molecule_window.drawVibration(-1)
        self.molecule_window.draw()

        # populate frequency_list
        for freq in vibData.GetFrequencies():
            item = QtGui.QListWidgetItem()
            item.setData(Qt.DisplayRole, freq)
            self.frequency_list.addItem(item)

        # window title
        obconv.SetOutFormat("smi")
        self.setWindowTitle(obconv.WriteString(mol))

        # SVG representation
        obconv.SetOutFormat("svg")
        obconv.AddOption("C", obconv.OUTOPTIONS)  # implicit carbons
        obconv.AddOption("d", obconv.OUTOPTIONS)  # no molecule name
        obconv.AddOption("d", obconv.GENOPTIONS)  # implicit hydrogens
        self.svgWidget.load(QtCore.QByteArray(obconv.WriteString(mol)))

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
        self.molecule_figure.save_molecule(imageFile.filePath())
        self.molecule_window.draw()

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
        self.spectrum_figure.save_spectrum(filename)

    def setWindowTitle(self, text=""):
        super(MainWindow, self).setWindowTitle(
            'QVibeplot' if not text else
            '%s - QVibeplot' % os.path.basename(text))

    def setBroadeningFunction(self, function_name):
        self.spectrum_figure.broadening = function_name
        self.spectrum_window.draw()

    def setFwhm(self, fwhm):
        self.spectrum_figure.width = fwhm
        self.spectrum_window.draw()

    def setMarker(self, frequency):
        try:
            frequency = float(frequency)
        except ValueError:
            # frequency == u""
            return
        self.spectrum_figure.mark_line(frequency)
        self.spectrum_window.draw()

    def setScalingFactor(self, scalingFactor):
        self.molecule_figure.scaling_factor = float(scalingFactor)
        self.molecule_window.draw()

    def setThreshold(self, threshold):
        self.molecule_figure.threshold = float(threshold)
        self.molecule_window.draw()

    def showAtomIndex(self, show=True):
        self.molecule_figure.show_atom_index = show
        self.molecule_window.draw()

    def setAllBlackAtomLabels(self, black=True):
        self.molecule_figure.black_and_white = black
        self.molecule_window.draw()

    def setFontSize(self, fontSize):
        self.molecule_figure.fontsize = int(fontSize)
        self.molecule_window.draw()

    def setLineWidth(self, lw):
        self.molecule_figure.linewidth = float(lw)
        self.molecule_window.draw()


def main():
    app = QtGui.QApplication(sys.argv)
    window = MainWindow()
    window.show()

    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
