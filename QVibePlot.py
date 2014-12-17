#!/usr/bin/env python
# Copyright (c) 2011-2014 Mathias Laurin, 3-clause BSD License

"""Qt graphical user interface (GUI) to vibeplot."""

from __future__ import print_function
import sys
import os.path
import platform
import csv
import six
from glob import glob
from functools import partial
from six.moves import zip_longest
from collections import defaultdict

# Import Qt and matplotlib modules
import matplotlib
try:
    from PyQt5.QtWidgets import *
    from PyQt5.QtGui import QKeySequence
    from PyQt5.QtCore import *
    from PyQt5 import uic
    Signal, Slot = pyqtSignal, pyqtSlot
    import rcc5
    matplotlib.use("Qt5Agg")
    from matplotlib.backends.backend_qt5agg import (NavigationToolbar2QT
                                                    as NavigationToolbar)
    _getOpenFileName = QFileDialog.getOpenFileName
except ImportError:
    import sip
    for qtype in "QString QTextStream QVariant".split():
        sip.setapi(qtype, 2)

    from PyQt4.QtGui import *
    from PyQt4.QtCore import *
    from PyQt4 import uic
    Signal, Slot = pyqtSignal, pyqtSlot
    import rcc4
    matplotlib.use("Qt4Agg")
    from matplotlib.backends.backend_qt4agg import (NavigationToolbar2QT
                                                    as NavigationToolbar)
    _getOpenFileName = QFileDialog.getOpenFileNameAndFilter
from matplotlib.backends.qt_compat import _getSaveFileName

if platform.system() == "Windows":
    from qvibeplot_ui import Ui_MainWindow
    class MainWindow(QMainWindow, Ui_MainWindow):

        def __init__(self):
            super(MainWindow, self).__init__()
            self.setupUi(self)

else:
    class MainWindow(QMainWindow):

        def __init__(self):
            super(MainWindow, self).__init__()
            uifile = QFile(":/qvibeplot.ui")
            uifile.open(QFile.ReadOnly)
            uic.loadUi(uifile, self)
            uifile.close()

import pybel
from vibeplot import __version__
import vibeplot.plotter as plotter


FORMATS = [
    "Generic Output file format",
    "ACES output format",
    "Aoforce output format",
    "GAMESS Output",
    "GAMESS-UK Output",
    "Gaussian Output",
    "Molden format",
    "Molpro output format",
    "MOPAC Output format",
    "NWChem output format",
    "Orca output format",
    "Q-Chem output format",
    "VASP format",
]


class FormatTextDelegate(QStyledItemDelegate):

    """QStyledItemDelegate formatting the text displayed."""

    def __init__(self, format="%0.2f", parent=None):
        super(FormatTextDelegate, self).__init__(parent)
        self._format = format

    def displayText(self, value, locale):
        return self._format % value


class ItemSelectionModel(QItemSelectionModel):

    def __init__(self, model, parent=None):
        super(ItemSelectionModel, self).__init__(model, parent)

    @Slot()
    def copy(self):
        if not self.hasSelection(): return
        previous = QModelIndex()
        fields = []
        for index in sorted(self.selectedIndexes()):
            if index.row() is previous.row():
                fields[-1].append(index.data())
            else:
                fields.append([index.data()])
            previous = index
        csvfile = six.StringIO()
        writer = csv.writer(csvfile)
        writer.writerows(fields)
        QApplication.clipboard().setText(csvfile.getvalue())


class QVibeplot(MainWindow):

    """vibeplot graphical user interface."""

    def __init__(self):
        super(QVibeplot, self).__init__()
        self.spectrumTable.setHorizontalHeaderLabels(
            [u"\u03bd\u0305 / cm\u207b\u00b9", u"Intensity / %"])
        self.spectrumTable.setSelectionModel(
            ItemSelectionModel(self.spectrumTable.model(), self.spectrumTable))
        self.spectrumTable.setItemDelegate(
            FormatTextDelegate("%0.0f", self.spectrumTable))

        self.toolbar = NavigationToolbar(self.moleculeCanvas, self, False)
        self.addToolBar(Qt.RightToolBarArea, self.toolbar)

        self.moleculeCanvas.setContextMenuPolicy(Qt.ActionsContextMenu)
        self.spectrumTable.setContextMenuPolicy(Qt.ActionsContextMenu)
        self.spectrumCanvas.setContextMenuPolicy(Qt.ActionsContextMenu)

        self.moleculePlotter = plotter.MoleculePlotter(
            self.moleculeCanvas.figure.add_subplot(111))
        self.spectrumPlotter = plotter.SpectrumPlotter(
            self.spectrumCanvas.figure.add_subplot(111))

        self._settings = QSettings("Mathias Laurin", "QVibePlot")
        for setting in "imagePath dataPath".split():
            if not self._settings.contains(setting):
                self._settings.setValue(setting, QDir.homePath())

        # Create actions
        self.saveSpectrumDataAction = QAction(
            u"Save Spectral Broadening Data...", self.spectrumCanvas,
            triggered=self.saveSpectrum, enabled=False)
        self.copySpectrumDataAction = QAction(
            u"Copy", self.spectrumTable,
            shortcut=QKeySequence.Copy,
            triggered=self.spectrumTable.selectionModel().copy)
        self.openFileAction = QAction(
            u"Open", self.fileMenu,
            shortcut=QKeySequence.Open,
            triggered=self.loadFile)
        self.saveImageAction = QAction(
            u"Save Image", self.fileMenu,
            shortcut=QKeySequence.Save,
            triggered=self.toolbar.save_figure)
        self.quitAction = QAction(
            u"Quit", self.fileMenu,
            shortcut=QKeySequence.Quit,
            triggered=self.close)
        self.showAtomIndexAction = QAction(
            u"Atom Index", self.viewMenu, checkable=True,
            triggered=self.moleculePlotter.show_atom_index)

        # Connect widgets
        self.fontSizeComboBox.currentIndexChanged[str].connect(
            self.moleculePlotter.set_fontsize)
        self.lineWidthComboBox.currentIndexChanged[str].connect(
            self.moleculePlotter.set_linewidth)
        self.colorLabelCheckBox.stateChanged.connect(
            lambda checked:
            self.moleculePlotter.set_black_labels(checked==Qt.Checked))
        self.bondlengthFilter.valueChanged.connect(
            lambda value: self._drawVibration())
        self.angleFilter.valueChanged.connect(
            lambda value: self._drawVibration())
        self.torsionFilter.valueChanged.connect(
            lambda value: self._drawVibration())
        self.broadeningComboBox.currentIndexChanged[str].connect(
            self.spectrumPlotter.set_broadening_function)
        self.broadeningComboBox.currentIndexChanged.connect(
            self.saveSpectrumDataAction.setEnabled)
        self.fwhmDoubleSpinBox.valueChanged.connect(
            self.spectrumPlotter.set_fwhm)
        self.spectrumTable.selectionModel().currentRowChanged.connect(
            lambda current, previous:
            self.spectrumPlotter.set_vibration(
                self.spectrumTable.item(current.row(), 0)
                .data(Qt.DisplayRole)))
        self.spectrumTable.selectionModel().currentRowChanged.connect(
            lambda current, previous: self._drawVibration())
        self.broadeningComboBox.setCurrentIndex(1)
        # Add actions
        self.spectrumTable.addAction(self.copySpectrumDataAction)
        self.spectrumCanvas.addAction(self.saveSpectrumDataAction)
        self.moleculeCanvas.addActions(
            (self.saveImageAction, self.showAtomIndexAction))
        # File menu
        self.fileMenu.addActions((self.openFileAction, self.saveImageAction,
                                  self.quitAction))
        # File > OrbiMol menu
        self.orbiMolDbMenu = QMenu(u"OrbiMol DB Molecules")
        self.orbiMolDbMenu.addActions([
            QAction(os.path.splitext(os.path.basename(filename))[0],  # text
                    self.orbiMolDbMenu,
                    triggered=partial(self.loadFile, filename, "g03"))
            for filename in glob("data/orbimol/*.freq")])
        self.orbiMolDbMenu.addSeparator()
        self.orbiMolDbMenu.addAction(
            QAction(u"About OrbiMol", self.orbiMolDbMenu,
                    triggered=partial(QMessageBox.about, self, u"About OrbiMol",
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
        self.viewMenu.addActions((self.showAtomIndexAction,))
        # Help menu
        self.helpMenu.addActions((
            QAction(u"About", self.helpMenu, triggered=partial(
                QMessageBox.about,
                self, "About QVibePlot", " ".join((
            u"""
            QVibePlot visualizes vibrational analysis performed by
            density functional theory calculations (DFT) in terms of
            changes of internal coordinates.
            <p>QVibePlot is written in Python and depends on matplotlib
            for the graphics and numpy for the maths. The GUI is
            written using PyQt4.</p>
            <p>Copyright (c) 2011-2014
            <a href="mailto:Mathias.Laurin+vibeplot@gmail.com"> Mathias
            Laurin</a></p>
            <p>QVibePlot %s is available under the modified BSD
            License.</p>
            <p>Support the program by citing: Laurin, M.  QVibeplot: A
            Program To Visualize Molecular Vibrations in Two
            Dimensions. <i>J. Chem. Educ.</i> <b>2013</b>
            <a href="http://dx.doi.org/10.1021/ed300554z">DOI:
            10.1021/ed300554z</a>.</p>

            """ % __version__).splitlines()))),
            QAction(u"About Qt", self.helpMenu,
                    triggered=partial(QMessageBox.aboutQt, self)),
            QAction(u"About Open Babel", self.helpMenu,
                    triggered=partial(QMessageBox.about, self,
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
            QAction(u"About Matplotlib", self.helpMenu,
                    triggered=partial(QMessageBox.about, self,
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

            """.format(matplotlib.__version__)).splitlines())))
        ))

    def _drawVibration(self):
        """Pass parameters from GUI to `moleculePlotter.draw_vibration`."""
        self.moleculePlotter.draw_vibration(
            self.spectrumTable.currentRow(),
            self.bondlengthFilter.value(),
            self.angleFilter.value(),
            self.torsionFilter.value())

    def _getFilename(self):
        items, values = six.iteritems, six.itervalues
        informats = defaultdict(str)
        for (fmt, ext) in ((fmt, ext) for (ext, fmt)
                           in items(pybel.informats)
                           if fmt in FORMATS):
            informats[fmt] += (" *.%s" if ext not in "CONTCAR POSCAR".split()
                               else " %s") % ext
        informats = {fmt: ext[1:] for (fmt, ext) in items(informats)}
        formats, extensions = items(informats), values(informats)
        filename, __ = _getOpenFileName(
            self,
            u"Open File",
            self._settings.value("dataPath"),
            ";;".join(
                ["Common formats (%s)" % " ".join(extensions)] +
                ["%s (%s)" % __ for __ in sorted(formats)] +
                ["all files (*)"]
            ))
        return filename

    def setWindowTitle(self, text=""):
        """Write QVibeplot in the window title."""
        super(QVibeplot, self).setWindowTitle(
            'QVibeplot' if not text else
            '%s - QVibeplot' % os.path.basename(text))

    def loadFile(self, filename=None, inFormat=None):
        """Use Open Babel to load a molecule from a file."""
        filename = filename if filename else self._getFilename()
        if filename:
            self._settings.setValue("dataPath", os.path.dirname(filename))
        else:
            return

        if inFormat is None:
            inFormat = str(os.path.splitext(filename)[1][1:])
        if not inFormat and os.path.basename(filename).lower() in (
                "poscar", "contcar"):
            inFormat = "vasp"

        # load data
        mol = next(pybel.readfile(inFormat, str(filename)))
        if not mol.atoms:
            self.statusBar().showMessage(
                "".join((
                    "Extension or file format '%s' unknown, ",
                    "see http://openbabel.org for the list of ",
                    "supported files.")) % inFormat)
        self.moleculePlotter.clear()
        self.spectrumPlotter.clear()
        self.moleculePlotter.set_molecule(mol)
        self.spectrumPlotter.set_molecule(mol)
        self.moleculePlotter.draw_molecule(
            lw=str(self.lineWidthComboBox.currentText()),
            fontsize=str(self.fontSizeComboBox.currentText())
        )
        self.spectrumPlotter.draw_spectrum()

        # populate spectrumTable
        self.spectrumTable.setRowCount(len(self.spectrumPlotter.frequencies))
        intensities = self.spectrumPlotter.intensities
        intensities *= 100.0 / intensities.max()
        for row, (frequency, intensity) in enumerate(
                zip_longest(self.spectrumPlotter.frequencies,
                            intensities,
                            fillvalue=100.0)):
            for column, data in enumerate((frequency, intensity)):
                item = QTableWidgetItem()
                item.setData(Qt.DisplayRole, data.item())
                self.spectrumTable.setItem(row, column, item)

        # window title
        mol.title = mol.write(
            "can", opt=dict(n=None,  # no molecule name
                            i=None,  # no chiral markings
                            ))
        self.setWindowTitle(mol.title)

        # SVG representation
        mol.removeh()      # implicit h
        self.svgWidget.load(QByteArray(
            mol.write("svg", opt=dict(
                d=None,    # no molecule name
                u=None,    # black and white
                C=None,    # implicit C
                b="none",  # transparent background
                # next line is an undocumented trick to write:
                # width="100%" height="100%" in the svg header instead of
                # the default fixed size.  This is required for the image
                # to resize properly in the widget.
                svgwritechemobject=None))))

    def saveSpectrum(self):
        """Save the broadened spectrum to file."""
        filename, __ = _getSaveFileName(
            self,
            u"Save Spectrum Values",
            self._settings.value("imagePath"),
            "plain text (*.txt)")
        if not filename:
            return
        if "." not in filename:
            filename += ".txt"
        self.spectrumPlotter.save_spectrum(filename)


def main():
    app = QApplication(sys.argv)
    app.lastWindowClosed.connect(app.quit)
    window = QVibeplot()
    window.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
