from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from PyQt4.QtGui import QPalette


class MplCanvas(FigureCanvas):

    def __init__(self, parent=None):
        fig = Figure()
        super(MplCanvas, self).__init__(fig)
        self.setParent(parent)

    def setParent(self, parent):
        super(MplCanvas, self).setParent(parent)
        if parent:
            color = parent.palette().brush(QPalette.Window).color()
            self.figure.set_facecolor("#%X%X%X" % (color.red(),
                                                   color.green(),
                                                   color.blue()))

    def draw(self):
        super(MplCanvas, self).draw()

