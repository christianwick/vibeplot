from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from PyQt4.QtGui import QPalette


class MplCanvas(FigureCanvas):

    def __init__(self, parent=None):
        fig = Figure()
        super(MplCanvas, self).__init__(fig)
        self.setParent(parent)
        self._oldSize = None
        self._background = None

    def setParent(self, parent):
        super(MplCanvas, self).setParent(parent)
        if parent:
            color = parent.palette().brush(QPalette.Window).color()
            self.figure.set_facecolor("#%X%X%X" % (color.red(),
                                                   color.green(),
                                                   color.blue()))

    def restore_background(self):
        self.figure.restore_region(self._background)

    def draw(self):
        super(MplCanvas, self).draw()
        self._oldSize = self.figure.bbox.width, self.figure.bbox.height
        self._background = self.copy_from_bbox(self.figure.bbox)

    def _handleResize(self):
        if self._oldSize != (self.figure.bbox.width, self.figure.bbox.height):
            self.draw()

