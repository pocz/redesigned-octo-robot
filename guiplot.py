from PyQt5 import  QtCore, QtGui, QtWidgets
from pyqtgraph import PlotWidget, plot
import pyqtgraph as pg
import sys
import os

import BoyerMoore
import FastQ

class RectItem(pg.GraphicsObject):
    def __init__(self, rect, parent=None):
        super().__init__(parent)
        self._rect = rect
        self.picture = QtGui.QPicture()
        self._generate_picture()

    @property
    def rect(self):
        return self._rect

    def _generate_picture(self):
        painter = QtGui.QPainter(self.picture)
        painter.setPen(pg.mkPen("w"))
        painter.setBrush(pg.mkBrush("g"))
        painter.drawRect(self.rect)
        painter.end()

    def paint(self, painter, option, widget=None):
        painter.drawPicture(0, 0, self.picture)

    def boundingRect(self):
        return QtCore.QRectF(self.picture.boundingRect())

class MainWindow(QtWidgets.QMainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)

        self.graphWidget = pg.PlotWidget()
        self.setCentralWidget(self.graphWidget)

        kmer = "AGGATTAC"
        seqpath = "data/timings/50m.fq"
        reads = FastQ.Sample(seqpath).reads
        first_key = list(reads.keys())[0]
        seq = FastQ.Sample(seqpath).reads[first_key].seq
        start_indexes = BoyerMoore.start_indexes(seq, kmer)
        nhits = len(start_indexes)

        print(f'{nhits} hits collected')

        ypos = 0
        width = len(kmer)
        height = 1
        for i in range(nhits):
            xpos = start_indexes[i]
            self.graphWidget.addItem(RectItem(QtCore.QRectF(xpos,ypos,width,height)))

def main():
    app = QtWidgets.QApplication(sys.argv)
    main = MainWindow()
    main.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
