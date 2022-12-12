#! /usr/bin/env python

from PyQt5.QtWidgets import QMainWindow

from PyQt5 import  QtCore, QtGui, QtWidgets
from pyqtgraph import PlotWidget, plot

import argparse
import pyqtgraph as pg
import sys
import os

import AhoCorasick
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

        def draw_graph(args):    
            def last_read(filepath):
                seq = ''
                for read in FastQ.Sample(filepath).reads:
                    seq = FastQ.Sample(filepath).reads[read].seq
                return seq

            kmers = []
            for file in args.kmer_files:
                kmer = last_read(file)
                kmers.append(kmer)

            seq = last_read(args.seq_file)
            
            if args.mode == 'ProcessPool':
                start_indexes = BoyerMoore.start_indexes_parallel(seq, kmers)
                start_indexes = [s for sublist in start_indexes for s in sublist]
            elif args.mode == "Threading":
                start_indexes = BoyerMoore.start_indexes_threading(seq, kmers)
                start_indexes = [s for sublist in start_indexes for s in sublist]
            elif args.mode == "AhoCorasick":
                automaton = AhoCorasick.Automaton(kmers,"DNA")
                start_index_dict = automaton.start_indexes(seq)
                start_indexes = []
                for indexes in start_index_dict.values():
                    start_indexes.append(indexes) 
                start_indexes = [s for sublist in start_indexes for s in sublist]
            print(f'{start_indexes}')

            nhits = len(start_indexes)
            print(f'{nhits} hit(s) collected')

            ypos = 0
            width = len(kmers[0])
            height = 1
            for i in range(nhits):
                xpos = start_indexes[i]
                self.graphWidget.addItem(RectItem(QtCore.QRectF(xpos,ypos,width,height)))

        def run():
            parser = argparse.ArgumentParser(description="Displays the position of each kmer hit in the sequence")
            parser.add_argument('--kmers',dest='kmer_files',  metavar='kmer', type=str, nargs='+', help="kmers to look for")
            parser.add_argument('--seq', dest='seq_file', type=str, required=True)
            parser.add_argument('--mode', dest='mode', metavar='ProcessPool|Threading|AhoCorasick', type=str, help="method used for string search", required=True)
            parser.set_defaults(func=draw_graph)
            args=parser.parse_args()
            args.func(args)

        super(MainWindow, self).__init__(*args, **kwargs)

        self.graphWidget = pg.PlotWidget()
        self.setCentralWidget(self.graphWidget)

        run()

if __name__ == '__main__':
#    appctxt = ApplicationContext()       # 1. Instantiate ApplicationContext
    appctxt = QtWidgets.QApplication(sys.argv)
    window = MainWindow()
    window.resize(250, 150)
    window.show()
    exit_code = appctxt.exec_()      # 2. Invoke appctxt.app.exec_()
    sys.exit(exit_code)

