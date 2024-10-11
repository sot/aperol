#!/usr/bin/env python

# from PyQt5 import QtCore as QtC, QtWidgets as QtW, QtGui as QtG
from PyQt5 import QtWidgets as QtW

from aperoll.widgets.main_window import MainWindow


def get_parser():
    import argparse

    parse = argparse.ArgumentParser()
    parse.add_argument("--date", default="2021-08-23 00:00:00.000")
    parse.add_argument("--ra", default="6.45483333")
    parse.add_argument("--dec", default="-26.03683611")
    return parse


def main():
    args = get_parser().parse_args()

    app = QtW.QApplication([])
    w = MainWindow(opts=vars(args))
    w.resize(1500, 1000)
    w.show()
    app.exec()


if __name__ == "__main__":
    main()
