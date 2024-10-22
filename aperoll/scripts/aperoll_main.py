#!/usr/bin/env python

# from PyQt5 import QtCore as QtC, QtWidgets as QtW, QtGui as QtG
from PyQt5 import QtWidgets as QtW

from aperoll.utils import AperollException, logger
from aperoll.widgets.main_window import MainWindow


def get_parser():
    import argparse

    parse = argparse.ArgumentParser()
    parse.add_argument("file", nargs="?", default=None)
    parse.add_argument("--obsid", help="Specify the OBSID", type=int)
    levels = ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]
    levels = [lvl.lower() for lvl in levels]
    parse.add_argument(
        "--log-level", help="Set the log level", default="INFO", choices=levels
    )
    return parse


def main():
    parser = get_parser()
    args = parser.parse_args()

    logger.setLevel(args.log_level)

    try:
        app = QtW.QApplication([])
        w = MainWindow(opts=vars(args))
        w.resize(1500, 1000)
        w.show()
        app.exec()
    except AperollException as e:
        logger.error(f"Error: {e}")
        parser.exit(1)
    except FileNotFoundError as e:
        logger.error(f"Error: {e}")
        parser.exit(1)


if __name__ == "__main__":
    main()
