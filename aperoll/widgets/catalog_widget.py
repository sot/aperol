import io

import numpy as np
from astropy.io import ascii
from astropy.table import Table
from cxotime import CxoTime
from PyQt5 import QtCore as QtC
from PyQt5 import QtGui as QtG
from PyQt5 import QtWidgets as QtW

from aperoll.utils import log_exception, logger
from aperoll.widgets.attitude_widget import (
    AttitudeWidget,
    QuatRepresentation,
    hstack,
    vstack,
)
from aperoll.widgets.line_edit import LineEdit


class CatalogWidget(QtW.QWidget):
    track = QtC.pyqtSignal(bool)

    def __init__(self):
        super().__init__()

        self.date_edit = LineEdit(self)
        self.date_edit.setReadOnly(True)
        self.attitude_widget = AttitudeWidget(
            self,
            columns={
                QuatRepresentation.EQUATORIAL: 0,
                QuatRepresentation.QUATERNION: 0,
                QuatRepresentation.SUN: 0,
            },
        )
        self.text_widget = QtW.QTextEdit()
        self.text_widget.setReadOnly(True)
        font = self.text_widget.document().defaultFont()
        font.setPixelSize(12)
        font.setFamily("Courier New")
        width = QtG.QFontMetrics(font).width("M" * 60)
        self.text_widget.document().setDefaultFont(font)
        self.text_widget.setMinimumSize(width, 200)

        self.track_button = QtW.QPushButton("Start tracking")
        self.track_button.setCheckable(True)
        self.track_button.clicked.connect(self._track)

        date_box = QtW.QGroupBox("Date")
        date_box.setLayout(
            vstack(
                self.date_edit,
            )
        )

        self.starcat_source_edit = QtW.QComboBox(self)
        self.starcat_source_edit.addItems(["Kadi", "Proseco", "Find Attitude"])

        source_box = QtW.QGroupBox("Source")
        source_box.setLayout(
            vstack(
                self.starcat_source_edit,
            )
        )

        att_box = QtW.QGroupBox("Attitude")
        att_box.setLayout(
            vstack(
                self.attitude_widget,
            )
        )

        table_box = QtW.QGroupBox("Catalog")
        table_box.setLayout(vstack(self.text_widget))

        layout = QtW.QHBoxLayout()
        layout.addStretch()
        layout.addLayout(
            vstack(
                hstack(date_box, source_box),
                att_box,
                table_box,
                hstack(self.track_button, stretch=True),
                stretch=True,
            )
        )
        layout.addStretch()
        self.setLayout(layout)
        self.reset()

    def _track(self, checked):
        text = "Start" if checked else "Stop"
        self.track_button.setText(f"{text} tracking")
        self.track.emit(checked)

    def set_attitude(self, attitude):
        self.attitude_widget.set_attitude(attitude)

    def get_attitude(self):
        self.attitude_widget.get_attitude()

    attitude = property(get_attitude, set_attitude)

    def set_date(self, date):
        date = CxoTime(date).date if date is not None else ""
        if date != self.date_edit.text():
            self.date_edit.setText(date)

    def get_date(self):
        text = self.date_edit.text()
        if text:
            return CxoTime(self.date_edit.text())

    date = property(get_date, set_date)

    def get_starcat_source(self):
        return self.starcat_source_edit.currentText()

    def set_starcat_source(self, source):
        self.starcat_source_edit.setCurrentText(source)

    starcat_source = property(get_starcat_source, set_starcat_source)

    def reset(self):
        self.set_date(None)
        self.attitude = None
        self.table = Table()
        self.table["slot"] = np.arange(8)
        self.table["id"] = np.ma.masked_all(8, dtype=int)
        self.table["yang"] = np.ma.masked_all(8, dtype=float)
        self.table["zang"] = np.ma.masked_all(8, dtype=float)
        self.table["mag"] = np.ma.masked_all(8, dtype=float)
        self.table["enabled"] = np.zeros(8, dtype=bool)
        self.table["yang"].format = "8.2f"
        self.table["zang"].format = "8.2f"
        self.table["mag"].format = "5.2f"
        self.display_text()

    def set_catalog(self, catalog):
        try:
            logger.debug("Setting catalog")
            self.reset()
            if np.any(self.table["slot"] != np.arange(8)):
                # this is a sanity check, it really should never happen
                raise Exception("Slots in table are not sorted")

            self.date = catalog.date
            self.attitude = catalog.att

            guides = catalog[np.in1d(catalog["type"], ["BOT", "GUI"])]
            idx = np.searchsorted(np.arange(8), guides["slot"])
            self.table["slot"][idx] = guides["slot"]
            self.table["id"][idx] = guides["id"]
            self.table["yang"][idx] = guides["yang"]
            self.table["zang"][idx] = guides["zang"]
            self.table["mag"][idx] = guides["mag"]
            self.table["enabled"][idx] = True
            self.table["yang"].format = "8.2f"
            self.table["zang"].format = "8.2f"
            self.table["mag"].format = "5.2f"
            self.table.pprint()
            self.display_text()
        except Exception as exc:
            log_exception("Error setting catalog", exc)

    def display_text(self):
        buf = io.StringIO()
        ascii.write(
            self.table,
            format="fixed_width_two_line",
            header_rows=["name", "dtype"],
            output=buf,
        )
        self.text_widget.setText(buf.getvalue())

    def read_text(self):
        buf = io.StringIO(self.text_widget.text())
        self.table = ascii.read(
            buf,
            format="fixed_width_two_line",
            header_rows=["name", "dtype"],
            output=buf,
        )

    def enable_slot(self, slot, enable):
        self.table["enabled"][slot] = enable
        self.display_text()
