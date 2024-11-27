import json
from enum import Enum

import numpy as np
from cxotime import CxoTime
from PyQt5 import QtCore as QtC
from PyQt5 import QtGui as QtG
from PyQt5 import QtWidgets as QtW
from Quaternion import Quat, normalize
from ska_sun import get_sun_pitch_yaw, off_nominal_roll

from aperoll.widgets.error_message import ErrorMessage


class QuatRepresentation(Enum):
    QUATERNION = "Quaternion"
    EQUATORIAL = "Equatorial"
    SUN = "Sun Position"


def stack(layout, *items, spacing=None, content_margins=None, stretch=False):
    if spacing is not None:
        layout.setSpacing(spacing)
    if content_margins is not None:
        layout.setContentsMargins(*content_margins)
    for item in items:
        if isinstance(item, QtW.QWidget):
            layout.addWidget(item)
        elif isinstance(item, QtW.QLayout):
            layout.addLayout(item)
        elif isinstance(item, QtW.QSpacerItem):
            layout.addItem(item)
        else:
            print(f"unknown type {type(item)}")
    if stretch:
        layout.addStretch()
    return layout


def hstack(*items, **kwargs):
    return stack(QtW.QHBoxLayout(), *items, **kwargs)


def vstack(*items, **kwargs):
    return stack(QtW.QVBoxLayout(), *items, **kwargs)


class TextEdit(QtW.QTextEdit):
    values_changed = QtC.pyqtSignal(list)

    def __init__(self, size=4, digits=12, width=None, parent=None):
        super().__init__(parent=parent)
        self.installEventFilter(self)
        self.setSizePolicy(
            QtW.QSizePolicy(
                QtW.QSizePolicy.MinimumExpanding,
                # QtW.QSizePolicy.Fixed
                QtW.QSizePolicy.Ignored,
            )
        )
        width = width or digits
        font = self.font()
        fm = QtG.QFontMetrics(font)
        font_size = fm.width("M")
        # font_height = fm.width("M")
        self.setMinimumWidth(width * font_size)
        self.setMinimumHeight(
            size * fm.lineSpacing() + fm.lineSpacing() // 2 + 2
            # size * fm.lineSpacing() + fm.lineSpacing() // 2 + 2 * self.frameWidth()
        )

        self.fmt = f"{{:.{digits}f}}"
        self.length = size
        self._vals = None

        self.reset()

    def sizeHint(self):
        return QtC.QSize(125, 20)

    def get_values(self):
        return self._vals

    def set_values(self, values):
        if values is None:
            self.reset()
            return
        if not hasattr(values, "__iter__"):
            raise ValueError("values must be an iterable")
        values = np.array(values)
        if len(values) != self.length:
            raise ValueError(f"expected {self.length} values, got {len(values)}")
        if np.all(values == self._vals):
            return
        self._vals = values
        self._display_values()

    values = property(get_values, set_values)

    def _parse_values(self, text):
        """
        Parse a string to get the values.

        The string usually comes from the text box, or from the clipboard.
        """
        # we expect a string of floats separated by commas or whitespace with length == self.length
        unknown = set(text) - set("-e0123456789., \n\t")
        if unknown:
            raise ValueError(f"invalid characters: {unknown}")
        vals = [float(s.strip()) for s in text.replace(",", " ").split()]
        if len(vals) != self.length:
            raise ValueError(f"expected {self.length} values, got {len(vals)}")
        return vals

    def _update_values(self):
        # take the text, parse it, and set the values
        try:
            vals = self._parse_values(self.toPlainText())
            self._vals = vals
            pos = self.textCursor().position()
            self._display_values()
            cursor = self.textCursor()
            cursor.setPosition(pos)
            self.setTextCursor(cursor)
            self.values_changed.emit(self._vals)
        except ValueError as exc:
            error_dialog = ErrorMessage(title="Value Error", message=str(exc))
            error_dialog.exec()

    def _display_values(self):
        """
        Display the values in the text box.
        """
        text = "\n".join(self.fmt.format(v) for v in self._vals)
        self.setPlainText(text)

    def reset(self):
        """
        Clear the contents of the text box and set the values to None.
        """
        self._vals = None
        self.setPlainText("\n".join("" for _ in range(self.length)))

    def focusOutEvent(self, event):
        super().focusOutEvent(event)
        # originally I had the following, but this causes a horrible error on exit which I still
        # need to investigate
        # self._update_values()

    def keyPressEvent(self, event):
        """
        Listen for Key_Return to save and escape to discard changes.
        """
        if event.key() == QtC.Qt.Key_Return:
            self._update_values()
        elif event.key() == QtC.Qt.Key_Escape:
            # discard any changes to the text box
            self._display_values()
        elif event.matches(QtG.QKeySequence.Copy):
            # copy the selected text (if it is selected) or all values to the clipboard
            # when copying all the values, they are converted to a json string
            cursor = self.textCursor()
            if cursor.hasSelection():
                text = cursor.selectedText()
                QtW.QApplication.clipboard().setText(text)
            else:
                vals = self._parse_values(self.toPlainText())
                text = json.dumps(vals)
                QtW.QApplication.clipboard().setText(text)
        else:
            return super().keyPressEvent(event)

    def insertFromMimeData(self, data):
        """
        Insert data from the clipboard.
        """
        try:
            # if this succeeds, presumably we are pasting the whole thing, so values are set
            vals = json.loads(data.text())
            self.set_values(vals)
        except ValueError:
            # if it fails, paste it and the user can edit it
            self.insertPlainText(data.text())


class UnpaddedLabel(QtW.QLabel):
    def __init__(self, text, parent=None):
        super().__init__(text, parent)
        self.setStyleSheet("QLabel { padding: 0px; }")


class AttitudeWidget(QtW.QWidget):
    attitude_changed = QtC.pyqtSignal(Quat)
    attitude_cleared = QtC.pyqtSignal()

    def __init__(self, parent=None, columns=None):
        super(AttitudeWidget, self).__init__(parent)

        if columns is None:
            columns = {
                QuatRepresentation.QUATERNION: 0,
                QuatRepresentation.EQUATORIAL: 1,
                QuatRepresentation.SUN: 2,
            }

        self._q = TextEdit()
        self._eq = TextEdit(size=3, digits=5, width=8)
        self._sun_pos = TextEdit(size=3, digits=5, width=8)

        self._q.values_changed.connect(self._set_attitude)
        self._eq.values_changed.connect(self._set_attitude)

        self._sun_pos.setReadOnly(True)

        self._set_layout(columns)

        self._attitude = None
        self._date = None

    def _set_layout(self, columns):
        layout = QtW.QHBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(layout)

        margins = (20, 0, 20, 0)
        layout_q = vstack(
            hstack(
                vstack(
                    QtW.QSpacerItem(0, 5),
                    UnpaddedLabel("Q1"),
                    UnpaddedLabel("Q2"),
                    UnpaddedLabel("Q3"),
                    UnpaddedLabel("Q4"),
                    spacing=0,
                    content_margins=margins,
                    stretch=True,
                ),
                vstack(
                    self._q,
                    spacing=0,
                    content_margins=margins,
                    # stretch=True,
                ),
            ),
            content_margins=(0, 0, 0, 0),
            stretch=True,
        )

        layout_eq = vstack(
            hstack(
                vstack(
                    QtW.QSpacerItem(0, 3),
                    UnpaddedLabel("ra  "),
                    UnpaddedLabel("dec "),
                    UnpaddedLabel("roll"),
                    spacing=0,
                    content_margins=margins,
                    stretch=True,
                ),
                vstack(
                    self._eq,
                    spacing=0,
                    content_margins=margins,
                    # stretch=True,
                ),
            ),
            content_margins=(0, 0, 0, 0),
            stretch=True,
        )

        layout_sun = vstack(
            hstack(
                vstack(
                    QtW.QSpacerItem(0, 3),
                    UnpaddedLabel("pitch"),
                    UnpaddedLabel("yaw"),
                    UnpaddedLabel("roll"),
                    spacing=0,
                    content_margins=margins,
                    stretch=True,
                ),
                vstack(
                    self._sun_pos,
                    spacing=0,
                    content_margins=margins,
                    # stretch=True,
                ),
            ),
            content_margins=(0, 0, 0, 0),
            stretch=True,
        )

        layouts = {
            QuatRepresentation.QUATERNION: layout_q,
            QuatRepresentation.EQUATORIAL: layout_eq,
            QuatRepresentation.SUN: layout_sun,
        }
        name = {
            QuatRepresentation.QUATERNION: "Quaternion",
            QuatRepresentation.EQUATORIAL: "Equatorial",
            QuatRepresentation.SUN: "Sun",
        }

        self.tab_widgets = {
            col: QtW.QTabWidget() for col in set(columns.values()) if col is not None
        }

        for representation, col in columns.items():
            if col is None:
                continue
            w = QtW.QWidget()
            w.setLayout(layouts[representation])
            self.tab_widgets[col].addTab(w, name[representation])

        for widget in self.tab_widgets.values():
            layout.addWidget(widget)
            widget.setCurrentIndex(0)

        self.update()

    def get_attitude(self):
        return self._attitude

    def set_attitude(self, attitude):
        self._set_attitude(attitude, emit=False)

    def _set_attitude(self, attitude, emit=True):
        # work around the requirement that q be normalized
        if attitude is None:
            self._clear()
            if emit:
                self.attitude_cleared.emit()
            return
        if (
            not isinstance(attitude, Quat)
            and len(attitude) == 4
        ):
            attitude = normalize(attitude)
        # this check is to break infinite recursion because in the connections
        q1 = None if attitude is None else Quat(attitude).q
        q2 = None if self._attitude is None else self._attitude.q
        if np.any(q1 != q2):
            self._attitude = Quat(attitude)
            self._display_attitude_at_date(self._attitude, self._date)
            if emit:
                self.attitude_changed.emit(self._attitude)

    attitude = property(get_attitude, set_attitude)

    def set_date(self, date):
        date = None if date is None else CxoTime(date)
        if self._date == date:
            return
        self._date = date
        self._display_attitude_at_date(self._attitude, self._date)

    def get_date(self):
        return self._date

    date = property(get_date, set_date)

    def _display_attitude_at_date(self, attitude, date):
        if attitude is None:
            self._clear()
            return
        self._q.set_values(attitude.q)
        self._eq.set_values(attitude.equatorial)

        if date is None:
            self._sun_pos.reset()
        else:
            pitch, yaw = get_sun_pitch_yaw(attitude.ra, attitude.dec, date)
            roll = off_nominal_roll(attitude, date)
            self._sun_pos.set_values([pitch, yaw, roll])

    def _clear(self):
        self._q.reset()
        self._eq.reset()
        self._sun_pos.reset()


if __name__ == "__main__":
    app = QtW.QApplication([])
    widget = AttitudeWidget()
    q = Quat([344.571937, 1.026897, 302.0])
    widget.set_attitude(q)
    widget.set_date("2021:001:00:00:00")
    widget.resize(1200, 200)
    widget.show()
    app.exec()
