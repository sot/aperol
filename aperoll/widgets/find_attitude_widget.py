import io
from pprint import pformat

import numpy as np
from astropy.io import ascii
from astropy.table import Table
from cxotime import CxoTime
from find_attitude import Constraints, find_attitude_solutions
from PyQt5 import QtCore as QtC
from PyQt5 import QtGui as QtG
from PyQt5 import QtWidgets as QtW
from ska_sun import get_sun_pitch_yaw

from aperoll.utils import AperollException, log_exception, logger
from aperoll.widgets.attitude_widget import (
    AttitudeWidget,
    QuatRepresentation,
    hstack,
    vstack,
)
from aperoll.widgets.error_message import ErrorMessage
from aperoll.widgets.line_edit import LineEdit


class NoSolutionError(AperollException):
    pass


class Catalog(Table):
    def __init__(self, date=None, att=None, **kwargs):
        super().__init__(**kwargs)
        self.date = date
        self.att = att


class OptionalParameter(QtW.QWidget):
    def __init__(
            self,
            label,
            widget=None,
            parent=None,
            enabled=False,
            layout="v",
            default="",
            formatter=lambda val: f"{val}",
            value_type=float
        ):
        super().__init__(parent)
        self.formatter = formatter
        self.value_type = value_type
        self.default = default
        self.label = QtW.QLabel(label)
        self.widget = LineEdit(self) if widget is None else widget
        self.widget.setText(f"{default}")
        self.widget.setEnabled(enabled)
        self.checkbox = QtW.QCheckBox()
        self.checkbox.setChecked(enabled)
        stack = hstack if layout == "h" else vstack
        hs = hstack(
            self.checkbox,
            self.label,
            stretch=True,
        )
        hs.setSpacing(10)
        hs.setContentsMargins(0, 0, 0, 0)
        lyt = stack(
            hs,
            self.widget,
        )
        s = 0
        lyt.setSpacing(0)
        lyt.setContentsMargins(s, 0, s, 0)

        self.setLayout(lyt)
        self.checkbox.stateChanged.connect(self.setEnabled)

    def setEnabled(self, enabled):
        self.widget.setEnabled(enabled)
        self.checkbox.setChecked(enabled)

    def isEnabled(self):
        return self.checkbox.isChecked()

    enabled = property(isEnabled, setEnabled)

    def setText(self, text):
        self.widget.setText(text)

    def setValue(self, value):
        self.widget.setText(self.formatter(value))

    def reset(self):
        self.setValue(self.default)

    def value(self):
        return self.value_type(self.widget.text()) if self.enabled else None

class FindAttitudeWidget(QtW.QWidget):

    found_attitude = QtC.pyqtSignal(Catalog)

    def __init__(self):
        super().__init__()

        self.date_edit = LineEdit(self)
        self.date_edit.setReadOnly(True)
        self.tolerance_edit = LineEdit(self)
        self.tolerance_edit.setText("3")

        self.rough_attitude_estimate = OptionalParameter(
            "Rough Attitude Estimate (q1, q2, q3, q4)",
            value_type=lambda val: np.array([float(x) for x in val.split(",")]),
        )
        self.sun_pitch = OptionalParameter("Sun Pitch (degrees)")
        self.sun_pitch_sigma = OptionalParameter(
            "Sun Pitch Uncertainty (degrees)", default=1.5
        )
        self.radial_uncertainty = OptionalParameter(
            "Radial Uncertainty (degrees)", default=4.0
        )
        self.max_off_nominal_roll = OptionalParameter(
            "Max Off-nominal roll (degrees)", default=2.0
        )
        self.min_n_stars = OptionalParameter("Min N stars", default=2)
        self.max_mag_diff = OptionalParameter("Max Mag diff", default=1.5)

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
        self.text_widget.setHorizontalScrollBarPolicy(QtC.Qt.ScrollBarAlwaysOff)

        font = self.text_widget.document().defaultFont()
        font.setPixelSize(12)
        font.setFamily("Courier New")
        fm = QtG.QFontMetrics(font)
        height = 12 * fm.lineSpacing()
        width = QtG.QFontMetrics(font).width("M" * 60)
        self.text_widget.document().setDefaultFont(font)
        self.text_widget.setMinimumSize(width, height)

        self.auto_update_params = QtW.QCheckBox()
        self.update_params_button = QtW.QPushButton("Update Params")
        self.find_attitude_button = QtW.QPushButton("Find Attitude")

        self._set_layout()

        self.update_params_button.clicked.connect(self.update_params)
        self.find_attitude_button.clicked.connect(self.find_attitude)

    def _set_layout(self):
        date_box = QtW.QGroupBox("Date")
        date_box.setLayout(
            vstack(
                self.date_edit,
            )
        )

        button_box = QtW.QGroupBox()
        button_box.setLayout(
            hstack(
                self.auto_update_params,
                QtW.QLabel("auto-update params"),
                QtW.QSpacerItem(
                    0,
                    0,
                    hPolicy=QtW.QSizePolicy.Expanding,
                    vPolicy=QtW.QSizePolicy.Minimum,
                ),
                self.update_params_button,
                self.find_attitude_button,
            )
        )

        tolerance_box = QtW.QGroupBox("Tolerance")
        tolerance_box.setLayout(
            vstack(
                self.tolerance_edit,
            )
        )

        att_box = QtW.QGroupBox("Latest Solution")
        att_box.setLayout(
            vstack(
                self.attitude_widget,
            )
        )

        table_box = QtW.QGroupBox("Star data")
        table_box.setLayout(vstack(self.text_widget))

        other_box = QtW.QGroupBox("Advanced Parameters")
        grid = QtW.QGridLayout()
        grid.setContentsMargins(10, 10, 10, 0)
        # grid.setSpacing(0)
        grid.setVerticalSpacing(10)
        grid.setHorizontalSpacing(0)

        grid.addWidget(self.rough_attitude_estimate, 0, 0, 1, 2)
        grid.addWidget(self.sun_pitch, 1, 0, 1, 1)
        grid.addWidget(self.sun_pitch_sigma, 1, 1, 1, 1)
        grid.addWidget(self.radial_uncertainty, 2, 0)
        grid.addWidget(self.max_off_nominal_roll, 2, 1)
        grid.addWidget(self.min_n_stars, 3, 0)
        grid.addWidget(self.max_mag_diff, 3, 1)

        other_box.setLayout(grid)

        layout = QtW.QHBoxLayout()
        layout.addStretch()
        layout.addLayout(
            vstack(
                hstack(date_box, tolerance_box),
                table_box,
                other_box,
                QtW.QSpacerItem(
                    0,
                    0,
                    hPolicy=QtW.QSizePolicy.Minimum,
                    vPolicy=QtW.QSizePolicy.Expanding,
                ),
                att_box,
                button_box,
                stretch=False,
            )
        )
        layout.addStretch()
        self.setLayout(layout)
        self.reset()

    def set_attitude(self, attitude):
        self.attitude_widget.set_attitude(attitude)
        if attitude is None:
            self.rough_attitude_estimate.enabled = False
            self.sun_pitch.enabled = False
            self.sun_pitch_sigma.enabled = False
            self.rough_attitude_estimate.reset()
            self.sun_pitch.reset()
        else:
            q = ", ".join([f"{x:.12f}" for x in attitude.q])
            self.rough_attitude_estimate.setText(q)
            pitch, _ = get_sun_pitch_yaw(attitude.ra, attitude.dec, self.date)
            self.sun_pitch.setText(f"{pitch:.2f}")

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

    def reset(self):
        self.set_date(None)
        self.attitude = None
        self.table = Table()
        self.table["slot"] = np.arange(8)
        self.table["yang"] = np.ma.masked_all(8, dtype=float)
        self.table["zang"] = np.ma.masked_all(8, dtype=float)
        self.table["mag"] = np.ma.masked_all(8, dtype=float)
        self.table["enabled"] = np.ones(8, dtype=bool)
        self.table["yang"].format = "8.2f"
        self.table["zang"].format = "8.2f"
        self.table["mag"].format = "5.2f"
        self.display_text()

    def update_params(self):
        self.set_date(None)
        self.attitude = None
        self.table = Table()
        self.table["slot"] = np.arange(8)
        self.table["yang"] = np.ma.masked_all(8, dtype=float)
        self.table["zang"] = np.ma.masked_all(8, dtype=float)
        self.table["mag"] = np.ma.masked_all(8, dtype=float)
        self.table["enabled"] = np.ones(8, dtype=bool)
        self.table["yang"].format = "8.2f"
        self.table["zang"].format = "8.2f"
        self.table["mag"].format = "5.2f"
        self.display_text()

    def set_centroids(self, centroids):
        if not np.all(centroids["IMGNUM"] == np.arange(8)):
            logger.debug(
                f"Got centroids that do not match: IMGNUM={centroids['IMGNUM']}"
            )
            return
        self.table["yang"][:] = centroids["YAGS"]
        self.table["zang"][:] = centroids["ZAGS"]
        self.table["mag"][:] = centroids["AOACMAG"]
        self.table["enabled"][centroids["IMGFID"]] = False
        self.display_text()

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

    def find_attitude(self):
        try:
            logger.debug("Find Attitude")
            print(logger.level)

            constraints = self.get_constraints()
            constraints = Constraints(**constraints)

            stars = self.table[self.table["enabled"]][["slot", "yang", "zang", "mag"]].copy()
            stars.rename_columns(
                ["yang", "zang", "mag"],
                ["YAG", "ZAG", "MAG_ACA"],
            )

            logger.debug("Star data")
            logger.debug("\n" + "\n".join(stars.pformat()))
            logger.debug("Constraints")
            logger.debug("\n" + pformat(constraints))

            solutions = find_attitude_solutions(
                np.asarray(stars),
                tolerance=float(self.tolerance_edit.text()),
                constraints=constraints
            )
            if not solutions:
                raise NoSolutionError("No solutions found")
            logger.debug(f"Found {len(solutions)} solutions")

            solution = solutions[0]
            catalog = Catalog(
                date=self.date,
                att=solution["att_fit"],
            )
            table = solution["summary"][~solution["summary"]["m_agasc_id"].mask]
            catalog["slot"] = table["slot"]
            catalog["id"] = table["m_agasc_id"]
            catalog["type"] = "GUI"
            catalog["yang"] = table["m_yag"]
            catalog["zang"] = table["m_zag"]
            catalog["mag"] = table["m_mag"]

            self.found_attitude.emit(catalog)
        except NoSolutionError as exc:
            log_exception("No solutions found", exc, level="ERROR")
            error_dialog = ErrorMessage(title="Error", message=str(exc))
            error_dialog.exec()
        except Exception as exc:
            log_exception("Error finding attitude", exc, level="ERROR")
            error_dialog = ErrorMessage(title="Error finding attitude", message=str(exc))
            error_dialog.exec()

    def get_constraints(self):
        parameters = {
            "att": self.rough_attitude_estimate.value(),
            "pitch": self.sun_pitch.value(),
            "pitch_err": self.sun_pitch_sigma.value(),
            "att_err": self.radial_uncertainty.value(),
            "off_nom_roll_max": self.max_off_nominal_roll.value(),
            "min_stars": self.min_n_stars.value(),
            "mag_err": self.max_mag_diff.value(),
        }
        return parameters
