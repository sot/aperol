import functools
import pprint

import numpy as np
from cxotime import CxoTime
from PyQt5 import QtCore as QtC
from PyQt5 import QtWidgets as QtW
from ska_sun import get_sun_pitch_yaw

from aperoll.proseco_data import ProsecoData
from aperoll.utils import logger
from aperoll.widgets.attitude_widget import (
    AttitudeWidget,
    QuatRepresentation,
    TextEdit,
    hstack,
    vstack,
)
from aperoll.widgets.catalog_widget import CatalogWidget
from aperoll.widgets.error_message import ErrorMessage
from aperoll.widgets.find_attitude_widget import FindAttitudeWidget
from aperoll.widgets.line_edit import LineEdit
from aperoll.widgets.proseco_params import ProsecoParams
from aperoll.widgets.star_plot import StarPlot
from aperoll.widgets.starcat_review import StarcatReview


class TelemetryView(QtW.QWidget):
    def __init__(self):  # noqa: PLR0915
        super().__init__()

        self.create_widgets()
        self.set_connections()
        self.set_layout()

    def create_widgets(self):
        """
        Creates all the widgets.
        """
        self.data = ProsecoData()

        self.star_plot = StarPlot(self)
        self.star_plot.scene.state = "Telemetry"
        self.star_plot.view.setBackgroundBrush(QtC.Qt.white)

        self.starcat_review = StarcatReview()
        self.catalog_widget = CatalogWidget()
        self.find_attitude_widget = FindAttitudeWidget()
        self.date_edit = LineEdit(self)
        self.date_edit.setReadOnly(True)
        self.obsid_edit = LineEdit(self)
        self.obsid_edit.setReadOnly(True)
        self.t_ccd_edit = LineEdit(self)
        self.t_ccd_edit.setReadOnly(True)
        self.onboard_attitude_widget = AttitudeWidget(
            self,
            columns={
                QuatRepresentation.EQUATORIAL: 0,
                QuatRepresentation.QUATERNION: 0,
                QuatRepresentation.SUN: 0,
            },
        )
        self.view_attitude_widget = AttitudeWidget(
            self,
            columns={
                QuatRepresentation.EQUATORIAL: 0,
                QuatRepresentation.QUATERNION: 0,
                QuatRepresentation.SUN: 0,
            },
        )
        self.delta_q_edit = TextEdit(size=4)
        self.delta_sun_pitch_edit = LineEdit(self)
        self.rasl_edit = LineEdit(self)

        self.delta_q_edit.setReadOnly(True)
        self.delta_sun_pitch_edit.setReadOnly(True)
        self.rasl_edit.setReadOnly(True)

        self.proseco_params_widget = ProsecoParams(show_buttons=True)
        self.instrument_edit = QtW.QComboBox(self)
        self.instrument_edit.addItems(["ACIS-S", "ACIS-I", "HRC-S", "HRC-I"])

        self.get_catalog_button = QtW.QPushButton("Get Catalog")
        self.find_attitude_button = QtW.QPushButton("Find Attitude")

        self.user_attitude_button = QtW.QRadioButton("User")
        self.telemetry_attitude_button = QtW.QRadioButton("Telemetry")
        self.track_attitude_button = QtW.QRadioButton("ACA Centroids")

        self.telemetry_attitude_button.setChecked(True)

        self.attitude_source_edit = QtW.QComboBox(self)
        self.attitude_source_edit.addItems(["Telemetry", "User", "Find Attitude"])

    def get_attitude_source(self):
        return self.attitude_source_edit.currentText()

    def set_attitude_source(self, source):
        self.attitude_source_edit.setCurrentText(source)
        if source == "Telemetry":
            self.star_plot.scene.set_attitude(self.onboard_attitude_widget.attitude)
            self.star_plot.scene.onboard_attitude_slot = "attitude"

    attitude_source = property(get_attitude_source, set_attitude_source)

    def get_starcat_source(self):
        return self.catalog_widget.starcat_source

    def set_starcat_source(self, source):
        self.catalog_widget.starcat_source = source

    starcat_source = property(get_starcat_source, set_starcat_source)

    def set_connections(self):
        """
        Connects signals to slots.
        """

        # connect widgets to proseco params
        self.date_edit.value_changed.connect(self.proseco_params_widget.set_date)
        self.onboard_attitude_widget.attitude_changed.connect(self._update_delta_quat)
        self.view_attitude_widget.attitude_changed.connect(self._update_delta_quat)
        self.star_plot.attitude_changed.connect(self._update_delta_quat)
        self.view_attitude_widget.attitude_changed.connect(self._set_user_attitude)
        self.star_plot.attitude_changed.connect(self._set_user_attitude)
        self.obsid_edit.value_changed.connect(
            functools.partial(self.proseco_params_widget.set_obsid)
        )
        self.instrument_edit.currentTextChanged.connect(
            functools.partial(self.proseco_params_widget.__setitem__, "detector")
        )

        # Special case when we drag the star field. We do not want to trigger a new proseco call
        # so we skip emitting the params_changed signal
        self.star_plot.attitude_changed.connect(
            lambda att: self.proseco_params_widget.set_attitude(att, emit=False)
        )
        # and because we skipped the connection, we explicitly connect the star plot to the attitude
        # widget
        self.star_plot.attitude_changed.connect(self.view_attitude_widget.set_attitude)

        # connect buttons
        self.get_catalog_button.clicked.connect(
            lambda _check: self._get_catalog(interactive=True)
        )

        self.attitude_source_edit.currentTextChanged.connect(self.set_attitude_source)

        self.star_plot.include_slot.connect(self.find_attitude_widget.enable_slot)
        self.star_plot.include_star.connect(self.proseco_params_widget.include_star)

        self.find_attitude_button.clicked.connect(self._find_attitude)
        self.find_attitude_widget.found_attitude.connect(self._attitude_found)

    def set_layout(self):
        """
        Arranges the widgets.
        """
        self.main_tab = QtW.QTabWidget()

        general_info_layout = QtW.QGridLayout()
        general_info_layout.addWidget(QtW.QLabel("OBSID"), 0, 0, 1, 1)
        general_info_layout.addWidget(self.obsid_edit, 0, 1, 1, 2)
        general_info_layout.addWidget(QtW.QLabel("date"), 1, 0, 1, 1)
        general_info_layout.addWidget(self.date_edit, 1, 1, 1, 2)
        general_info_layout.addWidget(QtW.QLabel("T_CCD"), 2, 0, 1, 1)
        general_info_layout.addWidget(self.t_ccd_edit, 2, 1, 1, 2)
        general_info_layout.addWidget(QtW.QLabel("instrument"), 3, 0, 1, 1)
        general_info_layout.addWidget(self.instrument_edit, 3, 1, 1, 2)

        controls_group_box = QtW.QGroupBox()
        controls_group_box.setLayout(
            vstack(
                hstack(
                    QtW.QLabel("attitude source"),
                    self.attitude_source_edit,
                ),
                hstack(self.get_catalog_button, self.find_attitude_button),
            )
        )

        v_layout_left = QtW.QVBoxLayout()
        v_layout_right = QtW.QVBoxLayout()

        info_box = QtW.QGroupBox("Info")
        info_box.setLayout(general_info_layout)
        v_layout_left.addWidget(info_box)

        onboard_att_box = QtW.QGroupBox("On-board Attitude")
        onboard_att_box.setLayout(hstack(self.onboard_attitude_widget))
        v_layout_right.addWidget(onboard_att_box)

        view_att_box = QtW.QGroupBox("View Attitude")
        view_att_box.setLayout(hstack(self.view_attitude_widget))
        v_layout_right.addWidget(view_att_box)

        delta_quat_group_box = QtW.QGroupBox("Delta Quaternion")
        delta_quat_group_box.setLayout(
            hstack(
                vstack(
                    QtW.QLabel("Sun Pitch"),
                    QtW.QLabel("Sun Yaw (RASL)"),
                    QtW.QLabel("Quaternion"),
                    stretch=True,
                    spacing=15,
                ),
                vstack(
                    self.delta_sun_pitch_edit,
                    self.rasl_edit,
                    self.delta_q_edit,
                    stretch=True,
                ),
            )
        )
        v_layout_right.addWidget(delta_quat_group_box)
        v_layout_right.addStretch(1)
        v_layout_left.addStretch(1)
        v_layout_left.addWidget(controls_group_box)

        self.main_tab.addTab(self.star_plot, "Star Plot")
        self.main_tab.addTab(self.find_attitude_widget, "Find Attitude")
        self.main_tab.addTab(self.proseco_params_widget, "Proseco Parameters")
        self.main_tab.addTab(self.starcat_review, "Catalog Review")
        self.main_tab.addTab(self.catalog_widget, "Catalog")

        layout = QtW.QHBoxLayout()
        layout.addLayout(v_layout_left, 1)
        layout.addWidget(self.main_tab, 3)
        layout.addLayout(v_layout_right, 1)

        self.setLayout(layout)

    def _set_user_attitude(self):
        self.attitude_source = "User"

    def _find_attitude(self):
        # catalog should be cleared before setting new catalog, in case there is an error
        self.catalog_widget.reset()
        self.find_attitude_widget.find_attitude()

    def _attitude_found(self, solution):
        logger.debug("Solution found")
        logger.debug(solution.att)
        logger.debug(solution.date)
        logger.debug(pprint.pformat(solution))
        self.catalog_widget.set_catalog(solution)


    def _include_star(self, star_id, action, include):
        print(f"Include star {star_id} {action} {include}")

    def _include_slot(self, slot_id, include):
        print(f"Include slot {slot_id} {include}")

    def _get_catalog(self, interactive=True):
        try:
            # catalog should be cleared before setting new catalog, in case there is an error
            self.data.parameters = self.proseco_params_widget._parameters
            if self.data.proseco:
                self.starcat_source = "Proseco"
                self.set_catalog(self.data.proseco["catalog"], "Proseco")
        except Exception as exc:
            if interactive:
                msg = ErrorMessage(title="Error", message=str(exc))
                msg.exec()

    def _reset(self):
        self.starcat_review.reset()

    def set_catalog(self, catalog, source):
        if self.starcat_source == source:
            self.starcat_review.set_catalog(catalog)
            self.catalog_widget.set_catalog(catalog)
            self.star_plot.set_catalog(catalog)

    def set_time(self, time):
        self.proseco_params_widget.set_date(time, emit=True)
        self.star_plot.set_time(time)
        self.onboard_attitude_widget.set_date(time)
        self.view_attitude_widget.set_date(time)
        self.date_edit.setText(CxoTime(time).date)
        self._update_delta_quat()
        self.find_attitude_widget.set_date(time)

    def set_centroids(self, centroids):
        self.star_plot.set_centroids(centroids)
        self.find_attitude_widget.set_centroids(centroids)

    def set_onboard_attitude(self, att):
        self.onboard_attitude_widget.set_attitude(att)
        self.star_plot.set_onboard_attitude(att)
        self.find_attitude_widget.set_attitude(att)
        self._update_delta_quat()

    def set_telemetry(self, telemetry):
        if not isinstance(telemetry["COBSRQID"], np.ma.core.MaskedConstant):
            self.proseco_params_widget.set_obsid(telemetry["COBSRQID"])
            self.obsid_edit.setText(str(telemetry["COBSRQID"]))
        if not isinstance(telemetry["AACCCDPT"], np.ma.core.MaskedConstant):
            t_ccd = telemetry["AACCCDPT"] - 273.15
            self.t_ccd_edit.setText(f"{t_ccd:.2f}")
            self.proseco_params_widget.set_t_ccd(t_ccd)

    def _update_delta_quat(self):
        if (
            self.onboard_attitude_widget.attitude is None
            or self.view_attitude_widget.attitude is None
        ):
            return
        # dq is the quaternion such that:
        # view_attitude = onboard_attitude * dq
        dq = self.onboard_attitude_widget.attitude.dq(
            self.view_attitude_widget.attitude
        )
        self.delta_q_edit.set_values(dq.q)

        date = self.date_edit.text().strip()
        if date:
            date = CxoTime(date)
            onboard_pitch, onboard_yaw = get_sun_pitch_yaw(
                self.onboard_attitude_widget.attitude.ra,
                self.onboard_attitude_widget.attitude.dec,
                date,
            )
            view_pitch, view_yaw = get_sun_pitch_yaw(
                self.view_attitude_widget.attitude.ra,
                self.view_attitude_widget.attitude.dec,
                date,
            )
            self.delta_sun_pitch_edit.setText(f"{view_pitch - onboard_pitch:.2f}")
            self.rasl_edit.setText(f"{view_yaw - onboard_yaw:.2f}")


def _test_data():
    # example data to test the display
    from Quaternion import Quat

    telem = np.array(
        (
            8.47854822e08,
            1.6959999,
            8.47854823e08,
            28843,
            4,
            99810,
            36,
            -0.76007281,
            -0.30605962,
            -0.1090998,
            0.56277355,
            "NPNT",
            "KALM",
            273.15,
            0,
            0,
            "ORIG",
            "DYNB",
            12775716,
            "2024292",
            "2024:318:02:40:22.946",
        ),
        dtype=[
            ("TIME", "<f8"),
            ("INTEG", "<f4"),
            ("END_INTEG_TIME", "<f8"),
            ("COBSRQID", "<u4"),
            ("AOKALSTR", "i1"),
            ("MJF", "<u4"),
            ("MNF", "u1"),
            ("AOATTQT1", "<f8"),
            ("AOATTQT2", "<f8"),
            ("AOATTQT3", "<f8"),
            ("AOATTQT4", "<f8"),
            ("AOPCADMD", "<U4"),
            ("AOACASEQ", "<U4"),
            ("AACCCDPT", "<f8"),
            ("COMMCNT", "u1"),
            ("COMMPROG", "u1"),
            ("AAPIXTLM", "<U4"),
            ("AABGDTYP", "<U4"),
            ("VCDUCTR", "<u4"),
            ("dark_cal_id", "<U7"),
            ("starcat_date", "<U21"),
        ],
    )
    time = "2024:318:03:12:32.530"
    q = Quat(
        [
            -0.7600728109882766,
            -0.30605961532546644,
            -0.10909980369319783,
            0.5627735485595622,
        ]
    )
    centroids = np.array(
        [
            (0, 155.0, -377.0, 155, -377, 7.0625, -764.975, -1876.875, True),
            (1, -430.0, 8.0, -430, 8, 7.1875, 2146.425, 31.1, True),
            (2, 368.0, 6.0, 368, 6, 7.1875, -1819.0, 26.6, True),
            (3, -168.0, -279.0, -168, -279, 7.875, 843.05, -1395.725, False),
            (4, 137.0, 152.0, 137, 152, 9.0625, -673.325, 755.15, False),
            (5, 141.0, 413.0, 141, 413, 9.375, -686.725, 2047.4, False),
            (6, 139.0, -137.0, 139, -137, 10.3125, -682.9, -688.475, False),
            (7, 81.0, 208.0, 81, 208, 9.3125, -390.85, 1030.925, False),
        ],
        dtype=[
            ("IMGNUM", "<i8"),
            ("IMGROW0_8X8", "<f8"),
            ("IMGCOL0_8X8", "<f8"),
            ("IMGROW0", "<i8"),
            ("IMGCOL0", "<i8"),
            ("AOACMAG", "<f8"),
            ("YAGS", "<f8"),
            ("ZAGS", "<f8"),
            ("IMGFID", "?"),
        ],
    )
    return {"telemetry": telem, "time": time, "q": q, "centroids": centroids}


if __name__ == "__main__":
    import logging
    import os
    import sys

    from aca_view import resources
    sys.path.insert(0, os.path.dirname(resources.__file__))

    import qdarkstyle

    # this is just a simple and quick way to test the widget, not intended for normal use.
    logger.setLevel("DEBUG")

    # find_attitude logger has level="INFO" by default
    # I would rather redirect that logger's info into debug.
    logging.getLogger("find_attitude").setLevel(logging.WARNING)

    app = QtW.QApplication([])
    app.setStyleSheet(qdarkstyle.load_stylesheet())

    widget = QtW.QMainWindow()
    telemetry_view = TelemetryView()
    widget.setCentralWidget(telemetry_view)
    widget.resize(1400, 800)

    test_data = _test_data()
    telemetry_view.set_time(test_data["time"])
    telemetry_view.set_onboard_attitude(test_data["q"])
    telemetry_view.set_centroids(test_data["centroids"])
    telemetry_view.set_telemetry(test_data["telemetry"])

    widget.show()
    app.exec()
