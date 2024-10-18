import gzip
import json
import pickle
from pprint import pprint

import maude
import Ska.Sun as sun
from astropy import units as u
from cxotime.cxotime import CxoTime
from kadi.commands.observations import get_detector_and_sim_offset
from PyQt5 import QtCore as QtC
from PyQt5 import QtWidgets as QtW
from Quaternion import Quat


class LineEdit(QtW.QLineEdit):
    """
    A QLineEdit with a signal emitted when pressing Enter, loosing focus or calling setText.

    The signal is emitted only if the text changed.
    """
    value_changed = QtC.pyqtSignal(str)

    def __init__(self, parent):
        super().__init__(parent)
        self._prev_text = self.text()

    def setText(self, text):
        super().setText(text)
        self._check_value()

    def keyPressEvent(self, event):
        super().keyPressEvent(event)
        if event.key() == QtC.Qt.Key_Return:
            self._check_value()

    def focusOutEvent(self, event):
        super().focusOutEvent(event)
        self._check_value()

    def _check_value(self):
        if self.text() != self._prev_text:
            self._prev_text = self.text()
            self.value_changed.emit(self.text())


def get_default_parameters():
    """
    Get default initial parameters from current telemetry.
    """

    msid_list = ["3TSCPOS", "AACCCDPT"] + [f"aoattqt{i}".upper() for i in range(1, 5)]
    msids = maude.get_msids(msid_list)
    data = {msid: msids["data"][i]["values"][-1] for i, msid in enumerate(msid_list)}
    q = Quat(q=[data[f"AOATTQT{i}"] for i in range(1, 5)])

    instrument, sim_offset = get_detector_and_sim_offset(data["3TSCPOS"])
    t_ccd = (data["AACCCDPT"] - 32) * 5 / 9

    result = {
        "date": CxoTime().date,
        "attitude": q,
        "ra": q.ra,
        "dec": q.dec,
        "roll": q.roll,
        "instrument": instrument,
        "sim_offset": sim_offset,
        "t_ccd": t_ccd,
        "obsid": 0,
        "man_angle": 0,
        "dither_acq_y": 16,
        "dither_acq_z": 16,
        "dither_guide_y": 16,
        "dither_guide_z": 16,
        "n_fid": 0,
        "n_guide": 8,
    }

    return result


def get_parameters_from_yoshi(filename, obsid=None):
    """
    Get initial parameters from a Yoshi JSON file.
    """

    with open(filename) as fh:
        contents = json.load(fh)
        if obsid is not None:
            contents = [obs for obs in contents if obs["obsid"] == obsid]
            if not contents:
                raise Exception(f"OBSID {obsid} not found in {filename}")

        if not contents:
            raise Exception(f"No entries found in {filename}")

        yoshi_params = contents[0]  # assuming there is only one entry

        yoshi_params["date"] = yoshi_params["obs_date"]
        yoshi_params["ra"] = yoshi_params["ra_targ"]
        yoshi_params["dec"] = yoshi_params["dec_targ"]
        yoshi_params["roll"] = yoshi_params["roll_targ"]
        yoshi_params["instrument"] = yoshi_params["detector"]
        for key in [
            "obs_date",
            "ra_targ",
            "dec_targ",
            "roll_targ",
            "detector",
        ]:
            del yoshi_params[key]

    if abs(yoshi_params.get("obsid", 0)) < 38000:
        yoshi_params["n_fid"] = "3"
        yoshi_params["n_guide"] = "5"
    else:
        yoshi_params["n_fid"] = "0"
        yoshi_params["n_guide"] = "8"

    att = Quat(
        equatorial=(yoshi_params["ra"], yoshi_params["dec"], yoshi_params["roll"])
    )

    default = get_default_parameters()
    parameters = {
        "obsid": yoshi_params.get("obsid", default.get("obsid", 0)),
        "man_angle": yoshi_params.get("man_angle", default.get("man_angle", 0)),
        "dither_acq_y": yoshi_params.get("dither_y", default.get("dither_acq_y", 16)),
        "dither_acq_z": yoshi_params.get("dither_z", default.get("dither_acq_z", 16)),
        "dither_guide_y": yoshi_params.get(
            "dither_y", default.get("dither_guide_y", 16)
        ),
        "dither_guide_z": yoshi_params.get(
            "dither_z", default.get("dither_guide_z", 16)
        ),
        "date": yoshi_params.get("date", default["date"]),
        "attitude": att,
        "ra": yoshi_params.get("ra", default["attitude"].ra),
        "dec": yoshi_params.get("dec", default["attitude"].dec),
        "roll": yoshi_params.get("roll", default["attitude"].roll),
        "t_ccd": yoshi_params.get("t_ccd", default["t_ccd"]),
        "instrument": yoshi_params.get("instrument", default["instrument"]),
        "n_guide": yoshi_params["n_guide"],
        "n_fid": yoshi_params["n_fid"],
    }
    return parameters


def get_parameters_from_pickle(filename, obsid=None):
    """
    Get initial parameters from a proseco pickle file.
    """
    open_fcn = open if filename.endswith(".pkl") else gzip.open
    with open_fcn(filename, "rb") as fh:
        catalogs = pickle.load(fh)

    if float(obsid) not in catalogs:
        raise Exception(f"OBSID {obsid} not found in {filename}")

    catalog = catalogs[float(obsid)]

    parameters = {
        "obsid": obsid,
        "man_angle": catalog.man_angle,
        "dither_acq_y": catalog.dither_acq.y,
        "dither_acq_z": catalog.dither_acq.z,
        "dither_guide_y": catalog.dither_guide.y,
        "dither_guide_z": catalog.dither_guide.z,
        "date": CxoTime(catalog.date).date,  # date is not guaranteed to be a fixed type in pickle
        "attitude": catalog.att,
        "ra": catalog.att.ra,
        "dec": catalog.att.dec,
        "roll": catalog.att.roll,
        "t_ccd": catalog.t_ccd,
        "instrument": catalog.detector,
        "n_guide": catalog.n_guide,
        "n_fid": catalog.n_fid,
    }
    return parameters


class Parameters(QtW.QWidget):
    do_it = QtC.pyqtSignal()
    run_sparkles = QtC.pyqtSignal()
    draw_test = QtC.pyqtSignal()
    reset = QtC.pyqtSignal()
    parameters_changed = QtC.pyqtSignal()

    def __init__(self, **kwargs):  # noqa: PLR0915
        super().__init__()
        self.date_edit = LineEdit(self)
        self.obsid_edit = LineEdit(self)
        self.ra_edit = LineEdit(self)
        self.dec_edit = LineEdit(self)
        self.roll_edit = LineEdit(self)
        self.n_guide_edit = LineEdit(self)
        self.n_fid_edit = LineEdit(self)
        self.n_t_ccd_edit = LineEdit(self)
        self.man_angle_edit = LineEdit(self)
        self.instrument_edit = QtW.QComboBox(self)
        self.instrument_edit.addItems(["ACIS-S", "ACIS-I", "HRC-S", "HRC-I"])
        self.do = QtW.QPushButton("Get Catalog")
        self.reset_button = QtW.QPushButton("Reset")
        self.run_sparkles_button = QtW.QPushButton("Run Sparkles")
        self.draw_test_button = QtW.QPushButton("Draw Test")
        self.include = {
            "acq": QtW.QListWidget(self),
            "guide": QtW.QListWidget(self),
        }
        self.exclude = {
            "acq": QtW.QListWidget(self),
            "guide": QtW.QListWidget(self),
        }
        self.dither_acq_y_edit = LineEdit(self)
        self.dither_acq_z_edit = LineEdit(self)
        self.dither_guide_y_edit = LineEdit(self)
        self.dither_guide_z_edit = LineEdit(self)

        self.date_edit.value_changed.connect(self._values_changed)
        self.obsid_edit.value_changed.connect(self._values_changed)
        self.ra_edit.value_changed.connect(self._values_changed)
        self.dec_edit.value_changed.connect(self._values_changed)
        self.roll_edit.value_changed.connect(self._values_changed)
        self.n_guide_edit.value_changed.connect(self._values_changed)
        self.n_fid_edit.value_changed.connect(self._values_changed)
        self.n_t_ccd_edit.value_changed.connect(self._values_changed)
        self.man_angle_edit.value_changed.connect(self._values_changed)
        self.instrument_edit.currentIndexChanged.connect(self._values_changed)
        self.include["acq"].model().rowsInserted.connect(self._values_changed)
        self.include["guide"].model().rowsInserted.connect(self._values_changed)
        self.exclude["acq"].model().rowsInserted.connect(self._values_changed)
        self.exclude["guide"].model().rowsInserted.connect(self._values_changed)
        self.include["acq"].model().rowsRemoved.connect(self._values_changed)
        self.include["guide"].model().rowsRemoved.connect(self._values_changed)
        self.exclude["acq"].model().rowsRemoved.connect(self._values_changed)
        self.exclude["guide"].model().rowsRemoved.connect(self._values_changed)
        self.dither_acq_y_edit.value_changed.connect(self._values_changed)
        self.dither_acq_z_edit.value_changed.connect(self._values_changed)
        self.dither_guide_y_edit.value_changed.connect(self._values_changed)
        self.dither_guide_z_edit.value_changed.connect(self._values_changed)

        layout = QtW.QHBoxLayout()

        info_group_box = QtW.QGroupBox()
        info_group_box_layout = QtW.QGridLayout()
        info_group_box_layout.addWidget(QtW.QLabel("OBSID"), 0, 0, 1, 1)
        info_group_box_layout.addWidget(self.obsid_edit, 0, 1, 1, 2)
        info_group_box_layout.addWidget(QtW.QLabel("date"), 1, 0, 1, 1)
        info_group_box_layout.addWidget(self.date_edit, 1, 1, 1, 2)
        info_group_box_layout.addWidget(QtW.QLabel("instrument"), 2, 0, 1, 1)
        info_group_box_layout.addWidget(self.instrument_edit, 2, 1, 1, 2)
        for i in range(3):
            info_group_box_layout.setColumnStretch(i, 2)

        info_group_box.setLayout(info_group_box_layout)

        layout.addWidget(info_group_box)

        attitude_group_box = QtW.QWidget()
        attitude_group_box_layout = QtW.QGridLayout()
        attitude_group_box_layout.addWidget(QtW.QLabel("ra"), 0, 0, 1, 1)
        attitude_group_box_layout.addWidget(QtW.QLabel("dec"), 1, 0, 1, 1)
        attitude_group_box_layout.addWidget(QtW.QLabel("roll"), 2, 0, 1, 1)
        attitude_group_box_layout.addWidget(self.ra_edit, 0, 1, 1, 1)
        attitude_group_box_layout.addWidget(self.dec_edit, 1, 1, 1, 1)
        attitude_group_box_layout.addWidget(self.roll_edit, 2, 1, 1, 1)
        attitude_group_box.setLayout(attitude_group_box_layout)
        for i in range(3):
            attitude_group_box_layout.setColumnStretch(i, 10)

        info_2_group_box = QtW.QWidget()
        info_2_group_box_layout = QtW.QGridLayout()
        info_2_group_box_layout.addWidget(QtW.QLabel("n_guide"), 0, 0, 1, 1)
        info_2_group_box_layout.addWidget(self.n_guide_edit, 0, 1, 1, 1)
        info_2_group_box_layout.addWidget(QtW.QLabel("n_fid"), 0, 2, 1, 1)
        info_2_group_box_layout.addWidget(self.n_fid_edit, 0, 3, 1, 1)
        info_2_group_box_layout.addWidget(QtW.QLabel("t_ccd"), 1, 0, 1, 1)
        info_2_group_box_layout.addWidget(self.n_t_ccd_edit, 1, 1, 1, 1)
        info_2_group_box_layout.addWidget(QtW.QLabel("Man. angle"), 1, 2, 1, 1)
        info_2_group_box_layout.addWidget(self.man_angle_edit, 1, 3, 1, 1)
        info_2_group_box.setLayout(info_2_group_box_layout)
        for i in range(4):
            info_2_group_box_layout.setColumnStretch(i, 1)

        dither_group_box = QtW.QWidget()
        dither_group_box_layout = QtW.QGridLayout()
        dither_group_box_layout.addWidget(QtW.QLabel(""), 0, 0, 1, 4)
        dither_group_box_layout.addWidget(QtW.QLabel("y"), 0, 4, 1, 4)
        dither_group_box_layout.addWidget(QtW.QLabel("z"), 0, 8, 1, 4)
        dither_group_box_layout.addWidget(QtW.QLabel("acq"), 1, 0, 1, 4)
        dither_group_box_layout.addWidget(self.dither_acq_y_edit, 1, 4, 1, 4)
        dither_group_box_layout.addWidget(self.dither_acq_z_edit, 1, 8, 1, 4)

        dither_group_box_layout.addWidget(QtW.QLabel("guide"), 2, 0, 1, 4)
        dither_group_box_layout.addWidget(self.dither_guide_y_edit, 2, 4, 1, 4)
        dither_group_box_layout.addWidget(self.dither_guide_z_edit, 2, 8, 1, 4)
        dither_group_box.setLayout(dither_group_box_layout)

        tab_2 = QtW.QTabWidget()
        tab_2.addTab(attitude_group_box, "Attitude")
        tab_2.addTab(dither_group_box, "Dither")
        tab_2.addTab(info_2_group_box, "Other")
        tab_2.setCurrentIndex(0)
        layout.addWidget(tab_2)

        tab = QtW.QTabWidget()
        tab.addTab(self.include["acq"], "Include Acq.")
        tab.addTab(self.exclude["acq"], "Exclude Acq.")
        tab.addTab(self.include["guide"], "Include Guide")
        tab.addTab(self.exclude["guide"], "Exclude Guide")
        tab.setCurrentIndex(0)
        layout.addWidget(tab)

        controls_group_box = QtW.QGroupBox()
        controls_group_box_layout = QtW.QVBoxLayout()
        controls_group_box_layout.addWidget(self.do)
        controls_group_box_layout.addWidget(self.run_sparkles_button)
        controls_group_box_layout.addWidget(self.reset_button)
        controls_group_box.setLayout(controls_group_box_layout)

        layout.addWidget(controls_group_box)

        self.setLayout(layout)

        self.do.clicked.connect(self._do_it)
        self.draw_test_button.clicked.connect(self._draw_test)
        self.run_sparkles_button.clicked.connect(self.run_sparkles)
        self.reset_button.clicked.connect(self.reset)

        self.set_parameters(**kwargs)

    def set_parameters(self, **kwargs):
        if "file" in kwargs and (
            kwargs["file"].endswith(".pkl") or kwargs["file"].endswith(".pkl.gz")
        ):
            params = get_parameters_from_pickle(
                kwargs["file"], obsid=kwargs.get("obsid", None)
            )
        elif "file" in kwargs and kwargs["file"].endswith(".json"):
            params = get_parameters_from_yoshi(
                kwargs["file"], obsid=kwargs.get("obsid", None)
            )
        else:
            params = get_default_parameters()

        pprint(params)
        self.obsid_edit.setText(f"{params['obsid']}")
        self.man_angle_edit.setText(f"{params['man_angle']}")
        self.dither_acq_y_edit.setText(f"{params['dither_acq_y']}")
        self.dither_acq_z_edit.setText(f"{params['dither_acq_z']}")
        self.dither_guide_y_edit.setText(f"{params['dither_guide_y']}")
        self.dither_guide_z_edit.setText(f"{params['dither_guide_z']}")
        self.date_edit.setText(kwargs.get("date", params["date"]))
        self.ra_edit.setText(f"{params['ra']:.5f}")
        self.dec_edit.setText(f"{params['dec']:.5f}")
        self.roll_edit.setText(f"{params['roll']:.5f}")
        self.n_guide_edit.setText(f"{params['n_guide']}")
        self.n_fid_edit.setText(f"{params['n_fid']}")
        self.n_t_ccd_edit.setText(f"{params['t_ccd']:.2f}")
        self.instrument_edit.setCurrentText(params["instrument"])

        self.values = self._validate()

    def _draw_test(self):
        self.values = self._validate()
        if self.values:
            self.draw_test.emit()

    def _values_changed(self):
        # values are empty if validation fails, but the signal is still emitted to notify anyone
        # that the values have changed
        self.values = self._validate(quiet=True)
        self.parameters_changed.emit()

    def _validate(self, quiet=False):
        try:
            n_fid = int(self.n_fid_edit.text())
            n_guide = int(self.n_guide_edit.text())
            obsid = int(self.obsid_edit.text())
            assert self.date_edit.text() != "", "No date"
            assert self.ra_edit.text() != "", "No RA"
            assert self.dec_edit.text() != "", "No dec"
            assert n_fid + n_guide == 8, "n_fid + n_guide != 8"
            ra = float(self.ra_edit.text()) * u.deg
            dec = float(self.dec_edit.text()) * u.deg
            time = CxoTime(self.date_edit.text())
            if self.roll_edit.text() == "":
                roll = sun.nominal_roll(ra, dec, time)
            else:
                roll = float(self.roll_edit.text())
            return {
                "date": self.date_edit.text(),
                "ra": ra,
                "dec": dec,
                "roll": roll,
                "n_guide": n_guide,
                "n_fid": n_fid,
                "t_ccd": float(self.n_t_ccd_edit.text()),
                "instrument": self.instrument_edit.currentText(),
                "obsid": obsid,
                "exclude_ids_acq": [
                    int(self.exclude["acq"].item(i).text())
                    for i in range(self.exclude["acq"].count())
                ],
                "include_ids_acq": [
                    int(self.include["acq"].item(i).text())
                    for i in range(self.include["acq"].count())
                ],
                "exclude_ids_guide": [
                    int(self.exclude["guide"].item(i).text())
                    for i in range(self.exclude["guide"].count())
                ],
                "include_ids_guide": [
                    int(self.include["guide"].item(i).text())
                    for i in range(self.include["guide"].count())
                ],
                "dither_acq": (
                    float(self.dither_acq_y_edit.text()),
                    float(self.dither_acq_z_edit.text()),
                ),
                "dither_guide": (
                    float(self.dither_guide_y_edit.text()),
                    float(self.dither_guide_z_edit.text()),
                ),
                "man_angle": float(self.man_angle_edit.text()),
            }
        except Exception as e:
            if not quiet:
                print(e)
            return {}

    def _do_it(self):
        self.values = self._validate()
        if self.values:
            self.do_it.emit()

    def set_ra_dec(self, ra, dec, roll):
        self.ra_edit.setText(f"{ra:.8f}")
        self.dec_edit.setText(f"{dec:.8f}")
        self.roll_edit.setText(f"{roll:.8f}")

    def include_star(self, star, type, include):
        if include is True:
            self._include_star(star, type, True)
            self._exclude_star(star, type, False)
        elif include is False:
            self._include_star(star, type, False)
            self._exclude_star(star, type, True)
        else:
            self._include_star(star, type, include=False)
            self._exclude_star(star, type, exclude=False)

    def _include_star(self, star, type, include):
        items = self.include[type].findItems(f"{star}", QtC.Qt.MatchExactly)
        if include:
            if not items:
                self.include[type].addItem(f"{star}")
        else:
            for it in items:
                self.include[type].takeItem(self.include[type].row(it))

    def _exclude_star(self, star, type, exclude):
        items = self.exclude[type].findItems(f"{star}", QtC.Qt.MatchExactly)
        if exclude:
            if not items:
                self.exclude[type].addItem(f"{star}")
        else:
            for it in items:
                self.exclude[type].takeItem(self.exclude[type].row(it))

    def proseco_args(self):
        obsid = self.values["obsid"]
        ra, dec = self.values["ra"], self.values["dec"]
        roll = self.values["roll"]
        time = CxoTime(self.values["date"])

        aca_attitude = Quat(equatorial=(float(ra / u.deg), float(dec / u.deg), roll))

        args = {
            "obsid": obsid,
            "att": aca_attitude,
            "date": time,
            "n_fid": self.values["n_fid"],
            "n_guide": self.values["n_guide"],
            "dither_acq": self.values["dither_acq"],
            "dither_guide": self.values["dither_guide"],
            "t_ccd_acq": self.values["t_ccd"],
            "t_ccd_guide": self.values["t_ccd"],
            "man_angle": self.values["man_angle"],
            "detector": self.values["instrument"],
            "sim_offset": 0,  # docs say this is optional, but it does not seem to be
            "focus_offset": 0,  # docs say this is optional, but it does not seem to be
        }

        for key in [
            "exclude_ids_guide",
            "include_ids_guide",
            "exclude_ids_acq",
            "include_ids_acq",
        ]:
            if self.values[key]:
                args[key] = self.values[key]

        return args
