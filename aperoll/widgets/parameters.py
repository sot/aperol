import json
from pprint import pprint

import Ska.Sun as sun
from astropy import units as u
from cxotime.cxotime import CxoTime
from PyQt5 import QtCore as QtC
from PyQt5 import QtWidgets as QtW


def get_parameters():
    import maude
    from Quaternion import Quat

    msid_list = ["3TSCPOS", "AACCCDPT"] + [f"aoattqt{i}".upper() for i in range(1, 5)]
    msids = maude.get_msids(msid_list)
    data = {msid: msids["data"][i]["values"][-1] for i, msid in enumerate(msid_list)}
    q = Quat(q=[data[f"AOATTQT{i}"] for i in range(1, 5)])
    from kadi.commands.observations import get_detector_and_sim_offset

    instrument, sim_offset = get_detector_and_sim_offset(data["3TSCPOS"])
    t_ccd = (data["AACCCDPT"] - 32) * 5 / 9

    result = {
        "date": CxoTime().date,
        "attitude": q,
        "instrument": instrument,
        "sim_offset": sim_offset,
        "t_ccd": t_ccd,
    }
    from pprint import pprint

    pprint(result)
    return result


class Parameters(QtW.QWidget):
    do_it = QtC.pyqtSignal()
    run_sparkles = QtC.pyqtSignal()
    draw_test = QtC.pyqtSignal()

    def __init__(self, **kwargs):  # noqa: PLR0915
        super().__init__()
        self.date_edit = QtW.QLineEdit(self)
        self.obsid_edit = QtW.QLineEdit(self)
        self.ra_edit = QtW.QLineEdit(self)
        self.dec_edit = QtW.QLineEdit(self)
        self.roll_edit = QtW.QLineEdit(self)
        self.n_guide_edit = QtW.QLineEdit(self)
        self.n_fid_edit = QtW.QLineEdit(self)
        self.n_t_ccd_edit = QtW.QLineEdit(self)
        self.man_angle_edit = QtW.QLineEdit(self)
        self.instrument_edit = QtW.QComboBox(self)
        self.instrument_edit.addItems(["ACIS-S", "ACIS-I", "HRC-S", "HRC-I"])
        self.do = QtW.QPushButton("Get Catalog")
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
        self.dither_acq_y_edit = QtW.QLineEdit(self)
        self.dither_acq_z_edit = QtW.QLineEdit(self)
        self.dither_guide_y_edit = QtW.QLineEdit(self)
        self.dither_guide_z_edit = QtW.QLineEdit(self)

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
        controls_group_box.setLayout(controls_group_box_layout)

        layout.addWidget(controls_group_box)

        file = kwargs.pop("file", None)
        if file:
            with open(file) as fh:
                file_kwargs = json.load(fh)[0]  # assuming there is only one entry

                file_kwargs["date"] = file_kwargs["obs_date"]
                file_kwargs["ra"] = file_kwargs["ra_targ"]
                file_kwargs["dec"] = file_kwargs["dec_targ"]
                file_kwargs["roll"] = file_kwargs["roll_targ"]
                file_kwargs["dither"] = (
                    file_kwargs["dither_y"],
                    file_kwargs["dither_z"],
                )
                file_kwargs["instrument"] = file_kwargs["detector"]
                for key in [
                    "obs_date",
                    "ra_targ",
                    "dec_targ",
                    "roll_targ",
                    "dither_y",
                    "dither_z",
                    "detector",
                ]:
                    del file_kwargs[key]
                kwargs.update(file_kwargs)

        params = get_parameters()

        if abs(kwargs.get("obsid", 0)) < 38000:
            kwargs["n_fid"] = "3"
            kwargs["n_guide"] = "5"
        else:
            kwargs["n_fid"] = "0"
            kwargs["n_guide"] = "8"

        pprint(kwargs)
        self.obsid_edit.setText(f"{kwargs.get('obsid', params.get('obsid', 0))}")
        self.man_angle_edit.setText(
            f"{kwargs.get('man_angle', params.get('man_angle', 0))}"
        )
        self.dither_acq_y_edit.setText(
            f"{kwargs.get('dither_acq_y', params.get('dither_acq_y', 16))}"
        )
        self.dither_acq_z_edit.setText(
            f"{kwargs.get('dither_acq_z', params.get('dither_acq_z', 16))}"
        )
        self.dither_guide_y_edit.setText(
            f"{kwargs.get('dither_guide_y', params.get('dither_guide_y', 16))}"
        )
        self.dither_guide_z_edit.setText(
            f"{kwargs.get('dither_guide_z', params.get('dither_guide_z', 16))}"
        )
        self.date_edit.setText(kwargs.get("date", params["date"]))
        self.ra_edit.setText(f"{kwargs.get('ra', params['attitude'].ra):.5f}")
        self.dec_edit.setText(f"{kwargs.get('dec', params['attitude'].dec):.5f}")
        self.roll_edit.setText(f"{kwargs.get('roll', params['attitude'].roll):.5f}")
        self.n_guide_edit.setText(f"{kwargs['n_guide']}")
        self.n_fid_edit.setText(f"{kwargs['n_fid']}")
        self.n_t_ccd_edit.setText(kwargs.get("t_ccd", f"{params['t_ccd']:.2f}"))
        self.instrument_edit.setCurrentText(
            kwargs.get("instrument", params["instrument"])
        )

        self.setLayout(layout)

        self.values = self._validate()
        self.do.clicked.connect(self._do_it)
        self.draw_test_button.clicked.connect(self._draw_test)
        self.run_sparkles_button.clicked.connect(self.run_sparkles)

    def _draw_test(self):
        self.values = self._validate()
        if self.values:
            self.draw_test.emit()

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
                    int(self.dither_acq_y_edit.text()),
                    int(self.dither_acq_z_edit.text()),
                ),
                "dither_guide": (
                    int(self.dither_guide_y_edit.text()),
                    int(self.dither_guide_z_edit.text()),
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
