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
    draw_test = QtC.pyqtSignal()

    def __init__(self, **kwargs):  # noqa: PLR0915
        kwargs = {}

        super().__init__()
        self.date_label = QtW.QLabel("date")
        self.date_edit = QtW.QLineEdit(self)
        self.ra_label = QtW.QLabel("ra")
        self.ra_edit = QtW.QLineEdit(self)
        self.dec_label = QtW.QLabel("dec")
        self.dec_edit = QtW.QLineEdit(self)
        self.roll_label = QtW.QLabel("roll")
        self.roll_edit = QtW.QLineEdit(self)
        self.n_guide_label = QtW.QLabel("n_guide")
        self.n_guide_edit = QtW.QLineEdit(self)
        self.n_fid_label = QtW.QLabel("n_fid")
        self.n_fid_edit = QtW.QLineEdit(self)
        self.n_t_ccd_label = QtW.QLabel("t_ccd")
        self.n_t_ccd_edit = QtW.QLineEdit(self)
        self.instrument_label = QtW.QLabel("instrument")
        self.instrument_edit = QtW.QComboBox(self)
        self.instrument_edit.addItems(["ACIS-S", "ACIS-I", "HRC-S", "HRC-I"])
        self.do = QtW.QPushButton("Starcheck")
        self.draw_test_button = QtW.QPushButton("Draw Test")
        v_layout = QtW.QVBoxLayout(self)
        layout = QtW.QGridLayout()
        layout.addWidget(self.date_label, 0, 0)
        layout.addWidget(self.date_edit, 0, 1)
        layout.addWidget(self.ra_label, 1, 0)
        layout.addWidget(self.ra_edit, 1, 1)
        layout.addWidget(self.dec_label, 2, 0)
        layout.addWidget(self.dec_edit, 2, 1)
        layout.addWidget(self.roll_label, 3, 0)
        layout.addWidget(self.roll_edit, 3, 1)
        layout.addWidget(self.n_guide_label, 4, 0)
        layout.addWidget(self.n_guide_edit, 4, 1)
        layout.addWidget(self.n_fid_label, 5, 0)
        layout.addWidget(self.n_fid_edit, 5, 1)
        layout.addWidget(self.n_t_ccd_label, 6, 0)
        layout.addWidget(self.n_t_ccd_edit, 6, 1)
        layout.addWidget(self.instrument_label, 7, 0)
        layout.addWidget(self.instrument_edit, 7, 1)
        layout.addWidget(self.do, 8, 1)
        layout.addWidget(self.draw_test_button, 9, 1)
        v_layout.addLayout(layout)
        v_layout.addStretch(0)

        params = get_parameters()

        self.date_edit.setText(kwargs.get("date", params["date"]))
        self.ra_edit.setText(kwargs.get("ra", f"{params['attitude'].ra:.8f}"))
        self.dec_edit.setText(kwargs.get("dec", f"{params['attitude'].dec:.8f}"))
        self.roll_edit.setText(kwargs.get("roll", f"{params['attitude'].roll:.8f}"))
        self.n_guide_edit.setText(kwargs.get("n_guide", "5"))
        self.n_fid_edit.setText(kwargs.get("n_fid", "3"))
        self.n_t_ccd_edit.setText(kwargs.get("t_ccd", f"{params['t_ccd']:.2f}"))
        self.instrument_edit.setCurrentText(
            kwargs.get("instrument", params["instrument"])
        )

        self.setLayout(v_layout)

        self.values = self._validate()
        self.do.clicked.connect(self._do_it)
        self.draw_test_button.clicked.connect(self._draw_test)

    def _draw_test(self):
        self.values = self._validate()
        if self.values:
            self.draw_test.emit()

    def _validate(self, quiet=False):
        try:
            n_fid = int(self.n_fid_edit.text())
            n_guide = int(self.n_guide_edit.text())
            obsid = 0 if n_fid else 38000
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
