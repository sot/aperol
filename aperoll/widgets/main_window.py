# from PyQt5 import QtCore as QtC, QtWidgets as QtW, QtGui as QtG
from pathlib import Path
from pprint import pprint
from tempfile import TemporaryDirectory

import PyQt5.QtGui as QtG
import PyQt5.QtWebEngineWidgets as QtWe
import PyQt5.QtWidgets as QtW
import sparkles
from astropy import units as u
from cxotime import CxoTime
from proseco import get_aca_catalog
from PyQt5 import QtCore as QtC
from Quaternion import Quat

from .parameters import Parameters
from .star_plot import StarPlot

STYLE = """
  <style>
h1,h2,h3,h4 {
  color: #990000;
}

table.table-striped {
        border-width: thin thin thin thin;
        border-spacing: 1px;
        border-style: outset outset outset outset;
        border-color: gray gray gray gray;
        border-collapse: separate;
        background-color: white;
}

table.table-striped th {
        border-width: 1px 1px 1px 1px;
        padding: 1px 3px 1px 3px;
        border-color: gray;
        border-style: inset;
}
table.table-striped td {
        border-width: 1px 1px 1px 1px;
        padding: 1px 3px 1px 3px;
        border-style: inset;
        border-color: gray;
        text-align: right;
}
span.critical {
  color:#ff0000;
  font-weight:bold;
}
span.warning {
  color:#ff6400;
}
span.caution {
  color:#009900;
}
span.info {
  color:#000099;
}

span.monospace {
  font-family:monospace;
}

.callargs {
  unicode-bidi: embed;
  font-family: monospace;
  white-space: pre;
}

.hidden {
  display: none;
}

.shown {
  display: block;
}

</style>
"""


class WebPage(QtWe.QWebEnginePage):
    def __init__(self, parent=None):
        super().__init__(parent)

    # trick from https://stackoverflow.com/questions/54920726/how-make-any-link-blank-open-in-same-window-using-qwebengine
    def createWindow(self, _type):
        page = WebPage(self)
        page.urlChanged.connect(self.on_url_changed)
        return page

    @QtC.pyqtSlot(QtC.QUrl)
    def on_url_changed(self, url):
        page = self.sender()
        self.setUrl(url)
        page.deleteLater()


class MainWindow(QtW.QWidget):
    def __init__(self, opts=None):
        super().__init__()
        opts = {} if opts is None else opts
        opts = {k: opts[k] for k in opts if opts[k] is not None}

        pprint(opts)

        self.web_page = None

        self._tmp_dir = TemporaryDirectory()
        self._dir = Path(self._tmp_dir.name)

        self.plot = StarPlot()
        self.parameters = Parameters(**opts)
        self.textEdit = QtW.QTextEdit()
        font = QtG.QFont("Courier New")  # setting a fixed-width font (close enough)
        font.setPixelSize(5)  # setting a pixel size so it can be changed later
        self.textEdit.setFont(font)

        layout = QtW.QHBoxLayout(self)
        layout_2 = QtW.QVBoxLayout(self)
        layout_2.addWidget(self.parameters)
        layout_2.addWidget(self.textEdit)
        layout.addLayout(layout_2)
        layout.addWidget(self.plot)

        layout.setStretch(0, 3)
        layout.setStretch(1, 4)
        self.setLayout(layout)

        self.parameters.do_it.connect(self._run_proseco)
        self.parameters.run_sparkles.connect(self._run_sparkles)
        self.parameters.draw_test.connect(self._draw_test)
        self.plot.attitude_changed.connect(self.parameters.set_ra_dec)

        self._init()

    def closeEvent(self, event):
        if self.web_page is not None:
            del self.web_page
            self.web_page = None
        event.accept()

    def _init(self):
        if self.parameters.values:
            # obsid = self.parameters.values["obsid"]
            ra, dec = self.parameters.values["ra"], self.parameters.values["dec"]
            roll = self.parameters.values["roll"]
            time = CxoTime(self.parameters.values["date"])

            # aca_attitude = calc_aca_from_targ(
            #     Quat(equatorial=(float(ra / u.deg), float(dec / u.deg), nominal_roll)),
            #     0,
            #     0
            # )
            aca_attitude = Quat(
                equatorial=(float(ra / u.deg), float(dec / u.deg), roll)
            )
            # print("ra, dec, roll =", (float(ra / u.deg), float(dec / u.deg), roll))
            self.plot.set_base_attitude(aca_attitude, update=False)
            self.plot.set_time(time, update=True)

    def _draw_test(self):
        if self.parameters.values:
            ra, dec = self.parameters.values["ra"], self.parameters.values["dec"]
            roll = self.parameters.values["roll"]
            aca_attitude = Quat(
                equatorial=(float(ra / u.deg), float(dec / u.deg), roll)
            )
            # self.plot.show_test_stars_q(aca_attitude)
            dq = self.plot._base_attitude.dq(aca_attitude)
            self.plot.show_test_stars(
                ra_offset=dq.ra, dec_offset=dq.dec, roll_offset=dq.roll
            )

    def _proseco_args(self):
        obsid = self.parameters.values["obsid"]
        ra, dec = self.parameters.values["ra"], self.parameters.values["dec"]
        roll = self.parameters.values["roll"]
        time = CxoTime(self.parameters.values["date"])

        aca_attitude = Quat(equatorial=(float(ra / u.deg), float(dec / u.deg), roll))

        args = {
            "obsid": obsid,
            "att": aca_attitude,
            "date": time,
            "n_fid": self.parameters.values["n_fid"],
            "n_guide": self.parameters.values["n_guide"],
            "dither_acq": (16, 16),  # standard dither with ACIS
            "dither_guide": (16, 16),  # standard dither with ACIS
            "t_ccd_acq": self.parameters.values["t_ccd"],
            "t_ccd_guide": self.parameters.values["t_ccd"],
            "man_angle": 0,  # what is a sensible number to use??
            "detector": self.parameters.values["instrument"],
            "sim_offset": 0,  # docs say this is optional, but it does not seem to be
            "focus_offset": 0,  # docs say this is optional, but it does not seem to be
        }
        return args

    def _run_proseco(self):
        print("parameters:", self.parameters.values)
        if self.parameters.values:
            args = self._proseco_args()
            pprint(args)
            catalog = get_aca_catalog(**args)
            self.plot.set_catalog(catalog, update=False)

            aca = catalog.get_review_table()

            sparkles.core.check_catalog(aca)

            # aca.messages

            # self.textEdit.setText(table_to_html(catalog))
            self.textEdit.setText(f"{STYLE}<pre>{aca.get_text_pre()}</pre>")

    def resizeEvent(self, _size):
        font = self.textEdit.font()
        header = (
            "idx slot    id    type  sz p_acq  mag  mag_err "
            "maxmag   yang    zang   row    col   dim res halfw"
        )
        n_lines = 35
        scale_x = float(0.9 * self.textEdit.width()) / QtG.QFontMetrics(font).width(
            header
        )
        scale_y = float(0.9 * self.textEdit.height()) / (
            n_lines * QtG.QFontMetrics(font).height()
        )
        pix_size = int(font.pixelSize() * min(scale_x, scale_y))
        if pix_size > 0:
            font.setPixelSize(pix_size)
        self.textEdit.setFont(font)

    def _run_sparkles(self):
        # print("parameters:", self.parameters.values)
        if self.parameters.values:
            args = self._proseco_args()
            pprint(args)
            catalog = get_aca_catalog(**args)

            sparkles.run_aca_review(
                "Exploration",
                acars=[catalog.get_review_table()],
                report_dir=self._dir / "sparkles",
                report_level="all",
                roll_level="none",
            )
            print(f"sparkles report at {self._dir / 'sparkles'}")
            try:
                w = QtW.QMainWindow(self)
                w.resize(1400, 1000)
                web = QtWe.QWebEngineView(w)
                w.setCentralWidget(web)
                self.web_page = WebPage()
                web.setPage(self.web_page)
                url = self._dir / "sparkles" / "index.html"
                web.load(QtC.QUrl(f"file://{url}"))
                web.show()
                w.show()
            except Exception as e:
                print(e)

            self.plot.set_catalog(catalog, update=False)
