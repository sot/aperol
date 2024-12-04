# from PyQt5 import QtCore as QtC, QtWidgets as QtW, QtGui as QtG
import gzip
import os
import pickle
import tarfile
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import PyQt5.QtWebEngineWidgets as QtWe
import PyQt5.QtWidgets as QtW
import sparkles
from astropy import units as u
from cxotime import CxoTime
from proseco import get_aca_catalog
from PyQt5 import QtCore as QtC
from Quaternion import Quat
from ska_helpers import utils

from aperoll.utils import logger

from .parameters import Parameters
from .star_plot import StarPlot
from .starcat_view import StarcatView

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


class MainWindow(QtW.QMainWindow):
    def __init__(self, opts=None):  # noqa: PLR0915
        super().__init__()
        opts = {} if opts is None else opts
        opts = {k: opts[k] for k in opts if opts[k] is not None}

        self.opts = opts

        self._main = QtW.QWidget()
        self.setCentralWidget(self._main)

        self.menu_bar = QtW.QMenuBar()
        self.setMenuBar(self.menu_bar)

        self.fileMenu = self.menu_bar.addMenu("&File")
        export_action = QtW.QAction("&Export Pickle", self)
        export_action.triggered.connect(self._export_proseco)
        self.fileMenu.addAction(export_action)
        export_action = QtW.QAction("&Export Sparkles", self)
        export_action.triggered.connect(self._export_sparkles)
        self.fileMenu.addAction(export_action)

        self.web_page = None

        self.plot = StarPlot()
        self.parameters = Parameters(**opts)
        self.starcat_view = StarcatView()

        layout = QtW.QVBoxLayout(self._main)
        layout_2 = QtW.QHBoxLayout()

        layout.addWidget(self.parameters)
        layout_2.addWidget(self.starcat_view)
        layout_2.addWidget(self.plot)
        layout.addLayout(layout_2)

        layout.setStretch(0, 1)  # the dialog on top should not stretch much
        layout.setStretch(1, 10)
        self._main.setLayout(layout)

        self.plot.include_star.connect(self.parameters.include_star)
        # self.plot.exclude_star.connect(self.parameters.exclude_star)

        self.parameters.do_it.connect(self._run_proseco)
        self.plot.update_proseco.connect(self._run_proseco)
        self.parameters.run_sparkles.connect(self._run_sparkles)
        self.parameters.reset.connect(self._reset)
        self.parameters.draw_test.connect(self._draw_test)
        self.parameters.parameters_changed.connect(self._parameters_changed)
        self.plot.attitude_changed_eq.connect(self.parameters.set_ra_dec)

        self._data = Data(self.parameters.proseco_args())
        self.outdir = Path(os.getcwd())

        self._init()

        starcat = None
        if "file" in opts:
            filename = opts.get("file")
            catalogs = {}
            if filename.endswith(".pkl"):
                with open(filename, "rb") as fh:
                    catalogs = pickle.load(fh)
            elif filename.endswith(".pkl.gz"):
                with gzip.open(filename, "rb") as fh:
                    catalogs = pickle.load(fh)
            if catalogs:
                obsids = [int(np.round(float(k))) for k in catalogs]
                if "obsid" not in opts or opts["obsid"] is None:
                    starcat = catalogs[obsids[0]]
                else:
                    starcat = catalogs[opts["obsid"]]
                aca = starcat.get_review_table()
                sparkles.core.check_catalog(aca)

        if starcat is not None:
            self.plot.set_catalog(starcat)
            self.starcat_view.set_catalog(aca)
            # make sure the catalog is not overwritten automatically
            self.plot.scene.state.auto_proseco = False

        if self.plot.scene.state.auto_proseco:
            self._run_proseco()

    def closeEvent(self, event):
        if self.web_page is not None:
            del self.web_page
            self.web_page = None
        event.accept()

    def _parameters_changed(self):
        proseco_args = self.parameters.proseco_args()
        self.plot.set_base_attitude(proseco_args["att"])
        self._data.reset(proseco_args)
        if self.plot.scene.state.auto_proseco and not self.plot.view.moving:
            self._run_proseco()

    def _init(self):
        if self.parameters.values:
            ra, dec = self.parameters.values["ra"], self.parameters.values["dec"]
            roll = self.parameters.values["roll"]
            time = CxoTime(self.parameters.values["date"])
            aca_attitude = Quat(
                equatorial=(float(ra / u.deg), float(dec / u.deg), roll)
            )
            self.plot.set_base_attitude(aca_attitude)
            self.plot.set_time(time)
            self.plot.scene.state = "Proseco"

    def _reset(self):
        self.parameters.set_parameters(**self.opts)
        self.starcat_view.reset()
        self._data.reset(self.parameters.proseco_args())
        self._init()

    def _draw_test(self):
        if self.parameters.values:
            ra, dec = self.parameters.values["ra"], self.parameters.values["dec"]
            roll = self.parameters.values["roll"]
            aca_attitude = Quat(
                equatorial=(float(ra / u.deg), float(dec / u.deg), roll)
            )
            dq = self.plot._base_attitude.dq(aca_attitude)
            self.plot.show_test_stars(
                ra_offset=dq.ra, dec_offset=dq.dec, roll_offset=dq.roll
            )

    def _run_proseco(self):
        """
        Display the star catalog.
        """
        if self._data.proseco:
            self.starcat_view.set_catalog(self._data.proseco["aca"])
            self.plot.set_catalog(self._data.proseco["catalog"])

    def _export_proseco(self):
        """
        Save the star catalog in a pickle file.
        """
        if self._data.proseco:
            catalog = self._data.proseco["catalog"]
            dialog = QtW.QFileDialog(
                self,
                "Export Pickle",
                str(self.outdir / f"aperoll-obsid_{catalog.obsid:.0f}.pkl"),
            )
            dialog.setAcceptMode(QtW.QFileDialog.AcceptSave)
            dialog.setDefaultSuffix("pkl")
            rc = dialog.exec()
            if rc:
                self._data.export_proseco(dialog.selectedFiles()[0])

    def _export_sparkles(self):
        """
        Save the sparkles report to a tarball.
        """
        if self._data.sparkles:
            catalog = self._data.proseco["catalog"]
            # for some reason, the extension hidden but it works
            dialog = QtW.QFileDialog(
                self,
                "Export Pickle",
                str(self.outdir / f"aperoll-obsid_{catalog.obsid:.0f}.tgz"),
            )
            dialog.setAcceptMode(QtW.QFileDialog.AcceptSave)
            dialog.setDefaultSuffix(".tgz")
            rc = dialog.exec()
            if rc:
                self._data.export_sparkles(dialog.selectedFiles()[0])

    def _run_sparkles(self):
        """
        Display the sparkles report in a web browser.
        """
        if self._data.sparkles:
            try:
                w = QtW.QMainWindow(self)
                w.resize(1400, 1000)
                web = QtWe.QWebEngineView(w)
                w.setCentralWidget(web)
                self.web_page = WebPage()
                web.setPage(self.web_page)
                url = self._data.sparkles / "index.html"
                web.load(QtC.QUrl(f"file://{url}"))
                web.show()
                w.show()
            except Exception as e:
                logger.warning(e)


class CachedVal:
    def __init__(self, func):
        self._func = func
        self.reset()

    def reset(self):
        self._value = utils.LazyVal(self._func)

    @property
    def val(self):
        return self._value.val


class Data:
    def __init__(self, parameters=None) -> None:
        self._proseco = CachedVal(self.run_proseco)
        self._sparkles = CachedVal(self.run_sparkles)
        self.parameters = parameters
        self._tmp_dir = TemporaryDirectory()
        self._dir = Path(self._tmp_dir.name)

    def reset(self, parameters):
        self.parameters = parameters
        self._proseco.reset()
        self._sparkles.reset()

    @property
    def proseco(self):
        return self._proseco.val

    @property
    def sparkles(self):
        return self._sparkles.val

    def export_proseco(self, outfile):
        if self.proseco:
            outfile = self._dir / outfile
            catalog = self.proseco["catalog"]
            if catalog:
                with open(outfile, "wb") as fh:
                    pickle.dump({catalog.obsid: catalog}, fh)

    def export_sparkles(self, outfile):
        if self.sparkles:
            outfile = self._dir / outfile
            if self.sparkles:
                dest = Path(outfile.name.replace(".tar", "").replace(".gz", ""))
                with tarfile.open(outfile, "w") as tar:
                    for name in self.sparkles.glob("**/*"):
                        tar.add(
                            name,
                            arcname=dest / name.relative_to(self._dir / "sparkles"),
                        )

    def run_proseco(self):
        if self.parameters:
            catalog = get_aca_catalog(**self.parameters)
            aca = catalog.get_review_table()
            sparkles.core.check_catalog(aca)

            return {
                "catalog": catalog,
                "aca": aca,
            }
        return {}

    def run_sparkles(self):
        if self.proseco and self.proseco["catalog"]:
            sparkles.run_aca_review(
                "Exploration",
                acars=[self.proseco["catalog"].get_review_table()],
                report_dir=self._dir / "sparkles",
                report_level="all",
                roll_level="none",
            )
            return self._dir / "sparkles"
