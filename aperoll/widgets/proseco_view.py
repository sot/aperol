import functools
from pprint import pformat

import PyQt5.QtWebEngineWidgets as QtWe
from PyQt5 import QtCore as QtC
from PyQt5 import QtWidgets as QtW
from Quaternion import Quat

from aperoll.proseco_data import ProsecoData
from aperoll.utils import logger
from aperoll.widgets.attitude_widget import (
    AttitudeWidget,
    QuatRepresentation,
)
from aperoll.widgets.error_message import ErrorMessage
from aperoll.widgets.parameters import (
    LineEdit,
    get_default_parameters,
    get_parameters_from_pickle,
    get_parameters_from_yoshi,
)
from aperoll.widgets.proseco_params import ProsecoParams
from aperoll.widgets.star_plot import StarPlot
from aperoll.widgets.starcat_review import StarcatReview


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


class ProsecoView(QtW.QWidget):
    def __init__(self, opts=None):  # noqa: PLR0915
        super().__init__()
        opts = {} if opts is None else opts
        opts = {k: opts[k] for k in opts if opts[k] is not None}

        self.create_widgets()
        self.set_connections()
        self.set_layout()
        self.add_menu()

        self.opts = opts
        self.set_parameters(**self.opts)

        self._auto_proseco()

    def add_menu(self):
        """
        Add menu actions to toolbar.
        """
        application = QtW.QApplication.instance()
        main_windows = [
            w for w in application.topLevelWidgets() if isinstance(w, QtW.QMainWindow)
        ]
        for window in main_windows:
            menu_bar = window.menuBar()
            actions = [
                action
                for action in menu_bar.actions()
                if action.text().replace("&", "") == "File"
            ]
            if actions:
                file_menu = actions[0].menu()
            else:
                file_menu = menu_bar.addMenu("&File")
            export_action = QtW.QAction("&Export Pickle", self)
            export_action.triggered.connect(self.data.export_proseco_dialog)
            file_menu.addAction(export_action)
            export_action = QtW.QAction("&Export Sparkles", self)
            export_action.triggered.connect(self.data.export_sparkles_dialog)
            file_menu.addAction(export_action)

    def create_widgets(self):
        """
        Creates all the widgets.
        """
        self.data = ProsecoData()

        self.star_plot = StarPlot(self)
        self.star_plot.scene.state = "Proseco"

        self.starcat_review = StarcatReview()
        self.date_edit = LineEdit(self)
        self.obsid_edit = LineEdit(self)
        self.attitude_widget = AttitudeWidget(
            self,
            columns={
                QuatRepresentation.EQUATORIAL: 0,
                QuatRepresentation.QUATERNION: 0,
                QuatRepresentation.SUN: 0,
            },
        )
        self.proseco_params_widget = ProsecoParams(show_buttons=True)
        self.instrument_edit = QtW.QComboBox(self)
        self.instrument_edit.addItems(["ACIS-S", "ACIS-I", "HRC-S", "HRC-I"])

        self.get_catalog_button = QtW.QPushButton("Get Catalog")
        self.reset_button = QtW.QPushButton("Reset")
        self.run_sparkles_button = QtW.QPushButton("Run Sparkles")

    def set_connections(self):
        """
        Connects signals to slots.
        """
        # connect proseco params to all other widgets
        # this is the main connection, when the proseco params change, all widgets update,
        # proseco can be called, and whatever else.
        self.proseco_params_widget.params_changed.connect(self._values_changed)

        # connect widgets to proseco params
        self.date_edit.value_changed.connect(self.proseco_params_widget.set_date)
        self.attitude_widget.attitude_changed.connect(
            self.proseco_params_widget.set_attitude
        )
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
        self.star_plot.attitude_changed.connect(self.attitude_widget.set_attitude)

        # connect buttons
        self.get_catalog_button.clicked.connect(self._get_catalog)
        self.reset_button.clicked.connect(self._reset)
        self.run_sparkles_button.clicked.connect(self._run_sparkles)

        # auto-update catalog whenever the attitude changes (other than in the proseco params)
        self.star_plot.update_proseco.connect(self._auto_proseco)
        self.attitude_widget.attitude_changed.connect(self._auto_proseco)

    def set_layout(self):
        """
        Arranges the widgets.
        """
        v_layout = QtW.QVBoxLayout()

        general_info_layout = QtW.QGridLayout()
        general_info_layout.addWidget(QtW.QLabel("OBSID"), 0, 0, 1, 1)
        general_info_layout.addWidget(self.obsid_edit, 0, 1, 1, 2)
        general_info_layout.addWidget(QtW.QLabel("date"), 1, 0, 1, 1)
        general_info_layout.addWidget(self.date_edit, 1, 1, 1, 2)
        general_info_layout.addWidget(QtW.QLabel("instrument"), 2, 0, 1, 1)
        general_info_layout.addWidget(self.instrument_edit, 2, 1, 1, 2)

        v_layout.addLayout(general_info_layout, 1)
        v_layout.addWidget(QtW.QLabel("Attitude"))
        v_layout.addWidget(self.attitude_widget, 0)
        v_layout.addWidget(QtW.QLabel("Proseco parameters"))
        v_layout.addWidget(self.proseco_params_widget, 1)

        controls_group_box = QtW.QGroupBox()
        controls_group_box_layout = QtW.QHBoxLayout()
        controls_group_box_layout.addWidget(self.get_catalog_button)
        controls_group_box_layout.addWidget(self.run_sparkles_button)
        controls_group_box_layout.addWidget(self.reset_button)
        controls_group_box.setLayout(controls_group_box_layout)

        v_layout.addWidget(controls_group_box)

        layout = QtW.QHBoxLayout()
        layout.addLayout(v_layout, 1)
        layout.addWidget(self.starcat_review, 4)
        layout.addWidget(self.star_plot, 4)

        self.setLayout(layout)

    def set_parameters(self, **kwargs):
        """
        Set the initial parameters.
        """
        if "file" in kwargs and (
            kwargs["file"].endswith(".pkl") or kwargs["file"].endswith(".pkl.gz")
        ):
            logger.debug(f"Loading parameters from {kwargs['file']}")
            params = get_parameters_from_pickle(
                kwargs["file"], obsid=kwargs.get("obsid", None)
            )
        elif "file" in kwargs and kwargs["file"].endswith(".json"):
            logger.debug(f"Loading parameters from {kwargs['file']}")
            params = get_parameters_from_yoshi(
                kwargs["file"], obsid=kwargs.get("obsid", None)
            )
        else:
            logger.debug(f"Using default parameters (file={kwargs.get('file', None)})")
            params = get_default_parameters()
            # obsid is a command-line argument, so I set it here
            if "obsid" in kwargs:
                params["obsid"] = kwargs["obsid"]

        logger.debug(pformat(params))
        self.proseco_params_widget.date = params["date"]
        self.proseco_params_widget.attitude = Quat(
            [params["ra"], params["dec"], params["roll"]]
        )
        self.proseco_params_widget["obsid"] = int(params["obsid"])
        self.proseco_params_widget["detector"] = params["instrument"]

        self.proseco_params_widget["t_ccd_acq"] = params.get(
            "t_ccd_acq", params.get("t_ccd", None)
        )
        self.proseco_params_widget["t_ccd_guide"] = params.get(
            "t_ccd_guide", params.get("t_ccd", None)
        )
        self.proseco_params_widget["man_angle"] = params["man_angle"]
        self.proseco_params_widget["dither_acq"] = (
            params["dither_acq_y"],
            params["dither_acq_z"],
        )
        self.proseco_params_widget["dither_guide"] = (
            params["dither_guide_y"],
            params["dither_guide_z"],
        )
        self.proseco_params_widget["n_acq"] = int(params.get("n_acq", 8))
        self.proseco_params_widget["n_fid"] = int(params.get("n_fid", 0))
        self.proseco_params_widget["n_guide"] = int(params.get("n_guide", 8))
        self.proseco_params_widget["sim_offset"] = int(params.get("sim_offset", 0))
        self.proseco_params_widget["focus_offset"] = int(params.get("focus_offset", 0))

    def _values_changed(self):
        self.star_plot.set_time(self.proseco_params_widget["date"])
        self.star_plot.set_base_attitude(self.proseco_params_widget.attitude)
        self.obsid_edit.setText(f"{self.proseco_params_widget['obsid']}")
        self.date_edit.setText(self.proseco_params_widget["date"])
        self.attitude_widget.date = self.proseco_params_widget.date
        self.attitude_widget.attitude = self.proseco_params_widget.attitude
        self.instrument_edit.setCurrentText(self.proseco_params_widget["detector"])

        if self.star_plot.scene.state.auto_proseco:
            self._auto_proseco()

    def _auto_proseco(self):
        self._get_catalog(interactive=False)

    def _get_catalog(self, interactive=True):
        try:
            self.data.parameters = self.proseco_params_widget._parameters
            if self.data.proseco:
                self.starcat_review.set_catalog(self.data.proseco["review_table"])
                self.star_plot.set_catalog(self.data.proseco["catalog"])
        except Exception as exc:
            if interactive:
                msg = ErrorMessage(title="Error", message=str(exc))
                msg.exec()

    def _reset(self):
        self.set_parameters(**self.opts)
        self.starcat_review.reset()

    def _run_sparkles(self):
        if self.data.sparkles:
            try:
                w = QtW.QMainWindow(self)
                w.resize(1400, 1000)
                web = QtWe.QWebEngineView(w)
                w.setCentralWidget(web)
                self.web_page = WebPage()
                web.setPage(self.web_page)
                url = self.data.sparkles / "index.html"
                web.load(QtC.QUrl(f"file://{url}"))
                web.show()
                w.show()
            except Exception as e:
                logger.warning(e)


if __name__ == "__main__":
    # this is just a simple and quick way to test the widget, not intended for normal use.
    import sys

    kwargs = {"file": sys.argv[1]} if sys.argv[1:] else {}
    logger.setLevel("INFO")
    app = QtW.QApplication([])
    widget = QtW.QMainWindow()
    widget.setCentralWidget(ProsecoView(**kwargs))
    widget.resize(1400, 800)
    widget.show()
    app.exec()
