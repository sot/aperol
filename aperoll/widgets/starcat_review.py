# from PyQt5 import QtCore as QtC, QtWidgets as QtW, QtGui as QtG
import traceback

import PyQt5.QtGui as QtG
import PyQt5.QtWebEngineWidgets as QtWe
import PyQt5.QtWidgets as QtW
from PyQt5 import QtCore as QtC
from sparkles.core import ACAReviewTable, check_catalog

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


class StarcatReview(QtW.QTextEdit):
    def __init__(self, catalog=None, parent=None):
        super().__init__(parent)
        font = QtG.QFont("Courier New")  # setting a fixed-width font (close enough)
        font.setPixelSize(12)  # setting a pixel size so it can be changed later
        self.setFont(font)

        self.set_catalog(catalog)

    def reset(self):
        self.set_catalog(None)

    def set_catalog(self, catalog):
        if catalog is None:
            self.setText("")
        else:
            try:
                self.the_cat = catalog
                if not (
                    catalog.acqs
                    and catalog.guides
                    and catalog.dither_acq
                    and catalog.dither_guide
                ):
                    lines = catalog.pformat()
                    lines += [
                        "\n\n<span class='critical'> No review performed "
                        "(acqs/guides are empty or dither_acq/guide are None) <span>"
                    ]
                    text = "\n".join(lines)
                else:
                    aca = ACAReviewTable(catalog)
                    check_catalog(aca)
                    text = aca.get_text_pre()
            except Exception as exc:
                lines = [
                    "<span class='critical'>A review could not be performed:</span>",
                    "",
                    f'<span class="critical">    {type(exc).__name__} {exc} </span>',
                ]
                trace = traceback.extract_tb(exc.__traceback__)
                for step in trace:
                    lines.append(f"    in {step.filename}:{step.lineno}/{step.name}:")
                    lines.append(f"        {step.line}")
                text = "\n".join(lines)

            self.setText(f"{STYLE}<pre>{text}</pre>")

    def resizeEvent(self, _size):
        super().resizeEvent(_size)
        font = self.font()
        header = (
            "idx slot    id    type  sz p_acq  mag  mag_err "
            "maxmag   yang    zang   row    col   dim res halfw"
        )
        n_lines = 35
        scale_x = float(0.9 * self.width()) / QtG.QFontMetrics(font).width(header)
        scale_y = float(0.9 * self.height()) / (
            n_lines * QtG.QFontMetrics(font).height()
        )
        pix_size = int(font.pixelSize() * min(scale_x, scale_y))
        if pix_size > 0:
            font.setPixelSize(pix_size)
        self.setFont(font)
