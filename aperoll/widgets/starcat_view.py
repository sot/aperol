# from PyQt5 import QtCore as QtC, QtWidgets as QtW, QtGui as QtG
import PyQt5.QtGui as QtG
import PyQt5.QtWebEngineWidgets as QtWe
import PyQt5.QtWidgets as QtW
from PyQt5 import QtCore as QtC

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


class StarcatView(QtW.QTextEdit):
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
            self.setText(f"{STYLE}<pre>{catalog.get_text_pre()}</pre>")

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
