from PyQt5 import QtCore as QtC
from PyQt5 import QtWidgets as QtW


class ErrorMessage(QtW.QDialog):
    """
    Dialog to configure data fetching.
    """

    def __init__(self, title="", message=""):
        QtW.QDialog.__init__(self)
        self.setLayout(QtW.QVBoxLayout())
        self.resize(QtC.QSize(400, 300))

        text = f"""<h1> {title} </h1>

        <p>{message}</p>
        """
        text_box = QtW.QTextBrowser()
        text_box.setText(text)

        button_box = QtW.QDialogButtonBox(QtW.QDialogButtonBox.Ok)
        button_box.accepted.connect(self.accept)
        self.layout().addWidget(text_box)
        self.layout().addWidget(button_box)


if __name__ == "__main__":
    from aca_view.tests.utils import qt

    with qt():
        app = ErrorMessage("This is the title", "This is the message")
        app.resize(1200, 800)
        app.show()
