from PyQt5 import QtCore as QtC
from PyQt5 import QtWidgets as QtW


class LineEdit(QtW.QLineEdit):
    """
    A QLineEdit with a signal emitted when pressing Enter, loosing focus or calling setText.

    The signal is emitted only if the text changed.
    """

    value_changed = QtC.pyqtSignal(str)

    def __init__(self, parent):
        super().__init__(parent)
        self._prev_text = self.text()

    def setText(self, text, emit=False):
        super().setText(text)
        # in the following, emit=False so the signal is not emitted
        self._check_value(emit=emit)

    def keyPressEvent(self, event):
        super().keyPressEvent(event)
        if event.key() == QtC.Qt.Key_Return:
            self._check_value()

    def focusOutEvent(self, event):
        super().focusOutEvent(event)
        self._check_value()

    def _check_value(self, emit=True):
        if self.text() != self._prev_text:
            self._prev_text = self.text()
            if emit:
                self.value_changed.emit(self.text())
