import json

from pygments import highlight
from pygments.formatters import HtmlFormatter
from pygments.lexers import get_lexer_by_name
from PyQt5 import QtCore as QtC
from PyQt5 import QtWidgets as QtW

from aperoll.utils import AperollException
from aperoll.widgets.error_message import ErrorMessage


class ValidationError(AperollException):
    pass


class JsonEditor(QtW.QWidget):
    """
    A widget to edit parameters.

    The parameters are stored as a dictionary and displayed as a JSON string in a text editor.
    The widget provides the text editor and two buttons to save and discard changes.

    Derived classes can override the `default_params` method to provide default parameters and
    the `validate` method to validate the parameters before saving. The `params_changed` signal
    is emitted when the parameters are saved. If there is an error in the JSON, an error dialog
    is shown.
    """

    params_changed = QtC.pyqtSignal(dict)

    def __init__(self, show_buttons=False):
        super().__init__()

        self.installEventFilter(self)

        self.text_widget = QtW.QTextEdit()
        self.discard_button = QtW.QPushButton("Discard")
        self.discard_button.clicked.connect(self.display_text)
        self.save_button = QtW.QPushButton("Save")
        self.save_button.clicked.connect(self.save)

        layout = QtW.QVBoxLayout()
        layout.addWidget(self.text_widget)
        if show_buttons:
            h_layout = QtW.QHBoxLayout()
            layout.addLayout(h_layout)
            h_layout.addWidget(self.discard_button)
            h_layout.addWidget(self.save_button)
        self.setLayout(layout)

        self._parameters = {}

        self.reset()

    def display_text(self):
        """
        Display the parameters in the text editor.
        """
        lexer = get_lexer_by_name("json")
        formatter = HtmlFormatter(full=False)
        code = json.dumps(self._parameters, indent=2)
        args_str = highlight(code, lexer, formatter)
        style = formatter.get_style_defs()
        self.text_widget.document().setDefaultStyleSheet(
            style
        )  # NOTE: not using self.styleSheet
        self.text_widget.setText(args_str)

    def reset(self):
        """
        Set the parameters to the default values.
        """
        self._parameters = self.default_params()
        self.display_text()

    def save(self):
        """
        Set the parameter values from the text editor.
        """
        try:
            params = json.loads(self.text_widget.toPlainText())
            self.validate(params)
        except json.JSONDecodeError as exc:
            error_dialog = ErrorMessage(title="JSON Error", message=str(exc))
            error_dialog.exec()
            return
        except ValidationError as exc:
            msg = str(exc)  # .replace(",", "</br>")
            error_dialog = ErrorMessage(title="Validation Error", message=msg)
            error_dialog.exec()
            return
        if params == self._parameters:
            return
        self._parameters = params
        self.params_changed.emit(self.default_params())

    def eventFilter(self, obj, event):
        """
        Listen for ctrl-S to save and escape to discard changes.
        """
        if obj == self and event.type() == QtC.QEvent.KeyPress:
            if (
                event.key() == QtC.Qt.Key_S
                and event.modifiers() == QtC.Qt.ControlModifier
            ):
                self.save()
                return True
            elif (
                event.key() == QtC.Qt.Key_Z
                and event.modifiers() == QtC.Qt.ControlModifier
            ):
                self.display_text()
                return True
            elif event.key() == QtC.Qt.Key_Escape:
                self.display_text()
                return True
        return super().eventFilter(obj, event)

    @staticmethod
    def validate(params):
        """
        Validate the parameters before saving. Raises an exception if the parameters are invalid.
        """

    @classmethod
    def default_params(cls):
        """
        Default parameters to show.
        """
        params = {}
        return params

    def get_parameters(self):
        return self._parameters

    def set_parameters(self, parameters):
        self._parameters = parameters
        self.params_changed.emit(self._parameters)

    parameters = property(get_parameters, set_parameters)

    def __getitem__(self, key):
        return self._parameters[key]

    def __setitem__(self, key, value):
        self.set_value(key, value, emit=True)

    def set_value(self, key, value, emit=True):
        if self._parameters[key] != value:
            self._parameters[key] = value
            self.display_text()
            if emit:
                self.params_changed.emit(self._parameters)


if __name__ == "__main__":
    from aca_view.tests.utils import qt

    with qt():
        app = JsonEditor()
        app.resize(1200, 800)
        app.show()
