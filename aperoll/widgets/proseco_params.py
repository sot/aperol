from cxotime import CxoTime
from Quaternion import Quat

from aperoll.widgets.json_editor import JsonEditor, ValidationError


class ProsecoParams(JsonEditor):
    """
    A widget to edit Proseco parameters.
    """

    def __init__(self, show_buttons=False):
        super().__init__(show_buttons)
        self.setMinimumHeight(400)

    @classmethod
    def default_params(cls):
        params = {
            "date": None,
            "att": None,
            "detector": None,
            "sim_offset": None,
            "focus_offset": None,
            "t_ccd_acq": None,
            "t_ccd_guide": None,
            "obsid": 0,
            "man_angle": 0,
            "dither_acq": (16, 16),
            "dither_guide": (16, 16),
            "n_acq": 8,
            "n_fid": 0,
            "n_guide": 8,
            "exclude_ids_acq": [],
            "include_ids_acq": [],
            "exclude_ids_guide": [],
            "include_ids_guide": [],
            # "monitors": None, # proseco chokes on this
        }
        return params

    @staticmethod
    def validate(params):
        errors = []
        if params["date"] is None:
            errors.append("No date")
        if params["att"] is None:
            errors.append("No attitude")
        if params["detector"] is None:
            errors.append("No detector")
        if errors:
            raise ValidationError(", ".join(errors))

    # the following are methods to set some parameters that can be set in other parts of the GUI:
    # - the user can drag the star view to set the attitude,
    # - the attitude, the date, t_ccd and maybe others can be set from telemetry
    # - excluded/included stars can be set from the star view
    #
    # They have the option to skip the params_changed signal.

    def set_date(self, date, emit=True):
        self.set_value("date", CxoTime(date).date, emit=emit)

    def get_date(self):
        return CxoTime(self["date"])

    date = property(get_date, set_date)

    def set_attitude(self, attitude, emit=True):
        # calling self.set_value so I can skip emitting the signal
        self.set_value("att", Quat(attitude).equatorial.tolist(), emit=emit)

    def get_attitude(self):
        if self["att"] is not None:
            return Quat(self["att"])

    attitude = property(get_attitude, set_attitude)

    def set_detector(self, instrument, emit=True):
        self.set_value("detector", instrument, emit=emit)

    def set_instrument(self, instrument, emit=True):
        self.set_value("detector", instrument, emit=emit)

    def set_t_ccd(self, t_ccd, emit=True):
        self.set_value("t_ccd_guide", float(t_ccd), emit=emit)
        self.set_value("t_ccd_acq", float(t_ccd), emit=emit)

    def set_t_ccd_acq(self, t_ccd, emit=True):
        self.set_value("t_ccd_acq", float(t_ccd), emit=emit)

    def set_t_ccd_guide(self, t_ccd, emit=True):
        self.set_value("t_ccd_guide", float(t_ccd), emit=emit)

    def set_obsid(self, obsid, emit=True):
        # convenience method that should work with 1234 or "1234" or "1234.0" or 1234.1
        self.set_value("obsid", int(float(obsid) // 1), emit=emit)

    def include_star(self, star, type, include):
        """
        Force-include/exclude stars from the acq or guide list.

        Note that the state for a particular star is True, False or None. The possibilities are:

        - `include` is True: star will be included.
        - `include` is False: star will be excluded.
        - `include` is None: star will be neither included nor excluded.

        Parameters
        ----------
        star : int
            The AGASC star ID.
        type : str
            Either "acq" or "guide".
        include : bool
            Whether to include or exclude the star. True, False or None.
        """
        if include is True:
            if star not in self._parameters[f"include_ids_{type}"]:
                self._parameters[f"include_ids_{type}"].append(star)
            if star in self._parameters[f"exclude_ids_{type}"]:
                self._parameters[f"exclude_ids_{type}"].remove(star)
        elif include is False:
            if star in self._parameters[f"include_ids_{type}"]:
                self._parameters[f"include_ids_{type}"].remove(star)
            if star not in self._parameters[f"exclude_ids_{type}"]:
                self._parameters[f"exclude_ids_{type}"].append(star)
        else:
            if star in self._parameters[f"include_ids_{type}"]:
                self._parameters[f"include_ids_{type}"].remove(star)
            if star in self._parameters[f"exclude_ids_{type}"]:
                self._parameters[f"exclude_ids_{type}"].remove(star)


if __name__ == "__main__":
    from aca_view.tests.utils import qt

    with qt():
        app = ProsecoParams()
        app.resize(1200, 800)
        app.show()
