import functools
import gzip
import json
import pickle
import traceback

import maude
import numpy as np
from chandra_aca.transform import (
    pixels_to_yagzag,
    yagzag_to_pixels,
)
from cxotime.cxotime import CxoTime
from kadi.commands.observations import get_detector_and_sim_offset
from Quaternion import Quat
from ska_helpers import logging

logger = logging.basic_logger("aperoll")


# The nominal origin of the CCD, in pixel coordinates (yagzag_to_pixels(0, 0))
CCD_ORIGIN = yagzag_to_pixels(
    0, 0
)  # (6.08840495576943, 4.92618563916467) as of this writing
# The (0,0) point of the CCD coordinates in (yag, zag)
YZ_ORIGIN = pixels_to_yagzag(0, 0)


class AperollException(RuntimeError):
    pass


def log_exception(msg, exc, level="DEBUG"):
    import logging

    trace = traceback.extract_tb(exc.__traceback__)
    level = logging.__getattribute__(level)
    logger.log(level, f"{msg}: {exc}")
    for step in trace:
        logger.log(level, f"    in {step.filename}:{step.lineno}/{step.name}:")
        logger.log(level, f"        {step.line}")


def single_entry(func):
    """
    Decorator to prevent a function from being called again before it finished (recursively).

    This can be used in a function that emits a signal. If the signal triggers a second function
    which in turn calls the first function, that can potentially create an infinite loop. This
    decorator just drops the second call.

    NOTES:

    - This decorator assumes that the decorated function has no return value.
    - This prevents any concurrent use of the funtion.

    This is not completely fail-safe, because the second function call might happen after the first
    function finished. This happens if the connection is queued. Most connections are direct
    (the slot is called immediately) but they can be queued (the slot is called later) if the caller
    and callee are in different threads, or if the signals is explicitly created as queued.

    Generally speaking, this is a hack and should not be needed, but it is useful.

    There are a few ways to break the recursion other than using this decorator:
    - make sure the setters do not enter the function if there will be no effect
      (i.e.: a set_value function checks the current value and returns if it is the same).
    - Have two signals, one that is emitted in response to "internal" changes and another that is
      emitted in response to external changes. One can then respond to internal changes only
      (e.g.: The LineEdit class as a textEdited signal emitted whenever the text is edited in the
      widget and a textChanged signal emitted whenever the text is changed using setText).
    - Have two slots: one is private and emits a signal, and the other is public and does not emit.
    """

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            if not wrapper.busy:
                wrapper.busy = True
                res = func(*args, **kwargs)
                wrapper.busy = False
                return res
        except Exception:
            wrapper.busy = False
            raise

    wrapper.busy = False
    return wrapper


def get_camera_fov_frame():
    """
    Paths that correspond ot the edges of the ACA CCD and the quadrant boundaries.
    """
    frame = {}
    N = 100
    edge_1 = np.array(
        [[-520, i] for i in np.linspace(-512, 512, N)]
        + [[i, 512] for i in np.linspace(-520, 520, N)]
        + [[520, i] for i in np.linspace(512, -512, N)]
        + [[i, -512] for i in np.linspace(520, -520, N)]
        + [[-520, 0]]
    ).T
    frame["edge_1"] = {
        "row": edge_1[0],
        "col": edge_1[1],
    }

    edge_2 = np.array(
        [[-512, i] for i in np.linspace(-512, 512, N)]
        + [[i, 512] for i in np.linspace(-512, 512, N)]
        + [[512, i] for i in np.linspace(512, -512, N)]
        + [[i, -512] for i in np.linspace(512, -512, N)]
        + [[-512, 0]]
    ).T
    frame["edge_2"] = {
        "row": edge_2[0],
        "col": edge_2[1],
    }

    cross_2 = np.array([[i, 0] for i in np.linspace(-511, 511, N)]).T
    frame["cross_2"] = {
        "row": cross_2[0],
        "col": cross_2[1],
    }

    cross_1 = np.array([[0, i] for i in np.linspace(-511, 511, N)]).T
    frame["cross_1"] = {
        "row": cross_1[0],
        "col": cross_1[1],
    }

    for value in frame.values():
        value["yag"], value["zag"] = pixels_to_yagzag(
            value["row"], value["col"], allow_bad=True
        )

    return frame


def get_default_parameters():
    """
    Get default initial parameters from current telemetry.
    """

    msid_list = ["3TSCPOS", "AACCCDPT"] + [f"aoattqt{i}".upper() for i in range(1, 5)]
    msids = maude.get_msids(msid_list)
    data = {msid: msids["data"][i]["values"][-1] for i, msid in enumerate(msid_list)}
    q = Quat(q=[data[f"AOATTQT{i}"] for i in range(1, 5)])

    instrument, sim_offset = get_detector_and_sim_offset(data["3TSCPOS"])
    t_ccd = (data["AACCCDPT"] - 32) * 5 / 9

    result = {
        "date": CxoTime().date,
        "attitude": q,
        "ra": q.ra,
        "dec": q.dec,
        "roll": q.roll,
        "instrument": instrument,
        "sim_offset": sim_offset,
        "t_ccd": t_ccd,
        "obsid": 0,
        "man_angle": 0,
        "dither_acq_y": 16,
        "dither_acq_z": 16,
        "dither_guide_y": 16,
        "dither_guide_z": 16,
        "n_fid": 0,
        "n_guide": 8,
    }

    return result


def get_parameters_from_yoshi(filename, obsid=None):
    """
    Get initial parameters from a Yoshi JSON file.
    """

    with open(filename) as fh:
        contents = json.load(fh)
        if obsid is not None:
            contents = [obs for obs in contents if obs["obsid"] == obsid]
            if not contents:
                raise AperollException(f"OBSID {obsid} not found in {filename}")

        if not contents:
            raise AperollException(f"No entries found in {filename}")

        yoshi_params = contents[0]  # assuming there is only one entry

        yoshi_params["date"] = yoshi_params["obs_date"]
        yoshi_params["ra"] = yoshi_params["ra_targ"]
        yoshi_params["dec"] = yoshi_params["dec_targ"]
        yoshi_params["roll"] = yoshi_params["roll_targ"]
        yoshi_params["instrument"] = yoshi_params["detector"]
        for key in [
            "obs_date",
            "ra_targ",
            "dec_targ",
            "roll_targ",
            "detector",
        ]:
            del yoshi_params[key]

    if abs(yoshi_params.get("obsid", 0)) < 38000:
        yoshi_params["n_fid"] = "3"
        yoshi_params["n_guide"] = "5"
    else:
        yoshi_params["n_fid"] = "0"
        yoshi_params["n_guide"] = "8"

    att = Quat(
        equatorial=(yoshi_params["ra"], yoshi_params["dec"], yoshi_params["roll"])
    )

    default = get_default_parameters()
    parameters = {
        "obsid": yoshi_params.get("obsid", default.get("obsid", 0)),
        "man_angle": yoshi_params.get("man_angle", default.get("man_angle", 0)),
        "dither_acq_y": yoshi_params.get("dither_y", default.get("dither_acq_y", 16)),
        "dither_acq_z": yoshi_params.get("dither_z", default.get("dither_acq_z", 16)),
        "dither_guide_y": yoshi_params.get(
            "dither_y", default.get("dither_guide_y", 16)
        ),
        "dither_guide_z": yoshi_params.get(
            "dither_z", default.get("dither_guide_z", 16)
        ),
        "date": yoshi_params.get("date", default["date"]),
        "attitude": att,
        "ra": yoshi_params.get("ra", default["attitude"].ra),
        "dec": yoshi_params.get("dec", default["attitude"].dec),
        "roll": yoshi_params.get("roll", default["attitude"].roll),
        "t_ccd": yoshi_params.get("t_ccd", default["t_ccd"]),
        "instrument": yoshi_params.get("instrument", default["instrument"]),
        "n_guide": yoshi_params["n_guide"],
        "n_fid": yoshi_params["n_fid"],
    }
    return parameters


def get_parameters_from_pickle(filename, obsid=None):
    """
    Get initial parameters from a proseco pickle file.
    """
    open_fcn = open if filename.endswith(".pkl") else gzip.open
    with open_fcn(filename, "rb") as fh:
        catalogs = pickle.load(fh)

    if not catalogs:
        raise AperollException(f"No entries found in {filename}")

    if obsid is None:
        # this is ugly but it works whether the keys are strings of floats or ints
        obsid = int(np.round(float(list(catalogs.keys())[0])))

    if float(obsid) not in catalogs:
        raise AperollException(f"OBSID {obsid} not found in {filename}")

    catalog = catalogs[float(obsid)]

    parameters = {
        "obsid": obsid,
        "man_angle": catalog.man_angle,
        "dither_acq_y": catalog.dither_acq.y,
        "dither_acq_z": catalog.dither_acq.z,
        "dither_guide_y": catalog.dither_guide.y,
        "dither_guide_z": catalog.dither_guide.z,
        "date": CxoTime(
            catalog.date
        ).date,  # date is not guaranteed to be a fixed type in pickle
        "attitude": catalog.att,
        "ra": catalog.att.ra,
        "dec": catalog.att.dec,
        "roll": catalog.att.roll,
        "t_ccd_acq": catalog.t_ccd_acq,
        "t_ccd_guide": catalog.t_ccd_guide,
        "instrument": catalog.detector,
        "n_guide": catalog.n_guide,
        "n_fid": catalog.n_fid,
    }
    return parameters
