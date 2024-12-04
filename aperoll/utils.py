import functools
import traceback

import numpy as np
from chandra_aca.transform import (
    pixels_to_yagzag,
    yagzag_to_pixels,
)
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
