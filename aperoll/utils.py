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
