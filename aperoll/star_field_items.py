"""
Collection of QGraphicsItem subclasses to represent star field items in the star field view.
"""

import numpy as np
from astropy.table import Table
from chandra_aca.transform import (
    radec_to_yagzag,
    yagzag_to_pixels,
    yagzag_to_radec,
)
from proseco.acq import get_acq_candidates_mask
from PyQt5 import QtCore as QtC
from PyQt5 import QtGui as QtG
from PyQt5 import QtWidgets as QtW

from aperoll import utils

__all__ = [
    "Star",
    "Catalog",
    "FidLight",
    "StarcatLabel",
    "GuideStar",
    "AcqStar",
    "MonBox",
    "Centroid",
    "CameraOutline",
    "FieldOfView",
]


def symsize(mag):
    # map mags to figsizes, defining
    # mag 6 as 40 and mag 11 as 3
    # interp should leave it at the bounding value outside
    # the range
    return np.interp(mag, [6.0, 11.0], [32.0, 8.0])


class Star(QtW.QGraphicsEllipseItem):
    def __init__(self, star, parent=None, highlight=False):
        s = symsize(star["MAG_ACA"])
        rect = QtC.QRectF(-s / 2, -s / 2, s, s)
        super().__init__(rect, parent)
        # self._stars = Table([star], names=star.colnames, dtype=star.dtype)
        self.star = star
        self.highlight = highlight
        color = self.color()
        self.setBrush(QtG.QBrush(color))
        self.setPen(QtG.QPen(color))
        self.included = {
            "acq": None,
            "guide": None,
        }
        # stars are stacked in z by magnitude, so small stars never hide behind big ones
        # the brightest entry in the catalog has MAG_ACA = -1.801
        # the faintest entry in the catalog has MAG_ACA ~ 21.5
        self.setZValue(20 + star["MAG_ACA"])

    def __repr__(self):
        return f"Star({self.star['AGASC_ID']})"

    def color(self):
        if self.highlight:
            return QtG.QColor("red")
        if self.star["MAG_ACA"] > 10.5:
            return QtG.QColor("lightGray")
        if self.bad():
            return QtG.QColor(255, 99, 71, 191)
        return QtG.QColor("black")

    def bad(self):
        return not get_acq_candidates_mask(self.star)

    def text(self):
        return (
            "<pre>"
            f"ID:      {self.star['AGASC_ID']}\n"
            f"mag:     {self.star['MAG_ACA']:.2f} +- {self.star['MAG_ACA_ERR']/100:.2}\n"
            f"color:   {self.star['COLOR1']:.2f}\n"
            f"ASPQ1:   {self.star['ASPQ1']}\n"
            f"ASPQ2:   {self.star['ASPQ2']}\n"
            f"class:   {self.star['CLASS']}\n"
            f"pos err: {self.star['POS_ERR']/1000} mas\n"
            f"VAR:     {self.star['VAR']}"
            "</pre>"
        )


class Catalog(QtW.QGraphicsItem):
    """
    Utility class to keep together all graphics item for a star catalog.

    Note that the position of the catalog is ALLWAYS (0,0) and the item positions need to be set
    separately for a given attitude.
    """

    def __init__(self, catalog=None, parent=None):
        super().__init__(parent)
        self.reset(catalog)
        self.setZValue(50)

    def reset(self, catalog):
        self.starcat = None
        for item in self.childItems():
            item.setParentItem(None)
            item.hide()

        if catalog is None:
            return

        self.starcat = catalog.copy()  # will add some columns

        cat = Table(self.starcat)
        # item positions are set from row/col
        cat["row"], cat["col"] = yagzag_to_pixels(
            cat["yang"], cat["zang"], allow_bad=True
        )
        # when attitude changes, the positions (row, col) are recalculated from (ra, dec)
        # so these items move with the corresponding star.
        cat["ra"], cat["dec"] = yagzag_to_radec(
            cat["yang"], cat["zang"], self.starcat.att
        )
        gui_stars = cat[(cat["type"] == "GUI") | (cat["type"] == "BOT")]
        acq_stars = cat[(cat["type"] == "ACQ") | (cat["type"] == "BOT")]
        fids = cat[cat["type"] == "FID"]
        mon_wins = cat[cat["type"] == "MON"]

        self.star_labels = [StarcatLabel(star, self) for star in cat]
        self.guide_stars = [GuideStar(gui_star, self) for gui_star in gui_stars]
        self.acq_stars = [AcqStar(acq_star, self) for acq_star in acq_stars]
        self.mon_stars = [MonBox(mon_box, self) for mon_box in mon_wins]
        self.fid_lights = [FidLight(fid, self) for fid in fids]

    def setPos(self, *_args, **_kwargs):
        # the position of the catalog is ALLWAYS (0,0)
        pass

    def set_pos_for_attitude(self, attitude):
        """
        Set the position of all items in the catalog for a given attitude.

        Calling QGraphicsItem.set_pos would not work. Children positions are relative to the
        parent, but in reality the relative distances between items changes with the attitude.
        One cannot change the position of a single item and then get the rest as a relative shift.
        Each item needs to be set individually.
        """
        # item positions are relative to the parent's position (self)
        # but the parent's position is (or should be) always (0, 0)
        for item in self.childItems():
            yag, zag = radec_to_yagzag(
                item.starcat_row["ra"], item.starcat_row["dec"], attitude
            )
            row, col = yagzag_to_pixels(yag, zag, allow_bad=True)
            # item.setPos(-yag, -zag)
            item.setPos(row, -col)

    def boundingRect(self):
        return QtC.QRectF(0, 0, 1, 1)

    def paint(self, _painter, _option, _widget):
        # this item draws nothing, it just holds children
        pass

    def __repr__(self):
        return repr(self.starcat)


class Centroid(QtW.QGraphicsEllipseItem):
    def __init__(self, imgnum, parent=None):
        self.imgnum = imgnum
        self.excluded = False
        s = 18
        w = 3
        rect = QtC.QRectF(-s, -s, 2 * s, 2 * s)
        super().__init__(rect, parent)
        color = QtG.QColor("blue")
        pen = QtG.QPen(color, w)
        self.setPen(pen)

        self._line_1 = QtW.QGraphicsLineItem(-s, 0, s, 0, self)
        self._line_1.setPen(pen)
        self._line_2 = QtW.QGraphicsLineItem(0, -s, 0, s, self)
        self._line_2.setPen(pen)

        self._label = QtW.QGraphicsTextItem(f"{imgnum}", self)
        self._label.setFont(QtG.QFont("Arial", 30))
        self._label.setDefaultTextColor(color)
        self._label.setPos(30, -30)

    def set_excluded(self, excluded):
        self.excluded = excluded
        pen = self.pen()
        color = pen.color()
        color.setAlpha(85 if excluded else 255)
        pen.setColor(color)
        self.setPen(pen)
        self._line_1.setPen(pen)
        self._line_2.setPen(pen)


class FidLight(QtW.QGraphicsEllipseItem):
    def __init__(self, fid, parent=None):
        self.starcat_row = fid
        s = 27
        w = 3
        rect = QtC.QRectF(-s, -s, 2 * s, 2 * s)
        super().__init__(rect, parent)
        self.fid = fid
        pen = QtG.QPen(QtG.QColor("red"), w)
        self.setPen(pen)
        # self.setPos(-fid["yang"], -fid["zang"])
        self.setPos(fid["row"], -fid["col"])

        line = QtW.QGraphicsLineItem(-s, 0, s, 0, self)
        line.setPen(pen)
        line = QtW.QGraphicsLineItem(0, -s, 0, s, self)
        line.setPen(pen)


class StarcatLabel(QtW.QGraphicsTextItem):
    def __init__(self, star, parent=None):
        self.starcat_row = star
        super().__init__(f"{star['idx']}", parent)
        self._offset = 30
        self.setFont(QtG.QFont("Arial", 30))
        self.setDefaultTextColor(QtG.QColor("red"))
        # self.setPos(-star["yang"], -star["zang"])
        self.setPos(star["row"], -star["col"])

    def setPos(self, x, y):
        rect = self.boundingRect()
        super().setPos(
            x + self._offset - rect.width() / 2, y - self._offset - rect.height() / 2
        )


class GuideStar(QtW.QGraphicsEllipseItem):
    def __init__(self, star, parent=None):
        self.starcat_row = star
        s = 27
        w = 5
        rect = QtC.QRectF(-s, -s, s * 2, s * 2)
        super().__init__(rect, parent)
        self.setPen(QtG.QPen(QtG.QColor("green"), w))
        # self.setPos(-star["yang"], -star["zang"])
        self.setPos(star["row"], -star["col"])


class AcqStar(QtW.QGraphicsRectItem):
    def __init__(self, star, parent=None):
        self.starcat_row = star
        hw = star["halfw"] / 5
        w = 5
        super().__init__(-hw, -hw, hw * 2, hw * 2, parent)
        self.setPen(QtG.QPen(QtG.QColor("blue"), w))
        # self.setPos(-star["yang"], -star["zang"])
        self.setPos(star["row"], -star["col"])


class MonBox(QtW.QGraphicsRectItem):
    def __init__(self, star, parent=None):
        self.starcat_row = star
        # starcheck convention was to plot monitor boxes at 2X halfw
        hw = star["halfw"] / 5
        w = 5
        super().__init__(-(hw * 2), -(hw * 2), hw * 4, hw * 4, parent)
        self.setPen(QtG.QPen(QtG.QColor(255, 165, 0), w))
        self.setPos(star["row"], -star["col"])


class FieldOfView(QtW.QGraphicsItem):
    def __init__(self, attitude=None, alternate_outline=False):
        super().__init__()
        self.camera_outline = None
        self.alternate_outline = alternate_outline
        self.attitude = attitude
        self.centroids = [Centroid(i, parent=self) for i in range(8)]
        self._centroids = np.array(
            [(0, 511, 511, 511, 511, 14, 0, 0, 0, 0) for _ in self.centroids],
            dtype=[
                ("IMGNUM", int),
                ("IMGROW0_8X8", float),
                ("IMGCOL0_8X8", float),
                ("IMGROW0", int),
                ("IMGCOL0", int),
                ("AOACMAG", float),
                ("YAGS", float),
                ("ZAGS", float),
                ("IMGFID", bool),
                ("excluded", bool),
            ],
        )
        self._centroids["IMGNUM"] = np.arange(8)
        self.show_centroids = True
        self.setZValue(80)

    def get_attitude(self):
        return self._attitude

    def set_attitude(self, attitude):
        if hasattr(self, "_attitude") and attitude == self._attitude:
            return
        self._attitude = attitude
        if self.camera_outline is None:
            self.camera_outline = CameraOutline(
                attitude, parent=self, simple=self.alternate_outline
            )
        else:
            self.camera_outline.attitude = attitude
        if self.scene() is not None and self.scene().attitude is not None:
            self.set_pos_for_attitude(self.scene().attitude)

    attitude = property(get_attitude, set_attitude)

    def set_pos_for_attitude(self, attitude):
        if self.camera_outline is not None:
            self.camera_outline.set_pos_for_attitude(attitude)
        self._set_centroid_pos_for_attitude(attitude)

    def _set_centroid_pos_for_attitude(self, attitude):
        yag, zag = self._centroids["YAGS"], self._centroids["ZAGS"]
        row0, col0 = yagzag_to_pixels(yag, zag, allow_bad=True)
        if attitude != self.attitude:
            # the first attitude is the attitude of this frame,
            # so we first get the ra/dec pointed to by the centroids
            # and then calculate the position in the scene coordinate system
            # for those ra/dec values
            ra, dec = yagzag_to_radec(yag, zag, self.attitude)
            yag, zag = radec_to_yagzag(ra, dec, attitude)
        row, col = yagzag_to_pixels(yag, zag, allow_bad=True)

        off_ccd = (
            (row0 < -511)
            | (row0 > 511)
            | (col0 < -511)
            | (col0 > 511)
            | (self._centroids["IMGFID"])
        )
        for i, centroid in enumerate(self.centroids):
            if off_ccd[i]:
                centroid.setVisible(False)
            else:
                centroid.setVisible(self._centroids_visible)
            centroid.setPos(row[i], -col[i])

    def boundingRect(self):
        return QtC.QRectF(0, 0, 1, 1)

    def paint(self, _painter, _option, _widget):
        # this item draws nothing, it just holds children
        pass

    def set_centroids(self, centroids):
        missing_cols = set(centroids.dtype.names) - set(self._centroids.dtype.names)
        if missing_cols:
            raise ValueError(f"Missing columns in centroids: {missing_cols}")

        cols = list(centroids.dtype.names)
        for col in cols:
            self._centroids[col] = centroids[col]
        self._set_centroid_pos_for_attitude(self.scene().attitude)

    def set_show_centroids(self, show=True):
        self._centroids_visible = show
        row, col = yagzag_to_pixels(
            self._centroids["YAGS"], self._centroids["ZAGS"], allow_bad=True
        )
        off_ccd = (
            (row < -511)
            | (row > 511)
            | (col < -511)
            | (col > 511)
            | (self._centroids["IMGFID"])
        )
        for i, centroid in enumerate(self.centroids):
            if off_ccd[i]:
                centroid.setVisible(False)
            else:
                centroid.setVisible(show)

    def get_show_centroids(self):
        return self._centroids_visible

    show_centroids = property(get_show_centroids, set_show_centroids)

    def get_centroid_table(self):
        self._centroids["MAG_ACA"] = [
            (centroid.excluded or not centroid.isVisible())
            for centroid in self.centroids
        ]

        table = Table()
        table["IMGNUM"] = self._centroids["IMGNUM"]
        table["MAG_ACA"] = self._centroids["MAG_ACA"]
        table["excluded"] = self._centroids["MAG_ACA"]
        table["YAG"] = self._centroids["YAGS"]
        table["ZAG"] = self._centroids["ZAGS"]

        return table


class CameraOutline(QtW.QGraphicsItem):
    def __init__(self, attitude, parent=None, simple=False):
        super().__init__(parent)
        self.simple = simple
        self._frame = utils.get_camera_fov_frame()
        self.attitude = attitude
        self.setZValue(100)

    def boundingRect(self):
        # roughly half of the CCD size in arcsec (including some margin)
        w = 2650
        return QtC.QRectF(-w, -w, 2 * w, 2 * w)

    def paint(self, painter, _option, _widget):
        color = "lightGray" if self.simple else "black"
        pen = QtG.QPen(QtG.QColor(color))
        pen.setWidth(1)
        pen.setCosmetic(True)

        # I want to use antialising for these lines regardless of what is set for the scene,
        # because they are large and otherwise look hideous. It will be reset at the end.
        anti_aliasing_set = painter.testRenderHint(QtG.QPainter.Antialiasing)
        painter.setRenderHint(QtG.QPainter.Antialiasing, True)

        painter.setPen(pen)
        for i in range(len(self._frame["edge_1"]["row"]) - 1):
            painter.drawLine(
                QtC.QPointF(
                    self._frame["edge_1"]["row"][i], -self._frame["edge_1"]["col"][i]
                ),
                QtC.QPointF(
                    self._frame["edge_1"]["row"][i + 1],
                    -self._frame["edge_1"]["col"][i + 1],
                ),
            )
        if self.simple:
            painter.setRenderHint(QtG.QPainter.Antialiasing, anti_aliasing_set)
            return

        for i in range(len(self._frame["edge_2"]["row"]) - 1):
            painter.drawLine(
                QtC.QPointF(
                    self._frame["edge_2"]["row"][i], -self._frame["edge_2"]["col"][i]
                ),
                QtC.QPointF(
                    self._frame["edge_2"]["row"][i + 1],
                    -self._frame["edge_2"]["col"][i + 1],
                ),
            )

        magenta_pen = QtG.QPen(QtG.QColor("magenta"))
        magenta_pen.setCosmetic(True)
        magenta_pen.setWidth(1)
        painter.setPen(magenta_pen)
        for i in range(len(self._frame["cross_2"]["row"]) - 1):
            painter.drawLine(
                QtC.QPointF(
                    self._frame["cross_2"]["row"][i], -self._frame["cross_2"]["col"][i]
                ),
                QtC.QPointF(
                    self._frame["cross_2"]["row"][i + 1],
                    -self._frame["cross_2"]["col"][i + 1],
                ),
            )
        for i in range(len(self._frame["cross_1"]["row"]) - 1):
            painter.drawLine(
                QtC.QPointF(
                    self._frame["cross_1"]["row"][i], -self._frame["cross_1"]["col"][i]
                ),
                QtC.QPointF(
                    self._frame["cross_1"]["row"][i + 1],
                    -self._frame["cross_1"]["col"][i + 1],
                ),
            )

        painter.setRenderHint(QtG.QPainter.Antialiasing, anti_aliasing_set)

    def set_pos_for_attitude(self, attitude):
        """
        Set the item position given the scene attitude.

        Note that the given attitude is NOT the attitude of the camera corresponding to this frame.
        It's the origin of the scene coordinate system.
        """
        if self._attitude is None:
            raise Exception("FieldOfView attitude is not set. Can't set position.")
        for key in self._frame:
            # self._frame[key]["yag"] must not be modified
            yag2, zag2 = radec_to_yagzag(
                self._frame[key]["ra"], self._frame[key]["dec"], attitude
            )
            self._frame[key]["row"], self._frame[key]["col"] = yagzag_to_pixels(
                yag2, zag2, allow_bad=True
            )
        self.update()

    def set_attitude(self, attitude):
        """
        Set the attitude of the camera corresponding to this frame.

        Note that this is not the attitude of the scene coordinate system.
        """
        if hasattr(self, "_attitude") and attitude == self._attitude:
            return
        self._attitude = attitude
        if self._attitude is None:
            for key in self._frame:
                self._frame[key]["ra"] = None
                self._frame[key]["dec"] = None
        else:
            for key in self._frame:
                self._frame[key]["ra"], self._frame[key]["dec"] = yagzag_to_radec(
                    self._frame[key]["yag"], self._frame[key]["zag"], self._attitude
                )

    def get_attitude(self):
        return self._attitude

    attitude = property(get_attitude, set_attitude)
