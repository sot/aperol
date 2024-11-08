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
from PyQt5 import QtCore as QtC
from PyQt5 import QtGui as QtG
from PyQt5 import QtWidgets as QtW

__all__ = [
    "Star",
    "Catalog",
    "FidLight",
    "StarcatLabel",
    "GuideStar",
    "AcqStar",
    "MonBox",
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
        ok = (
            (self.star["CLASS"] == 0)
            & (self.star["MAG_ACA"] > 5.3)
            & (self.star["MAG_ACA"] < 11.0)
            & (~np.isclose(self.star["COLOR1"], 0.7))
            & (self.star["MAG_ACA_ERR"] < 100.0)  # mag_err is in 0.01 mag
            & (self.star["ASPQ1"] < 40)
            & (  # Less than 2 arcsec centroid offset due to nearby spoiler
                self.star["ASPQ2"] == 0
            )
            & (self.star["POS_ERR"] < 3000)  # Position error < 3.0 arcsec
            & (
                (self.star["VAR"] == -9999) | (self.star["VAR"] == 5)
            )  # Not known to vary > 0.2 mag
        )
        return not ok

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

    def __init__(self, catalog, parent=None):
        super().__init__(parent)
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
            row, col = yagzag_to_pixels(
                yag, zag, allow_bad=True
            )
            # item.setPos(-yag, -zag)
            item.setPos(row, -col)

    def boundingRect(self):
        return QtC.QRectF(0, 0, 1, 1)

    def paint(self, _painter, _option, _widget):
        # this item draws nothing, it just holds children
        pass

    def __repr__(self):
        return repr(self.starcat)


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

