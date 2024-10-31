import agasc
import numpy as np
from astropy.table import Table
from chandra_aca.transform import (
    pixels_to_yagzag,
    radec_to_yagzag,
    yagzag_to_pixels,
    yagzag_to_radec,
)
from cxotime import CxoTime
from PyQt5 import QtCore as QtC
from PyQt5 import QtGui as QtG
from PyQt5 import QtWidgets as QtW
from Quaternion import Quat

# The nominal origin of the CCD, in pixel coordinates (yagzag_to_pixels(0, 0))
CCD_ORIGIN = yagzag_to_pixels(
    0, 0
)  # (6.08840495576943, 4.92618563916467) as of this writing


def symsize(mag):
    # map mags to figsizes, defining
    # mag 6 as 40 and mag 11 as 3
    # interp should leave it at the bounding value outside
    # the range
    return np.interp(mag, [6.0, 11.0], [16.0, 4.0])


def get_stars(starcat_time, quaternion, radius=2):
    stars = agasc.get_agasc_cone(
        quaternion.ra, quaternion.dec, radius=radius, date=starcat_time
    )

    if "yang" not in stars.colnames or "zang" not in stars.colnames:
        stars["yang"], stars["zang"] = radec_to_yagzag(
            stars["RA_PMCORR"], stars["DEC_PMCORR"], quaternion
        )

    return stars


class Star(QtW.QGraphicsEllipseItem):
    def __init__(self, star, parent=None, highlight=False):
        s = symsize(star["MAG_ACA"]) * 10
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
        # the brightest entry in the catalog is has MAG_ACA = -1.801
        # the faintest entry in the catalog is has MAG_ACA ~ 21.5
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
        cat["row"], cat["col"] = yagzag_to_pixels(
            cat["yang"], cat["zang"], allow_bad=True
        )
        cat["ra"], cat["dec"] = yagzag_to_radec(
            cat["yang"], cat["zang"], self.starcat.att
        )
        cat["angle_halfw"] = cat["halfw"]
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
            item.setPos(-yag, -zag)

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
        s = 25
        w = 3
        s *= 5
        w *= 5
        rect = QtC.QRectF(-s, -s, 2 * s, 2 * s)
        super().__init__(rect, parent)
        self.fid = fid
        pen = QtG.QPen(QtG.QColor("red"), w)
        self.setPen(pen)
        self.setPos(-fid["yang"], -fid["zang"])

        line = QtW.QGraphicsLineItem(-s, 0, s, 0, self)
        line.setPen(pen)
        line = QtW.QGraphicsLineItem(0, -s, 0, s, self)
        line.setPen(pen)


class StarcatLabel(QtW.QGraphicsTextItem):
    def __init__(self, star, parent=None):
        self.starcat_row = star
        super().__init__(f"{star['idx']}", parent)
        self._offset = 150
        self.setFont(QtG.QFont("Arial", 150))
        self.setDefaultTextColor(QtG.QColor("red"))
        self.setPos(-star["yang"], -star["zang"])

    def setPos(self, x, y):
        rect = self.boundingRect()
        super().setPos(
            x + self._offset - rect.width() / 2, y - self._offset - rect.height() / 2
        )


class GuideStar(QtW.QGraphicsEllipseItem):
    def __init__(self, star, parent=None):
        self.starcat_row = star
        s = 90
        w = 15
        rect = QtC.QRectF(-s, -s, s * 2, s * 2)
        super().__init__(rect, parent)
        self.setPen(QtG.QPen(QtG.QColor("green"), w))
        self.setPos(-star["yang"], -star["zang"])


class AcqStar(QtW.QGraphicsRectItem):
    def __init__(self, star, parent=None):
        self.starcat_row = star
        hw = star["angle_halfw"]
        w = 15
        super().__init__(-hw, -hw, hw * 2, hw * 2, parent)
        self.setPen(QtG.QPen(QtG.QColor("blue"), w))
        self.setPos(-star["yang"], -star["zang"])


class MonBox(QtW.QGraphicsRectItem):
    def __init__(self, star, parent=None):
        self.starcat_row = star
        # starcheck convention was to plot monitor boxes at 2X halfw
        hw = star["angle_halfw"]
        w = 15
        super().__init__(-(hw * 2), -(hw * 2), hw * 4, hw * 4, parent)
        self.setPen(QtG.QPen(QtG.QColor(255, 165, 0), w))
        self.setPos(-star["yang"], -star["zang"])


class StarView(QtW.QGraphicsView):
    include_star = QtC.pyqtSignal(int, str, object)

    def __init__(self, scene=None):
        super().__init__(scene)
        # mouseTracking is set so we can show tooltips
        self.setMouseTracking(True)
        # Antialiasing could be disabled if it affects performance
        self.setRenderHint(QtG.QPainter.Antialiasing)

        self._start = None
        self._rotating = False
        self._moving = False

    def mouseMoveEvent(self, event):
        pos = event.pos()

        items = [item for item in self.items(event.pos()) if isinstance(item, Star)]
        if items:
            global_pos = event.globalPos()
            # supposedly, the following should cause the tooltip to stay for a long time
            # but it is the same
            # QtW.QToolTip.showText(global_pos, items[0].text(), self, QtC.QRect(), 1000000000)
            QtW.QToolTip.showText(global_pos, items[0].text())

        if self._start is None:
            return

        if pos != self._start:
            if event.modifiers() == QtC.Qt.ShiftModifier:
                self._rotating = True
            else:
                self._moving = True

        if self._moving or self._rotating:
            end_pos = self.mapToScene(pos)
            start_pos = self.mapToScene(self._start)
            if self._moving:
                dx, dy = end_pos.x() - start_pos.x(), end_pos.y() - start_pos.y()
                self.scene().shift_scene(dx, dy)
            elif self._rotating:
                center = self.mapToScene(self.viewport().rect().center())
                x1 = start_pos.x() - center.x()
                y1 = start_pos.y() - center.y()
                x2 = end_pos.x() - center.x()
                y2 = end_pos.y() - center.y()
                r1 = np.sqrt(x1**2 + y1**2)
                r2 = np.sqrt(x2**2 + y2**2)
                angle = np.rad2deg(np.arcsin((x1 * y2 - x2 * y1) / (r1 * r2)))
                self.scene().rotate_scene(angle, center)

            self._start = pos

    def mouseReleaseEvent(self, event):
        if event.button() == QtC.Qt.LeftButton:
            self._start = None

    def mousePressEvent(self, event):
        if event.button() == QtC.Qt.LeftButton:
            self._moving = False
            self._rotating = False
            self._start = event.pos()

    def wheelEvent(self, event):
        scale = 1 + 0.5 * event.angleDelta().y() / 360
        self.scale(scale, scale)

    def drawForeground(self, painter, _rect):
        # I want to use antialising for these lines regardless of what is set for the scene,
        # because they are large and otherwise look hideous. It will be reset at the end.
        anti_aliasing_set = painter.testRenderHint(QtG.QPainter.Antialiasing)
        painter.setRenderHint(QtG.QPainter.Antialiasing, True)

        black_pen = QtG.QPen()
        black_pen.setWidth(10)
        center = QtC.QPoint(self.viewport().width() // 2, self.viewport().height() // 2)
        center = self.mapToScene(center)

        # The following draws the edges of the CCD
        frame = _get_camera_fov_frame()

        row, col = "yag", "zag"
        painter.setPen(black_pen)
        for i in range(len(frame["edge_1"][row]) - 1):
            painter.drawLine(
                QtC.QPointF(frame["edge_1"][row][i], frame["edge_1"][col][i]),
                QtC.QPointF(frame["edge_1"][row][i + 1], frame["edge_1"][col][i + 1]),
            )
        for i in range(len(frame["edge_2"][row]) - 1):
            painter.drawLine(
                QtC.QPointF(frame["edge_2"][row][i], frame["edge_2"][col][i]),
                QtC.QPointF(frame["edge_2"][row][i + 1], frame["edge_2"][col][i + 1]),
            )

        magenta_pen = QtG.QPen(QtG.QColor("magenta"))
        magenta_pen.setWidth(10)
        painter.setPen(magenta_pen)
        for i in range(len(frame["cross_2"][row]) - 1):
            painter.drawLine(
                QtC.QPointF(frame["cross_2"][row][i], frame["cross_2"][col][i]),
                QtC.QPointF(frame["cross_2"][row][i + 1], frame["cross_2"][col][i + 1]),
            )
        for i in range(len(frame["cross_1"][row]) - 1):
            painter.drawLine(
                QtC.QPointF(frame["cross_1"][row][i], frame["cross_1"][col][i]),
                QtC.QPointF(frame["cross_1"][row][i + 1], frame["cross_1"][col][i + 1]),
            )
        painter.setRenderHint(QtG.QPainter.Antialiasing, anti_aliasing_set)

    def contextMenuEvent(self, event):
        items = [item for item in self.items(event.pos()) if isinstance(item, Star)]
        if not items:
            return
        item = items[0]

        menu = QtW.QMenu()

        incl_action = QtW.QAction("include acq", menu, checkable=True)
        incl_action.setChecked(item.included["acq"] is True)
        menu.addAction(incl_action)

        excl_action = QtW.QAction("exclude acq", menu, checkable=True)
        excl_action.setChecked(item.included["acq"] is False)
        menu.addAction(excl_action)

        incl_action = QtW.QAction("include guide", menu, checkable=True)
        incl_action.setChecked(item.included["guide"] is True)
        menu.addAction(incl_action)

        excl_action = QtW.QAction("exclude guide", menu, checkable=True)
        excl_action.setChecked(item.included["guide"] is False)
        menu.addAction(excl_action)

        result = menu.exec_(event.globalPos())
        if result is not None:
            action, action_type = result.text().split()
            if items:
                if action == "include":
                    item.included[action_type] = True if result.isChecked() else None
                elif action == "exclude":
                    item.included[action_type] = False if result.isChecked() else None
                self.include_star.emit(
                    item.star["AGASC_ID"], action_type, item.included[action_type]
                )
        event.accept()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        if event.oldSize().width() == -1 and event.oldSize().height() == -1:
            # this fits the viewport to a circle of radius 7200 plus some margine
            # (this assumes the scene is in arcsec, where the diagonal of the CCD is ~7200 arcsec)
            scale = min(event.size().height(), event.size().width()) / 10000
            self.scale(scale, scale)

class StarField(QtW.QGraphicsScene):
    attitude_changed = QtC.pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)

        self.attitude = None
        self.stars = []
        self._catalog = None

    def add_stars(self, stars):
        self.stars = stars
        self.stars["row"], self.stars["col"] = yagzag_to_pixels(
            self.stars["yang"], self.stars["zang"], allow_bad=True
        )
        black_pen = QtG.QPen()
        black_pen.setWidth(2)
        self.stars = [Star(star, highlight=False) for star in self.stars]
        for item in self.stars:
            # note that the coordinate system is (row, -col), which is (-yag, -zag)
            yag, zag = radec_to_yagzag(
                item.star["RA_PMCORR"], item.star["DEC_PMCORR"], self.attitude
            )
            item.setPos(-yag, -zag)
            self.addItem(item)

    def shift_scene(self, dx, dy):
        """
        Apply an active transformation on the scene, shifting the items.
        """
        if self.attitude is None:
            return
        dq = Quat(equatorial=[dx / 3600, dy / 3600, 0])
        self.set_attitude(self.attitude * dq)

    def rotate_scene(self, angle, around=None):
        """
        Apply an active transformation on the scene, rotating the items around the given point.
        """
        if self.attitude is None:
            return

        dq = Quat(equatorial=[0, 0, -angle])
        self.set_attitude(self.attitude * dq)

    def set_attitude(self, q):
        """
        Set the attitude of the scene, rotating the items to the given attitude.
        """
        self.attitude = q
        for item in self.stars:
            yag, zag = radec_to_yagzag(
                item.star["RA_PMCORR"], item.star["DEC_PMCORR"], self.attitude
            )
            item.setPos(-yag, -zag)

        if self._catalog is not None:
            self._catalog.set_pos_for_attitude(self.attitude)
        self.attitude_changed.emit()

    def add_catalog(self, starcat):
        if self._catalog is not None:
            self.removeItem(self._catalog)

        self._catalog = Catalog(starcat)
        self.addItem(self._catalog)

class StarPlot(QtW.QWidget):
    attitude_changed = QtC.pyqtSignal(float, float, float)
    include_star = QtC.pyqtSignal(int, str, object)

    def __init__(self, parent=None):
        super().__init__(parent)

        layout = QtW.QVBoxLayout(self)
        self.setLayout(layout)

        self._origin = [6.08840495576943, 4.92618563916467]

        self.scene = StarField(self)
        self.scene.setSceneRect(-100, -100, 200, 200)

        self.view = StarView(self.scene)

        self.layout().addWidget(self.view)

        self.stars = None
        self._time = None
        self._highlight = []
        self._catalog = None

        self.scene.attitude_changed.connect(self._attitude_changed)

        self.view.include_star.connect(self.include_star)

    def _attitude_changed(self):
        if self.scene.attitude is not None:
            self.attitude_changed.emit(
                self.scene.attitude.ra,
                self.scene.attitude.dec,
                self.scene.attitude.roll,
            )

    def set_base_attitude(self, q, update=True):
        """
        Sets the base attitude

        The base attitude is the attitude corresponding to the origin of the scene.
        When the base attitude changes, the star positions must be updated. Not doing so will
        leave the display in an inconsistent state. The "update" argument is there as a convenience
        to delay the update in case one wants to call several setters.
        """
        self.scene.set_attitude(q)
        if update:
            self.show_stars()

    def set_time(self, t, update=True):
        self._time = CxoTime(t)
        if update:
            self.show_stars()

    def highlight(self, agasc_ids, update=True):
        self._highlight = agasc_ids
        if update:
            self.show_stars()

    def set_catalog(self, catalog, update=True):
        self.set_time(catalog.date, update=False)
        self._catalog = catalog
        self.show_catalog()
        if update:
            self.show_stars()

    def show_stars(self):
        self.scene.clear()
        if self.scene.attitude is None or self._time is None:
            return
        self.stars = get_stars(self._time, self.scene.attitude)
        self.scene.add_stars(self.stars)

    def show_catalog(self):
        if self._catalog is not None:
            self.scene.add_catalog(self._catalog)


def _get_camera_fov_frame():
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

    for key in frame:
        frame[key]["yag"], frame[key]["zag"] = pixels_to_yagzag(
            frame[key]["row"], frame[key]["col"], allow_bad=True
        )

    return frame


def main():
    from aperoll.widgets.parameters import get_default_parameters

    params = get_default_parameters()

    app = QtW.QApplication([])
    w = StarPlot()
    w.set_base_attitude(params["attitude"], update=False)
    w.set_time(params["date"], update=True)
    w.resize(1500, 1000)
    w.show()
    app.exec()


if __name__ == "__main__":
    main()
