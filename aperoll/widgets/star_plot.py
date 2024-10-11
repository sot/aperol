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
    import agasc
    from Ska.quatutil import radec2yagzag

    stars = agasc.get_agasc_cone(
        quaternion.ra, quaternion.dec, radius=radius, date=starcat_time
    )

    if "yang" not in stars.colnames or "zang" not in stars.colnames:
        # Add star Y angle and Z angle in arcsec to the stars table.
        # radec2yagzag returns degrees.
        yags, zags = radec2yagzag(stars["RA_PMCORR"], stars["DEC_PMCORR"], quaternion)
        stars["yang"] = yags * 3600
        stars["zang"] = zags * 3600

    return stars


class StarView(QtW.QGraphicsView):
    roll_changed = QtC.pyqtSignal(float)

    def __init__(self, scene=None):
        super().__init__(scene)

        self._start = None
        self._rotating = False
        self._moving = False

    def mouseMoveEvent(self, event):
        pos = event.pos()
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

                scene_rect = self.scene().sceneRect()
                new_scene_rect = QtC.QRectF(
                    scene_rect.x() - dx,
                    scene_rect.y() - dy,
                    scene_rect.width(),
                    scene_rect.height(),
                )
                self.scene().setSceneRect(new_scene_rect)
            elif self._rotating:
                center = self.mapToScene(self.viewport().rect().center())
                x1 = start_pos.x() - center.x()
                y1 = start_pos.y() - center.y()
                x2 = end_pos.x() - center.x()
                y2 = end_pos.y() - center.y()
                r1 = np.sqrt(x1**2 + y1**2)
                r2 = np.sqrt(x2**2 + y2**2)
                angle = np.rad2deg(np.arcsin((x1 * y2 - x2 * y1) / (r1 * r2)))
                transform = self.viewportTransform().rotate(angle)
                self.setTransform(transform)
                self.roll_changed.emit(self.get_roll_offset())

            self._start = pos

    def mouseReleaseEvent(self, event):
        self._start = None

    def mousePressEvent(self, event):
        self._moving = False
        self._rotating = False
        self._start = event.pos()

    def wheelEvent(self, event):
        scale = 1 + 0.5 * event.angleDelta().y() / 360
        self.scale(scale, scale)

    def drawForeground(self, painter, rect):
        black_pen = QtG.QPen()
        black_pen.setWidth(2)
        center = QtC.QPoint(self.viewport().width() // 2, self.viewport().height() // 2)
        center = self.mapToScene(center)

        # this rectangle is adapted to the "current" field of view
        # and it moves together with the stars in the scene
        # b1hw = 512
        # painter.drawRect(-b1hw, -b1hw, 2 * b1hw, 2 * b1hw)

        # these was to fix the frame to the view by rotating the viewport
        # everything drawn after this will not rotate when the canvas is rotated
        # I don't remember why it was needed. Leaving it here for now.
        # since then, I chenged the transforms
        # transform = self.viewportTransform()
        # t11 = transform.m11()
        # t12 = transform.m12()
        # angle = np.rad2deg(np.arctan2(t12, t11))
        # painter.translate(center)
        # painter.rotate(-angle)
        # painter.translate(-center)

        # The following draws the edges of the CCD
        att_offset = self.get_attitude_offset(no_roll=True).inv()

        painter.setPen(black_pen)
        N = 100
        edge_1 = np.array(
            [[-520, i] for i in np.linspace(-512, 512, N)]
            + [[i, 512] for i in np.linspace(-520, 520, N)]
            + [[520, i] for i in np.linspace(512, -512, N)]
            + [[i, -512] for i in np.linspace(520, -520, N)]
            + [[-520, 0]]
        ).T
        row2, col2 = self.transform_pixels_to_attitude(edge_1[0], edge_1[1], att_offset)
        for i in range(len(row2) - 1):
            painter.drawLine(
                QtC.QPointF(row2[i], col2[i]), QtC.QPointF(row2[i + 1], col2[i + 1])
            )

        edge_2 = np.array(
            [[-512, i] for i in np.linspace(-512, 512, N)]
            + [[i, 512] for i in np.linspace(-512, 512, N)]
            + [[512, i] for i in np.linspace(512, -512, N)]
            + [[i, -512] for i in np.linspace(512, -512, N)]
            + [[-512, 0]]
        ).T
        row2, col2 = self.transform_pixels_to_attitude(edge_2[0], edge_2[1], att_offset)
        for i in range(len(row2) - 1):
            painter.drawLine(
                QtC.QPointF(row2[i], col2[i]), QtC.QPointF(row2[i + 1], col2[i + 1])
            )

        painter.setPen(QtG.QPen(QtG.QColor("magenta")))
        cross_2 = np.array([[i, 0] for i in np.linspace(-511, 511, N)]).T
        row2, col2 = self.transform_pixels_to_attitude(
            cross_2[0], cross_2[1], att_offset
        )
        for i in range(len(row2) - 1):
            painter.drawLine(
                QtC.QPointF(row2[i], col2[i]), QtC.QPointF(row2[i + 1], col2[i + 1])
            )

        cross_1 = np.array([[0, i] for i in np.linspace(-511, 511, N)]).T
        row2, col2 = self.transform_pixels_to_attitude(
            cross_1[0], cross_1[1], att_offset
        )
        for i in range(len(row2) - 1):
            painter.drawLine(
                QtC.QPointF(row2[i], col2[i]), QtC.QPointF(row2[i + 1], col2[i + 1])
            )

    def transform_pixels_to_attitude(self, row, col, q):
        row, col = self.transform_pixels_to_attitude_2(row, col, q)
        row, col = self.transform_pixels_to_attitude_1(row, col, q)
        return row, col

    def transform_pixels_to_attitude_1(self, row, col, q):
        q0 = Quat(q=[0, 0, 0, 1])
        y, z = pixels_to_yagzag(row, col, allow_bad=True)
        ra, dec = yagzag_to_radec(y, z, q0)
        y2, z2 = radec_to_yagzag(ra, dec, q)
        row, col = yagzag_to_pixels(y2, z2, allow_bad=True)
        return row, col

    def transform_pixels_to_attitude_2(self, row, col, q):
        transform = self.viewportTransform()
        row_col = np.array([row, col])
        transform = np.array(
            [[transform.m11(), transform.m12()], [transform.m21(), transform.m22()]]
        )
        # all scaling transforms we apply are isotropic, therefore the two eigenvalues are equal,
        # so the following turns the transform into a rotation matrix
        transform /= np.sqrt(np.linalg.det(transform))
        # the following is a sanity check (that it is actually a rotation)
        assert np.allclose(transform[0, 1], -transform[1, 0])
        assert np.allclose(np.linalg.det(transform), 1)
        row, col = transform @ row_col
        return row, col

    def get_origin_offset(self):
        """
        Get the translation offset (in pixels) of the current view from (511, 511)
        """
        center = QtC.QPoint(self.viewport().width() // 2, self.viewport().height() // 2)
        center = self.mapToScene(center)
        return center.x(), center.y()

    def get_roll_offset(self):
        transform = self.viewportTransform()
        return np.rad2deg(np.arctan2(transform.m12(), transform.m11()))

    def get_attitude_offset(self, no_roll=False):
        x, y = self.get_origin_offset()
        yag, zag = pixels_to_yagzag(
            CCD_ORIGIN[0] + x, CCD_ORIGIN[1] - y, allow_bad=True
        )
        if no_roll:
            roll = 0
        else:
            roll = self.get_roll_offset()
        q = Quat(equatorial=[yag / 3600, -zag / 3600, roll])
        return q

    def re_center(self):
        scene_rect = self.scene().sceneRect()
        w, h = scene_rect.width(), scene_rect.height()
        new_scene_rect = QtC.QRectF(-w / 2, -h / 2, w, h)
        # print(f'recentering {w}, {h}')
        self.scene().setSceneRect(new_scene_rect)
        transform = self.viewportTransform().rotate(-self.get_roll_offset())
        self.setTransform(transform)


class StarPlot(QtW.QWidget):
    attitude_changed = QtC.pyqtSignal(float, float, float)

    def __init__(self, parent=None):
        super().__init__(parent)

        layout = QtW.QVBoxLayout(self)
        self.setLayout(layout)

        self._origin = [6.08840495576943, 4.92618563916467]

        self.scene = QtW.QGraphicsScene(self)
        self.scene.setSceneRect(-100, -100, 200, 200)

        self.view = StarView(self.scene)
        scale = 1
        self.view.scale(scale, scale)  # I should not need this but...

        self.layout().addWidget(self.view)

        self.stars = None
        # "base attitude" refers to the attitude when the viewport is at the origin and not rotated
        # the actual attitude takes the base attitude and applies a displacement and a rotation
        self._base_attitude = None
        # "current attitude" refers to the attitude taking into account the viewport's position
        self._current_attitude = None
        self._time = None
        self._highlight = None
        self._catalog = None

        self._catalog_items = []
        self._test_stars_items = []

        self.scene.sceneRectChanged.connect(self._radec_changed)
        self.view.roll_changed.connect(self._roll_changed)

    def _radec_changed(self):
        # RA/dec change when the scene rectangle changes, and its given by the rectangle's center
        # the base attitude corresponds to RA/dec at the origin, se we take the displacement
        # of the view offset, apply it from the origin, and get ra/dec for the offset origin
        if self._base_attitude is None:
            return
        x, y = self.view.get_origin_offset()
        yag, zag = pixels_to_yagzag(
            self._origin[0] + x, self._origin[1] - y, allow_bad=True
        )
        ra, dec = yagzag_to_radec(yag, zag, self._base_attitude)
        # print('RA/dec changed', ra, dec)
        self._current_attitude = Quat(equatorial=[ra, dec, self._current_attitude.roll])
        # print(f'Attitude changed. RA: {ra}, dec: {dec}, roll: {roll} ({x}, {y})')
        self.attitude_changed.emit(
            self._current_attitude.ra,
            self._current_attitude.dec,
            self._current_attitude.roll,
        )

    def _roll_changed(self, roll_offset):
        if self._current_attitude is None:
            return
        # roll changes when the viewport is rotated.
        # the view class keeps track of this.
        # print('roll changed', roll_offset)
        self._current_attitude = Quat(
            equatorial=[
                self._current_attitude.ra,
                self._current_attitude.dec,
                self._base_attitude.roll - roll_offset,
            ]
        )
        self.attitude_changed.emit(
            self._current_attitude.ra,
            self._current_attitude.dec,
            self._current_attitude.roll,
        )

    def set_base_attitude(self, q, update=True):
        """
        Sets the base attitude

        The base attitude is the attitude corresponding to the origin of the scene.
        When the base attitude changes, the star positions must be updated. Not doing so will
        leave the display in an inconsistent state. The "update" argument is there as a convenience
        to delay the update in case one wants to call several setters.
        """
        self._base_attitude = Quat(q)
        self._current_attitude = Quat(self._base_attitude)
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

    def clear_test_stars(self):
        for item in self._test_stars_items:
            self.scene.removeItem(item)
        self._test_stars_items = []

    def show_test_stars(
        self, ra_offset=0 / 3600, dec_offset=3000 / 3600, roll_offset=15, N=100, mag=10
    ):
        dq = Quat(equatorial=[ra_offset, dec_offset, roll_offset])
        q = self._base_attitude * dq
        self.show_test_stars_q(q, N, mag)

    def show_test_stars_q(self, q, N=100, mag=10):
        self.clear_test_stars()
        if self._base_attitude is None or self._time is None:
            return
        red_pen = QtG.QPen(QtG.QColor("red"))
        red_brush = QtG.QBrush(QtG.QColor("red"))

        # The test stars will be placed so they appar at the edge of the camera
        # if the attitude is the one given.
        # These are the positions at the edge, in pixel coordinates
        rows, cols = np.array(
            [[-511, i] for i in np.linspace(-511, 511, N)]
            + [[i, 511] for i in np.linspace(-511, 511, N)]
            + [[511, i] for i in np.linspace(511, -511, N)]
            + [[i, -511] for i in np.linspace(511, -511, N)]
        ).T
        yags, zags = pixels_to_yagzag(rows, cols, allow_bad=True)
        # and these are the positions in RA/dec, assuming the given attitude
        ra, dec = yagzag_to_radec(yags, zags, q)
        # now transform the other way, assuming the current attitude
        yags2, zags2 = radec_to_yagzag(ra, dec, self._base_attitude)
        rows2, cols2 = yagzag_to_pixels(yags2, zags2, allow_bad=True)
        # and plot them
        # note that the coordinate system is (row, -col)
        for row, col in zip(rows2, cols2, strict=False):
            s = symsize(mag)
            rect = QtC.QRectF(row - s / 2, -col - s / 2, s, s)
            self._test_stars_items.append(
                self.scene.addEllipse(rect, red_pen, red_brush)
            )

    def show_stars(self):
        self.scene.clear()
        if self._base_attitude is None or self._time is None:
            return
        self.stars = get_stars(self._time, self._base_attitude)
        # Update table to include row/col values corresponding to yag/zag
        self.stars["row"], self.stars["col"] = yagzag_to_pixels(
            self.stars["yang"], self.stars["zang"], allow_bad=True
        )
        black_pen = QtG.QPen()
        black_pen.setWidth(2)
        black_brush = QtG.QBrush(QtG.QColor("black"))
        red_pen = QtG.QPen(QtG.QColor("red"))
        red_brush = QtG.QBrush(QtG.QColor("red"))
        for star in self.stars:
            s = symsize(star["MAG_ACA"])
            # note that the coordinate system is (row, -col)
            rect = QtC.QRectF(star["row"] - s / 2, -star["col"] - s / 2, s, s)
            if self._highlight is not None and star["AGASC_ID"] in self._highlight:
                self.scene.addEllipse(rect, red_pen, red_brush)
            else:
                self.scene.addEllipse(rect, black_pen, black_brush)

        # self.view.centerOn(QtC.QPointF(self._origin[0], self._origin[1]))
        self.view.re_center()

    def clear_catalog(self):
        for item in self._catalog_items:
            self.scene.removeItem(item)
        self._catalog_items = []

    def show_catalog(self):
        self.clear_catalog()
        if self._catalog is not None:
            cat = Table(self._catalog)
            # the catalog was made using self._current_attitude, but the stars are plotted using
            # self._base_attitude
            ra, dec = yagzag_to_radec(cat["yang"], cat["zang"], self._current_attitude)
            yag, zag = radec_to_yagzag(ra, dec, self._base_attitude)
            cat["row"], cat["col"] = yagzag_to_pixels(yag, zag, allow_bad=True)
            gui_stars = cat[(cat["type"] == "GUI") | (cat["type"] == "BOT")]
            acq_stars = cat[(cat["type"] == "ACQ") | (cat["type"] == "BOT")]
            fids = cat[cat["type"] == "FID"]
            mon_wins = cat[cat["type"] == "MON"]

            for gui_star in gui_stars:
                w = 20
                # note that the coordinate system is (row, -col)
                rect = QtC.QRectF(
                    gui_star["row"] - w, -gui_star["col"] - w, w * 2, w * 2
                )
                self._catalog_items.append(
                    self.scene.addEllipse(rect, QtG.QPen(QtG.QColor("green"), 3))
                )
            for acq_star in acq_stars:
                self._catalog_items.append(
                    self.scene.addRect(
                        acq_star["row"] - acq_star["halfw"] / 5,
                        -acq_star["col"] - acq_star["halfw"] / 5,
                        acq_star["halfw"] * 2 / 5,
                        acq_star["halfw"] * 2 / 5,
                        QtG.QPen(QtG.QColor("blue"), 3),
                    )
                )
            for mon_box in mon_wins:
                # starcheck convention was to plot monitor boxes at 2X halfw
                self._catalog_items.append(
                    self.scene.addRect(
                        mon_box["row"] - (mon_box["halfw"] * 2 / 5),
                        -mon_box["col"] - (mon_box["halfw"] * 2 / 5),
                        mon_box["halfw"] * 4 / 5,
                        mon_box["halfw"] * 4 / 5,
                        QtG.QPen(QtG.QColor(255, 165, 0), 3),
                    )
                )
            for fid in fids:
                w = 25
                rect = QtC.QRectF(fid["row"] - w, -fid["col"] - w, w * 2, w * 2)
                self._catalog_items.append(
                    self.scene.addEllipse(rect, QtG.QPen(QtG.QColor("red"), 3))
                )
                self._catalog_items.append(
                    self.scene.addLine(
                        fid["row"] - w,
                        -fid["col"],
                        fid["row"] + w,
                        -fid["col"],
                        QtG.QPen(QtG.QColor("red"), 3),
                    )
                )
                self._catalog_items.append(
                    self.scene.addLine(
                        fid["row"],
                        -fid["col"] - w,
                        fid["row"],
                        -fid["col"] + w,
                        QtG.QPen(QtG.QColor("red"), 3),
                    )
                )
