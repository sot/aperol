
import numpy as np
from astropy.table import Table

from PyQt5 import QtCore as QtC, QtWidgets as QtW, QtGui as QtG

from Quaternion import Quat
from cxotime import CxoTime

from chandra_aca.transform import yagzag_to_pixels, pixels_to_yagzag, yagzag_to_radec


def symsize(mag):
    # map mags to figsizes, defining
    # mag 6 as 40 and mag 11 as 3
    # interp should leave it at the bounding value outside
    # the range
    return np.interp(mag, [6.0, 11.0], [16.0, 4.0])


def get_stars(starcat_time, quaternion, radius=3):
    import agasc
    from Ska.quatutil import radec2yagzag
    stars = agasc.get_agasc_cone(quaternion.ra, quaternion.dec,
                                 radius=radius,
                                 date=starcat_time)

    if 'yang' not in stars.colnames or 'zang' not in stars.colnames:
        # Add star Y angle and Z angle in arcsec to the stars table.
        # radec2yagzag returns degrees.
        yags, zags = radec2yagzag(stars['RA_PMCORR'], stars['DEC_PMCORR'], quaternion)
        stars['yang'] = yags * 3600
        stars['zang'] = zags * 3600

    return stars


class StarView(QtW.QGraphicsView):
    roll_changed = QtC.pyqtSignal(float)

    def __init__(self, scene=None):
        super().__init__(scene)

        self._start = None
        self._moving = False
        b1hw = 512.
        self.fov = self.scene().addRect(-b1hw, -b1hw, 2 * b1hw, 2 * b1hw)

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
                new_scene_rect = QtC.QRectF(scene_rect.x() - dx, scene_rect.y() - dy,
                                            scene_rect.width(), scene_rect.height())
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
        b1hw = 512.
        center = QtC.QPoint(self.viewport().width() / 2, self.viewport().height() / 2)
        center = self.mapToScene(center)

        transform = self.viewportTransform()
        t11 = transform.m11()
        t12 = transform.m12()
        angle = np.rad2deg(np.arctan2(t12, t11))
        painter.translate(center)
        painter.rotate(-angle)
        painter.translate(-center)

        painter.drawRect(center.x() - b1hw, center.y() - b1hw, 2 * b1hw, 2 * b1hw)
        b2w = 520
        painter.drawRect(center.x() - b2w, center.y() - b1hw, 2 * b2w, 2 * b1hw)

        painter.setPen(QtG.QPen(QtG.QColor('magenta')))
        painter.drawLine(center.x() - 511, center.y(), center.x() + 511, center.y())
        painter.drawLine(center.x(), center.y() - 511, center.x(), center.y() + 511)

    def get_origin_offset(self):
        """
        Get the translation offset (in pixels) of the current view from (511, 511)
        """
        center = QtC.QPoint(self.viewport().width() / 2, self.viewport().height() / 2)
        center = self.mapToScene(center)
        return center.x(), center.y()

    def get_roll_offset(self):
        transform = self.viewportTransform()
        return np.rad2deg(np.arctan2(transform.m12(), transform.m11()))

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

        self.scene.sceneRectChanged.connect(self._radec_changed)
        self.view.roll_changed.connect(self._roll_changed)

    def _radec_changed(self):
        # RA/dec change when the scene rectangle changes, and its given by the rectangle's center
        # the base attitude corresponds to RA/dec at the origin, se we take the displacement
        # of the view offset, apply it from the origin, and get ra/dec for the offset origin
        if self._base_attitude is None:
            return
        x, y = self.view.get_origin_offset()
        yag, zag = pixels_to_yagzag(self._origin[0] + x, self._origin[1] - y, allow_bad=True)
        ra, dec = yagzag_to_radec(yag, zag, self._base_attitude)
        # print('RA/dec changed', ra, dec)
        self._current_attitude = Quat(
            equatorial=[ra, dec, self._current_attitude.roll]
        )
        # print(f'Attitude changed. RA: {ra}, dec: {dec}, roll: {roll} ({x}, {y})')
        self.attitude_changed.emit(
            self._current_attitude.ra, self._current_attitude.dec, self._current_attitude.roll
        )

    def _roll_changed(self, roll_offset):
        if self._current_attitude is None:
            return
        # roll changes when the viewport is rotated.
        # the view class keeps track of this.
        # print('roll changed', roll_offset)
        self._current_attitude = Quat(
            equatorial=[self._current_attitude.ra,
                        self._current_attitude.dec,
                        self._base_attitude.roll - roll_offset]
        )
        self.attitude_changed.emit(
            self._current_attitude.ra, self._current_attitude.dec, self._current_attitude.roll
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
        self.set_base_attitude(catalog.att, update=False)
        self.set_time(catalog.date, update=False)
        self._catalog = catalog
        if update:
            self.show_stars()

    def show_stars(self):
        self.scene.clear()
        if self._base_attitude is None or self._time is None:
            return
        self.stars = get_stars(self._time, self._base_attitude)
        # Update table to include row/col values corresponding to yag/zag
        self.stars['row'], self.stars['col'] = yagzag_to_pixels(
            self.stars['yang'], self.stars['zang'], allow_bad=True
        )
        b1hw = 512.
        self.fov = self.scene.addRect(-b1hw, -b1hw, 2 * b1hw, 2 * b1hw)
        black_pen = QtG.QPen()
        black_pen.setWidth(2)
        black_brush = QtG.QBrush(QtG.QColor("black"))
        red_pen = QtG.QPen(QtG.QColor("red"))
        red_brush = QtG.QBrush(QtG.QColor("red"))
        for star in self.stars:
            s = symsize(star['MAG'])
            rect = QtC.QRectF(star['row'] - s / 2, -star['col'] - s / 2, s, s)
            if self._highlight is not None and star['AGASC_ID'] in self._highlight:
                self.scene.addEllipse(rect, red_pen, red_brush)
            else:
                self.scene.addEllipse(rect, black_pen, black_brush)

        if self._catalog is not None:
            cat = Table(self._catalog)
            cat['row'], cat['col'] = yagzag_to_pixels(cat['yang'], cat['zang'], allow_bad=True)
            gui_stars = cat[(cat['type'] == 'GUI') | (cat['type'] == 'BOT')]
            acq_stars = cat[(cat['type'] == 'ACQ') | (cat['type'] == 'BOT')]
            fids = cat[cat['type'] == 'FID']
            mon_wins = cat[cat['type'] == 'MON']

            for gui_star in gui_stars:
                w = 20
                rect = QtC.QRectF(
                    gui_star['row'] - w,
                    -gui_star['col'] - w,
                    w * 2,
                    w * 2
                )
                self.scene.addEllipse(
                    rect,
                    QtG.QPen(QtG.QColor("green"), 3)
                )
            for acq_star in acq_stars:
                self.scene.addRect(
                    acq_star['row'] - acq_star['halfw'] / 5,
                    -acq_star['col'] - acq_star['halfw'] / 5,
                    acq_star['halfw'] * 2 / 5,
                    acq_star['halfw'] * 2 / 5,
                    QtG.QPen(QtG.QColor("blue"), 3)
                )
            for mon_box in mon_wins:
                # starcheck convention was to plot monitor boxes at 2X halfw
                self.scene.addRect(
                    mon_box['row'] - (mon_box['halfw'] * 2 / 5),
                    -mon_box['col'] - (mon_box['halfw'] * 2 / 5),
                    mon_box['halfw'] * 4 / 5,
                    mon_box['halfw'] * 4 / 5,
                    QtG.QPen(QtG.QColor(255, 165, 0), 3)
                )

            for fid in fids:
                w = 25
                rect = QtC.QRectF(
                    fid['row'] - w,
                    -fid['col'] - w,
                    w * 2,
                    w * 2
                )
                self.scene.addEllipse(
                    rect,
                    QtG.QPen(QtG.QColor("red"), 3)
                )
                self.scene.addLine(
                    fid['row'] - w, -fid['col'],
                    fid['row'] + w, -fid['col'],
                    QtG.QPen(QtG.QColor("red"), 3)
                )
                self.scene.addLine(
                    fid['row'], -fid['col'] - w,
                    fid['row'], -fid['col'] + w,
                    QtG.QPen(QtG.QColor("red"), 3)
                )

        # self.view.centerOn(QtC.QPointF(self._origin[0], self._origin[1]))
        self.view.re_center()
