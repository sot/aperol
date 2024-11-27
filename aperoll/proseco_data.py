# from PyQt5 import QtCore as QtC, QtWidgets as QtW, QtGui as QtG
import os
import pickle
import tarfile
import traceback
from pathlib import Path
from tempfile import TemporaryDirectory

import PyQt5.QtWidgets as QtW
import sparkles
from proseco import get_aca_catalog
from ska_helpers import utils

from aperoll.utils import logger


class CachedVal:
    def __init__(self, func):
        self._func = func
        self.reset()

    def reset(self):
        self._value = utils.LazyVal(self._func)

    @property
    def val(self):
        return self._value.val


class ProsecoData:
    """
    Class to deal with calling Proseco/Sparkles, temporary directories and exporting the results.

    Parameters
    ----------
    parameters : dict
        The parameters to pass to Proseco. Optional.
    """

    def __init__(self, parameters=None) -> None:
        self._proseco = CachedVal(self.run_proseco)
        self._sparkles = CachedVal(self.run_sparkles)
        self._parameters = parameters
        self._tmp_dir = TemporaryDirectory()
        self._dir = Path(self._tmp_dir.name)

    def reset(self, parameters):
        self._parameters = parameters
        self._proseco.reset()
        self._sparkles.reset()

    @property
    def proseco(self):
        return self._proseco.val

    @property
    def sparkles(self):
        return self._sparkles.val

    def set_parameters(self, parameters):
        self.reset(parameters.copy())

    def get_parameters(self):
        return self._parameters

    parameters = property(get_parameters, set_parameters)

    def export_proseco(self, outfile=None):
        if self.proseco and self.proseco["catalog"]:
            catalog = self.proseco["catalog"]
            if outfile is None:
                outfile = f"aperoll-proseco-obsid_{catalog.obsid:.0f}.pkl"
            with open(outfile, "wb") as fh:
                pickle.dump({catalog.obsid: catalog}, fh)

    def export_sparkles(self, outfile=None):
        if self.sparkles:
            if outfile is None:
                catalog = self.proseco["catalog"]
                outfile = Path(f"aperoll-sparkles-obsid_{catalog.obsid:.0f}.tar.gz")
            dest = Path(str(outfile).replace(".tar", "").replace(".gz", ""))
            with tarfile.open(outfile, "w") as tar:
                for name in self.sparkles.glob("**/*"):
                    tar.add(
                        name,
                        arcname=dest / name.relative_to(self._dir / "sparkles"),
                    )

    def export_proseco_dialog(self):
        """
        Save the star catalog in a pickle file.
        """
        if self.proseco:
            catalog = self.proseco["catalog"]
            dialog = QtW.QFileDialog(
                caption="Export Pickle",
                directory=str(
                    Path(os.getcwd()) / f"aperoll-proseco-obsid_{catalog.obsid:.0f}.pkl"
                ),
            )
            dialog.setAcceptMode(QtW.QFileDialog.AcceptSave)
            dialog.setDefaultSuffix("pkl")
            rc = dialog.exec()
            if rc:
                self.export_proseco(dialog.selectedFiles()[0])

    def export_sparkles_dialog(self):
        """
        Save the sparkles report to a tarball.
        """
        if self.sparkles:
            catalog = self.proseco["catalog"]
            # for some reason, the extension hidden but it works
            dialog = QtW.QFileDialog(
                caption="Export Pickle",
                directory=str(
                    Path(os.getcwd())
                    / f"aperoll-sparkles-obsid_{catalog.obsid:.0f}.tar.gz"
                ),
            )
            dialog.setAcceptMode(QtW.QFileDialog.AcceptSave)
            dialog.setDefaultSuffix(".tgz")
            rc = dialog.exec()
            if rc:
                self.export_sparkles(dialog.selectedFiles()[0])

    def run_proseco(self):
        if self._parameters:
            try:
                params = self._parameters.copy()
                # remove some optional arguments and let proseco deal with it.
                keys = [
                    "exclude_ids_acq",
                    "include_ids_acq",
                    "exclude_ids_guide",
                    "include_ids_guide",
                ]
                for key in keys:
                    if not params[key]:
                        del params[key]
                catalog = get_aca_catalog(**params)
                aca_review = catalog.get_review_table()
                sparkles.core.check_catalog(aca_review)

                return {
                    "catalog": catalog,
                    "review_table": aca_review,
                }
            except Exception as exc:
                logger.debug(f"ProsecoData failed calling proseco ({type(exc).__name__}): {exc}")
                trace = traceback.extract_tb(exc.__traceback__)
                for step in trace:
                    logger.debug(f"    in {step.filename}:{step.lineno}/{step.name}:")
                    logger.debug(f"        {step.line}")
                raise Exception(f"ProsecoData failed calling proseco: {exc}") from None
        return {}

    def run_sparkles(self):
        if self.proseco and self.proseco["catalog"]:
            try:
                sparkles.run_aca_review(
                    "Exploration",
                    acars=[self.proseco["catalog"].get_review_table()],
                    report_dir=self._dir / "sparkles",
                    report_level="all",
                    roll_level="none",
                )
                return self._dir / "sparkles"
            except Exception as exc:
                logger.debug(f"ProsecoData failed calling sparkles ({type(exc).__name__}): {exc}")
                trace = traceback.extract_tb(exc.__traceback__)
                for step in trace:
                    logger.debug(f"    in {step.filename}:{step.lineno}/{step.name}:")
                    logger.debug(f"        {step.line}")
                raise Exception(f"ProsecoData failed calling sparkles: {exc}") from None


    def open_export_proseco_dialog(self):
        """
        Save the star catalog in a pickle file.
        """
        if self.proseco:
            catalog = self.proseco["catalog"]
            dialog = QtW.QFileDialog(
                self,
                "Export Pickle",
                str(self.outdir / f"aperoll-obsid_{catalog.obsid:.0f}.pkl"),
            )
            dialog.setAcceptMode(QtW.QFileDialog.AcceptSave)
            dialog.setDefaultSuffix("pkl")
            rc = dialog.exec()
            if rc:
                self.export_proseco(dialog.selectedFiles()[0])

    def open_export_sparkles_dialog(self):
        """
        Save the sparkles report to a tarball.
        """
        if self.sparkles:
            catalog = self.proseco["catalog"]
            # for some reason, the extension hidden but it works
            dialog = QtW.QFileDialog(
                self,
                "Export Pickle",
                str(self.outdir / f"aperoll-obsid_{catalog.obsid:.0f}.tgz"),
            )
            dialog.setAcceptMode(QtW.QFileDialog.AcceptSave)
            dialog.setDefaultSuffix(".tgz")
            rc = dialog.exec()
            if rc:
                self.export_sparkles(dialog.selectedFiles()[0])
