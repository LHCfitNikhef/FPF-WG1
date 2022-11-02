# -*- coding: utf-8 -*-
"""Compute the Structure Function predictions and dump them into LHAPDF."""
import logging
import pathlib
import tarfile
import tempfile
from typing import Optional

import lhapdf
import numpy as np
import numpy.typing as npt
import pineappl

from . import utils
from .theory import defs

_logger = logging.getLogger(__name__)

Prediction = tuple[npt.NDArray[np.float_], int, slice, str]
ROUNDING = 6
# The labels and IDs in utils.LHAPDF_ID have to match exactly
SFS_LABEL = ["F2nu", "FLnu", "xF3nu", "F2nub", "FLnub", "xF3nub"]


def theory_error(
    grid: pineappl.grid.Grid,
    pdf: str,
    prescription: list[tuple[float, float]],
    xgrid: npt.NDArray[np.float_],
    reshape: Optional[bool] = True,
) -> Prediction:
    """Compute theory uncertainties using scale variation methods."""
    pdfset = lhapdf.mkPDF(pdf)
    pred = grid.convolute_with_one(
        2212, pdfset.xfxQ2, pdfset.alphasQ2, xi=prescription
    )
    if reshape:
        pred = np.array(pred).T.reshape((*xgrid.shape, len(pred)))
    else:
        pred = np.array(pred).T
    return pred, 4, slice(0, -1), "9 pts."


def pdf_error(
    grid: pineappl.grid.Grid,
    pdf: str,
    xgrid: npt.NDArray[np.float_],
    reshape: Optional[bool] = True,
) -> Prediction:
    """Compute PDF uncertainties for a given observables."""
    pred = []
    for pdfset in lhapdf.mkPDFs(pdf):
        member_pred = grid.convolute_with_one(
            2212, pdfset.xfxQ2, pdfset.alphasQ2
        )
        pred.append(member_pred)

    if reshape:
        pred = np.array(pred).T.reshape((*xgrid.shape, len(pred)))
    else:
        pred = np.array(pred).T
    return pred, 0, slice(1, -1), "PDF replicas"


def construct_stacked_predictions(predictions_dict: dict):
    """Stack the predictions in the same order as the LHAPDF IDs."""
    stacked_list = [predictions_dict[k] for k in SFS_LABEL]
    # Stack the predictions such that the shape now becomes
    # (n_x, n_q2, n_rep, n_sfs)
    stacked_predictions = np.stack(stacked_list, axis=-1)
    # Move axis to get shape (n_rep, n_q2, n_x, n_sfs)
    return np.moveaxis(stacked_predictions, [0, 2], [2, 0])


def parse_yadism_pred(
    predictions: np.ndarray,
    x_specs: np.ndarray,
    q2_dic_specs: np.ndarray,
):
    prediction_infoq2 = [round(q, ROUNDING) for q in q2_dic_specs]

    # Append the average to the array block
    copied_pred = np.copy(predictions)
    for i in range(copied_pred.shape[-1] // 2):
        avg = (copied_pred[:, :, :, i] + copied_pred[:, :, :, i + 3]) / 2
        average = np.expand_dims(avg, axis=-1)
        predictions = np.concatenate([predictions, average], axis=-1)

    # Parse the array blocks as a Dictionary
    combined_replica = []
    for replica in predictions:  # loop over the replica
        q2rep_dic = {}
        # Loop over the results for all Q2 values
        for idq, q2rep in enumerate(replica):
            sfs_q2rep = {}
            q2rep_idx = prediction_infoq2[idq]
            # loop over the Structure Functions
            for idx in range(q2rep.shape[-1]):
                sfs_q2rep[utils.LHAPDF_ID[idx]] = q2rep[:, idx]
            q2rep_dic[round(q2rep_idx, ROUNDING)] = sfs_q2rep
        combined_replica.append(q2rep_dic)

    grids_info_specs = {
        "x_grids": x_specs.tolist(),
        "q2_grids": prediction_infoq2,
        "nrep": len(combined_replica),
    }
    return grids_info_specs, combined_replica


def main(
    grids: pathlib.Path,
    pdf: str,
    destination: pathlib.Path,
    err: str = "pdf",
):
    """Run predictions computation.

    Parameters
    ----------
    grids: pathlib.Path
        path to grids archive
    pdf: str
        LHAPDF name of the PDF to be used
    err: str
        type of error to be used
    """
    utils.mkdest(destination)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = pathlib.Path(tmpdir).absolute()

        # extract tar content
        if grids.suffix == ".tar":
            utils.extract_tar(grids, tmpdir)
            grids = tmpdir / "grids"

        preds_dest = tmpdir / "predictions"
        preds_dest.mkdir()

        # This has to be exactly the same as in the runcard generation
        # TODO: Read the details about the kinematics from PineAPPL grid
        x = dict(min=1e-5, max=1.0, num=25)
        q2 = dict(min=5, max=1e5, num=30)
        q2_grid = np.geomspace(q2["min"], q2["max"], num=int(q2["num"]))
        lognx = int(x["num"] / 3)
        linnx = int(x["num"] - lognx)
        xgrid_log = np.logspace(np.log10(x["min"]), -1, lognx + 1)
        xgrid_lin = np.linspace(0.1, 1, linnx)
        x_grid = np.concatenate([xgrid_log[:-1], xgrid_lin])
        xgrid, _ = np.meshgrid(*tuple((q2_grid, x_grid)))

        predictions_dictionary = {}
        nucleus_predictions = []
        for gpath in grids.iterdir():
            if "pineappl" not in gpath.name:
                continue
            obs = gpath.stem.split(".")[0]
            kind = obs.split("_")[0]

            grid = pineappl.grid.Grid.read(gpath)

            if err == "theory":
                pred, central, bulk, err_source = theory_error(
                    grid,
                    pdf,
                    defs.nine_points,
                    xgrid,
                )
            elif err == "pdf":
                pred, central, bulk, err_source = pdf_error(
                    grid,
                    pdf,
                    xgrid,
                )
            else:
                raise ValueError(f"Invalid error type '{err}'")

            gridname = gpath.stem.split(".")[0]
            nuc = gridname.split("-")[0].split("_")[-1]
            sf_type = gridname.split("-")[-1]
            neutrino_type = gridname.split("_")[0]
            if sf_type == "F3":
                obsname = "x" + sf_type + neutrino_type
            else:
                obsname = sf_type + neutrino_type
            nucleus_predictions.append(nuc)
            predictions_dictionary[obsname] = pred

        tardest = destination / "predictions.tar"
        with tarfile.open(tardest, "w") as tar:
            for path in preds_dest.iterdir():
                tar.add(path.absolute(), path.relative_to(tmpdir))

        _logger.info(f"Preedictions saved in '{tardest}'.")

        nuc_info = [i for i in list(set(nucleus_predictions))]
        assert len(nuc_info) == 1
        stacked_pred = construct_stacked_predictions(predictions_dictionary)
        grid_info, cpred = parse_yadism_pred(stacked_pred, x_grid, q2_grid)
        utils.dump_lhapdf(f"YADISM_{nuc_info[0]}", nuc_info[0], cpred, grid_info)
