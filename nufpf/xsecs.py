# -*- coding: utf-8 -*-
"""
Compute the differential cross section given a PineAPPL grid and
dump the results into a .txt file.
"""
import logging
import pathlib
import tarfile
import tempfile
from typing import Optional, Union

import lhapdf
import numpy as np
import numpy.typing as npt
import pandas as pd
import pineappl

from . import utils
from .theory import defs

_logger = logging.getLogger(__name__)

Prediction = tuple[npt.NDArray[np.float_], int, slice, str]
ROUNDING = 6
PDF_MEMBER = 0
GEV_CM2_CONV = 3.893793e10
MW = 80.398
G_FERMI = 1.663787e-5 # GeV^-2
INVGEV2_TO_PB = GEV_CM2_CONV / 100.0  # pb
SFS_LABEL = ["F2nu", "FLnu", "xF3nu", "F2nub", "FLnub", "xF3nub"]


def theory_error(
    grid: pineappl.grid.Grid,
    pdf: str,
    prescription: list[tuple[float, float]],
    xgrid: Optional[Union[npt.NDArray[np.float_], None]] = None,
    reshape: Optional[bool] = False,
) -> Prediction:
    """Compute theory uncertainties using scale variation methods."""
    pdfset = lhapdf.mkPDF(pdf)
    pred = grid.convolute_with_one(
        2212, pdfset.xfxQ2, pdfset.alphasQ2, xi=prescription
    )
    if reshape and xgrid is not None:
        pred = np.array(pred).T.reshape((*xgrid.shape, len(pred)))
    else:
        pred = np.array(pred).T
    return pred, 4, slice(0, -1), "9 pts."


def pdf_error(
    grid: pineappl.grid.Grid,
    pdf: str,
    xgrid: Optional[Union[npt.NDArray[np.float_], None]] = None,
    reshape: Optional[bool] = False,
) -> Prediction:
    """Compute PDF uncertainties for a given observables."""
    pred = []
    for pdfset in lhapdf.mkPDFs(pdf):
        member_pred = grid.convolute_with_one(
            2212, pdfset.xfxQ2, pdfset.alphasQ2
        )
        pred.append(member_pred)

    if reshape and xgrid is not None:
        pred = np.array(pred).T.reshape((*xgrid.shape, len(pred)))
        # Make sure to not include Member 0
        pred = pred[:, :, 1:]
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

    # Compute the central replicas and prepend to array
    central = np.expand_dims(np.mean(predictions, axis=0), axis=0)
    predictions = np.concatenate([central, predictions], axis=0)

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


def dump_sfs_as_lahpdf_grids(
    grids: pathlib.Path,
    input_grids: npt.NDArray,
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

        x_grid = np.unique(input_grids[:, 0])
        q2_grid = np.unique(input_grids[:, -1])
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
                    reshape=True,
                )
            elif err == "pdf":
                pred, central, bulk, err_source = pdf_error(
                    grid,
                    pdf,
                    xgrid,
                    reshape=True,
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


def dump_xsecs_as_datfile(
    grids: pathlib.Path,
    input_grids: npt.NDArray,
    pdf: str,
    destination: pathlib.Path,
    file_name: str,
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

        output_name = grids.stem.replace("grids", "diffxsec")
        output_name += f"_{str(pdf)}"
        # extract tar content
        if grids.suffix == ".tar":
            utils.extract_tar(grids, tmpdir)
            grids = tmpdir / "grids"

        preds_dest = tmpdir / "predictions"
        preds_dest.mkdir()

        predictions_dict = {}
        for gpath in grids.iterdir():
            if "pineappl" not in gpath.name:
                continue
            kind = gpath.stem.split(".")[0].split("_")[0] # nu or nub
            a_nb = gpath.stem.split("-")[0].split("_")[-1] # A value

            grid = pineappl.grid.Grid.read(gpath)

            if err == "theory":
                pred, central, bulk, err_source = theory_error(
                    grid,
                    pdf,
                    defs.nine_points,
                )
            elif err == "pdf":
                pred, central, bulk, err_source = pdf_error(
                    grid,
                    pdf,
                )
            else:
                raise ValueError(f"Invalid error type '{err}'")
            predictions_dict[kind] = np.expand_dims(pred.T[0], axis=-1)

        combined_predictions = np.concatenate(
            [input_grids, predictions_dict["nu"], predictions_dict["nub"]],
            axis=-1,
        )
        np.savetxt(
            f"{destination}/{output_name}.txt",
            combined_predictions,
            header=f"x y Q2 xsec_nu xsec_nub",
            fmt="%e %e %e %e %e",
        )
        _logger.info(
            f"The .txt file containing the xsec predictions is stored in "
            f"'{destination.relative_to(pathlib.Path.cwd())}/{output_name}.txt'"
        )


def multiply_sfs_xsecs(
    input_grids: npt.NDArray,
    pdf: str,
    destination: pathlib.Path,
):
    mkpdf = lhapdf.mkPDF(pdf, PDF_MEMBER)

    def xsec(x, y, Q2, type):
        """Given x, y and Q2 as input it returns the differential
        cross section according to the XSHERACC definition of
        https://yadism.readthedocs.io/en/latest/theory/intro.html.
        Use l=0 for neutrinos and l=1 for antineutrinos.
        """

        yp, ym = 1 + (1 - y)**2, 1 - (1 - y)**2
        # Extract the individual cross structure functions
        f2  = mkpdf.xfxQ2(1001 + 1000 * type, x, Q2)
        fl  = mkpdf.xfxQ2(1002 + 1000 * type, x, Q2)
        xf3 = mkpdf.xfxQ2(1003 + 1000 * type, x, Q2)

        result = (yp * f2 - (y**2) * fl + (-1) ** type * ym * xf3)

        # Normalization as defined in https://arxiv.org/pdf/1808.02034.pdf
        prefac_const = (G_FERMI**2 * INVGEV2_TO_PB) / (4.0 * np.pi)
        # TODO: Check where the factor of 2 is coming from
        prefac_funct =  1 / (2 * x * (1 + (Q2 / MW**2))**2) / 2
        return prefac_funct * prefac_const * result

    results = []
    for kins in input_grids:
        results.append([*kins, xsec(*kins, 0), xsec(*kins, 1)])

    df = pd.DataFrame(results, columns=["x", "y", "Q2", "xsec_nu", "xsec_nub"])
    # df.to_csv(f"{destination}/res_from_lhapdf.csv", index=False)
    combined_predictions = np.concatenate(
            [
                input_grids,
                np.expand_dims(df["xsec_nu"].values, axis=-1),
                np.expand_dims(df["xsec_nub"].values, axis=-1)
            ],
            axis=-1,
        )
    np.savetxt(
        f"{destination}/diff_xsecs_sfs.txt",
        combined_predictions,
        header=f"x y Q2 xsec_nu xsec_nub",
        fmt="%e %e %e %e %e",
    )
    _logger.info(
        f"The .txt file containing the xsec predictions is stored in "
        f"'{destination.relative_to(pathlib.Path.cwd())}/diff_xsecs_sfs.txt'"
    )
