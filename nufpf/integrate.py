# -*- coding: utf-8 -*-
"""
Compute the total inclusive cross section given a LHAPDF set
of structure functions.
"""
import lhapdf
import pathlib
import numpy as np

import numpy.typing as npt
from rich.console import Console
from scipy import integrate, interpolate
from typing import Optional

console = Console()

M_NEUTRON = 939.565346 * 0.001
M_PROTON = 938.272013 * 0.001
# Dictionary that maps A value to Z
MAP_AZ = {1: [1, 1], 56: [26, 52.019], 208: [82, 193.729]}


def xs(diffxs: npt.NDArray, xgrid: npt.NDArray, ygrid: npt.NDArray):
    """
    Given a two 1D inputs (x, y) and a 2D resulting functions `diffxs(x,y)`
    compute the integral of the function wrt dx and dy.
    """
    interp_func = interpolate.RectBivariateSpline(xgrid, ygrid, diffxs)

    integral = integrate.dblquad(
        interp_func,
        ygrid.min(),
        ygrid.max(),
        lambda y_: xgrid.min(),
        lambda y_: xgrid.max(),
    )

    return integral


def xsec_as_enu_unc(
    mkpdfs: list,
    a_value: int,
    x: npt.NDArray,
    q2: npt.NDArray,
    enu: float,
    type: int,
) -> npt.NDArray[np.float64]:
    """Compute the differential cross section for a fixed $E_\nu$."""

    def xsec_as_enu(pdf, x, q2, enu, type):
        # Compute the nucleus mass
        A, Z = a_value, MAP_AZ[a_value][0]
        m_nucleus = MAP_AZ[a_value][1] / (Z * M_PROTON + (A - Z) * M_NEUTRON)

        # Compute the value of `y` for a given Enu
        y = q2 / (2. * x * m_nucleus * enu)
        yp, ym = 1 + (1 - y)**2, 1 - (1 - y)**2

        # Extract the individual cross structure functions
        f2  = pdf.xfxQ2(1001 + 1000 * type, x, q2)
        fl  = pdf.xfxQ2(1002 + 1000 * type, x, q2)
        xf3 = pdf.xfxQ2(1003 + 1000 * type, x, q2)

        result = (f2 - (y**2) / yp * fl + (-1) ** type * ym / yp * xf3)
        return (1.0 / 4.0 * yp) * result

    results = np.zeros((len(mkpdfs), x.size, q2.size))
    for p in range(len(mkpdfs)):
        for xi in range(x.size):
            for qi in range(q2.size):
                results[p][xi][qi] += xsec_as_enu(
                    mkpdfs[p],
                    x[xi],
                    q2[qi],
                    enu=enu,
                    type=type,
                )
    return results


def integrate_diffxs(
    pdfname: str,
    a_value: int,
    enu: float,
    type: int,
    x: Optional[npt.NDArray] = None,
    q2: Optional[npt.NDArray] = None,
) -> npt.NDArray:
    mkpdfs = lhapdf.mkPDFs(pdfname)

    # Extract the (x, q2) extremum from the central PDF
    pdf = mkpdfs[0]
    if x is None:
        x = np.geomspace(pdf.xMin, pdf.xMax, 50)
    if q2 is None:
        q2 = np.geomspace(pdf.q2Min, pdf.q2Max, 50)

    results = xsec_as_enu_unc(mkpdfs, a_value, x, q2, enu, type)

    integrated_results = []
    # Loop over the replica results
    for res in results:
        integrated_results.append(xs(res, x, q2))

    return np.asarray(integrated_results)


def main(
    pdfname: str,
    a_value: int,
    enu: float,
    type: int,
    x: Optional[npt.NDArray] = None,
    q2: Optional[npt.NDArray] = None,
) -> None:
    results = integrate_diffxs(
        pdfname,
        a_value,
        enu,
        type,
        x,
        q2,
    )
    console.print(f"mean={results.mean()} +/- {results.std()}", style="bold cyan")
