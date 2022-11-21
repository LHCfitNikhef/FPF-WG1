# -*- coding: utf-8 -*-
"""
Compute the total inclusive cross section given a LHAPDF set
of structure functions.
"""
import lhapdf
import numpy as np

from rich.console import Console
from scipy import integrate
from typing import Optional

console = Console()

M_NEUTRON = 939.565346 * 0.001
M_PROTON = 938.272013 * 0.001
# Dictionary that maps A value to Z
MAP_AZ = {1: [1, 1], 56: [26, 52.019], 208: [82, 193.729]}

# Various constants as defined in https://arxiv.org/pdf/1808.02034.pdf
MW = 80.385 # GeV
G_FERMI = 1.663787e-5 # GeV^-2
INVGEV2_TO_PB = 3.8937966e8 # pb


def xs(
    integrand,
    x_range: list,
    q2_range: list,
    nucleon_mass: float,
    e_nu: float,
):
    """
    Given a two 1D inputs (x, y) and a 2D resulting functions `diffxs(x,y)`
    compute the integral of the function wrt dx and dy.
    """

    threshold = 2.0 * nucleon_mass * e_nu
    integral = integrate.dblquad(
        integrand,
        max(q2_range[0], 1.65), # Q2min
        threshold, # Q2max
        lambda x: max(x / threshold, x_range[0]), #xmin
        1, # xmax
    )

    return integral


def xsec_as_enu_unc(
    mkpdfs,
    a_value: int,
    enu: float,
    type: int,
):
    """Compute the differential cross section for a fixed $E_\nu$."""

    # Compute the nucleus mass
    A, Z = a_value, MAP_AZ[a_value][0]
    if A != 1:
        m_nucleus = MAP_AZ[a_value][1] / (Z * M_PROTON + (A - Z) * M_NEUTRON)
    else:
        m_nucleus = M_PROTON
    console.print(f"Mass of the target Mn={m_nucleus} GeV.")

    def xsec_as_enu(x, q2):

        # Compute the value of `y` for a given Enu
        y = q2 / (2. * x * m_nucleus * enu)
        yp, ym = 1 + (1 - y)**2, 1 - (1 - y)**2

        # Extract the individual cross structure functions
        f2  = mkpdfs.xfxQ2(1001 + 1000 * type, x, q2)
        fl  = mkpdfs.xfxQ2(1002 + 1000 * type, x, q2)
        xf3 = mkpdfs.xfxQ2(1003 + 1000 * type, x, q2)

        result = (yp * f2 - (y**2) * fl + (-1) ** type * ym * xf3)

        # Normalization as defined in https://arxiv.org/pdf/1808.02034.pdf
        prefac_const = (G_FERMI**2 * INVGEV2_TO_PB) / (4.0 * np.pi)
        prefac_funct =  1 / (x * (1 + (q2 / MW**2))**2)
        return prefac_funct * prefac_const * result

    return xsec_as_enu, m_nucleus, enu


def integrate_diffxs(
    pdfname: str,
    a_value: int,
    enu: float,
    type: int,
    x_range: Optional[list] = None,
    q2_range: Optional[list] = None,
):
    mkpdfs = lhapdf.mkPDF(pdfname, 0) # Consider only replica 0

    # Extract the (x, q2) extremum from the central PDF
    if x_range is None:
        x_range = [mkpdfs.xMin, mkpdfs.xMax]
    if q2_range is None:
        q2_range = [mkpdfs.q2Min, mkpdfs.q2Max]

    integrand, m_nuc, enu = xsec_as_enu_unc(mkpdfs, a_value, enu, type)

    return xs(integrand, x_range, q2_range, m_nuc, enu)


def main(
    pdfname: str,
    a_value: int,
    enu: float,
    type: int,
    x_range: Optional[list] = None,
    q2_range: Optional[list] = None,
) -> None:
    results = integrate_diffxs(
        pdfname,
        a_value,
        enu,
        type,
        x_range,
        q2_range,
    )
    console.print(f"mean={results}.", style="bold cyan")
