"""script to compute the differential cross section given the structure functions as LHAPDF grids"""

import lhapdf
import numpy as np
import pandas as pd


grid = "YADISM_A1"
pdf_member = 0
F = lhapdf.mkPDF(grid, pdf_member)


def xsec(x, y, Q2, l):
    """Given x, y and Q2 as input it returns the differential cross section according to the
    XSHERACC definition of https://yadism.readthedocs.io/en/latest/theory/intro.html.
    Use l=0 for neutrinos and l=1 for antineutrinos"""

    F2 = F.xfxQ2(1001 + 1000 * l, x, Q2)
    FL = F.xfxQ2(1002 + 1000 * l, x, Q2)
    xF3 = F.xfxQ2(1003 + 1000 * l, x, Q2)

    yp = 1 + (1 - y) ** 2
    ym = 1 - (1 - y) ** 2
    yL = y**2

    N = 1.0 / 4.0 * yp

    res = N * (F2 - yL / yp * FL + (-1) ** l * ym / yp * xF3)
    return res


def load_grid(path_to_grid):
    """Load grid containing kinematic points"""
    logx, y, logQ2 = np.loadtxt(
        path_to_grid, delimiter=",", skiprows=1, usecols=(0, 1, 2), unpack=True
    )
    return 10**logx, y, 10**logQ2


if __name__ == "__main__":

    xlist, ylist, Q2list = load_grid("./../igrids/xyQ2.csv")

    res = []
    for x, y, Q2 in zip(xlist, ylist, Q2list):
        res.append([x, y, Q2, xsec(x, y, Q2, 0), xsec(x, y, Q2, 1)])

    df = pd.DataFrame(res, columns=["x", "y", "Q2", "xsec_nu", "xsec_nub"])
    df.to_csv("res_from_lhapdf.csv", index=False)
