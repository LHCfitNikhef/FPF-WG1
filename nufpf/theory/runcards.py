# -*- coding: utf-8 -*-
"""Generate runcards with x and A fixed."""
import copy

import numpy as np
from yadmark.data.observables import default_card

from . import defs


def observables(x: dict, q2: dict, A: int) -> dict:
    """Collect all yadism runcards."""
    run_nu = copy.deepcopy(default_card)
    run_nu["prDIS"] = "CC"
    run_nu["ProjectileDIS"] = "neutrino"

    # Construct the (x, Q2) knots for the predictions
    q2_grid = np.geomspace(q2["min"], q2["max"], num=int(q2["num"]))
    lognx = int(x["num"] / 3)
    linnx = int(x["num"] - lognx)
    xgrid_log = np.logspace(np.log10(x["min"]), -1, lognx + 1)
    xgrid_lin = np.linspace(0.1, 1, linnx)
    x_grid = np.concatenate([xgrid_log[:-1], xgrid_lin])

    kins = []
    for xv in x_grid:
        for q2v in q2_grid:
            kins.append({"x": float(xv), "Q2": float(q2v), "y": 0.0})

    run_nu["observables"] = {"F2": kins, "F3": kins, "FL": kins}
    run_nu["TargetDIS"] = defs.targets[A]
    run_nb = copy.deepcopy(run_nu)
    run_nb["ProjectileDIS"] = "antineutrino"

    return {f"nu_A_{A}": run_nu, f"nub_A_{A}": run_nb}
