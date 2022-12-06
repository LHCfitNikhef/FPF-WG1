# -*- coding: utf-8 -*-
"""Generate runcards with x and A fixed."""
import copy

import numpy.typing as npt
from yadmark.data.observables import default_card

from . import defs


def observables(input_grid: npt.NDArray, A: int, obs: str) -> dict:
    """Collect all yadism runcards."""
    run_nu = copy.deepcopy(default_card)
    run_nu["prDIS"] = "CC"
    run_nu["ProjectileDIS"] = "neutrino"

    # Construct the input kinematics as a dictionary
    if obs == "XSEC":
        kins = [
            dict(zip(['x', 'y', 'Q2'], [float(k) for k in kin]))
            for kin in input_grid
        ]
        run_nu["observables"] = {"XSFPFCC": kins}
    elif obs == "SF":
        input_grid[:, 1] = 0
        kins = [
            dict(zip(['x', 'y', 'Q2'], [float(k) for k in kin]))
            for kin in input_grid
        ]
        run_nu["observables"] = {"F2": kins, "F3": kins, "FL": kins}
    else:
        raise ValueError("Observable type non-recognised!")

    run_nu["TargetDIS"] = defs.targets[A]
    run_nb = copy.deepcopy(run_nu)
    run_nb["ProjectileDIS"] = "antineutrino"

    return {f"nu_A_{A}": run_nu, f"nub_A_{A}": run_nb}
