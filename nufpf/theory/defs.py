# -*- coding: utf-8 -*-
"""Physics constants and definitions."""

import itertools

sfmap = dict(F2="F2_total", FL="FL_total", F3="F3_total")
projectiles = dict(NU="neutrino", NB="antineutrino")
targets = {1: "proton", 20: "neon", 100: "marble", 56: "iron", 208: "lead"}

three_points = [0.5, 1.0, 2.0]

nine_points = list(itertools.product(three_points, three_points))

seven_points = [
    (0.5, 0.5),
    (1.0, 0.5),
    (2.0, 1.0),
    (1.0, 1.0),
    (0.5, 1.0),
    (1.0, 2.0),
    (2.0, 2.0),
]
