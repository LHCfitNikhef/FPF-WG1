# -*- coding: utf-8 -*-
"""Provide extra scripts subcommand."""
import click
import numpy as np

from .. import integrate
from . import base


@base.command.group("extra")
def subcommand():
    """Extra commands for various scripts."""


@subcommand.command("integrate")
@click.argument("pdfset", type=str)
@click.option(
    "-a",
    "--a_value",
    type=int,
    default=1,
    help="""Value of A for the nucleus. Default: A=1""",
)
@click.option(
    "-e",
    "--e_nu",
    type=float,
    default=5e3,
    help="""Value of the (anti-)neutrino energy. Default: Enu=5e3""",
)
@click.option(
    "--type",
    type=click.Choice(["neutrino", "antineutrino"], case_sensitive=False),
    default="neutrino",
)
@click.option(
    "-x",
    "--x_grids",
    default=None,
    help="""Stringified dictionary containing specs for x-grid"""
    """" e.g. '{"min": 0.01, "max": 1.0, "num": 100}'.""",
)
@click.option(
    "-q",
    "--q2_grids",
    default=None,
    help="""Stringified dictionary containing specs for Q2-grid"""
    """" e.g. '{"min": 0.001, "max": 100000, "num": 200}'.""",
)
def sub_integrate(pdfset, a_value, e_nu, type, x_grids, q2_grids):
    """Compute the integral of a Structure Function given a LHAPDF ID.

    nnu extra integrate <pdfset_name>
    """
    if x_grids is not None:
        x_grids = eval(x_grids)
        xgrid = np.geomspace(x_grids["min"], x_grids["max"], num=x_grids["num"])
    else:
        xgrid = None

    if q2_grids is not None:
        q2_grids = eval(q2_grids)
        q2grid = np.geomspace(
            q2_grids["min"], q2_grids["max"], num=q2_grids["num"]
        )
    else:
        q2grid = None

    if type == "neutrino":
        lepton_type = 0
    elif type == "antineutrino":
        lepton_type = 1
    else:
        raise ValueError("Type not determined")

    integrate.main(pdfset, a_value, e_nu, lepton_type, xgrid, q2grid)
