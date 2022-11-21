# -*- coding: utf-8 -*-
"""Provide extra scripts subcommand."""
import click

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
    help="""Stringified list containing specs for x-range"""
    """" e.g. '[1e-5, 1]'.""",
)
@click.option(
    "-q",
    "--q2_grids",
    default=None,
    help="""Stringified dictionary containing specs for Q2-range"""
    """" e.g. '[1.65, 1e4]'.""",
)
def sub_integrate(pdfset, a_value, e_nu, type, x_grids, q2_grids):
    """Compute the integral of a Structure Function given a LHAPDF ID.

    nnu extra integrate <pdfset_name>
    """
    if x_grids is not None:
        x_grids = eval(x_grids)
    else:
        x_grids = None

    if q2_grids is not None:
        q2_grids = eval(q2_grids)
    else:
        q2_grids = None

    if type == "neutrino":
        lepton_type = 0
    elif type == "antineutrino":
        lepton_type = 1
    else:
        raise ValueError("Type not determined")

    integrate.main(pdfset, a_value, e_nu, lepton_type, x_grids, q2_grids)
