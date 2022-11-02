# -*- coding: utf-8 -*-
"""Proviede commands and sub-commands for the xsecs modules"""

import click
import pathlib
import numpy as np

from . import base
from ..grids import grids
from ..runcards import runcards
from .. import xsecs

@base.command.group("xsecs")
def subcommand():
    """Commands for generating cross-sections"""


@subcommand.command("runcards")
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
@click.option(
    "-a",
    "--a_value",
    type=int,
    default=1,
    help="""Atomic mass number value. Default: 1""",
)
@click.option(
    "-d",
    "--destination",
    type=click.Path(path_type=pathlib.Path),
    default=pathlib.Path.cwd().absolute() / "theory",
    help="Destination to store the PineAPPL grid (default: $PWD/theory)",
)
def sub_runcards(x_grids, q2_grids, a_value, destination):
    """Generate the Yadism cards for a fixed proton/nucleus A."""
    if x_grids is not None:
        x_grids = eval(x_grids)
    else:
        x_grids = dict(min=1e-5, max=1.0, num=25)

    if q2_grids is not None:
        q2_grids = eval(q2_grids)
    else:
        q2_grids = dict(min=5, max=1e5, num=30)

    runcards.generate_cards(x_grids, q2_grids, a_value, destination)


@subcommand.command("grids")
@click.argument(
    "runcards",
    type=click.Path(exists=True, path_type=pathlib.Path),
)
@click.option(
    "-d",
    "--destination",
    type=click.Path(path_type=pathlib.Path),
    default=pathlib.Path.cwd().absolute(),
    help="Destination to store the PineAPPL grid (default: $PWD/theory)",
)
def sub_grids(runcards, destination):
    """Generate grids with yadism provided with a theory card.

    RUNCARDS is a path to folder (or tar folder) containing the runcards:
    - only one theory card is expected, whose name has to be `theory.yaml`
    - several observable cards might be provided

    The exact name of the observable cards files are mostly ignored but for
    prefix and suffix: it has to start with `obs`, and extension has to be
    `.yaml`.
    The internal `name` key is used for the generated grids.
    """
    grids.main(runcards.absolute(), destination.absolute())


@subcommand.command("generate_lhapdf")
@click.argument("grids", type=click.Path(exists=True, path_type=pathlib.Path))
@click.argument("pdf")
@click.option(
    "--err",
    type=click.Choice(["pdf", "theory"], case_sensitive=False),
    default="theory",
)
@click.option(
    "-d",
    "--destination",
    type=click.Path(path_type=pathlib.Path),
    default=pathlib.Path.cwd().absolute(),
    help="Destination to store the LHAPDF set (default: $PWD)",
)
def sub_predictions(grids, pdf, err, destination):
    """Generate predictions from yadism grids.

    GRIDS is a path to folder (or tar folder) containing the grids, one per
    observable.
    PDF is the pdf to be convoluted with the grids, in order to obtain the
    structure functions predictions.

    """
    xsecs.main(
        grids.absolute(),
        pdf,
        err=err,
        destination=destination.absolute(),
    )
