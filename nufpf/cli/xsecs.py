# -*- coding: utf-8 -*-
"""
Provide commands and sub-commands for computing the
differential cross sections using the xsecs.py modules.
"""

import click
import pathlib
import numpy as np
import pandas as pd

from . import base
from ..grids import grids
from ..runcards import runcards
from .. import xsecs

def parse_inputgrid(csv_grid: pd.DataFrame):
    """Parse the input grid knots given read from the CSV."""
    parsed_inputs = []
    igrid = csv_grid[["Log10(x)", "y", "Log10(Q^2 / GeV^2)"]]
    for index, kinematic in enumerate(igrid.values.T):
        if index != 1:
            parsed_inputs.append(10**kinematic)
        else:
            parsed_inputs.append(kinematic)
    return np.asarray(parsed_inputs).T


@base.command.group("xsecs")
def subcommand():
    """Commands for generating cross-sections"""


@subcommand.command("runcards")
@click.argument(
    "input_grids",
    type=click.Path(exists=True, path_type=pathlib.Path),
)
@click.option(
    "-a",
    "--a_value",
    type=int,
    default=1,
    help="""Atomic mass number value. Default: 1""",
)
@click.option(
    "--obs",
    type=click.Choice(["SF", "XSEC"], case_sensitive=False),
    default="XSEC",
    help="Observable type: SF [Struct Func], XSEC [Cross Sec]",
)
@click.option(
    "-d",
    "--destination",
    type=click.Path(path_type=pathlib.Path),
    default=pathlib.Path.cwd().absolute() / "theory",
    help="Destination to store the PineAPPL grid (default: $PWD/theory)",
)
def sub_runcards(input_grids, a_value, obs, destination):
    """Generate the Yadism cards for a fixed proton/nucleus A."""
    igrid = parse_inputgrid(pd.read_csv(input_grids))
    runcards.generate_cards(igrid, a_value, obs, destination)


@subcommand.command("grids")
@click.argument(
    "runcards",
    type=click.Path(exists=True, path_type=pathlib.Path),
)
@click.option(
    "-d",
    "--destination",
    type=click.Path(path_type=pathlib.Path),
    default=pathlib.Path.cwd().absolute() / "theory",
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


@subcommand.command("generate_sfs_lhapdf")
@click.argument("grids", type=click.Path(exists=True, path_type=pathlib.Path))
@click.argument(
    "input_grids",
    type=click.Path(exists=True, path_type=pathlib.Path),
)
@click.argument("pdf")
@click.option(
    "--err",
    type=click.Choice(["pdf", "theory"], case_sensitive=False),
    default="pdf",
)
@click.option(
    "-d",
    "--destination",
    type=click.Path(path_type=pathlib.Path),
    default=pathlib.Path.cwd().absolute(),
    help="Destination to store the LHAPDF set (default: $PWD)",
)
def sub_generate_sfs_lhapdf(grids, input_grids, pdf, err, destination):
    """Compute the Structure Function predictions and dump them into
    LHAPDF grid format.
    """
    input_grids = parse_inputgrid(pd.read_csv(input_grids))
    xsecs.dump_sfs_as_lahpdf_grids(
        grids.absolute(),
        input_grids,
        pdf,
        err=err,
        destination=destination.absolute(),
    )


@subcommand.command("generate_xsecs_datfile")
@click.argument("grids", type=click.Path(exists=True, path_type=pathlib.Path))
@click.argument(
    "input_grids",
    type=click.Path(exists=True, path_type=pathlib.Path),
)
@click.argument("pdf")
@click.option(
    "--err",
    type=click.Choice(["pdf", "theory"], case_sensitive=False),
    default="pdf",
)
@click.option(
    "-d",
    "--destination",
    type=click.Path(path_type=pathlib.Path),
    default=pathlib.Path.cwd().absolute(),
    help="Destination to store the LHAPDF set (default: $PWD)",
)
def sub_generate_xsecs_datfile(grids, input_grids, pdf, err, destination):
    """Compute the differential cross section and dump the results as
    a .dat file.
    """
    input_grids = parse_inputgrid(pd.read_csv(input_grids))
    xsecs.dump_xsecs_as_datfile(
        grids.absolute(),
        input_grids,
        pdf,
        err=err,
        destination=destination.absolute(),
    )


@subcommand.command("multiply_sfs")
@click.argument(
    "input_grids",
    type=click.Path(exists=True, path_type=pathlib.Path),
)
@click.argument("pdf")
@click.option(
    "-d",
    "--destination",
    type=click.Path(path_type=pathlib.Path),
    default=pathlib.Path.cwd().absolute(),
    help="Destination to store the LHAPDF set (default: $PWD)",
)
def sub_multiply_sfs(input_grids, pdf, destination):
    """
    Takes a LHAPDF set of structure functions and multiply with the
    corresponding factors to compute the differential cross sections.
    """
    input_grids = parse_inputgrid(pd.read_csv(input_grids))
    xsecs.multiply_sfs_xsecs(
        input_grids,
        pdf,
        destination=destination.absolute(),
    )
