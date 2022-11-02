# -*- coding: utf-8 -*-
"""
Module that computes the Structure Function predictions from the NN and
dump them as a LHAPDF-like grid.

Some of the functions below are taken and adapted from the eko library:

    https://github.com/NNPDF/eko
"""

import io
import logging
import math
import pathlib
import re
import shutil

import numpy as np
import scipy.special as sp
import yaml

_logger = logging.getLogger(__name__)


def generate_block(xfxQ2, xgrid, Q2grid, pids):
    """
    Generate an LHAPDF data block from a callable

    Parameters
    ----------
        xfxQ2 : callable
            LHAPDF like callable
        Q2grid : list(float)
            Q grid
        pids : list(int)
            Flavors list
        xgrid : list(float)
            x grid

    Returns
    -------
        dict :
            PDF block
    """
    block = dict(Q2grid=Q2grid, pids=pids, xgrid=xgrid)
    data = []
    for x in xgrid:
        for Q2 in Q2grid:
            data.append(np.array([xfxQ2(pid, x, Q2) for pid in pids]))
    block["data"] = np.array(data)
    return block


def create_info_file(sf_flavors, a_value, x_grids, q2_grids, nrep):
    """
    Generate a lhapdf info file from theory and operators card

    Parameters
    ----------
        inputs: type
            description

    Returns
    -------
        : dict
        info file in lhapdf format
    """
    template_info = {}
    template_info["SetDesc"] = f"Structure Function PDFs for A={a_value}"
    template_info["Authors"] = "NvSF"
    template_info["FlavorScheme"] = ""
    template_info["NumFlavors"] = len(sf_flavors) + 1
    template_info["Flavors"] = sf_flavors
    template_info["XMin"] = x_grids[0]
    template_info["XMax"] = x_grids[-1]
    template_info["NumMembers"] = nrep
    template_info["OrderQCD"] = ""
    template_info["QMin"] = round(math.sqrt(q2_grids[0]), 4)
    template_info["QMax"] = round(math.sqrt(q2_grids[-1]), 4)
    template_info["MZ"] = ""
    template_info["MUp"] = 0.0
    template_info["MDown"] = 0.0
    template_info["MStrange"] = 0.0
    template_info["MCharm"] = ""
    template_info["MBottom"] = ""
    template_info["MTop"] = ""
    template_info["AlphaS_MZ"] = 0.118000
    # The following Alphas entries are required by LHAPDF at loading
    template_info["AlphaS_OrderQCD"] = 0
    template_info["AlphaS_Type"] = "ipol"
    template_info["AlphaS_Qs"] = [
        1.295000e00,
        1.475010e00,
        1.691180e00,
        1.952560e00,
        2.270860e00,
        2.661400e00,
        3.144360e00,
        3.746580e00,
        4.504060e00,
        5.465600e00,
        6.697950e00,
        8.293410e00,
        1.038100e01,
        1.314320e01,
        1.684110e01,
        2.185310e01,
        2.873440e01,
        3.831200e01,
        5.183410e01,
        7.121450e01,
        9.943320e01,
        1.412090e02,
        2.041420e02,
        3.007030e02,
        4.517430e02,
        6.928310e02,
        1.085930e03,
        1.741380e03,
        2.860290e03,
        4.818130e03,
        8.334020e03,
        1.482260e04,
        2.714580e04,
        5.126670e04,
        1.000000e05,
    ]
    template_info["AlphaS_Vals"] = [
        4.200710e-01,
        3.895640e-01,
        3.619420e-01,
        3.368440e-01,
        3.139670e-01,
        2.930520e-01,
        2.738810e-01,
        2.562640e-01,
        2.400390e-01,
        2.250660e-01,
        2.112210e-01,
        1.983960e-01,
        1.864960e-01,
        1.754390e-01,
        1.651500e-01,
        1.555630e-01,
        1.466190e-01,
        1.382660e-01,
        1.304560e-01,
        1.231460e-01,
        1.162990e-01,
        1.098790e-01,
        1.038540e-01,
        9.819640e-02,
        9.287940e-02,
        8.787930e-02,
        8.317410e-02,
        7.874380e-02,
        7.457000e-02,
        7.063560e-02,
        6.692520e-02,
        6.342420e-02,
        6.011950e-02,
        5.699860e-02,
        5.405020e-02,
    ]

    return template_info


def list_to_str(ls, fmt="%.6e"):
    """
    Convert a list of numbers to a string

    Parameters
    ----------
        ls : list(float)
            list
        fmt : str
            format string

    Returns
    -------
        str :
            final string
    """
    return " ".join([fmt % x for x in ls])


def array_to_str(ar):
    """
    Convert an array of numbers to a string

    Parameters
    ----------
        ar : list(list(float))
            array

    Returns
    -------
        str :
            final string
    """
    table = ""
    for line in ar:
        table += f"{line[0]:.8e} " + list_to_str(line[1:], fmt="%.8e") + "\n"
    return table


def dump_blocks(name, member, blocks, pdf_type=None):
    """
    Write LHAPDF data file.

    Parameters
    ----------
        name : str or os.PathLike
            target name or path
        member : int
            PDF member
        blocks : list(dict)
            pdf blocks of data
        inherit : str
            str to be copied in the head of member files
    """
    path_name = pathlib.Path(name)
    target = path_name / f"{path_name.stem}_{member:04d}.dat"
    target.parent.mkdir(exist_ok=True)
    with open(target, "w", encoding="utf-8") as o:
        if pdf_type is None:
            if member == 0:
                o.write("PdfType: central\n")
            else:
                o.write("PdfType: replica\n")
        else:
            o.write(pdf_type)
        o.write("Format: lhagrid1\n---\n")
        for b in blocks:
            o.write(list_to_str(b["xgrid"]) + "\n")
            o.write(list_to_str(list(np.sqrt(b["Q2grid"]))) + "\n")
            o.write(list_to_str(b["pids"], "%d") + "\n")
            o.write(array_to_str(b["data"]))
            o.write("---\n")


def dump_info(name, info):
    """
    Write LHAPDF info file.

    NOTE: Since LHAPDF info files are not truely yaml files,
    we have to use a slightly more complicated function to
    dump the info file.

    Parameters
    ----------
        name : str or os.Pathlike
            target name or path
        info : dict
            info dictionary
    """
    path_name = pathlib.Path(name)
    target = path_name / f"{path_name.stem}.info"
    target.parent.mkdir(exist_ok=True)
    # write on string stream to capture output
    stream = io.StringIO()
    yaml.safe_dump(
        info,
        stream,
        sort_keys=False,
        default_flow_style=True,
        width=100000,
        line_break="\n",
    )
    cnt = stream.getvalue()
    # now insert some newlines for each key
    new_cnt = re.sub(r", ([A-Za-z_]+\d?):", r"\n\1:", cnt.strip()[1:-1])
    with open(target, "w", encoding="utf-8") as o:
        o.write(new_cnt)


def dump_set(name, info, member_blocks, pdf_type_list=None):
    """
    Dump a whole set.

    Parameters
    ----------
        name : str
            target name
        info : dict
            info dictionary
        member_blocks : list(list(dict))
            blocks for all members
        pdf_type : list(str)
            list of strings to be copied in the head of member files
    """
    dump_info(name, info)
    for mem, blocks in enumerate(member_blocks):
        if not isinstance(pdf_type_list, list):
            dump_blocks(name, mem, blocks)
        elif len(pdf_type_list) == 0:
            dump_blocks(name, mem, blocks)
        else:
            dump_blocks(name, mem, blocks, pdf_type=pdf_type_list[mem])


def install_pdf(name):
    """
    Install set into LHAPDF.

    The set to be installed has to be in the current directory.

    Parameters
    ----------
        name : str
            source pdf name
    """
    import lhapdf  # pylint: disable=import-error, import-outside-toplevel

    target = pathlib.Path(lhapdf.paths()[0]).joinpath(name)
    src = pathlib.Path(name)
    if not src.exists():
        raise FileExistsError(src)
    if target.exists():
        _logger.warning("Set already exists in LHAPDF and will be overwritten.")
        shutil.rmtree(target, ignore_errors=False)
    shutil.copytree(str(src), str(target))


def make_lambert_grid(n_pts, x_min=1e-7, x_max=1.0):
    r"""
    The following function has been copied from `eko.interpolation`. It creates
    a smoothly spaced grid that is linear near 1 and logarithmic near 0.

    It is generated by the relation:

    .. math::
        x(y) = \frac{1}{5} W_{0}(5 \exp(5-y))

    where :math:`W_{0}` is the principle branch of the
    :func:`Lambert W function <scipy.special.lambertw>` and
    :math:`y` is a variable which extremes are given as function of :math:`x`
    by the direct relation:

    .. math::
        y(x) = 5(1-x)-\log(x)

    This method is implemented in `PineAPPL`, :cite:`Carrazza_2020` eq 2.11 and
    relative paragraph.

    Parameters
    ----------
        n_pts : int
            number of points
        x_min : float
            minimum x (included)
        x_max : float
            maximum x (included)

    Returns
    -------
        xgrid : numpy.ndarray
            generated grid
    """

    def direct_relation(x):
        return 5 * (1 - x) - np.log(x)

    def inverse_relation(y):
        return np.real(1 / 5 * sp.lambertw(5 * np.exp(5 - y)))

    y_min = direct_relation(x_min)
    y_max = direct_relation(x_max)

    return np.array(
        [inverse_relation(y) for y in np.linspace(y_min, y_max, n_pts)]
    )
