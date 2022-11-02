# -*- coding: utf-8 -*-
"""Provide various utilities to various modules."""

import logging
import pathlib
import yaml
import tarfile

from rich.progress import track
from typing import Optional
from .export_lhapdf.utils import (
    create_info_file,
    dump_set,
    generate_block,
)

_logger = logging.getLogger(__name__)

pkg = pathlib.Path(__file__).parent.absolute()

LHAPDF_ID = [1001, 1002, 1003, 2001, 2002, 2003, 3001, 3002, 3003]


def mkdest(destination: pathlib.Path):
    """Make sure destination exists.

    Create it if does not exist, else make sure it is a directory.

    Parameters
    ----------
    destination: pathlib.Path
        path to check

    Raises
    ------
    NotADirectoryError
        if it exists but it is not a directory

    """
    if destination.exists():
        if not destination.is_dir():
            raise NotADirectoryError(
                f"The given destination exists, but is not a"
                f" directory - '{destination}'"
            )
    else:
        destination.mkdir(parents=True)


def write(content: dict, path: pathlib.Path, what=None):
    """Write a file, the suitable way."""
    # autodetect
    if what is None:
        what = path.suffix[1:]

    if what == "yaml":
        path.write_text(yaml.dump(content), encoding="utf-8")
    else:
        raise ValueError(f"Format to be read undetected {what}")


def theory(update: Optional[dict] = None) -> dict:
    """Load and return internal theory runcard."""
    assets = pkg / "theory" / "assets" / "theory_200.yaml" 
    runcard = yaml.safe_load(assets.read_text(encoding="utf-8"))

    if update is not None:
        runcard.update(update)
        _logger.info(f"Base theory updated with {update}")

    return runcard


def extract_tar(
    path: pathlib.Path, dest: pathlib.Path, subdirs: Optional[int] = None
):
    """Extract a tar archive to given directory."""
    with tarfile.open(path) as tar:
        tar.extractall(dest)

    if subdirs is not None:
        count = len(list(dest.iterdir()))
        if count == subdirs:
            return

        expected = (
            f"{subdirs} folders are" if subdirs > 1 else "A single folder is"
        )
        found = f"{count} files" if count > 1 else "a single file"
        raise ValueError(
            f"{expected} supposed to be contained by the tar file,"
            f" but more {found} have been detected"
        )


def read(path: pathlib.Path, what=None) -> dict:
    """Read a file, the suitable way."""
    # autodetect
    if what is None:
        what = path.suffix[1:]

    if what == "yaml":
        return yaml.safe_load(path.read_text(encoding="utf-8"))
    else:
        raise ValueError(
            f"Format to be read undetected (attempted for '{what}')"
        )


def dump_lhapdf(
    name: str, am: int, all_replicas: list, grids_info_specs: dict
):
    # Generate the dictionary containing the info file
    info_file = create_info_file(
        sf_flavors=LHAPDF_ID,
        a_value=int(am),
        x_grids=grids_info_specs["x_grids"],
        q2_grids=grids_info_specs["q2_grids"],
        nrep=grids_info_specs["nrep"],
    )

    all_blocks = []
    xgrid = grids_info_specs["x_grids"]
    for pred in track(all_replicas, description="Looping over Replicas:"):
        all_singular_blocks = []
        block = generate_block(
            lambda pid, x, q2, pred=pred: pred[q2][pid][xgrid.index(x)],
            xgrid=grids_info_specs["x_grids"],
            Q2grid=grids_info_specs["q2_grids"],
            pids=LHAPDF_ID,
        )
        all_singular_blocks.append(block)
        all_blocks.append(all_singular_blocks)

    dump_set(name, info_file, all_blocks)
