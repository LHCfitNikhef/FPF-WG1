# -*- coding: utf-8 -*-
"""Compute DIS grids out of given runcards."""
import logging
import pathlib
import tarfile
import tempfile

import yadism
from yadbox.export import dump_pineappl_to_file

from .. import utils

_logger = logging.getLogger(__name__)


def main(cards: pathlib.Path, destination: pathlib.Path):
    """Run grids computation.

    Parameters
    ----------
    cards: pathlib.Path
        path to runcards archive
    destination: pathlib.Path

    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = pathlib.Path(tmpdir).absolute()
        a_value = cards.stem.split("-")[-1]
        obs_name = cards.stem.split("-")[-2]

        # extract tar content
        if cards.suffix == ".tar":
            utils.extract_tar(cards, tmpdir, subdirs=1)
            cards = tmpdir / "runcards"

        grids_dest = tmpdir / "grids"
        grids_dest.mkdir()

        theory = utils.read(cards / "theory.yaml")
        for cpath in cards.iterdir():
            _logger.info(f"Attempt computing '{cpath.name}'")
            if cpath.name.startswith("obs"):
                observables = utils.read(cpath, what="yaml")
                output = yadism.run_yadism(theory, observables)

                data_name = cpath.name.split("-")[1][:-5]
                for obs in observables["observables"]:
                    file_name = (
                        f"{data_name}-{obs}.pineappl.lz4"
                        if data_name != ""
                        else f"{obs}.pineappl.lz4"
                    )
                    res_path = grids_dest / file_name
                    dump_pineappl_to_file(output, res_path, obs)
                    _logger.info(f"Dumped {res_path.name}")

        file_name = (
            f"grids-{obs_name}-{a_value}.tar"
            if a_value != ""
            else f"grids.tar"
        )
        with tarfile.open(destination / file_name, "w") as tar:
            for path in grids_dest.iterdir():
                tar.add(path.absolute(), path.relative_to(tmpdir))
