# -*- coding: utf-8 -*-
"""Utilities and runners to generate yadism runcards."""

import logging
import pathlib
import tarfile
import tempfile

import numpy.typing as npt
from .. import theory, utils

_logger = logging.getLogger(__name__)


def generate_cards(igrid: npt.NDArray, A: int, obs: str, destination: pathlib.Path):
    """Generate yadism runcards with the given kinametics"""
    utils.mkdest(destination)
    ocards = theory.runcards.observables(igrid, A, obs)
    theory_card = utils.theory()

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = pathlib.Path(tmpdir)
        utils.write(theory_card, tmpdir / "theory.yaml")
        for name, observable_card in ocards.items():
            utils.write(observable_card, tmpdir / f"obs-{name}.yaml")
            tarpath = destination / f"runcards-{obs.lower()}-a{A}.tar"

            with tarfile.open(tarpath, "w") as tar:
                for tmppath in tmpdir.iterdir():
                    arcname = "runcards/" + tmppath.name
                    tar.add(tmppath.absolute(), arcname=arcname)

            store = tarpath.relative_to(pathlib.Path.cwd())
            _logger.info(f"Runcards have been dumped to '{store}'")
