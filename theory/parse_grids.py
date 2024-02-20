from os import path
import pathlib
import tarfile
import tempfile
import shutil

from typing import Optional

CURR_PATH = pathlib.Path(__file__).parent
TARGET_DIR = CURR_PATH.joinpath("parsed_grids")
TARGET_DIR.mkdir(exist_ok=True)

# Collect all the grids in the directory
path_grids = CURR_PATH.joinpath("grids")
all_pgrids = path_grids.glob("*grids-FASERv2FCC*.tar")

# Map grid to the correct FK table names
MAP_DATASET_NAMES = {
    "FASERv2FCC_wide_inclusive_nu": 'FASERV2NU_FCC_WIDE_INCLUSIVE',
    "FASERv2FCC_wide_inclusive_nub": 'FASERV2NB_FCC_WIDE_INCLUSIVE',
    "FASERv2FCC_inclusive_nu": 'FASERV2NU_FCC_INCLUSIVE',
    "FASERv2FCC_inclusive_nub": 'FASERV2NB_FCC_INCLUSIVE',
    "FASERv2FCC_deep_inclusive_nu": 'FASERV2NU_FCC_DEEP_INCLUSIVE',
    "FASERv2FCC_deep_inclusive_nub": 'FASERV2NB_FCC_DEEP_INCLUSIVE',
}


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

for grid in all_pgrids:
    raw_gridname = grid.stem
    gridname = raw_gridname.split("-")[1]
    projectile = gridname.split("_")[-1]

    # Name of the subgrid to copy from
    new_gridname = gridname.removesuffix(f"_{projectile}")
    correct_name = f"{projectile}_A_1-{new_gridname}-XSFPFCC.pineappl.lz4"
    new_correct_name = f"{MAP_DATASET_NAMES[gridname]}.pineappl.lz4"

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = pathlib.Path(tmpdir).absolute()

        # Extract all the grids into a temporary directory
        extract_tar(grid, tmpdir, subdirs=1)
        source_dir = f"{tmpdir}/grids/{correct_name}"
        destin_dir = TARGET_DIR.joinpath(new_correct_name)

        if not pathlib.Path(source_dir).is_file():
            raise FileNotFoundError(f"{source_dir} does not exist!")

        shutil.move(src=source_dir, dst=destin_dir)
