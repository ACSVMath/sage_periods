r"""A SageMath package for computing D-finite equations of rational diagonals and period integrals"""

import importlib.metadata
# __version__ = importlib.metadata.version(__name__) # TODO: Uncomment this (and delete hardcoded version) after package is uploaded to pip.
__version__ = "0.1.0"

# Checks if ore_algebra is installed
import importlib.util


# Check if ore_algebra is installed
spec = importlib.util.find_spec("ore_algebra")
_is_ore_algebra_installed = spec is not None
# _is_ore_algebra_installed = False # DELETE LATER! For testing purposes.
# We will only use this value in the called modules.


# Functions to be exported for users
from .picard_fuchs import (
    compute_homogenization,
    compute_prepared_fraction,
    compute_period_annihilator,
    compute_diagonal_annihilator
)

__all__ = [
    "compute_homogenization",
    "compute_prepared_fraction",
    "compute_period_annihilator",
    "compute_diagonal_annihilator"
]