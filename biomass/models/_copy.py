"""
`biomass.models` contains sample models and other scripts if you would prefer to learn how to use biomass by example.
For details, please refer to https://biomass-core.readthedocs.io/en/latest/models.html.
"""

import os
import shutil
from typing import Literal


def copy_to_current(
    model_name: Literal[
        "circadian_clock",
        "insulin_signaling",
        "mapk_cascade",
        "Nakakuki_Cell_2010",
        "nfkb_pathway",
        "pan_rtk",
        "tgfb_smad",
        "g1s_transition",
        "prolif_quies",
    ]
) -> None:
    """Copy an example model to the current working directory.

    To execute example models in your machine, please follow these steps:

    >>> from biomass import create_model
    >>> from biomass.models import copy_to_current
    >>> copy_to_current('XXX')
    >>> model = create_model('XXX')

    Parameters
    ----------
    model_name : str
        Name of the example model.
    """
    if model_name not in (
        available_models := [
            "circadian_clock",
            "insulin_signaling",
            "mapk_cascade",
            "Nakakuki_Cell_2010",
            "nfkb_pathway",
            "pan_rtk",
            "tgfb_smad",
            "g1s_transition",
            "prolif_quies",
        ]
    ):
        raise ValueError(f"model_name must be one of [{', '.join(available_models)}].")
    this_directory = os.path.abspath(os.path.dirname(__file__))
    shutil.copytree(
        os.path.join(this_directory, model_name), os.path.join(os.getcwd(), model_name)
    )
