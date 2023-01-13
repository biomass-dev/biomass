import os
import shutil
from pathlib import Path


def copy_to_current(model_name: str) -> None:
    """Copy an example model to working directory.

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
    this_directory = os.path.abspath(os.path.dirname(__file__))
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
    shutil.copytree(os.path.join(this_directory, model_name), Path.cwd())
