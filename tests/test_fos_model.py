import os
import shutil
from distutils.dir_util import copy_tree
import warnings
from biomass import run_simulation, run_analysis

import pytest


os.makedirs("biomass/models/Nakakuki_Cell_2010/out", exist_ok=True)
copy_tree("tests/out", "biomass/models/Nakakuki_Cell_2010/out")
from biomass.models import Nakakuki_Cell_2010

MODEL_PATH = Nakakuki_Cell_2010.__path__[0]

for dir in ["/figure", "/simulation_data", "/sensitivity_coefficients"]:
    if os.path.isdir(MODEL_PATH + dir):
        shutil.rmtree(MODEL_PATH + dir)


def test_run_simulation():
    run_simulation(Nakakuki_Cell_2010, viz_type="average", stdev=True)
    assert os.path.isfile(MODEL_PATH + "/simulation_data/simulations_all.npy")


"""
def test_sensitivity_analysis():
    run_analysis(
        Nakakuki_Cell_2010, target='initial_condition', metric='integral'
    )
    assert os.path.isfile(
        MODEL_PATH
        + '/sensitivity_coefficients/initial_condition/integral/sc.npy'
    )
"""


def test_cleanup():
    for dir in ["/figure", "/simulation_data", "/out", "/sensitivity_coefficients"]:
        if os.path.isdir(MODEL_PATH + dir):
            shutil.rmtree(MODEL_PATH + dir)