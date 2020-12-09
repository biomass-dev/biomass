import os
import shutil
from biomass.models import tgfb_smad
from biomass.exec_model import ExecModel
from biomass import run_simulation

import pytest


for dir in ["/figure", "/simulation_data"]:
    if os.path.isdir(tgfb_smad.__path__[0] + dir):
        shutil.rmtree(tgfb_smad.__path__[0] + dir)


def test_simulate_successful():
    model = ExecModel(tgfb_smad)
    x = model.pval()
    y0 = model.ival()
    assert model.sim.simulate(x, y0) is None


def test_run_simulation():
    run_simulation(tgfb_smad, viz_type="original")
    assert os.path.isdir(tgfb_smad.__path__[0] + "/figure/simulation/original")


def test_cleanup():
    for dir in ["/figure", "/simulation_data"]:
        if os.path.isdir(tgfb_smad.__path__[0] + dir):
            shutil.rmtree(tgfb_smad.__path__[0] + dir)