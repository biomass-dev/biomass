import os
import shutil
import numpy as np

from biomass.models import circadian_clock
from biomass.exec_model import ExecModel
from biomass import run_simulation


MODEL_PATH = circadian_clock.__path__[0]

for dir in ["/figure", "/simulation_data"]:
    if os.path.isdir(MODEL_PATH + dir):
        shutil.rmtree(MODEL_PATH + dir)


def test_simulate_successful():
    model = ExecModel(circadian_clock)
    x = model.pval()
    y0 = model.ival()
    assert model.sim.simulate(x, y0) is None


def test_run_simulation():
    run_simulation(circadian_clock, viz_type="original")
    simulated_value = np.load(MODEL_PATH + "/simulation_data/simulations_original.npy")
    assert np.isfinite(simulated_value).all()


def test_cleanup():
    for dir in ["/figure", "/simulation_data"]:
        if os.path.isdir(MODEL_PATH + dir):
            shutil.rmtree(MODEL_PATH + dir)