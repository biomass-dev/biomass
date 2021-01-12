import os
import shutil
import numpy as np

from biomass.models import tgfb_smad
from biomass import run_simulation


model = tgfb_smad.create()


for dir in ["figure", "simulation_data"]:
    if os.path.isdir(os.path.join(model.path, dir)):
        shutil.rmtree(os.path.join(model.path, dir))


def test_simulate_successful():
    x = model.pval()
    y0 = model.ival()
    assert model.sim.simulate(x, y0) is None


def test_run_simulation():
    run_simulation(model, viz_type="original")
    simulated_value = np.load(
        os.path.join(
            model.path,
            "simulation_data",
            "simulations_original.npy",
        )
    )
    assert np.isfinite(simulated_value).all()


def test_cleanup():
    for dir in ["figure", "simulation_data"]:
        if os.path.isdir(os.path.join(model.path, dir)):
            shutil.rmtree(os.path.join(model.path, dir))