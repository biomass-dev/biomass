import os
import shutil
from biomass.models import circadian_clock
from biomass.exec_model import ExecModel
from biomass import run_simulation


for dir in ["/figure", "/simulation_data"]:
    if os.path.isdir(circadian_clock.__path__[0] + dir):
        shutil.rmtree(circadian_clock.__path__[0] + dir)


def test_simulate_successful():
    model = ExecModel(circadian_clock)
    x = model.pval()
    y0 = model.ival()
    assert model.sim.simulate(x, y0) is None


def test_run_simulation():
    run_simulation(circadian_clock, viz_type="original")
    assert os.path.isfile(circadian_clock.__path__[0] + "/simulation_data/simulations_original.npy")


def test_cleanup():
    for dir in ["/figure", "/simulation_data"]:
        if os.path.isdir(circadian_clock.__path__[0] + dir):
            shutil.rmtree(circadian_clock.__path__[0] + dir)