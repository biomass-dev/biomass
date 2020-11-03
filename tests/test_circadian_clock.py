from biomass.models import circadian_clock
from biomass.exec_model import ExecModel
from biomass import run_simulation

import pytest

def test_simulate_successful():
    model = ExecModel(circadian_clock)
    x = model.pval()
    y0 = model.ival()
    assert model.sim.simulate(x, y0) is None

def test_run_simulation():
    run_simulation(circadian_clock, viz_type='original')