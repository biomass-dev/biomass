import warnings
from biomass.models import Nakakuki_Cell_2010
from biomass.exec_model import ExecModel

import pytest

def test_simulate_failed():
    warnings.simplefilter('ignore')
    model = ExecModel(Nakakuki_Cell_2010)
    x = model.pval()
    y0 = model.ival()
    assert not model.sim.simulate(x, y0)