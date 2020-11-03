import os

from biomass.models import mapk_cascade
from biomass.exec_model import ExecModel
from biomass.result import OptimizationResults
from biomass import optimize, optimize_continue, run_simulation

import pytest

MODEL_PATH = mapk_cascade.__path__[0]


def test_simulate_successful():
    model = ExecModel(mapk_cascade)
    x = model.pval()
    y0 = model.ival()
    assert model.sim.simulate(x, y0) is None

def test_optimize():
    optimize(mapk_cascade, 1, 3, max_generation=10, overwrite=True)
    for i in range(1, 4):
        with open(MODEL_PATH + f'/out/{i:d}/optimization.log') as f:
            logs = f.readlines()
        assert logs[-1][:14] == 'Generation10: '

    optimize_continue(mapk_cascade, 1, 3, max_generation=20)
    for i in range(1, 4):
        with open(MODEL_PATH + f'/out/{i:d}/optimization.log') as f:
            logs = f.readlines()
        assert logs[-1][:14] == 'Generation20: '


def test_run_simulation():
    run_simulation(mapk_cascade, viz_type='original')
    run_simulation(mapk_cascade, viz_type='average', stdev=True)
    run_simulation(mapk_cascade, viz_type='best', show_all=True)


def test_save_result():
    res = OptimizationResults(mapk_cascade)
    res.get()
    assert os.path.isfile(
        MODEL_PATH + '/optimization_results/optimized_params.csv'
    )
    res.dynamic_assessment(include_original=True)
    assert os.path.isfile(
        MODEL_PATH + '/optimization_results/fitness_assessment.csv'
    )