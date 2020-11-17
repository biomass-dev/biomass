import os
import shutil

from biomass.models import mapk_cascade
from biomass.exec_model import ExecModel
from biomass.result import OptimizationResults
from biomass import optimize, optimize_continue, run_simulation, run_analysis

import pytest

MODEL_PATH = mapk_cascade.__path__[0]

if os.path.isdir(MODEL_PATH + '/figure'):
    shutil.rmtree(MODEL_PATH + '/figure')


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
    model = ExecModel(mapk_cascade)

    run_simulation(mapk_cascade, viz_type='original')
    run_simulation(mapk_cascade, viz_type='average', stdev=True)
    run_simulation(mapk_cascade, viz_type='best', show_all=True)
    for dir in ['/original', '/average', '/best']:
        files = os.listdir(MODEL_PATH + '/figure/simulation' + dir)
        n_pdf = 0
        for file in files:
            if '.pdf' in file:
                n_pdf += 1
        assert n_pdf == len(model.obs) + 1 # multiplot_observables


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


def test_sensitiity_analysis():
    for target in ['parameter', 'initial_condition']:
        run_analysis(
            mapk_cascade, target=target, metric='integral'
        )
        assert os.path.isfile(
            MODEL_PATH
            + '/sensitivity_coefficients/' + target + '/integral/sc.npy'
        )