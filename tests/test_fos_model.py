import os
import shutil
from distutils.dir_util import copy_tree

import numpy as np
import pytest
from scipy.optimize import OptimizeResult, differential_evolution

from biomass import Model, OptimizationResults, optimize, run_analysis, run_simulation
from biomass.estimation import ExternalOptimizer
from biomass.models import Nakakuki_Cell_2010

# from biomass import run_analysis

model = Model(Nakakuki_Cell_2010.__package__).create()


def test_initialization():
    for dir in ["figure", "out", "simulation_data", "sensitivity_coefficients"]:
        if os.path.isdir(os.path.join(model.path, dir)):
            shutil.rmtree(os.path.join(model.path, dir))
    os.mkdir(os.path.join("biomass", "models", "Nakakuki_Cell_2010", "out"))
    copy_tree(
        os.path.join("tests", "out"),
        os.path.join("biomass", "models", "Nakakuki_Cell_2010", "out"),
    )


def test_run_simulation():
    with pytest.warns(RuntimeWarning) as record:
        run_simulation(model, viz_type="original")
    # check that the message matches
    assert record[-1].message.args[0] in [
        "Simulation failed. #original",
        "invalid value encountered in true_divide",
    ]
    run_simulation(model, viz_type="average", stdev=True)
    assert os.path.isfile(
        os.path.join(
            model.path,
            "simulation_data",
            "simulations_all.npy",
        )
    )
    for obs_name in model.observables:
        assert os.path.isfile(
            os.path.join(model.path, "figure", "simulation", "average", f"{obs_name}.png")
        )


def test_save_resuts():
    res = OptimizationResults(model)
    res.savefig(figsize=(16, 5), boxplot_kws={"orient": "v"})
    assert os.path.isfile(
        os.path.join(
            model.path,
            "optimization_results",
            "estimated_parameter_sets.pdf",
        )
    )


def test_run_analysis():
    target = "reaction"
    metric = "integral"
    d_ln_vi = 0.009950330853168092
    run_analysis(model, target=target, metric=metric, options={"overwrite": True})
    layer = os.path.join("sensitivity_coefficients", f"{target}", f"{metric}.npy")
    loc_actual = os.path.join(model.path, layer)
    assert os.path.isfile(loc_actual)
    loc_desired = os.path.join(os.path.dirname(__file__), layer)
    actual = np.load(loc_actual)
    desired = np.load(loc_desired)
    assert np.shape(actual) == np.shape(desired)
    for i in range(actual.shape[0]):
        for j in range(actual.shape[1]):
            for k, obs_name in enumerate(model.observables):
                if obs_name != "Phosphorylated_MEKc":
                    np.allclose(
                        np.exp(actual[i, j, k] * d_ln_vi),
                        np.exp(desired[i, j, k] * d_ln_vi),
                        rtol=1e-3
                    )


def test_param_estim():
    optimize(
        model,
        x_id=11,
        options={
            "popsize": 3,
            "max_generation": 5,
            "allowable_error": 0.0,
            "local_search_method": "DE",
            "maxiter": 3,
            "workers": -1,
        },
    )
    with open(os.path.join(model.path, "out", "11", "optimization.log")) as f:
        logs = f.readlines()
    assert logs[-1][:13] == "Generation5: "


def test_external_optimizer():
    optimizer = ExternalOptimizer(model, differential_evolution)
    res = optimizer.run(
        model.problem.objective,
        model.problem.bounds,
        strategy="best2bin",
        maxiter=5,
        popsize=3,
        tol=1e-4,
        mutation=0.1,
        recombination=0.5,
        polish=False,
        workers=-1,
    )
    assert isinstance(res, OptimizeResult)
    optimizer.import_solution(res.x, x_id=0)
    n_iter: int = 0
    with open(
        os.path.join(model.path, "out", "0", "optimization.log"),
        mode="r",
        encoding="utf-8",
    ) as f:
        log_file = f.readlines()
    for message in log_file:
        if len(message.strip()) > 0:
            n_iter += 1
    for fname in [
        "count_num.npy",
        "generation.npy",
        f"fit_param{n_iter}.npy",
    ]:
        assert os.path.isfile(os.path.join(model.path, "out", "0", f"{fname}"))
    assert run_simulation(model, viz_type="0") is None
    assert os.path.isdir(os.path.join(model.path, "figure", "simulation", "0"))


def test_cleanup():
    for dir in ["figure", "simulation_data", "out", "sensitivity_coefficients"]:
        if os.path.isdir(os.path.join(model.path, dir)):
            shutil.rmtree(os.path.join(model.path, dir))
