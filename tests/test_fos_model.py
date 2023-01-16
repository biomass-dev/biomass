import multiprocessing
import os
import shutil
from distutils.dir_util import copy_tree
from typing import Final, Optional

import numpy as np
import pytest

from biomass import OptimizationResults, create_model, optimize, run_analysis, run_simulation
from biomass.estimation import InitialPopulation
from biomass.models import copy_to_current

MODEL_NAME: Final[str] = "Nakakuki_Cell_2010"

copy_to_current(MODEL_NAME)
assert os.path.exists(MODEL_NAME)
model = create_model(MODEL_NAME)


def test_initialization():
    os.mkdir(os.path.join(MODEL_NAME, "out"))
    copy_tree(
        os.path.join("tests", "out"),
        os.path.join(MODEL_NAME, "out"),
    )


def test_run_simulation():
    assert run_simulation(model, viz_type="experiment") is None
    with pytest.warns(RuntimeWarning) as record:
        run_simulation(model, viz_type="original")
    # check that the message matches
    assert record[-1].message.args[0] in [
        "Simulation failed. #original",
        "invalid value encountered in divide",
        "invalid value encountered in true_divide",
    ], f"got {record[-1].message.args[0]}"
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
                        rtol=1e-3,
                    )


def test_param_estim(n_proc: Optional[int] = None):
    if n_proc is None:
        n_proc = max(1, multiprocessing.cpu_count())
    initpop = InitialPopulation(model).generate(n_proc=n_proc)
    optimize(model, x_id=11, optimizer_options={"maxiter": 5, "init": initpop, "workers": -1})
    with open(os.path.join(model.path, "out", "11", "optimization.log")) as f:
        logs = f.readlines()
    assert logs[-1].startswith("differential_evolution step 5: ")
    assert run_simulation(model, viz_type="11") is None
    assert os.path.isdir(os.path.join(model.path, "figure", "simulation", "11"))


def test_cleanup():
    shutil.rmtree(MODEL_NAME)
