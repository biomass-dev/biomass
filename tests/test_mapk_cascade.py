import os
import shutil

import numpy as np

from biomass import (
    Model,
    OptimizationResults,
    optimize,
    optimize_continue,
    run_analysis,
    run_simulation,
)
from biomass.models import mapk_cascade

model = Model(mapk_cascade.__package__).create()

for dir in [
    "figure",
    "simulation_data",
    "sensitivity_coefficients",
    "optimization_results",
]:
    if os.path.isdir(os.path.join(model.path, dir)):
        shutil.rmtree(os.path.join(model.path, dir))


def test_simulate_successful():
    x = model.pval()
    y0 = model.ival()
    assert model.sim.simulate(x, y0) is None


def test_optimize():
    optimize(
        model=model,
        start=1,
        end=3,
        options={
            "popsize": 3,
            "max_generation": 3,
            "local_search_method": "mutation",
            "n_children": 15,
            "overwrite": True,
        },
    )
    for paramset in range(1, 4):
        with open(
            os.path.join(
                model.path,
                "out",
                f"{paramset:d}",
                "optimization.log",
            )
        ) as f:
            logs = f.readlines()
        assert logs[-1][:13] == "Generation3: "

    optimize_continue(
        model=model,
        start=1,
        end=3,
        options={
            "popsize": 3,
            "max_generation": 6,
            "local_search_method": "Powell",
        },
    )
    for paramset in range(1, 4):
        with open(
            os.path.join(
                model.path,
                "out",
                f"{paramset:d}",
                "optimization.log",
            )
        ) as f:
            logs = f.readlines()
        assert logs[-1][:13] == "Generation6: "

    optimize_continue(
        model=model,
        start=1,
        end=3,
        options={
            "popsize": 3,
            "max_generation": 9,
            "local_search_method": "DE",
            "workers": 1,
        },
    )
    for paramset in range(1, 4):
        with open(
            os.path.join(
                model.path,
                "out",
                f"{paramset:d}",
                "optimization.log",
            )
        ) as f:
            logs = f.readlines()
        assert logs[-1][:13] == "Generation9: "


def test_run_simulation():
    run_simulation(model, viz_type="original")
    run_simulation(model, viz_type="1")
    run_simulation(model, viz_type="2")
    run_simulation(model, viz_type="3")
    run_simulation(model, viz_type="average", stdev=True)
    run_simulation(model, viz_type="best", show_all=True)
    for npy_file in [
        ["simulation_data", "simulations_original.npy"],
        ["simulation_data", "simulations_1.npy"],
        ["simulation_data", "simulations_2.npy"],
        ["simulation_data", "simulations_3.npy"],
        ["simulation_data", "simulations_all.npy"],
        ["simulation_data", "simulations_best.npy"],
    ]:
        simulated_value = np.load(os.path.join(model.path, *npy_file))
        assert np.isfinite(simulated_value).all()
    for viz_type in ["original", "1", "2", "3", "average", "best"]:
        for obs_name in model.obs:
            assert os.path.isfile(
                os.path.join(model.path, "figure", "simulation", f"{viz_type}", f"{obs_name}.pdf")
            )


def test_save_result():
    res = OptimizationResults(model)
    res.to_csv()
    assert os.path.isfile(
        os.path.join(
            model.path,
            "optimization_results",
            "optimized_params.csv",
        )
    )
    res.savefig()
    assert os.path.isfile(
        os.path.join(
            model.path,
            "optimization_results",
            "estimated_parameter_sets.pdf",
        )
    )
    res.dynamic_assessment(include_original=True)
    assert os.path.isfile(
        os.path.join(
            model.path,
            "optimization_results",
            "fitness_assessment.csv",
        )
    )
    res.trace_obj()
    assert os.path.isfile(
        os.path.join(
            model.path,
            "optimization_results",
            "obj_func_traces.pdf",
        )
    )


def test_sensitivity_analysis():
    for target in [
        "parameter",
        "initial_condition",
        "reaction",
    ]:
        for metric in [
            "maximum",
            "minimum",
            "argmax",
            "argmin",
            "timepoint",
            "duration",
            "integral",
        ]:
            run_analysis(model, target=target, metric=metric)
            sensitivity_coefficients = np.load(
                os.path.join(
                    model.path,
                    "sensitivity_coefficients",
                    target,
                    f"{metric}.npy",
                )
            )
            assert np.isfinite(sensitivity_coefficients).all()


def test_cleanup():
    for dir in [
        "figure",
        "simulation_data",
        "out",
        "sensitivity_coefficients",
        "optimization_results",
    ]:
        if os.path.isdir(os.path.join(model.path, dir)):
            shutil.rmtree(os.path.join(model.path, dir))
