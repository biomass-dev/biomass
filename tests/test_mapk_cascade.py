import os
import shutil
import numpy as np

from biomass.models import mapk_cascade
from biomass.result import OptimizationResults
from biomass import optimize, optimize_continue, run_simulation, run_analysis


model = mapk_cascade.create()

for dir in [
    "/figure",
    "/simulation_data",
    "/sensitivity_coefficients",
    "/optimization_results",
]:
    if os.path.isdir(model.path + dir):
        shutil.rmtree(model.path + dir)


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
    for i in range(1, 4):
        with open(model.path + f"/out/{i:d}/optimization.log") as f:
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
    for i in range(1, 4):
        with open(model.path + f"/out/{i:d}/optimization.log") as f:
            logs = f.readlines()
        assert logs[-1][:13] == "Generation6: "

    optimize_continue(
        model=model,
        start=1,
        end=3,
        options={"popsize": 3, "max_generation": 9, "local_search_method": "DE", "workers": 1},
    )
    for i in range(1, 4):
        with open(model.path + f"/out/{i:d}/optimization.log") as f:
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
        "/simulation_data/simulations_original.npy",
        "/simulation_data/simulations_1.npy",
        "/simulation_data/simulations_2.npy",
        "/simulation_data/simulations_3.npy",
        "/simulation_data/simulations_all.npy",
        "/simulation_data/simulations_best.npy",
    ]:
        simulated_value = np.load(model.path + npy_file)
        assert np.isfinite(simulated_value).all()
    # assert os.path.isfile(model.path + "/simulation_data/simulations_original.npy")
    # assert os.path.isfile(model.path + "/simulation_data/simulations_all.npy")


def test_save_result():
    res = OptimizationResults(model)
    res.to_csv()
    assert os.path.isfile(model.path + "/optimization_results/optimized_params.csv")
    res.dynamic_assessment(include_original=True)
    assert os.path.isfile(model.path + "/optimization_results/fitness_assessment.csv")


def test_sensitivity_analysis():
    for target in ["parameter", "initial_condition"]:
        run_analysis(model, target=target, metric="integral")
        assert os.path.isfile(model.path + "/sensitivity_coefficients/" + target + "/integral/sc.npy")


def test_cleanup():
    for dir in [
        "/figure",
        "/simulation_data",
        "/out",
        "/sensitivity_coefficients",
        "/optimization_results",
    ]:
        if os.path.isdir(model.path + dir):
            shutil.rmtree(model.path + dir)