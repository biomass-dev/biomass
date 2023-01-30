import os
import shutil

import numpy as np

from biomass import OptimizationResults, create_model, optimize, run_analysis, run_simulation
from biomass.models import copy_to_current

MODEL_NAME: str = "mapk_cascade"

copy_to_current(MODEL_NAME)
assert os.path.exists(MODEL_NAME)
model = create_model(MODEL_NAME)


def test_simulate_successful():
    x = model.pval()
    y0 = model.ival()
    assert model.problem.simulate(x, y0) is None


def test_optimize():
    for x_id in range(1, 4):
        optimize(model, x_id=x_id, optimizer_options={"maxiter": 5, "workers": -1})
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
        assert logs[-1].startswith("differential_evolution step 5: ")


def test_run_simulation():
    run_simulation(model, viz_type="original")
    run_simulation(model, viz_type="1")
    run_simulation(model, viz_type="2")
    run_simulation(model, viz_type="3")
    run_simulation(model, viz_type="average", stdev=True)
    run_simulation(model, viz_type="best", show_all=True)
    for npy_file in [
        "simulations_original.npy",
        "simulations_1.npy",
        "simulations_2.npy",
        "simulations_3.npy",
        "simulations_all.npy",
        "simulations_best.npy",
    ]:
        simulated_value = np.load(os.path.join(model.path, "simulation_data", npy_file))
        assert np.isfinite(simulated_value).all()
    for viz_type in ["original", "1", "2", "3", "average", "best"]:
        for obs_name in model.observables:
            assert os.path.isfile(
                os.path.join(model.path, "figure", "simulation", f"{viz_type}", f"{obs_name}.png")
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


def _get_duration(
    time_course: np.ndarray,
    below_threshold: float = 0.5,
) -> int:
    """
    Calculation of the duration as the time it takes to decline below the threshold.

    Parameters
    ----------
    timecourse : array
        Simulated time course data.

    below_threshold : float (from 0.0 to 1.0)
        0.5 for 50% of its maximum.

    Returns
    -------
    duration : int

    """
    if not 0.0 < below_threshold < 1.0:
        raise ValueError("below_threshold must lie within (0.0, 1.0).")
    maximum_value = np.max(time_course)
    t_max = np.argmax(time_course)
    time_course = time_course - below_threshold * maximum_value
    time_course[time_course > 0.0] = -np.inf
    duration = np.argmax(time_course[t_max:]) + t_max

    return duration


def test_sensitivity_analysis():
    for target in [
        # "parameter",
        "initial_condition",
        # "reaction",
    ]:
        for metric in [
            "maximum",
            "minimum",
            "integral",
            "argmax",
            "timepoint",
            "duration",
        ]:
            run_analysis(
                model,
                target=target,
                metric=metric,
                create_metrics={
                    "argmax": np.argmax,
                    "timepoint": lambda time_course: time_course[model.problem.t[-1]],
                    "duration": _get_duration,
                },
            )
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
    shutil.rmtree(MODEL_NAME)
