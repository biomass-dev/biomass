import os
import shutil

import numpy as np

from biomass import Model, Text2Model, optimize, run_simulation

from .C import BOUNDS, CONDITIONS, EXPERIMENTAL_DATA, INDENT, NORMALIZATION, PPMEK

PATH_TO_MODEL = os.path.join(
    os.path.dirname(__file__),
    "text_files",
    "fos_model",
)


def test_preprocessing():
    if os.path.isdir(PATH_TO_MODEL):
        shutil.rmtree(PATH_TO_MODEL)


def test_model_construction():
    cfos = Text2Model(
        os.path.join(
            os.path.dirname(__file__),
            "text_files",
            "fos_model.txt",
        )
    )
    cfos.convert()
    assert os.path.isdir(PATH_TO_MODEL)

    with open(os.path.join(PATH_TO_MODEL, "ode.py"), mode="r") as f1:
        lines = f1.readlines()
    for line_num, line in enumerate(lines):
        if line.startswith(f"{INDENT}def diffeq(self, t, y, *x):"):
            lines[line_num - 1] = PPMEK
        elif line.startswith(f"{2*INDENT}dydt = [0] * V.NUM"):
            lines[line_num] += CONDITIONS
            break
    else:
        assert False
    with open(os.path.join(PATH_TO_MODEL, "ode.py"), mode="w") as f1:
        f1.writelines(lines)

    with open(os.path.join(PATH_TO_MODEL, "search_param.py"), mode="r") as f2:
        lines = f2.readlines()
    for line_num, line in enumerate(lines):
        if line.startswith(f"{2*INDENT}search_rgn = convert_scale("):
            lines[line_num - 1] = BOUNDS
            break
    else:
        assert False
    with open(os.path.join(PATH_TO_MODEL, "search_param.py"), mode="w") as f2:
        f2.writelines(lines)

    with open(os.path.join(PATH_TO_MODEL, "observable.py"), mode="r") as f3:
        lines = f3.readlines()
    for line_num, line in enumerate(lines):
        if line.startswith(f"{2*INDENT}self.normalization: dict = {{}}"):
            lines[line_num] += NORMALIZATION
        elif "y0 = get_steady_state(self.diffeq, y0, tuple(x))" in line:
            lines[line_num] = line.replace(
                "y0 = get_steady_state(self.diffeq, y0, tuple(x))",
                "y0 = get_steady_state(self.diffeq, y0, tuple(x), integrator='vode')",
            )
        elif "sol = solve_ode(self.diffeq, y0, self.t, tuple(x))" in line:
            lines[line_num] = line.replace(
                "sol = solve_ode(self.diffeq, y0, self.t, tuple(x))",
                "sol = solve_ode(self.diffeq, y0, self.t, tuple(x), method='BDF')",
            )
        elif line.startswith(f"{INDENT}def set_data(self)"):
            lines = lines[: line_num + 1]
            lines[line_num] = EXPERIMENTAL_DATA
            break
    else:
        assert False
    with open(os.path.join(PATH_TO_MODEL, "observable.py"), mode="w") as f3:
        f3.writelines(lines)


def test_parameter_estimation(workers: int = -1):
    model = Model("tests.test_text2model.text_files.fos_model").create()
    assert len(model.species) == 36
    assert len(model.observables) == 8
    assert len(model.problem.bounds) == 75

    param_idx = 1
    optimize(model, param_idx, optimizer_options={"maxiter": 10, "workers": workers})
    assert os.path.isdir(os.path.join(model.path, "out", f"{param_idx}"))
    assert run_simulation(model, viz_type=str(param_idx)) is None
    simulation_results = np.load(os.path.join(model.path, "simulation_data", "simulations_1.npy"))
    assert np.isfinite(simulation_results).all()
    assert len(
        [
            file
            for file in os.listdir(os.path.join(model.path, "figure", "simulation", "1"))
            if file.endswith(".pdf") or file.endswith(".png")
        ]
    ) == len(model.observables), ", ".join(
        os.listdir(os.path.join(model.path, "figure", "simulation", "1"))
    )


def test_cleanup():
    assert os.path.isdir(PATH_TO_MODEL)
    shutil.rmtree(PATH_TO_MODEL)
