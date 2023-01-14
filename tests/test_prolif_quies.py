import os
import shutil

import matplotlib.pyplot as plt
import numpy as np

from biomass import create_model, run_simulation
from biomass.models import copy_to_current

MODEL_NAME: str = "prolif_quies"

copy_to_current(MODEL_NAME)
assert os.path.exists(MODEL_NAME)
model = create_model(MODEL_NAME)


def test_simulate_successful():
    x = model.pval()
    y0 = model.ival()
    assert model.problem.simulate(x, y0) is None


def test_run_simulation():
    assert run_simulation(model) is None
    res = np.load(os.path.join(model.path, "simulation_data", "simulations_original.npy"))

    plt.figure(figsize=(10, 5))
    plt.rcParams["font.size"] = 32
    plt.rcParams["axes.linewidth"] = 1.5
    plt.rcParams["lines.linewidth"] = 6

    plt.plot(
        [t_ / 60 for t_ in model.problem.t],
        res[model.observables.index("CycA"), 0],
        color="slateblue",
        label="CycA",
    )
    plt.plot(
        [t_ / 60 for t_ in model.problem.t],
        res[model.observables.index("CycE"), 0],
        color="skyblue",
        label="CycE",
    )
    plt.plot(
        [t_ / 60 for t_ in model.problem.t],
        res[model.observables.index("active_RC"), 0],
        color="firebrick",
        label="aRC",
    )
    plt.plot(
        [t_ / 60 for t_ in model.problem.t],
        res[model.observables.index("p21_tot"), 0],
        color="limegreen",
        label="p21",
    )

    plt.xlim(0, 20)
    plt.xlabel("time from cytokinesis (h)")
    plt.ylim(0, 2)
    plt.ylabel("rel. level (AU)")
    plt.xticks([0, 5, 10, 15, 20])
    plt.yticks([0, 1, 2])
    plt.legend(loc="upper left", frameon=False, fontsize=18)

    plt.show()


def test_cleanup():
    shutil.rmtree(MODEL_NAME)
