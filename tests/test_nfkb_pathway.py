import os
import shutil

import matplotlib.pyplot as plt
import numpy as np

from biomass import create_model, run_simulation
from biomass.models import copy_to_current

MODEL_NAME: str = "nfkb_pathway"

copy_to_current(MODEL_NAME)
assert os.path.exists(MODEL_NAME)
model = create_model(MODEL_NAME)


def test_simulate_successful():
    x = model.pval()
    y0 = model.ival()
    assert model.problem.simulate(x, y0) is None


def test_example_plot():
    assert run_simulation(model, viz_type="original") is None
    res = np.load(os.path.join(model.path, "simulation_data", "simulations_original.npy"))

    fig = plt.figure(figsize=(9, 9))
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.size"] = 18
    plt.rcParams["axes.linewidth"] = 2
    plt.rcParams["lines.linewidth"] = 3
    plt.rcParams["lines.markersize"] = 16
    plt.subplots_adjust(wspace=0, hspace=0.3)

    for i, obs_name in enumerate(model.observables):
        plt.subplot(2, 1, i + 1)
        for j, (color, label) in enumerate(zip(["k", "r"], ["TNFα", "TNFα + DCF"])):
            plt.plot(model.problem.t, res[i, j], color=color, label=label)
        plt.title(f"{obs_name}".replace("_", " "))
        plt.xticks([0, 50, 100, 150, 200])
        if i == 0:
            plt.yticks([0, 0.05, 0.10, 0.15])
            plt.legend(loc="upper left", frameon=False)
    fig.text(0.5, 0.05, "time [min]", ha="center")
    fig.text(0.0, 0.5, "concentration [a.u.]", va="center", rotation="vertical")
    plt.show()


def test_cleanup():
    shutil.rmtree(MODEL_NAME)
