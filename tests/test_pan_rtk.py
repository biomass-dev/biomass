import os
import shutil

import matplotlib.pyplot as plt
import numpy as np

from biomass import create_model, run_simulation
from biomass.models import copy_to_current

MODEL_NAME: str = "pan_rtk"

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

    colors = ["k", "b", "c", "r", "y"]
    yticks = [
        np.arange(0.5, 3, 0.5),
        np.arange(-0.8, 1.2, 0.2),
        np.arange(-0.8, 1, 0.2),
        np.arange(-1.2, 0.2, 0.2),
        np.arange(-0.8, 0.2, 0.2),
        np.arange(-0.5, 0.3, 0.1),
    ]
    sd = 1.0e-1

    plt.rcParams["font.size"] = 6
    plt.rcParams["font.family"] = "Arial"
    # plt.rcParams['axes.linewidth'] = 1
    plt.rcParams["lines.linewidth"] = 0.8

    plt.subplots_adjust(wspace=0.5, hspace=0.4)

    for i in range(6):
        if i < 3:
            plt.subplot(2, 4, i + 1)
        else:
            plt.subplot(2, 4, i + 2)
        plt.figsize = (7, 4)
        plt.gca().spines["right"].set_visible(False)
        plt.gca().spines["top"].set_visible(False)

        for j, color in enumerate(colors):
            plt.plot(
                model.problem.t,
                res[i, j],
                color,
                label=model.problem.conditions[j].replace("_", "."),
            )
            plt.fill_between(
                model.problem.t, res[i, j] - sd, res[i, j] + sd, facecolor=color, lw=0, alpha=0.1
            )
        plt.title(model.observables[i][:-3], fontweight="bold")
        plt.xlabel("time [min]")
        plt.xticks([0, 60, 120, 180, 240])
        plt.ylabel("(conc.) [au]")
        plt.yticks(yticks[i])
        if i == 2:
            plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0, frameon=False)
    plt.show()


def test_cleanup():
    shutil.rmtree(MODEL_NAME)
