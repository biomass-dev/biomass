import os
import shutil

import matplotlib.pyplot as plt
import numpy as np

from biomass import create_model, run_simulation
from biomass.models import copy_to_current

MODEL_NAME: str = "insulin_signaling"

copy_to_current(MODEL_NAME)
assert os.path.exists(MODEL_NAME)
model = create_model(MODEL_NAME)


def test_simulate_successful():
    x = model.pval()
    y0 = model.ival()
    assert model.problem.simulate(x, y0) is None


def test_example_plot():
    assert run_simulation(model) is None
    res = np.load(os.path.join(model.path, "simulation_data", "simulations_original.npy"))

    plt.figure(figsize=(8, 8))
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["xtick.direction"] = "in"
    plt.rcParams["ytick.direction"] = "in"
    plt.rcParams["font.size"] = 18
    plt.rcParams["lines.linewidth"] = 4
    plt.rcParams["lines.markersize"] = 16

    plt.subplots_adjust(wspace=0.25, hspace=0.3)

    for i, obs_name in enumerate(model.observables):
        plt.subplot(2, 2, i + 1)
        for j, color in enumerate(
            ["darkblue", "cornflowerblue", "yellowgreen", "goldenrod", "brown"]
        ):
            plt.plot(model.problem.t, res[i, j], color=color)
        plt.xlim(0, 480)
        plt.xticks([0, 120, 240, 360, 480], fontweight="bold")
        plt.ylim(0, 1)
        plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0], fontweight="bold")
        plt.title(f"{obs_name}", fontweight="bold")
        if i % 2 == 0:
            plt.ylabel("Intensity (A.U.)", fontweight="bold")
        if i > 1:
            plt.xlabel("                           time (min)", fontweight="bold")
    plt.show()


def test_cleanup():
    shutil.rmtree(MODEL_NAME)
