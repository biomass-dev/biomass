import os
import shutil

import matplotlib.pyplot as plt
import numpy as np

from biomass import create_model, run_simulation
from biomass.models import copy_to_current

MODEL_NAME: str = "circadian_clock"

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

    plt.rcParams["font.size"] = 16
    fig, ax1 = plt.subplots(figsize=(6, 4))
    ax2 = ax1.twinx()

    ax1.plot(model.problem.t, res[model.observables.index("Cry_mRNA"), 0], "c")
    ax1.plot(model.problem.t, res[model.observables.index("Per_mRNA"), 0], "m")
    ax1.set_xlim([0, 72])
    ax1.set_xticks([0, 12, 24, 36, 48, 60, 72])
    ax1.set_xlabel("Time (h)")
    ax1.set_ylim([0, 5])
    ax1.set_ylabel(
        r"$\it{Per}$"
        + " "
        + r"$\sf{(M_P)}$"
        + " and "
        + r"$\it{Cry}$"
        + " "
        + r"$\sf{(M_C)}$"
        + "\nmRNAs, nM"
    )
    ax2.plot(model.problem.t, res[model.observables.index("Bmal1_mRNA"), 0], "y")
    ax2.set_ylim([7, 10])
    ax2.set_yticks([7, 8, 9, 10])
    ax2.set_ylabel(r"$\it{Bmal1}$" + " mRNA " + r"$\sf{(M_B)}$" + ", nM")
    plt.show()


def test_cleanup():
    shutil.rmtree(MODEL_NAME)
