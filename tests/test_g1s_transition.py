import os
import shutil

import matplotlib.pyplot as plt
import numpy as np

from biomass import Model, run_simulation
from biomass.models import g1s_transition

model = Model(g1s_transition.__package__).create()


for dir in ["figure", "simulation_data"]:
    if os.path.isdir(os.path.join(model.path, dir)):
        shutil.rmtree(os.path.join(model.path, dir))


def test_simulate_successful():
    x = model.pval()
    y0 = model.ival()
    assert model.problem.simulate(x, y0) is None


def test_run_simulation():
    assert run_simulation(model) is None
    res = np.load(os.path.join(model.path, "simulation_data", "simulations_original.npy"))

    plt.figure(figsize=(9,6))
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 16
    plt.rcParams['axes.linewidth'] = 1.5
    plt.rcParams['lines.linewidth'] = 2
    plt.rcParams['lines.markersize'] = 20

    plt.plot(model.problem.t, res[model.observables.index('p27_tot'), 0],'r-')
    plt.plot(model.problem.t, res[model.observables.index('CycE_tot'), 0],'b-')
    plt.plot(model.problem.t, res[model.observables.index('CycE'), 0],'b--')
    plt.plot(model.problem.t, res[model.observables.index('CycA_tot'), 0],'g-')
    plt.plot(model.problem.t, res[model.observables.index('CycA'), 0],'g--')

    plt.xticks([0,300,600,900],[-5,0,5,10])
    plt.xlim([0,900])
    plt.xlabel('time relative to G1/S transition(h)')
    plt.yticks([0,0.3,0.6,0.9])
    plt.ylabel('relative levels')

    plt.text(3.2*60,0.75,'p27',ha='center',va='bottom',color='r')
    plt.text(6.4*60,0.82,'CycE level',ha='center',va='bottom',color='b')
    plt.text(6.4*60,0.55,'free\nCycE:Cdk2',ha='center',va='bottom',color='b')
    plt.text(12.4*60,0.6,'CycA level',ha='center',va='bottom',color='g')
    plt.text(12.4*60,0.35,'free\nCycA:Cdk2',ha='center',va='bottom',color='g')

    plt.show()


def test_cleanup():
    for dir in ["figure", "simulation_data"]:
        if os.path.isdir(os.path.join(model.path, dir)):
            shutil.rmtree(os.path.join(model.path, dir))
