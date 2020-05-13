import os
import sys
import re
import numpy as np

from biomass.model import C, V, set_model, f_params, initial_values
from biomass.observable import observables, NumericalSimulation
from biomass.param_estim import update_param
from biomass.analysis import get_signaling_metric, dlnyi_dlnxj


def calc_sensitivity_coefficients(metric, num_reaction):
    """ Calculating Sensitivity Coefficients

    Parameters
    ----------
    metric: str
        - 'amplitude': The maximum value.
        - 'duration': The time it takes to decline below 10% of its maximum.
        - 'integral': The integral of concentration over the observation time.
    num_reaction: int
        len(v) in model/differential_equation.py

    Returns
    -------
    sensitivity_coefficients: numpy array
    
    """
    sim = NumericalSimulation()

    rate = 1.01  # 1% change

    x = f_params()
    y0 = initial_values()

    n_file = []
    fitparam_files = os.listdir('./out')
    for file in fitparam_files:
        if re.match(r'\d', file):
            n_file.append(int(file))

    signaling_metric = np.full(
        (len(n_file), num_reaction, len(observables), len(sim.conditions)),
        np.nan
    )
    for i, nth_paramset in enumerate(n_file):
        if os.path.isfile('./out/{:d}/generation.npy'.format(nth_paramset)):
            (x, y0) = update_param(nth_paramset, x, y0)
            for j in range(num_reaction):
                set_model.perturbation = [1] * num_reaction
                set_model.perturbation[j] = rate
                if sim.simulate(x, y0) is None:
                    for k, _ in enumerate(observables):
                        for l, _ in enumerate(sim.conditions):
                            signaling_metric[i, j, k, l] = get_signaling_metric(
                                metric, sim.simulations[k, :, l]
                            )
                sys.stdout.write(
                    '\r{:d} / {:d}'.format(
                        i*num_reaction+j+1, len(n_file)*num_reaction
                    )
                )
    sensitivity_coefficients = dlnyi_dlnxj(
        signaling_metric, n_file, range(num_reaction),
        observables, sim.conditions, rate, metric_idx=0
    )

    return sensitivity_coefficients