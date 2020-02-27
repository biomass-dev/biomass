import os
import sys
import re
import numpy as np
from scipy.integrate import simps

from biomass.model.name2idx import parameters as C
from biomass.model.name2idx import variables as V
from biomass.model.param_const import f_params
from biomass.model.initial_condition import initial_values
from biomass.observable import observables, NumericalSimulation
from biomass.param_estim.dynamics import update_param
from biomass.analysis.signaling_metric import get_duration, compute_sensitivity_coefficients


def analyze_sensitivity(metric, nonzero_idx):
    """Compute sensitivity coefficients

    Parameters
    ----------
    metric: str
        - 'amplitude': The maximum value.
        - 'duration': The time it takes to decline below 10% of its maximum.
        - 'integral': The integral of concentration over the observation time.
    nonzero_idx: list
        for i in nonzero_idx:
            y0[i] != 0.0

    Returns
    -------
    sensitivity_coefficients: numpy array
    
    """
    sim = NumericalSimulation()

    rate = 1.01  # 1% change
    epsilon = 1e-9  # If |M - M*| < epsilon, sensitivity_coefficient = 0
    

    x = f_params()
    y0 = initial_values()

    nonzero_idx = []
    for i, val in enumerate(y0):
        if val != 0.0:
            nonzero_idx.append(i)
    n_file = []
    fitparam_files = os.listdir('./out')
    for file in fitparam_files:
        if re.match(r'\d', file):
            n_file.append(int(file))

    signaling_metric = np.full(
        (len(n_file), len(nonzero_idx)+1, len(observables), len(sim.conditions)),
        np.nan
    )
    for i, nth_paramset in enumerate(n_file):
        if os.path.isfile('./out/%d/generation.npy' % (nth_paramset)):
            (x, y0) = update_param(nth_paramset, x, y0)
            copy_y0 = y0[:]
            for j, idx in enumerate(nonzero_idx):
                y0 = copy_y0[:]
                y0[idx] = copy_y0[idx] * rate
                if sim.simulate(x, y0) is None:
                    for k, _ in enumerate(observables):
                        for l, _ in enumerate(sim.conditions):
                            if metric == 'amplitude':
                                signaling_metric[i, j, k, l] = np.max(
                                    sim.simulations[k, :, l]
                                )
                            elif metric == 'duration':
                                signaling_metric[i, j, k, l] = get_duration(
                                    sim.simulations[k, :, l]
                                )
                            elif metric == 'integral':
                                signaling_metric[i, j, k, l] = simps(
                                    sim.simulations[k, :, l]
                                )
                            else:
                                raise ValueError(
                                    "metric ∈ {'amplitude', 'duration', 'integral'}"
                                )
                sys.stdout.write(
                    '\r%d / %d' % (
                        i*len(nonzero_idx)+j+1, len(n_file)*len(nonzero_idx)
                    )
                )
            y0 = copy_y0[:]
            if sim.simulate(x, y0) is None:
                for k, _ in enumerate(observables):
                    for l, _ in enumerate(sim.conditions):
                        if metric == 'amplitude':
                            signaling_metric[i, -1, k, l] = np.max(
                                sim.simulations[k, :, l]
                            )
                        elif metric == 'duration':
                            signaling_metric[i, -1, k, l] = get_duration(
                                sim.simulations[k, :, l]
                            )
                        elif metric == 'integral':
                            signaling_metric[i, -1, k, l] = simps(
                                sim.simulations[k, :, l]
                            )
                        else:
                            sys.exit()
    sensitivity_coefficients = compute_sensitivity_coefficients(
        signaling_metric, n_file, nonzero_idx,
        observables, sim.conditions, rate, metric_idx=-1
    )

    return sensitivity_coefficients
