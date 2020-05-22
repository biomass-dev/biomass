import os
import sys
import re
import numpy as np

from biomass.model import initial_values
from biomass.observable import observables, NumericalSimulation
from biomass.param_estim import load_param
from biomass.analysis import get_signaling_metric, dlnyi_dlnxj


def calc_sensitivity_coefficients(metric, nonzero_idx):
    """ Calculating Sensitivity Coefficients

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
        if os.path.isfile('./out/{:d}/generation.npy'.format(nth_paramset)):
            (x, y0) = load_param(nth_paramset)
            copy_y0 = y0[:]
            for j, idx in enumerate(nonzero_idx):
                y0 = copy_y0[:]
                y0[idx] = copy_y0[idx] * rate
                if sim.simulate(x, y0) is None:
                    for k, _ in enumerate(observables):
                        for l, _ in enumerate(sim.conditions):
                            signaling_metric[i, j, k, l] = get_signaling_metric(
                                metric, sim.simulations[k, :, l]
                            )
                sys.stdout.write(
                    '\r{:d} / {:d}'.format(
                        i*len(nonzero_idx)+j+1, len(n_file)*len(nonzero_idx)
                    )
                )
            # Signaling metric without perturbation (j=-1)
            y0 = copy_y0[:]
            if sim.simulate(x, y0) is None:
                for k, _ in enumerate(observables):
                    for l, _ in enumerate(sim.conditions):
                        signaling_metric[i, -1, k, l] = get_signaling_metric(
                            metric, sim.simulations[k, :, l]
                        )
    sensitivity_coefficients = dlnyi_dlnxj(
        signaling_metric, n_file, nonzero_idx,
        observables, sim.conditions, rate, metric_idx=-1
    )

    return sensitivity_coefficients
