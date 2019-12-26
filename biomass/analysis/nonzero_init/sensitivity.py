import os
import sys
import re
import numpy as np
from math import fabs, log
from scipy.integrate import simps

from biomass.model.name2idx import parameters as C
from biomass.model.name2idx import variables as V
from biomass.model.param_const import f_params
from biomass.model.initial_condition import initial_values
from biomass.observable import observables, NumericalSimulation
from biomass.param_estim.search_parameter import search_parameter_index


def get_duration(time_course_vector):
    """
    Calculation of the duration as the time it takes to decline below 10% of its maximum
    """
    maximum_value = np.max(time_course_vector)
    t_max = np.argmax(time_course_vector)
    time_course_vector = time_course_vector - 0.1*maximum_value
    time_course_vector[time_course_vector > 0.0] = -np.inf
    duration = np.argmax(time_course_vector[t_max:]) + t_max

    return duration


def analyze_sensitivity(metric,nonzero_idx):
    """Compute sensitivity coefficients

    Parameters
    ----------
    metric: str
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

    rate = 1.01 # 1% change

    x = f_params()
    y0 = initial_values()

    nonzero_idx = []
    for i,val in enumerate(y0):
        if val != 0.0:
            nonzero_idx.append(i)
    n_file = 0
    fitparam_files = os.listdir('./out')
    for file in fitparam_files:
        if re.match(r'\d',file):
            n_file += 1

    signaling_metric = np.full(
        (n_file,len(nonzero_idx)+1,len(observables),len(sim.conditions)),np.nan
    )
    search_idx = search_parameter_index()
    for i in range(n_file):
        if os.path.isfile('./out/%d/generation.npy'%(i+1)):
            best_generation = \
                np.load('./out/%d/generation.npy'%(i+1))
            best_indiv = \
                np.load('./out/%d/fit_param%d.npy'%(i+1,int(best_generation)))
            for m,n in enumerate(search_idx[0]):
                x[n] = best_indiv[m]
            for m,n in enumerate(search_idx[1]):
                y0[n] = best_indiv[m+len(search_idx[0])]
            copy_y0 = y0[:]
            for j,idx in enumerate(nonzero_idx):
                y0[idx] = copy_y0[idx]*rate
                if sim.simulate(x,y0) is None:
                    for k,_ in enumerate(observables):
                        for l,_ in enumerate(sim.conditions):
                            if metric == 'duration':
                                signaling_metric[i,j,k,l] = get_duration(
                                    sim.simulations[k,:,l]
                                )
                            elif metric == 'integral':
                                signaling_metric[i,j,k,l] = simps(
                                    sim.simulations[k,:,l]
                                )
                            else:
                                raise ValueError(
                                    "metric ∈ {'duration','integral'}"
                                )
                sys.stdout.write(
                    '\r%d/%d'%(i*len(nonzero_idx)+j+1,n_file*len(nonzero_idx))
                )
            y0 = copy_y0[:]
            if sim.simulate(x,y0) is None:
                for k,_ in enumerate(observables):
                    for l,_ in enumerate(sim.conditions):
                        if metric == 'duration':
                            signaling_metric[i,-1,k,l] = get_duration(
                                sim.simulations[k,:,l]
                            )
                        elif metric == 'integral':
                            signaling_metric[i,-1,k,l] = simps(
                                sim.simulations[k,:,l]
                            )
                        else:
                            raise ValueError(
                                "metric ∈ {'duration','integral'}"
                            )
    sensitivity_coefficients = np.empty(
        (n_file,len(nonzero_idx),len(observables),len(sim.conditions))
    )
    for i in range(n_file):
        for j,_ in enumerate(nonzero_idx):
            for k,_ in enumerate(observables):
                for l,_ in enumerate(sim.conditions):
                    if signaling_metric[i,j,k,l] < sys.float_info.epsilon or \
                        (signaling_metric[i,j,k,l]/signaling_metric[i,-1,k,l]) < 0:
                        sensitivity_coefficients[i,j,k,l] = np.nan
                    elif fabs(1 - signaling_metric[i,j,k,l]/signaling_metric[i,-1,k,l]) < sys.float_info.epsilon:
                        sensitivity_coefficients[i,j,k,l] = 0.0
                    else:
                        sensitivity_coefficients[i,j,k,l] = (
                            log(signaling_metric[i,j,k,l]/signaling_metric[i,-1,k,l])/log(rate)
                        )
    return sensitivity_coefficients