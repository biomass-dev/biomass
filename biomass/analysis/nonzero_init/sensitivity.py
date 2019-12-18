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
from biomass.param_estim.search_parameter import search_parameter_index


"""
Calculation of the duration as the time it takes to decline below 10% of its maximum
"""
def get_duration(time_course_vector):
    maximum_value = np.max(time_course_vector)
    t_max = np.argmax(time_course_vector)

    time_course_vector = time_course_vector - 0.1*maximum_value
    time_course_vector[time_course_vector > 0.0] = -np.inf

    duration = np.argmax(time_course_vector[t_max:]) + t_max

    return duration

def analyze_sensitivity(nonzero_idx):
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

    signaling_metric_cfos_mRNA = np.full((n_file,len(nonzero_idx)+1,len(sim.conditions)),np.nan)
    signaling_metric_PcFos     = np.full((n_file,len(nonzero_idx)+1,len(sim.conditions)),np.nan)

    search_idx = search_parameter_index()

    for i in range(n_file):
        if os.path.isfile('./out/%d/generation.npy'%(i+1)):
            generation = np.load('./out/%d/generation.npy'%(i+1))
            best_indiv = np.load('./out/%d/fit_param%d.npy'%(i+1,int(generation)))

            for m,n in enumerate(search_idx[0]):
                x[n] = best_indiv[m]
            for m,n in enumerate(search_idx[1]):
                y0[n] = best_indiv[m+len(search_idx[0])]

            copy_y0 = y0[:]

            for j,idx in enumerate(nonzero_idx):
                y0[idx] = copy_y0[idx]*rate
                
                if sim.simulate(x,y0) is None:
                    for k,_ in enumerate(sim.conditions):
                        signaling_metric_cfos_mRNA[i,j,k] = \
                            get_duration(
                                sim.simulations[observables.index('cfos_mRNA'),:,k]
                            )
                        signaling_metric_PcFos[i,j,k] = \
                            simps(
                                sim.simulations[observables.index('Phosphorylated_cFos'),:,k]
                            )

                sys.stdout.write('\r%d/%d'%(i*len(nonzero_idx)+j+1,n_file*len(nonzero_idx)))
            
            y0 = copy_y0[:]
            if sim.simulate(x,y0) is None:
                for k,_ in enumerate(sim.conditions):
                    signaling_metric_cfos_mRNA[i,len(nonzero_idx),k] = \
                        get_duration(
                            sim.simulations[observables.index('cfos_mRNA'),:,k]
                        )
                    signaling_metric_PcFos[i,len(nonzero_idx),k] = \
                        simps(
                            sim.simulations[observables.index('Phosphorylated_cFos'),:,k]
                        )
            

    sensitivity_coefficients_cfos_mRNA = np.empty_like(signaling_metric_cfos_mRNA)
    sensitivity_coefficients_PcFos = np.empty_like(signaling_metric_PcFos)

    for l in range(n_file):
        sensitivity_coefficients_cfos_mRNA[l,:,:] = \
            np.log(signaling_metric_cfos_mRNA[l,:,:]/signaling_metric_cfos_mRNA[l,len(nonzero_idx),:])/np.log(rate/1)

        sensitivity_coefficients_PcFos[l,:,:] = \
            np.log(signaling_metric_PcFos[l,:,:]/signaling_metric_PcFos[l,len(nonzero_idx),:])/np.log(rate/1)

    return sensitivity_coefficients_cfos_mRNA, sensitivity_coefficients_PcFos