import os
import re
import csv
import numpy as np

from biomass.model.name2idx import C, V
from biomass.model.set_search_param import get_search_index


def get_param():
    """
    Get optimized parameters as CSV file format

    Output
    ------
    optimized_params.csv
    optimized_initial_values.csv
    """
    n_file = []
    fitparam_files = os.listdir('./out')
    for file in fitparam_files:
        if re.match(r'\d', file):
            n_file.append(int(file))

    search_idx = get_search_index()

    if len(search_idx[0]) > 0:
        optimized_params = np.empty(
            (len(search_idx[0])+2, len(n_file)+1), dtype='<U21'
        )
        for i, param_index in enumerate(search_idx[0]):
            for j, nth_paramset in enumerate(n_file):
                best_generation = np.load(
                    './out/{:d}/generation.npy'.format(
                        nth_paramset
                    )
                )
                best_indiv = np.load(
                    './out/{:d}/fit_param{:d}.npy'.format(
                        nth_paramset, int(best_generation)
                    )
                )
                error = np.load(
                    './out/{:d}/best_fitness.npy'.format(
                        nth_paramset
                    )
                )
                optimized_params[0, 0] = ''
                optimized_params[1, 0] = '*Error*'
                optimized_params[i+2, 0] = C.param_names[param_index]
                optimized_params[0, nth_paramset] = str(nth_paramset)
                optimized_params[1, nth_paramset] = '{:8.3e}'.format(error)
                optimized_params[i+2, nth_paramset] = '{:8.3e}'.format(best_indiv[i])
        with open('optimized_params.csv', 'w') as f:
            writer = csv.writer(f, lineterminator='\n')
            writer.writerows(optimized_params)

    if len(search_idx[1]) > 0:
        optimized_initvars = np.empty(
            (len(search_idx[1])+2, len(n_file)+1), dtype='<U21'
        )
        for i, var_index in enumerate(search_idx[1]):
            for j, nth_paramset in enumerate(n_file):
                best_generation = np.load(
                    './out/{:d}/generation.npy'.format(
                        nth_paramset
                    )
                )
                best_indiv = np.load(
                    './out/{:d}/fit_param{:d}.npy'.format(
                        nth_paramset, int(best_generation)
                    )
                )
                error = np.load(
                    './out/{:d}/best_fitness.npy'.format(
                        nth_paramset
                    )
                )
                optimized_initvars[0, 0] = ''
                optimized_initvars[1, 0] = '*Error*'
                optimized_initvars[i+2, 0] = V.var_names[var_index]
                optimized_initvars[0, nth_paramset] = str(nth_paramset)
                optimized_initvars[1, nth_paramset] = '{:8.3e}'.format(error)
                optimized_initvars[i+2, nth_paramset] = \
                    '{:8.3e}'.format(best_indiv[i+len(search_idx[0])])
        with open('optimized_inital_varlues.csv', 'w') as f:
            writer = csv.writer(f, lineterminator='\n')
            writer.writerows(optimized_initvars)


if __name__ == '__main__':
    get_param()
