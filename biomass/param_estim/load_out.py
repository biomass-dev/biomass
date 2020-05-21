import os
import numpy as np

from biomass.model import C, V, f_params, initial_values
from .search_parameter import update_param


def load_best_param(paramset):
    best_generation = np.load(
        './out/{:d}/generation.npy'.format(
            paramset
        )
    )
    best_indiv = np.load(
        './out/{:d}/fit_param{:d}.npy'.format(
            paramset, int(best_generation)
        )
    )
    
    (x, y0) = update_param(best_indiv)

    return x, y0


def write_best_fit_param(best_paramset, search_idx):
    x = f_params()
    y0 = initial_values()

    best_generation = np.load(
        './out/{:d}/generation.npy'.format(
            best_paramset
        )
    )
    best_indiv = np.load(
        './out/{:d}/fit_param{:d}.npy'.format(
            best_paramset, int(best_generation)
        )
    )
    for i, j in enumerate(search_idx[0]):
        x[j] = best_indiv[i]
    for i, j in enumerate(search_idx[1]):
        y0[j] = best_indiv[i+len(search_idx[0])]
    
    with open('./out/best_fit_param.txt', mode='w') as f:
        f.write(
            '# param set: {:d}\n'.format(
                best_paramset
            )
        )
        f.write(
            '\n### f_params\n'
        )
        for i in range(C.len_f_params):
            f.write(
                'x[C.{}] = {:8.3e}\n'.format(
                    C.param_names[i], x[i]
                )
            )
        f.write(
            '\n### initial_values\n'
        )
        for i in range(V.len_f_vars):
            if y0[i] != 0:
                f.write(
                    'y0[V.{}] = {:8.3e}\n'.format(
                        V.var_names[i], y0[i]
                    )
                )


def get_optimized_param(n_file, search_idx):
    popt = np.empty(
        (len(n_file), len(search_idx[0]) + len(search_idx[1]))
    )
    for k, nth_paramset in enumerate(n_file):
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

        popt[k, :] = best_indiv

    return popt
