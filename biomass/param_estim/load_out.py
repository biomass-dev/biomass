import os
import numpy as np

from biomass.model import C, V
from biomass.param_estim import search_parameter_index


def update_param(paramset, x, y0):
    search_idx = search_parameter_index()

    if os.path.isfile('./out/{:d}/generation.npy'.format(paramset)):
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
        for i, j in enumerate(search_idx[0]):
            x[j] = best_indiv[i]
        for i, j in enumerate(search_idx[1]):
            y0[j] = best_indiv[i+len(search_idx[0])]

    return x, y0


def write_best_fit_param(best_paramset, x, y0):
    (x, y0) = update_param(best_paramset, x, y0)
    with open('./out/best_fit_param.txt', mode='w') as f:
        f.write(
            '# param set: {:d}\n'.format(
                best_paramset
            )
        )
        f.write(
            '\n### Param. const\n'
        )
        for i in range(C.len_f_params):
            f.write(
                'x[C.{}] = {:8.3e}\n'.format(
                    C.param_names[i], x[i]
                )
            )
        f.write(
            '\n### Non-zero initial conditions\n'
        )
        for i in range(V.len_f_vars):
            if y0[i] != 0:
                f.write(
                    'y0[V.{}] = {:8.3e}\n'.format(
                        V.var_names[i], y0[i]
                    )
                )


def get_search_param_matrix(n_file, x, y0):
    search_idx = search_parameter_index()
    search_param_matrix = np.empty(
        (len(n_file), len(search_idx[0]) + len(search_idx[1]))
    )
    for k, nth_paramset in enumerate(n_file):
        if os.path.isfile('./out/{:d}/generation.npy'.format(nth_paramset)):
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
        else:
            best_indiv = np.empty(
                len(search_idx[0]) + len(search_idx[1])
            )
            for i, j in enumerate(search_idx[0]):
                best_indiv[i] = x[j]
            for i, j in enumerate(search_idx[1]):
                best_indiv[i+len(search_idx[0])] = y0[j]

        search_param_matrix[k, :] = best_indiv

    return search_idx, search_param_matrix
