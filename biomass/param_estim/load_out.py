import os
import numpy as np

from biomass.model import C, V
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


def write_best_fit_param(best_paramset):
    (x, y0) = load_best_param(best_paramset)
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
