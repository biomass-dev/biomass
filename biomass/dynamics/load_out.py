import os
import re
import numpy as np


def _get_indiv(paramset):
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

    return best_indiv


def load_param(paramset, update):
    best_indiv = _get_indiv(paramset)
    (x, y0) = update(best_indiv)

    return x, y0


def write_best_fit_param(best_paramset, parameters, species, update):
    (x, y0) = load_param(best_paramset, update)
    
    with open('./out/best_fit_param.txt', mode='w') as f:
        f.write(
            '# param set: {:d}\n'.format(
                best_paramset
            )
        )
        f.write(
            '\n### f_params\n'
        )
        for i, param in enumerate(parameters):
            f.write('x[C.{}] = {:8.3e}\n'.format(param, x[i]))
        f.write(
            '\n### initial_values\n'
        )
        for i, specie in enumerate(species):
            if y0[i] != 0:
                f.write('y0[V.{}] = {:8.3e}\n'.format(specie, y0[i]))


def get_optimized_param(n_file, search_idx):
    popt = np.empty(
        (len(n_file), len(search_idx[0]) + len(search_idx[1]))
    )
    for k, nth_paramset in enumerate(n_file):
        best_indiv = _get_indiv(nth_paramset)
        popt[k, :] = best_indiv

    return popt


def get_executable():
    n_file = []
    fitparam_files = os.listdir('./out')
    for file in fitparam_files:
        if re.match(r'\d', file):
            n_file.append(int(file))
    empty_folder = []
    for i, nth_paramset in enumerate(n_file):
        if not os.path.isfile('./out/{:d}/generation.npy'.format(nth_paramset)):
            empty_folder.append(i)
    for i in sorted(empty_folder, reverse=True):
        n_file.pop(i)
    
    return n_file