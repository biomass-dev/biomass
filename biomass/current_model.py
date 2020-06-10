import os
import re
import csv
import numpy as np

from biomass.models.Nakakuki_Cell_2010 import *


def get_model_properties():
    """ Get number of model reactions, species and parameters
    """
    rxn = ReactionNetwork()
    biological_processes = rxn.group()
    model_reactions = 0
    for reactions_in_process in biological_processes:
        model_reactions += len(reactions_in_process)
    model_species = V.n_species
    model_parameters = C.n_parameters

    print(
        "{:d} reactions\n{:d} species\n{:d} parameters".format(
            model_reactions, model_species, model_parameters
        )
    )


def get_optimization_results():
    """ Get optimized parameters as CSV file format

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
        optimized_initials = np.empty(
            (len(search_idx[1])+2, len(n_file)+1), dtype='<U21'
        )
        for i, specie_index in enumerate(search_idx[1]):
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
                optimized_initials[0, 0] = ''
                optimized_initials[1, 0] = '*Error*'
                optimized_initials[i+2, 0] = V.species[specie_index]
                optimized_initials[0, nth_paramset] = str(nth_paramset)
                optimized_initials[1, nth_paramset] = '{:8.3e}'.format(error)
                optimized_initials[i+2, nth_paramset] = \
                    '{:8.3e}'.format(best_indiv[i+len(search_idx[0])])
        with open('optimized_inital_varlues.csv', 'w') as f:
            writer = csv.writer(f, lineterminator='\n')
            writer.writerows(optimized_initials)