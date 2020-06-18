import os
import re
import csv
import numpy as np

from biomass.models.Nakakuki_Cell_2010 import *
from biomass.dynamics import get_executable

class CurrentModel(object):
    def __init__(self,
                 parameters=C.NAMES,
                 species=V.NAMES,
                 reaction_system=set_model,
                 pval=param_values,
                 ival=initial_values,
                 obs=observables,
                 sim=NumericalSimulation(),
                 exp=ExperimentalData(),
                 sp=SearchParam(),
                 obj_func=objective,
                 rxn=ReactionNetwork()):
        self.parameters = parameters
        self.species = species
        self.reaction_system = reaction_system
        self.pval = pval
        self.ival = ival
        self.obs = obs
        self.obj_func = obj_func
        self.sim = sim
        self.exp = exp
        self.sp = sp
        self.rxn = rxn

    def get_properties(self):
        """ Get number of model reactions, species and parameters
        """
        biological_processes = self.rxn.group()
        model_reactions = 0
        for reactions_in_process in biological_processes:
            model_reactions += len(reactions_in_process)
        model_species = len(self.species)
        model_parameters = len(self.parameters)

        print(
            "{:d} reactions\n{:d} species\n{:d} parameters".format(
                model_reactions, model_species, model_parameters
            )
        )


    def get_optimization_results(self):
        """ Get optimized parameters as CSV file format

        Output
        ------
        optimized_params.csv
        optimized_initial_values.csv
        """
        n_file = get_executable()

        if len(self.sp.idx_params) > 0:
            optimized_params = np.empty(
                (len(self.sp.idx_params)+2, len(n_file)+1), dtype='<U21'
            )
            for i, param_index in enumerate(self.sp.idx_params):
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
                    optimized_params[i+2, 0] = self.parameters[param_index]
                    optimized_params[0, nth_paramset-1] = str(nth_paramset)
                    optimized_params[1, nth_paramset-1] = '{:8.3e}'.format(error)
                    optimized_params[i+2, nth_paramset-1] = '{:8.3e}'.format(best_indiv[i])
            with open('optimized_params.csv', 'w') as f:
                writer = csv.writer(f, lineterminator='\n')
                writer.writerows(optimized_params)

        if len(self.sp.idx_initials) > 0:
            optimized_initials = np.empty(
                (len(self.sp.idx_initials)+2, len(n_file)+1), dtype='<U21'
            )
            for i, specie_index in enumerate(self.sp.idx_initials):
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
                    optimized_initials[i+2, 0] = self.species[specie_index]
                    optimized_initials[0, nth_paramset-1] = str(nth_paramset)
                    optimized_initials[1, nth_paramset-1] = '{:8.3e}'.format(error)
                    optimized_initials[i+2, nth_paramset-1] = \
                        '{:8.3e}'.format(best_indiv[i+len(self.sp.idx_params)])
            with open('optimized_inital_varlues.csv', 'w') as f:
                writer = csv.writer(f, lineterminator='\n')
                writer.writerows(optimized_initials)