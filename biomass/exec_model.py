import os
import numpy as np
import csv
import multiprocessing
import warnings

from biomass.dynamics import SignalingSystems, get_executable
from biomass.ga import GeneticAlgorithmInit
from biomass.ga import GeneticAlgorithmContinue
from biomass.analysis.reaction import ReactionSensitivity
from biomass.analysis.nonzero_init import NonZeroInitSensitivity


class ExecModel(object):
    def __init__(self, model):
        self.model_path = model.__path__[0]
        self.parameters = model.C.NAMES
        self.species = model.V.NAMES
        self.reaction_system = model.set_model
        self.pval = model.param_values
        self.ival = model.initial_values
        self.obs = model.observables
        self.obj_func = model.objective
        self.sim = model.NumericalSimulation()
        self.exp = model.ExperimentalData()
        self.sp = model.SearchParam()
        self.rxn = model.ReactionNetwork()

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
        optimization_results/optimized_params.csv
        optimization_results/optimized_initials.csv
        """
        os.makedirs(self.model_path + '/optimization_results/', exist_ok=True)
        n_file = get_executable(self.model_path)

        if len(self.sp.idx_params) > 0:
            optimized_params = np.empty(
                (len(self.sp.idx_params)+2, len(n_file)+1), dtype='<U21'
            )
            for i, param_index in enumerate(self.sp.idx_params):
                for j, nth_paramset in enumerate(n_file):
                    best_generation = np.load(
                        self.model_path + '/out/{:d}/generation.npy'.format(
                            nth_paramset
                        )
                    )
                    best_indiv = np.load(
                        self.model_path + '/out/{:d}/fit_param{:d}.npy'.format(
                            nth_paramset, int(best_generation)
                        )
                    )
                    error = np.load(
                        self.model_path + '/out/{:d}/best_fitness.npy'.format(
                            nth_paramset
                        )
                    )
                    optimized_params[0, 0] = ''
                    optimized_params[1, 0] = '*Error*'
                    optimized_params[i+2, 0] = self.parameters[param_index]
                    optimized_params[0, nth_paramset-1] = str(nth_paramset)
                    optimized_params[1, nth_paramset-1] = \
                        '{:8.3e}'.format(error)
                    optimized_params[i+2, nth_paramset-1] = \
                        '{:8.3e}'.format(best_indiv[i])
            with open(
                    self.model_path
                    + '/optimization_results/optimized_params.csv', 'w') as f:
                writer = csv.writer(f, lineterminator='\n')
                writer.writerows(optimized_params)

        if len(self.sp.idx_initials) > 0:
            optimized_initials = np.empty(
                (len(self.sp.idx_initials)+2, len(n_file)+1), dtype='<U21'
            )
            for i, specie_index in enumerate(self.sp.idx_initials):
                for j, nth_paramset in enumerate(n_file):
                    best_generation = np.load(
                        self.model_path + '/out/{:d}/generation.npy'.format(
                            nth_paramset
                        )
                    )
                    best_indiv = np.load(
                        self.model_path + '/out/{:d}/fit_param{:d}.npy'.format(
                            nth_paramset, int(best_generation)
                        )
                    )
                    error = np.load(
                        self.model_path + '/out/{:d}/best_fitness.npy'.format(
                            nth_paramset
                        )
                    )
                    optimized_initials[0, 0] = ''
                    optimized_initials[1, 0] = '*Error*'
                    optimized_initials[i+2, 0] = self.species[specie_index]
                    optimized_initials[0, nth_paramset-1] = str(nth_paramset)
                    optimized_initials[1, nth_paramset-1] = \
                        '{:8.3e}'.format(error)
                    optimized_initials[i+2, nth_paramset-1] = \
                        '{:8.3e}'.format(best_indiv[i+len(self.sp.idx_params)])
            with open(
                    self.model_path
                    + '/optimization_results/optimized_initals.csv', 'w') as f:
                writer = csv.writer(f, lineterminator='\n')
                writer.writerows(optimized_initials)
    
    def run_simulation(self, viz_type='average', show_all=False, stdev=True):
        warnings.filterwarnings('ignore')
        SignalingSystems(
            model_path=self.model_path,
            parameters=self.parameters, 
            species=self.species,
            pval=self.pval, 
            ival=self.ival, 
            obs=self.obs,
            sim=self.sim, 
            exp=self.exp, 
            sp=self.sp
        ).simulate_all(viz_type=viz_type, show_all=show_all, stdev=stdev)

    def optimize(self, *args):
        warnings.filterwarnings('ignore')
        ga_init = GeneticAlgorithmInit(
            model_path=self.model_path,
            sp=self.sp,
            obj_func=self.obj_func
        )
        if len(args) == 1:
            ga_init.run(int(args[0]))
        elif len(args) == 2:
            n_proc = max(1, multiprocessing.cpu_count() - 1)
            p = multiprocessing.Pool(processes=n_proc)
            p.map(ga_init.run, range(int(args[0]), int(args[1]) + 1))
            p.close()
        else:
            raise TypeError(
                "optimze() takes 1 or 2 arguments ({:d} given)".format(
                    len(args)
                )
            )

    def optimize_continue(self, *args):
        warnings.filterwarnings('ignore')
        ga_continue = GeneticAlgorithmContinue(
            model_path=self.model_path,
            sp=self.sp,
            obj_func=self.obj_func
        )
        if len(args) == 1:
            ga_continue.run(int(args[0]))
        elif len(args) == 2:
            n_proc = max(1, multiprocessing.cpu_count() - 1)
            p = multiprocessing.Pool(processes=n_proc)
            p.map(ga_continue.run, range(int(args[0]), int(args[1]) + 1))
            p.close()
        else:
            raise TypeError(
                "optimze_continue() takes 1 or 2 arguments ({:d} given)".format(
                    len(args)
                )
            )
    
    def analyze(self, target, metric='integral', style='barplot'):
        warnings.filterwarnings('ignore')
        if target == 'raction':
            reaction = ReactionSensitivity(
                model_path=self.model_path,
                reaction_system=self.reaction_system,
                obs=self.obs,
                sim=self.sim,
                sp=self.sp,
                rxn=self.rxn
            )
            reaction.analyze(metric=metric, style=metric)
        elif target == 'initial_condition':
            nonzero_init = NonZeroInitSensitivity(
                model_path=self.model_path,
                species=self.species,
                ival=self.ival,
                obs=self.obs,
                sim=self.sim,
                sp=self.sp
            )
            nonzero_init.analyze(metric=metric, style=metric)
        else:
            raise ValueError(
                "Available targets are: 'reaction', 'initial_condition'"
            )
    