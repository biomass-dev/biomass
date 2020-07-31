import os
import numpy as np
import csv
import multiprocessing
import warnings

from biomass.exec_model import ExecModel
from biomass.dynamics import SignalingSystems
from biomass.ga import GeneticAlgorithmInit, GeneticAlgorithmContinue
from biomass.analysis import (ReactionSensitivity,
                              InitialConditionSensitivity,
                              ParameterSensitivity)


class OptimizationResults(ExecModel):
    def __init__(self, model):
        super().__init__(model)

    def get(self):
        """ Get optimized parameters as CSV file format

        Output
        ------
        optimization_results/optimized_params.csv
        optimization_results/optimized_initials.csv
        """
        os.makedirs(self.model_path + '/optimization_results/', exist_ok=True)
        n_file = self.get_executable()

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

    def dynamic_assessment(self, include_original=False):
        with open(self.model_path + '/fitness_assessment.csv', mode='w') as f:
            writer = csv.writer(f, lineterminator='\n')
            writer.writerow(['parameter set', 'Objective value'])
            if include_original:
                x = self.pval()
                y0 = self.ival()
                obj_val = self.obj_func(None, x, y0)
                writer.writerow(['original', '{:8.3e}'.format(obj_val)])
            n_file = self.get_executable()
            for paramset in sorted(n_file):
                (x, y0) = self.load_param(paramset)
                obj_val = self.obj_func(None, x, y0)
                writer.writerow(
                    ['{:d}'.format(paramset), '{:8.3e}'.format(obj_val)]
                )

def run_simulation(model, viz_type, show_all=False, stdev=False):
    warnings.filterwarnings('ignore')
    if not viz_type in ['best', 'average', 'original', 'experiment'] \
            and not viz_type.isdecimal():
        raise ValueError(
            "Avairable viz_type are: " \
            "'best','average','original','experiment','n(=1, 2, ...)'"
        )
    SignalingSystems(model).simulate_all(
        viz_type=viz_type, show_all=show_all, stdev=stdev
    )


def optimize(model, *args):
    warnings.filterwarnings('ignore')
    ga_init = GeneticAlgorithmInit(
        model,
        max_generation=10000,
        allowable_error=0.5
    )
    if len(args) == 1:
        ga_init.run(int(args[0]))
    elif len(args) == 2:
        n_proc = max(1, multiprocessing.cpu_count() - 1)
        p = multiprocessing.Pool(processes=n_proc)
        p.map(ga_init.run, range(int(args[0]), int(args[1]) + 1))
        p.close()
    else:
        raise ValueError('too many values to unpack (expected 2)')


def optimize_continue(model, *args):
    warnings.filterwarnings('ignore')
    ga_continue = GeneticAlgorithmContinue(
        model,
        max_generation=10000,
        allowable_error=0.5,
        p0_bounds=[0.1, 10.]  # [lower_bound, upper_bound]
    )
    if len(args) == 1:
        ga_continue.run(int(args[0]))
    elif len(args) == 2:
        n_proc = max(1, multiprocessing.cpu_count() - 1)
        p = multiprocessing.Pool(processes=n_proc)
        p.map(ga_continue.run, range(int(args[0]), int(args[1]) + 1))
        p.close()
    else:
        raise ValueError('too many values to unpack (expected 2)')


def run_analysis(
        model,
        target,
        metric='integral',
        style='barplot',
        excluded_params=[],
):
    warnings.filterwarnings('ignore')
    if target == 'reaction':
        ReactionSensitivity(model).analyze(metric=metric, style=style)
    elif target == 'initial_condition':
        InitialConditionSensitivity(model).analyze(metric=metric, style=style)
    elif target == 'parameter':
        ParameterSensitivity(
            model, excluded_params
        ).analyze(metric=metric, style=style)
        """
        >>> from biomass.models import Nakakuki_Cell_2010
        >>> from biomass import run_analysis
        >>> run_analysis(
                Nakakuki_Cell_2010,
                target='parameter',
                excluded_params=[
                    'a', 'Vn', 'Vc', 'Ligand', 'EGF', 'HRG', 'no_ligand'
                ]
            )
        """
    else:
        raise ValueError(
            "Available targets are: 'reaction', 'initial_condition' , 'parameter'"
        )

