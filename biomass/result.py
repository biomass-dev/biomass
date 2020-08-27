import os
import numpy as np
import csv

from biomass.exec_model import ExecModel


class OptimizationResults(ExecModel):
    def __init__(self, model):
        super().__init__(model)

    def get(self):
        """
        Get optimized parameters as CSV file format.

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
                for nth_paramset in n_file:
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
                    optimized_params[0, nth_paramset] = str(nth_paramset)
                    optimized_params[1, nth_paramset] = \
                        '{:8.3e}'.format(error)
                    optimized_params[i+2, nth_paramset] = \
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
                for nth_paramset in n_file:
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
                    optimized_initials[0, nth_paramset] = str(nth_paramset)
                    optimized_initials[1, nth_paramset] = \
                        '{:8.3e}'.format(error)
                    optimized_initials[i+2, nth_paramset] = \
                        '{:8.3e}'.format(best_indiv[i+len(self.sp.idx_params)])
            with open(
                    self.model_path
                    + '/optimization_results/optimized_initals.csv', 'w') as f:
                writer = csv.writer(f, lineterminator='\n')
                writer.writerows(optimized_initials)

    def dynamic_assessment(self, include_original=False):
        """ 
        Get objective values using estimated parameters.
        """
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