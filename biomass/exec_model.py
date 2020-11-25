import os
import re
import numpy as np
from typing import Tuple

class ExecModel(object):
    def __init__(self, model):
        self.model_path = model.__path__[0]
        self.parameters = model.C.NAMES
        self.species = model.V.NAMES
        self.pval = model.param_values
        self.ival = model.initial_values
        self.obs = model.observables
        self.sim = model.NumericalSimulation()
        self.exp = model.ExperimentalData()
        self.viz = model.Visualization()
        self.rxn = model.ReactionNetwork()
        self.sp = model.SearchParam()
        self.obj_func = model.objective

        if self.sim.normalization:
            for obs_name in (self.obs):
                if isinstance(
                        self.sim.normalization[obs_name]['timepoint'], int
                        ) and self.sim.t[-1] < self.sim.normalization[obs_name]['timepoint']:
                    raise ValueError(
                        f'Normalization timepoint must be smaller than {self.sim.t[-1]:d}.'
                    )
                if not self.sim.normalization[obs_name]['condition']:
                    self.sim.normalization[obs_name]['condition'] = self.sim.conditions
                else:
                    for c in self.sim.normalization[obs_name]['condition']:
                        if not c in self.sim.conditions:
                            raise ValueError(
                                f'Normalization condition `{c}`'
                                ' is not defined in sim.conditions.'
                            )

    def show_properties(self):
        print(
            'Model properties\n'
            '----------------\n'
            f'{len(self.species):d} species\n'
            f'{len(self.parameters):d} parameters, '
            f'of which {len(self.sp.idx_params):d} to be estimated'
        )

    def get_indiv(self, paramset: int) -> np.ndarray:
        best_generation = np.load(
            self.model_path
            + f'/out/{paramset:d}/generation.npy'
        )
        best_indiv = np.load(
            self.model_path
            + f'/out/{paramset:d}/fit_param{int(best_generation):d}.npy'
        )
        return best_indiv

    def load_param(self, paramset: int) -> Tuple[list, list]:
        best_indiv = self.get_indiv(paramset)
        (x, y0) = self.sp.update(best_indiv)
        return x, y0

    def get_executable(self) -> list:
        n_file = []
        try:
            fitparam_files = os.listdir(self.model_path + '/out')
            for file in fitparam_files:
                if re.match(r'\d', file):
                    n_file.append(int(file))
            empty_folder = []
            for i, nth_paramset in enumerate(n_file):
                if not os.path.isfile(
                        self.model_path 
                        + f'/out/{nth_paramset:d}/generation.npy'):
                    empty_folder.append(i)
            for i in sorted(empty_folder, reverse=True):
                n_file.pop(i)
        except FileNotFoundError:
            pass
        
        return n_file