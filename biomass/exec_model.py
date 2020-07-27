import os
import re
import numpy as np

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

    def get_indiv(self, paramset):
        best_generation = np.load(
            self.model_path + '/out/{:d}/generation.npy'.format(
                paramset
            )
        )
        best_indiv = np.load(
            self.model_path + '/out/{:d}/fit_param{:d}.npy'.format(
                paramset, int(best_generation)
            )
        )
        return best_indiv

    def load_param(self, paramset):
        best_indiv = self.get_indiv(paramset)
        (x, y0) = self.sp.update(best_indiv)
        return x, y0

    def get_executable(self):
        n_file = []
        fitparam_files = os.listdir(self.model_path + '/out')
        for file in fitparam_files:
            if re.match(r'\d', file):
                n_file.append(int(file))
        empty_folder = []
        for i, nth_paramset in enumerate(n_file):
            if not os.path.isfile(
                    self.model_path 
                    + '/out/{:d}/generation.npy'.format(nth_paramset)):
                empty_folder.append(i)
        for i in sorted(empty_folder, reverse=True):
            n_file.pop(i)
        
        return n_file