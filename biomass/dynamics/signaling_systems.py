import os
import re
import numpy as np

from . import plot_func
from .load_out import *


class SignalingSystems(object):
    def __init__(self, parameters, species, pval, ival, obs, sim, exp, sp):
        self.parameters = parameters
        self.species = species
        self.pval = pval
        self.ival = ival
        self.obs = obs
        self.sim = sim
        self.exp = exp
        self.sp = sp
    
    def _validate(self, nth_paramset):
        """
        Validates the dynamical viability of a set of estimated parameter values.
        """
        (x, y0) = load_param(nth_paramset, self.sp.update)
        if self.sim.simulate(x, y0) is None:
            return self.sim, True
        else:
            print(
                'Simulation failed. #{:d}\n'.format(nth_paramset)
            )
            return self.sim, False


    def simulate_all(self, viz_type, show_all, stdev):
        """Simulate ODE model with estimated parameter values.

        Parameters
        ----------
        viz_type : str
            - 'average':
                The average of simulation results with parameter sets in "out/".
            - 'best': 
                The best simulation result in "out/", simulation with 
                "best_fit_param".
            - 'original': 
                Simulation with the default parameters and initial values 
                defined in "set_model.py".
            - 'n(=1,2,...)':
                Use the parameter set in "out/n/".

        show_all : bool
            Whether to show all simulation results.
            
        stdev : bool
            If True, the standard deviation of simulated values will be shown
            (only available for 'average' visualization type).
            
        """

        search_idx = (self.sp.idx_params, self.sp.idx_initials)
        n_file = [] if viz_type == 'original' else get_executable()
        simulations_all = np.full(
            (len(self.obs), len(n_file), len(self.sim.t), len(self.sim.conditions)),
            np.nan
        )
        if viz_type != 'experiment':
            if len(n_file) > 0:
                if len(n_file) == 1 and viz_type == 'average':
                    viz_type = 'best'
                for i, nth_paramset in enumerate(n_file):
                    (sim, successful) = self._validate(nth_paramset)
                    if successful:
                        for j, _ in enumerate(self.obs):
                            simulations_all[j, i, :, :] = sim.simulations[j, :, :]
                best_fitness_all = np.full(len(n_file), np.inf)
                for i, nth_paramset in enumerate(n_file):
                    if os.path.isfile('./out/{:d}/best_fitness.npy'.format(nth_paramset)):
                        best_fitness_all[i] = np.load(
                            './out/{:d}/best_fitness.npy'.format(nth_paramset)
                        )
                best_paramset = n_file[np.argmin(best_fitness_all)]
                write_best_fit_param(
                    best_paramset,
                    self.parameters, self.species, self.sp.update
                )

                if viz_type == 'average':
                    pass
                elif viz_type == 'best':
                    sim, _ = self._validate(int(best_paramset))
                else:
                    sim, _ = self._validate(int(viz_type))

                if 2 <= len(n_file):
                    popt = get_optimized_param(n_file, search_idx)
                    plot_func.param_range(
                        search_idx, popt, 
                        self.parameters, self.species, self.sp, portrait=True
                    )
            else:
                x = self.pval()
                y0 = self.ival()
                if self.sim.simulate(x, y0) is not None:
                    print('Simulation failed.')
        plot_func.timecourse(
            sim, n_file, viz_type, show_all, stdev,
            simulations_all, self.obs, self.exp
        )