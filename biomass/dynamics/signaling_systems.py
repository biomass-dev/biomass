import os
import re
import numpy as np

from . import plot_func


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
        if not viz_type in ['best', 'average', 'original', 'experiment'] and \
                not viz_type.isdecimal():
            raise ValueError(
                "Avairable viz_type are: " \
                "'best','average','original','experiment','n(=1, 2, ...)'"
            )
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
                    (dynamic, successful) = self._validate(nth_paramset)
                    if successful:
                        for j, _ in enumerate(self.obs):
                            simulations_all[j, i, :, :] = dynamic.simulations[j, :, :]
                best_fitness_all = np.full(len(n_file), np.inf)
                for i, nth_paramset in enumerate(n_file):
                    if os.path.isfile('./out/{:d}/best_fitness.npy'.format(nth_paramset)):
                        best_fitness_all[i] = np.load(
                            './out/{:d}/best_fitness.npy'.format(nth_paramset)
                        )
                best_paramset = n_file[np.argmin(best_fitness_all)]
                self._write_best_fit_param(best_paramset)
                if viz_type == 'average':
                    pass
                elif viz_type == 'best':
                    dynamic, _ = self._validate(int(best_paramset))
                else:
                    dynamic, _ = self._validate(int(viz_type))

                if 2 <= len(n_file):
                    popt = np.empty(
                        (len(n_file), len(search_idx[0]) + len(search_idx[1]))
                    )
                    for i, nth_paramset in enumerate(n_file):
                        popt[i, :] = _get_indiv(nth_paramset)
                    plot_func.param_range(
                        search_idx, popt, 
                        self.parameters, self.species, self.sp, portrait=True
                    )
            else:
                x = self.pval()
                y0 = self.ival()
                if self.sim.simulate(x, y0) is not None:
                    print('Simulation failed.')
                else:
                    dynamic = self.sim
        plot_func.timecourse(
            dynamic, n_file, viz_type, show_all, stdev,
            simulations_all, self.obs, self.exp
        )
    
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

    def _write_best_fit_param(self, best_paramset):
        (x, y0) = load_param(best_paramset, self.sp.update)
        
        with open('./out/best_fit_param.txt', mode='w') as f:
            f.write(
                '# param set: {:d}\n'.format(
                    best_paramset
                )
            )
            f.write(
                '\n### f_params\n'
            )
            for i, param in enumerate(self.parameters):
                f.write('x[C.{}] = {:8.3e}\n'.format(param, x[i]))
            f.write(
                '\n### initial_values\n'
            )
            for i, specie in enumerate(self.species):
                if y0[i] != 0:
                    f.write('y0[V.{}] = {:8.3e}\n'.format(specie, y0[i]))


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