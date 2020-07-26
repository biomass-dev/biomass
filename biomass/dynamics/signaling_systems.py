import os
import numpy as np

from biomass.exec_model import ExecModel
from .plot_func import PlotFunc


class SignalingSystems(ExecModel):
    def __init__(self, model):
        super().__init__(model)
        self.model = model
        self.model_path = model.__path__[0]
        self.parameters = model.C.NAMES
        self.species = model.V.NAMES
        self.pval = model.param_values
        self.ival = model.initial_values
        self.obs = model.observables
        self.sim = model.NumericalSimulation()
        self.sp = model.SearchParam()

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
        n_file = [] if viz_type == 'original' else self.get_executable()
        simulations_all = np.full(
            (len(self.obs), len(n_file), len(self.sim.t), len(self.sim.conditions)),
            np.nan
        )
        if viz_type == 'experiment':
            dynamic = self.sim
        else:
            if len(n_file) > 0:
                if len(n_file) == 1 and viz_type == 'average':
                    viz_type = 'best'
                for i, nth_paramset in enumerate(n_file):
                    (dynamic, successful) = self._validate(nth_paramset)
                    if successful:
                        for j, _ in enumerate(self.obs):
                            simulations_all[j, i, :, :] = \
                                dynamic.simulations[j, :, :]
                best_fitness_all = np.full(len(n_file), np.inf)
                for i, nth_paramset in enumerate(n_file):
                    if os.path.isfile(
                            self.model_path
                            + '/out/{:d}/best_fitness.npy'.format(nth_paramset)):
                        best_fitness_all[i] = np.load(
                            self.model_path
                            + '/out/{:d}/best_fitness.npy'.format(nth_paramset)
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
                    popt = np.empty((len(n_file), 
                        len(self.sp.idx_params) + len(self.sp.idx_initials))
                    )
                    for i, nth_paramset in enumerate(n_file):
                        popt[i, :] = self.get_indiv(nth_paramset)
                    PlotFunc(self.model).param_range(popt, portrait=True)
            else:
                x = self.pval()
                y0 = self.ival()
                if self.sim.simulate(x, y0) is not None:
                    print('Simulation failed.\n')
                dynamic = self.sim            
        PlotFunc(self.model).timecourse(
            dynamic, n_file, viz_type, show_all, stdev, simulations_all
        )
    
    def _validate(self, nth_paramset):
        """
        Validates the dynamical viability of a set of estimated parameter values.
        """
        (x, y0) = self.load_param(nth_paramset)
        if self.sim.simulate(x, y0) is None:
            return self.sim, True
        else:
            print(
                'Simulation failed. #{:d}\n'.format(nth_paramset)
            )
            return self.sim, False

    def _write_best_fit_param(self, best_paramset):
        (x, y0) = self.load_param(best_paramset)
        
        with open(self.model_path + '/out/best_fit_param.txt', mode='w') as f:
            f.write('# param set: {:d}\n'.format(best_paramset))
            f.write(
                '\n### param_values\n'
            )
            for i, param in enumerate(self.parameters):
                f.write('x[C.{}] = {:8.3e}\n'.format(param, x[i]))
            f.write(
                '\n### initial_values\n'
            )
            for i, specie in enumerate(self.species):
                if y0[i] != 0:
                    f.write('y0[V.{}] = {:8.3e}\n'.format(specie, y0[i]))