import os
import numpy as np

from .temporal_dynamics import TemporalDynamics


class SignalingSystems(TemporalDynamics):
    def __init__(self, model):
        super().__init__(model)

    def simulate_all(self, viz_type, show_all, stdev):
        n_file = [] if viz_type in ['original', 'experiment'] \
            else self.get_executable()
        simulations_all = np.full(
            (len(self.obs), len(n_file), len(self.sim.t), len(self.sim.conditions)),
            np.nan
        )
        if viz_type == 'experiment':
            dynamic = self.sim
        else:
            if len(n_file) > 0:
                if len(n_file) == 1 and viz_type == 'average':
                    raise ValueError(
                        "viz_type should be 'best', not '{}'.".format(viz_type)
                    )
                    # viz_type = 'best'
                for i, nth_paramset in enumerate(n_file):
                    (dynamic, is_successful) = self._validate(nth_paramset)
                    if is_successful:
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
                    self.plot_param_range(popt, portrait=True)
            else:
                x = self.pval()
                y0 = self.ival()
                if self.sim.simulate(x, y0) is not None:
                    print('Simulation failed.\n')
                dynamic = self.sim            
        self.plot_timecourse(
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