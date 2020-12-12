import os
import numpy as np

from .temporal_dynamics import TemporalDynamics


class SignalingSystems(TemporalDynamics):
    def __init__(self, model):
        super().__init__(model)

    def simulate_all(self, viz_type: str, show_all: bool, stdev: bool):
        n_file = [] if viz_type in ["original", "experiment"] else self.get_executable()
        simulations_all = np.full((len(self.obs), len(n_file), len(self.sim.t), len(self.sim.conditions)), np.nan)
        if viz_type == "experiment":
            dynamic = self.sim
        else:
            os.makedirs(self.model_path + "/simulation_data", exist_ok=True)
            if len(n_file) > 0:
                if len(n_file) == 1 and viz_type == "average":
                    raise ValueError("viz_type should be 'best', not '{}'.".format(viz_type))
                for i, nth_paramset in enumerate(n_file):
                    (dynamic, is_successful) = self._validate(nth_paramset)
                    if is_successful:
                        for j, _ in enumerate(self.obs):
                            simulations_all[j, i, :, :] = dynamic.simulations[j, :, :]
                """
                simulations_all : numpy array
                    All simulated values with estimated parameter sets.
                """
                np.save(self.model_path + "/simulation_data/simulations_all.npy", simulations_all)
                best_fitness_all = np.full(len(n_file), np.inf)
                for i, nth_paramset in enumerate(n_file):
                    if os.path.isfile(self.model_path + "/out/{:d}/best_fitness.npy".format(nth_paramset)):
                        best_fitness_all[i] = np.load(
                            self.model_path + "/out/{:d}/best_fitness.npy".format(nth_paramset)
                        )
                best_paramset = n_file[np.argmin(best_fitness_all)]
                self._write_best_fit_param(best_paramset)
                if viz_type == "average":
                    pass
                elif viz_type == "best":
                    dynamic, _ = self._validate(int(best_paramset))
                else:  # viz_type == 'n(=1,2,...)'
                    dynamic, _ = self._validate(int(viz_type))
                """Visualization of estimated parameter values"""
                if 2 <= len(n_file):
                    popt = np.empty((len(n_file), len(self.sp.idx_params) + len(self.sp.idx_initials)))
                    for i, nth_paramset in enumerate(n_file):
                        popt[i, :] = self.get_individual(nth_paramset)
                    self.plot_param_range(popt, portrait=True)
            else:  # viz_type == 'original'
                x = self.pval()
                y0 = self.ival()
                if self.sim.simulate(x, y0) is not None:
                    print("Simulation failed.\n")
                dynamic = self.sim
                """
                simulations_original : numpy array
                    Simulated values with original parameter values.
                """
                np.save(self.model_path + "/simulation_data/simulations_original.npy", dynamic.simulations)
        self.plot_timecourse(dynamic, n_file, viz_type, show_all, stdev, simulations_all)

    def _validate(self, nth_paramset: int):
        """
        Validates the dynamical viability of a set of estimated parameter values.

        Parameters
        ----------
        nth_paramset : int
            Index of a parameter set.
        """
        (x, y0) = self.load_param(nth_paramset)
        if self.sim.simulate(x, y0) is None:
            return self.sim, True
        else:
            print("Simulation failed. #{:d}\n".format(nth_paramset))
            return self.sim, False

    def _write_best_fit_param(self, best_paramset: int):
        """
        Create best_fit_param.txt in out/.

        Parameters
        ----------
        best_paramset : int
            Index of parameter set with the best objective values.
        """
        (x, y0) = self.load_param(best_paramset)

        with open(self.model_path + "/out/best_fit_param.txt", mode="w") as f:
            f.write("# param set: {:d}\n".format(best_paramset))
            f.write("\n### param_values\n")
            for i, param in enumerate(self.parameters):
                f.write("x[C.{}] = {:8.3e}\n".format(param, x[i]))
            f.write("\n### initial_values\n")
            for i, specie in enumerate(self.species):
                if y0[i] != 0:
                    f.write("y0[V.{}] = {:8.3e}\n".format(specie, y0[i]))
