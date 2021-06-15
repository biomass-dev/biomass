import os
import sys
import warnings
from dataclasses import dataclass
from typing import List

import numpy as np

from ..exec_model import ModelObject
from .temporal_dynamics import TemporalDynamics


@dataclass
class SignalingSystems(TemporalDynamics):
    model: ModelObject

    def simulate_all(
        self,
        *,
        viz_type: str,
        show_all: bool,
        stdev: bool,
        save_format: str,
    ) -> None:
        """
        Run simulation and save figures.
        """
        n_file: List[int] = [] if viz_type in ["original", "experiment"] else self.get_executable()
        simulations_all = np.full(
            (
                len(self.model.observables),
                len(n_file),
                len(self.model.problem.t),
                len(self.model.problem.conditions),
            ),
            np.nan,
        )
        if viz_type != "experiment":
            os.makedirs(
                os.path.join(
                    self.model.path,
                    "simulation_data",
                ),
                exist_ok=True,
            )
            if len(n_file) > 0:
                if len(n_file) == 1 and viz_type == "average":
                    raise ValueError(f"viz_type should be 'best', not '{viz_type}'.")
                for j, nth_paramset in enumerate(n_file):
                    if self._validate(nth_paramset):
                        for i, _ in enumerate(self.model.observables):
                            simulations_all[i, j, :, :] = self.model.problem.simulations[i, :, :]
                # simulations_all : numpy array
                # All simulated values with estimated parameter sets.
                self._save_simulations(viz_type, simulations_all)
                best_fitness_all = self._get_best_objval(n_file)
                best_paramset = n_file[np.argmin(best_fitness_all)]
                self._write_best_fit_param(best_paramset)
                if viz_type == "average":
                    pass
                elif viz_type == "best":
                    is_successful = self._validate(int(best_paramset))
                    if is_successful:
                        self._save_simulations(viz_type, self.model.problem.simulations)
                else:
                    # viz_type == 'n(=1,2,...)'
                    is_successful = self._validate(int(viz_type))
                    if is_successful:
                        self._save_simulations(viz_type, self.model.problem.simulations)
            else:
                # viz_type == 'original'
                x = self.model.pval()
                y0 = self.model.ival()
                if self.model.problem.simulate(x, y0) is not None:
                    warnings.warn("Simulation failed. #original", RuntimeWarning)
                # simulations_original : numpy array
                # Simulated values with original parameter values.
                self._save_simulations(viz_type, self.model.problem.simulations)

        self.plot_timecourse(n_file, viz_type, show_all, stdev, save_format, simulations_all)

    def _preprocessing(self, simulated_values: np.ndarray) -> np.ndarray:
        """
        Replace small value in time-course simulated values to zero
        when all(abs(time_course) < sys.float_info.epsilon).

        Parameters
        ----------
        simulated_values : matrix (len(t) × len(self.model.problem.t))
        """
        if simulated_values.ndim == 2:
            for k, time_course in enumerate(simulated_values.T):
                if np.all(np.abs(time_course) < sys.float_info.epsilon):
                    simulated_values[:, k] = np.zeros(len(self.model.problem.t))
        elif simulated_values.ndim == 3:
            for i, _ in enumerate(self.model.observables):
                simulated_values[i] = self._preprocessing(simulated_values[i])
        elif simulated_values.ndim == 4:
            for j in range(simulated_values.shape[1]):
                for i, _ in enumerate(self.model.observables):
                    simulated_values[i, j] = self._preprocessing(simulated_values[i, j])
        return simulated_values

    def _save_simulations(self, viz_type: str, simulated_values: np.ndarray) -> None:
        """
        Save time course simulated values to simulation_data/.
        """
        np.save(
            os.path.join(
                self.model.path,
                "simulation_data",
                "simulations_{}.npy".format(viz_type if simulated_values.ndim < 4 else "all"),
            ),
            self._preprocessing(simulated_values),
        )

    def _get_best_objval(self, n_file: List[int]) -> np.ndarray:
        """
        Get best objective function values from parameter estimation results (out/).
        """
        best_fitness_all = np.full(len(n_file), np.inf)
        for i, nth_paramset in enumerate(n_file):
            if os.path.isfile(
                os.path.join(
                    self.model.path,
                    "out",
                    f"{nth_paramset:d}",
                    "best_fitness.npy",
                )
            ):
                best_fitness_all[i] = np.load(
                    os.path.join(
                        self.model.path,
                        "out",
                        f"{nth_paramset:d}",
                        "best_fitness.npy",
                    )
                )
        return best_fitness_all

    def _validate(self, nth_paramset: int) -> bool:
        """
        Validates the dynamical viability of a set of estimated parameter values.

        Parameters
        ----------
        nth_paramset : int
            Index of a parameter set.

        Returns
        -------
        is_successful : bool
            True if integration was successful.

        """
        optimized = self.load_param(nth_paramset)
        if self.model.problem.simulate(*optimized) is None:
            is_successful = True
        else:
            warnings.warn(f"Simulation failed. #{nth_paramset:d}", RuntimeWarning)
            is_successful = False

        return is_successful

    def _write_best_fit_param(
        self,
        best_paramset: int,
        indentation: str = " " * 4,
    ) -> None:
        """
        Create best_fit_param.txt in out/.

        Parameters
        ----------
        best_paramset : int
            Index of parameter set with the best objective value.

        """
        optimized = self.load_param(best_paramset)
        with open(
            os.path.join(
                self.model.path,
                "out",
                "best_fit_param.txt",
            ),
            mode="w",
        ) as f:
            f.write(f"# parameter set: {best_paramset:d}\n")
            f.write(f"\n\ndef param_values():\n{indentation}x = [0] * C.NUM\n")
            for i, name in enumerate(self.model.parameters):
                f.write(f"{indentation}x[C.{name}] = {optimized.params[i]:8.3e}\n")
            f.write(f"\n\ndef initial_values():\n{indentation}y0 = [0] * V.NUM\n")
            for i, name in enumerate(self.model.species):
                if optimized.initials[i] != 0:
                    f.write(f"{indentation}y0[V.{name}] = {optimized.initials[i]:8.3e}\n")
