import os
import numpy as np

from ..exec_model import BioMassModel
from .temporal_dynamics import TemporalDynamics


class SignalingSystems(TemporalDynamics):
    def __init__(self, model: BioMassModel) -> None:
        super().__init__(model)

    def simulate_all(
        self,
        viz_type: str,
        show_all: bool,
        stdev: bool,
        save_format: str,
    ) -> None:
        n_file = [] if viz_type in ["original", "experiment"] else self.get_executable()
        simulations_all = np.full(
            (
                len(self.model.obs),
                len(n_file),
                len(self.model.sim.t),
                len(self.model.sim.conditions),
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
                for i, nth_paramset in enumerate(n_file):
                    if self._validate(nth_paramset):
                        for j, _ in enumerate(self.model.obs):
                            simulations_all[j, i, :, :] = self.model.sim.simulations[j, :, :]
                """
                simulations_all : numpy array
                    All simulated values with estimated parameter sets.
                """
                np.save(
                    os.path.join(
                        self.model.path,
                        "simulation_data",
                        "simulations_all.npy",
                    ),
                    simulations_all,
                )
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
                best_paramset = n_file[np.argmin(best_fitness_all)]
                self._write_best_fit_param(best_paramset)
                if viz_type == "average":
                    pass
                elif viz_type == "best":
                    is_successful = self._validate(int(best_paramset))
                    if is_successful:
                        np.save(
                            os.path.join(
                                self.model.path,
                                "simulation_data",
                                "simulations_best.npy",
                            ),
                            self.model.sim.simulations,
                        )
                else:  # viz_type == 'n(=1,2,...)'
                    is_successful = self._validate(int(viz_type))
                    if is_successful:
                        np.save(
                            os.path.join(
                                self.model.path,
                                "simulation_data",
                                f"simulations_{int(viz_type):d}.npy",
                            ),
                            self.model.sim.simulations,
                        )
                """Visualization of estimated parameter values"""
                if 2 <= len(n_file):
                    popt = np.empty(
                        (
                            len(n_file),
                            len(self.model.sp.idx_params) + len(self.model.sp.idx_initials),
                        )
                    )
                    for i, nth_paramset in enumerate(n_file):
                        popt[i, :] = self.get_individual(nth_paramset)
                    self.plot_param_range(popt, save_format, portrait=True)
            else:  # viz_type == 'original'
                x = self.model.pval()
                y0 = self.model.ival()
                if self.model.sim.simulate(x, y0) is not None:
                    print("Simulation failed.\n")
                # dynamic = self.model.sim
                """
                simulations_original : numpy array
                    Simulated values with original parameter values.
                """
                np.save(
                    os.path.join(
                        self.model.path,
                        "simulation_data",
                        "simulations_original.npy",
                    ),
                    self.model.sim.simulations,
                )

        self.plot_timecourse(n_file, viz_type, show_all, stdev, save_format, simulations_all)

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
        (x, y0) = self.load_param(nth_paramset)
        if self.model.sim.simulate(x, y0) is None:
            is_successful = True
        else:
            print(f"Simulation failed. #{nth_paramset:d}\n")
            is_successful = False

        return is_successful

    def _write_best_fit_param(self, best_paramset: int) -> None:
        """
        Create best_fit_param.txt in out/.

        Parameters
        ----------
        best_paramset : int
            Index of parameter set with the best objective value.

        """
        (x, y0) = self.load_param(best_paramset)

        with open(
            os.path.join(
                self.model.path,
                "out",
                "best_fit_param.txt",
            ),
            mode="w",
        ) as f:
            f.write(f"# param set: {best_paramset:d}\n")
            f.write("\n### param_values\n")
            for i, param in enumerate(self.model.parameters):
                f.write(f"x[C.{param}] = {x[i]:8.3e}\n")
            f.write("\n### initial_values\n")
            for i, specie in enumerate(self.model.species):
                if y0[i] != 0:
                    f.write(f"y0[V.{specie}] = {y0[i]:8.3e}\n")
