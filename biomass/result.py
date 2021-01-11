import os
import numpy as np
import csv

from .exec_model import BioMassModel, ExecModel


class OptimizationResults(ExecModel):
    def __init__(self, model: BioMassModel) -> None:
        super().__init__(model)

    def to_csv(self) -> None:
        """
        Save optimized parameters as CSV file format.

        Output
        ------
        optimization_results/optimized_params.csv
        optimization_results/optimized_initials.csv

        Example
        -------
        >>> from biomass.models import Nakakuki_Cell_2010
        >>> from biomass.result import OptimizationResults
        >>> model = Nakakuki_Cell_2010.create()
        >>> res = OptimizationResults(model)
        >>> res.to_csv()

        """
        os.makedirs(
            os.path.join(
                self.model.path,
                "optimization_results",
            ),
            exist_ok=True,
        )
        n_file = self.get_executable()

        if len(self.model.sp.idx_params) > 0:
            optimized_params = np.empty((len(self.model.sp.idx_params) + 2, len(n_file) + 1), dtype="<U21")
            for i, param_index in enumerate(self.model.sp.idx_params):
                for j, nth_paramset in enumerate(sorted(n_file), start=1):
                    best_generation = np.load(
                        os.path.join(
                            self.model.path,
                            "out",
                            f"{nth_paramset:d}",
                            "generation.npy",
                        )
                    )
                    best_individual = np.load(
                        os.path.join(
                            self.model.path,
                            "out",
                            f"{nth_paramset:d}",
                            f"fit_param{int(best_generation):d}.npy",
                        )
                    )
                    error = np.load(
                        os.path.join(
                            self.model.path,
                            "out",
                            f"{nth_paramset:d}",
                            "best_fitness.npy",
                        )
                    )
                    optimized_params[0, 0] = ""
                    optimized_params[1, 0] = "*Error*"
                    optimized_params[i + 2, 0] = self.model.parameters[param_index]
                    optimized_params[0, j] = str(nth_paramset)
                    optimized_params[1, j] = f"{error:8.3e}"
                    optimized_params[i + 2, j] = f"{best_individual[i]:8.3e}"
            with open(
                os.path.join(
                    self.model.path,
                    "optimization_results",
                    "optimized_params.csv",
                ),
                mode="w",
            ) as f:
                writer = csv.writer(f, lineterminator="\n")
                writer.writerows(optimized_params)
        if len(self.model.sp.idx_initials) > 0:
            optimized_initials = np.empty((len(self.model.sp.idx_initials) + 2, len(n_file) + 1), dtype="<U21")
            for i, specie_index in enumerate(self.model.sp.idx_initials):
                for j, nth_paramset in enumerate(sorted(n_file), start=1):
                    best_generation = np.load(
                        os.path.join(
                            self.model.path,
                            "out",
                            f"{nth_paramset:d}",
                            "generation.npy",
                        )
                    )
                    best_individual = np.load(
                        os.path.join(
                            self.model.path,
                            "out",
                            f"{nth_paramset:d}",
                            f"fit_param{int(best_generation):d}.npy",
                        )
                    )
                    error = np.load(
                        os.path.join(
                            self.model.path,
                            "out",
                            f"{nth_paramset:d}",
                            "best_fitness.npy",
                        )
                    )
                    optimized_initials[0, 0] = ""
                    optimized_initials[1, 0] = "*Error*"
                    optimized_initials[i + 2, 0] = self.model.species[specie_index]
                    optimized_initials[0, j] = str(nth_paramset)
                    optimized_initials[1, j] = f"{error:8.3e}"
                    optimized_initials[i + 2, j] = f"{best_individual[i+len(self.model.sp.idx_params)]:8.3e}"
            with open(
                os.path.join(
                    self.model.path,
                    "optimization_results",
                    "optimized_initials.csv",
                ),
                mode="w",
            ) as f:
                writer = csv.writer(f, lineterminator="\n")
                writer.writerows(optimized_initials)

    def dynamic_assessment(self, include_original: bool = False) -> None:
        """
        Compute objective values using estimated parameters.

        Parameters
        ----------
        include_original : bool (default: False)
            If True, an objective value simulated with original parameters
            will also be shown.

        Output
        ------
        fitness_assessment.csv

        Example
        -------
        >>> from biomass.models import Nakakuki_Cell_2010
        >>> from biomass.result import OptimizationResults
        >>> model = Nakakuki_Cell_2010.create()
        >>> res = OptimizationResults(model)
        >>> res.dynamic_assessment()

        """
        os.makedirs(
            os.path.join(
                self.model.path,
                "optimization_results",
            ),
            exist_ok=True,
        )
        with open(
            os.path.join(
                self.model.path,
                "optimization_results",
                "fitness_assessment.csv",
            ),
            mode="w",
        ) as f:
            writer = csv.writer(f, lineterminator="\n")
            writer.writerow(["parameter set", "Objective value"])
            if include_original:
                x = self.model.pval()
                y0 = self.model.ival()
                obj_val = self.model.obj_func(None, x, y0)
                writer.writerow(["original", f"{obj_val:8.3e}"])
            n_file = self.get_executable()
            for paramset in sorted(n_file):
                (x, y0) = self.load_param(paramset)
                obj_val = self.model.obj_func(None, x, y0)
                writer.writerow([f"{paramset:d}", f"{obj_val:8.3e}"])
