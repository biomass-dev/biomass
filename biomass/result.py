import csv
import os
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np

from .exec_model import ExecModel, ModelObject


@dataclass
class OptimizationResults(ExecModel):
    model: ModelObject

    def __post_init__(self) -> None:
        os.makedirs(
            os.path.join(
                self.model.path,
                "optimization_results",
            ),
            exist_ok=True,
        )

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
        n_file = self.get_executable()

        if len(self.model.sp.idx_params) > 0:
            optimized_params = np.empty(
                (len(self.model.sp.idx_params) + 2, len(n_file) + 1), dtype="<U21"
            )
            for i, parameter_index in enumerate(self.model.sp.idx_params):
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
                    optimized_params[i + 2, 0] = self.model.parameters[parameter_index]
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
            optimized_initials = np.empty(
                (len(self.model.sp.idx_initials) + 2, len(n_file) + 1), dtype="<U21"
            )
            for i, species_index in enumerate(self.model.sp.idx_initials):
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
                    optimized_initials[i + 2, 0] = self.model.species[species_index]
                    optimized_initials[0, j] = str(nth_paramset)
                    optimized_initials[1, j] = f"{error:8.3e}"
                    optimized_initials[
                        i + 2, j
                    ] = f"{best_individual[i+len(self.model.sp.idx_params)]:8.3e}"
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
                optimized = self.load_param(paramset)
                obj_val = self.model.obj_func(None, *optimized)
                writer.writerow([f"{paramset:d}", f"{obj_val:8.3e}"])

    def trace_obj(self) -> None:
        """
        Visualize objective function traces for different optimization runs.

        Output
        ------
        obj_func_traces.pdf

        Example
        -------
        >>> from biomass.models import Nakakuki_Cell_2010
        >>> from biomass.result import OptimizationResults
        >>> model = Nakakuki_Cell_2010.create()
        >>> res = OptimizationResults(model)
        >>> res.trace_obj()

        """
        n_file = self.get_executable()
        # matplotlib
        plt.figure(figsize=(4, 3))
        plt.gca().spines["right"].set_visible(False)
        plt.gca().spines["top"].set_visible(False)
        plt.rcParams["font.size"] = 18
        plt.rcParams["lines.linewidth"] = 1.0
        # ---
        for paramset in n_file:
            with open(
                os.path.join(
                    self.model.path,
                    "out",
                    f"{paramset:d}",
                    "optimization.log",
                ),
                mode="r",
            ) as f:
                traces = f.readlines()
            iters = []
            obj_val = []
            for line in traces:
                if line.startswith("Generation"):
                    iters.append(line.lstrip("Generation").split(":")[0])
                    obj_val.append(line.split("=")[-1].strip())
            plt.plot(
                [int(num) for num in iters],
                [float(val) for val in obj_val],
            )
        plt.xlabel("Iteration")
        plt.ylabel("Objective function value")
        plt.savefig(
            os.path.join(
                self.model.path,
                "optimization_results",
                "obj_func_traces.pdf",
            ),
            bbox_inches="tight",
        )
