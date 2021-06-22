import csv
import os
from dataclasses import dataclass
from typing import Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .exec_model import ExecModel, ModelObject


@dataclass
class OptimizationResults(ExecModel):
    model: ModelObject

    def __post_init__(self) -> None:
        """
        Create optimization_results/ in the model folder.
        """
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

        Examples
        --------
        >>> from biomass import Model, OptimizationResults
        >>> from biomass.models import Nakakuki_Cell_2010
        >>> model = Model(Nakakuki_Cell_2010.__package__).create()
        >>> res = OptimizationResults(model)
        >>> res.to_csv()

        Notes
        -----
        Output:

        * optimization_results/optimized_params.csv
        * optimization_results/optimized_initials.csv
        """
        n_file = self.get_executable()

        if len(self.model.problem.idx_params) + len(self.model.problem.idx_initials) > 0:
            optimized_params = np.empty(
                (
                    len(self.model.problem.idx_params) + len(self.model.problem.idx_initials) + 2,
                    len(n_file) + 1,
                ),
                dtype="<U21",
            )
            for j, nth_paramset in enumerate(sorted(n_file), start=1):
                for i, parameter_index in enumerate(self.model.problem.idx_params):
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
                    optimized_params[0, j] = str(nth_paramset)
                    optimized_params[1, j] = f"{error:8.3e}"
                    optimized_params[i + 2, 0] = self.model.parameters[parameter_index]
                    optimized_params[i + 2, j] = f"{best_individual[i]:8.3e}"
                for i, species_index in enumerate(self.model.problem.idx_initials):
                    optimized_params[i + len(self.model.problem.idx_params) + 2, 0] = (
                        "init_" + self.model.species[species_index]
                    )
                    optimized_params[
                        i + len(self.model.problem.idx_params) + 2, j
                    ] = f"{best_individual[i+len(self.model.problem.idx_params)]:8.3e}"
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

    def savefig(
        self,
        *,
        figsize: Optional[Tuple[float, float]] = None,
        boxplot_kws: Optional[dict] = None,
    ) -> None:
        """
        Visualize estimated parameter sets using `seaborn.boxplot`.

        Parameters
        ----------
        figsize : Tuple[float, float], optional
            Width, height in inches.
        boxplot_kws : dict, optional
            Keyword arguments to pass to `seaborn.boxplot`.

        Examples
        --------
        >>> from biomass import Model, OptimizationResults
        >>> from biomass.models import Nakakuki_Cell_2010
        >>> model = Model(Nakakuki_Cell_2010.__package__).create()
        >>> res = OptimizationResults(model)
        >>> res.savefig(figsize=(16,5), boxplot_kws={"orient": "v"})

        Notes
        -----
        Output:

        * optimization_results/estimated_parameter_sets.pdf
        """
        if not os.path.isfile(
            os.path.join(
                self.model.path,
                "optimization_results",
                "optimized_params.csv",
            )
        ):
            self.to_csv()
        if boxplot_kws is None:
            boxplot_kws = {}
        boxplot_kws.setdefault("orient", "v")

        df = pd.read_csv(
            os.path.join(
                self.model.path,
                "optimization_results",
                "optimized_params.csv",
            ),
            index_col=0,
        )
        df.drop("*Error*", inplace=True)

        if isinstance(figsize, tuple) and len(figsize) == 2:
            plt.figure(figsize=figsize)
        ax = sns.boxplot(data=df.T, **boxplot_kws)
        if boxplot_kws["orient"] == "h":
            ax.set_xscale("log")
            ax.set_xlabel("Parameter value")
        elif boxplot_kws["orient"] == "v":
            ax.set_yscale("log")
            ax.set_ylabel("Parameter value")
            ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        else:
            raise ValueError("boxplot_kws['orient'] must be either 'v' or 'h'.")
        sns.despine()
        plt.savefig(
            os.path.join(
                self.model.path,
                "optimization_results",
                "estimated_parameter_sets.pdf",
            ),
            bbox_inches="tight",
        )

    def dynamic_assessment(self, include_original: bool = False) -> None:
        """
        Compute objective values using estimated parameters.

        Parameters
        ----------
        include_original : bool (default: False)
            If True, an objective value simulated with original parameters
            will also be shown.

        Examples
        --------
        >>> from biomass import Model, OptimizationResults
        >>> from biomass.models import Nakakuki_Cell_2010
        >>> model = Model(Nakakuki_Cell_2010.__package__).create()
        >>> res = OptimizationResults(model)
        >>> res.dynamic_assessment()

        Notes
        -----
        Output:

        * optimization_results/fitness_assessment.csv
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
                obj_val = self.model.problem.objective(None, x, y0)
                writer.writerow(["original", f"{obj_val:8.3e}"])
            n_file = self.get_executable()
            for paramset in sorted(n_file):
                optimized = self.load_param(paramset)
                obj_val = self.model.problem.objective(None, *optimized)
                writer.writerow([f"{paramset:d}", f"{obj_val:8.3e}"])

    def trace_obj(
        self,
        *,
        config: Optional[dict] = None,
        xlabel: str = "Iteration",
        ylabel: str = "Objective function value",
        xticks: Optional[list] = None,
        yticks: Optional[list] = None,
    ) -> None:
        """
        Visualize objective function traces for different optimization runs.

        Parameters
        ----------
        config : dict, optional
            A dictionary object for setting `matplotlib.rcParams`.
        xlabel: str (default: "Iteration")
            The label for the x-axis.
        ylabel: str (default: "Objective function value")
            The label for the x-axis.
        xticks: list, optional
            The list of xtick locations.
        yticks: list, optional
            The list of ytick locations.

        Examples
        --------
        >>> from biomass import Model, OptimizationResults
        >>> from biomass.models import Nakakuki_Cell_2010
        >>> model = Model(Nakakuki_Cell_2010.__package__).create()
        >>> res = OptimizationResults(model)
        >>> res.trace_obj()

        Notes
        -----
        Output:

        * optimization_results/obj_func_traces.pdf

        """
        n_file = self.get_executable()
        # matplotlib
        if config is None:
            config = {}
        config.setdefault("font.size", 15)
        config.setdefault("axes.linewidth", 1.5)
        config.setdefault("xtick.major.width", 1.5)
        config.setdefault("ytick.major.width", 1.5)
        config.setdefault("lines.linewidth", 1.5)
        plt.rcParams.update(config)
        plt.gca().spines["right"].set_visible(False)
        plt.gca().spines["top"].set_visible(False)
        # ---
        for paramset in n_file:
            with open(
                os.path.join(self.model.path, "out", f"{paramset:d}", "optimization.log"),
                mode="r",
            ) as f:
                traces = f.readlines()
            iters = []
            obj_val = []
            for line in traces:
                if line.startswith("Generation"):
                    iters.append(line.lstrip("Generation").split(":")[0])
                    obj_val.append(line.split("=")[-1].strip())
            plt.plot([int(num) - 1 for num in iters], [float(val) for val in obj_val])
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.xticks(xticks)
        plt.yticks(yticks)
        plt.savefig(
            os.path.join(self.model.path, "optimization_results", "obj_func_traces.pdf"),
            bbox_inches="tight",
        )
        plt.close()
