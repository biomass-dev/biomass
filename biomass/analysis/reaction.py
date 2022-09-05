import os
import sys
from dataclasses import dataclass
from typing import Callable, Dict, Final, List, Optional, Union

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from ..model_object import ModelObject
from ..plotting import SensitivityOptions
from .util import SignalingMetric, dlnyi_dlnxj, remove_nan


@dataclass
class ReactionSensitivity(SignalingMetric):
    """Sensitivity for rate equations"""

    model: ModelObject
    create_metrics: Optional[Dict[str, Callable[[np.ndarray], Union[int, float]]]]

    def __post_init__(self) -> None:
        self._plotting: SensitivityOptions = self.model.viz.get_sensitivity_options()
        self._coefficients: Callable[[str], str] = lambda metric: os.path.join(
            self.model.path,
            "sensitivity_coefficients",
            "reaction",
            f"{metric}.npy",
        )
        self._path_to_figs: Callable[[str], str] = lambda metric: os.path.join(
            self.model.path,
            "figure",
            "sensitivity",
            "reaction",
            f"{metric}",
        )
        if self.create_metrics is not None:
            for name, function in self.create_metrics.items():
                self.quantification[name] = function

    def _calc_sensitivity_coefficients(
        self,
        metric: str,
        reaction_indices: List[int],
    ) -> np.ndarray:
        """Calculating Sensitivity Coefficients
        Parameters
        ----------
        metric : str
            The signaling metric used for sensitivity analysis.
        reaction_indices : list of int
            List of reaction indices.
        Returns
        -------
        sensitivity_coefficients : numpy array
        """
        rate: Final[float] = 1.01  # 1% change
        n_file = self.model.get_executable()
        signaling_metric = np.full(
            (
                len(n_file),
                len(reaction_indices) + 1,
                len(self.model.observables),
                len(self.model.problem.conditions),
            ),
            np.nan,
        )
        for i, nth_paramset in enumerate(n_file):
            optimized = self.model.load_param(nth_paramset)
            for j, rxn_idx in enumerate(reaction_indices):
                perturbation: Dict[int, float] = {}
                for idx in reaction_indices:
                    perturbation[idx] = 1.0
                perturbation[rxn_idx] = rate
                if (
                    self.model.problem.simulate(optimized.params, optimized.initials, perturbation)
                    is None
                ):
                    for k, _ in enumerate(self.model.observables):
                        for l, _ in enumerate(self.model.problem.conditions):
                            signaling_metric[i, j, k, l] = self.quantification[metric](
                                self.model.problem.simulations[k, l]
                            )
                sys.stdout.write(
                    "\r{:d} / {:d}".format(
                        i * len(reaction_indices) + j + 1,
                        len(n_file) * len(reaction_indices),
                    )
                )
            if self.model.problem.simulate(optimized.params, optimized.initials) is None:
                for k, _ in enumerate(self.model.observables):
                    for l, _ in enumerate(self.model.problem.conditions):
                        signaling_metric[i, -1, k, l] = self.quantification[metric](
                            self.model.problem.simulations[k, l]
                        )
        sensitivity_coefficients = dlnyi_dlnxj(
            signaling_metric,
            len(n_file),
            len(reaction_indices),
            len(self.model.observables),
            len(self.model.problem.conditions),
            rate,
        )

        return sensitivity_coefficients

    def _load_sc(
        self,
        metric: str,
        reaction_indices: List[int],
    ) -> np.ndarray:
        """
        Load (or calculate) sensitivity coefficients.
        """
        if not os.path.isfile(self._coefficients(metric)):
            os.makedirs(
                os.path.join(
                    self.model.path,
                    "sensitivity_coefficients",
                    "reaction",
                ),
                exist_ok=True,
            )
            sensitivity_coefficients = self._calc_sensitivity_coefficients(
                metric, reaction_indices
            )
            np.save(self._coefficients(metric), sensitivity_coefficients)
        else:
            sensitivity_coefficients = np.load(self._coefficients(metric))

        return sensitivity_coefficients

    @staticmethod
    def _draw_vertical_span(
        biological_processes: List[List[int]],
        width: float,
    ) -> None:
        """
        Draw vertical span separating biological processes.
        """
        if len(biological_processes) > 1:
            left_end = 0
            for i, proc in enumerate(biological_processes):
                if i % 2 == 0:
                    plt.axvspan(
                        left_end - width,
                        left_end - width + len(proc),
                        facecolor="k",
                        alpha=0.1,
                    )
                left_end += len(proc)

    def _write_reaction_indices(
        self,
        reaction_indices: List[int],
        average: np.ndarray,
        stdev: np.ndarray,
        width: float,
    ) -> None:
        """
        Put reaction index on each bar.
        """
        distance = np.max(average) * 0.05
        for i, j in enumerate(reaction_indices):
            xp = i + width * 0.5 * (len(self.model.problem.conditions) - 1)
            yp = average[i, np.argmax(np.abs(average[i, :]))]
            yerr = stdev[i, np.argmax(stdev[i, :])]
            if yp > 0:
                plt.text(
                    xp,
                    yp + yerr + distance,
                    str(j),
                    ha="center",
                    va="bottom",
                    fontsize=10,
                    rotation=90,
                )
            else:
                plt.text(
                    xp,
                    yp - yerr - distance,
                    str(j),
                    ha="center",
                    va="top",
                    fontsize=10,
                    rotation=90,
                )

    def _barplot_sensitivity(
        self,
        metric: str,
        sensitivity_coefficients: np.ndarray,
        biological_processes: List[List[int]],
        reaction_indices: List[int],
        show_indices: bool,
    ) -> None:
        """
        Visualize sensitivity coefficients using barplot.
        """
        os.makedirs(os.path.join(self._path_to_figs(metric), "barplot"), exist_ok=True)
        # rcParams
        self.model.viz.set_sensitivity_rcParams()

        for k, obs_name in enumerate(self.model.observables):
            plt.figure(figsize=self._plotting.figsize)
            self._draw_vertical_span(biological_processes, self._plotting.width)

            sensitivity_array = sensitivity_coefficients[:, :, k, :]
            # Remove NaN
            nan_idx = []
            for i in range(sensitivity_array.shape[0]):
                for j in range(sensitivity_array.shape[1]):
                    if np.isnan(sensitivity_array[i, j, :]).any():
                        nan_idx.append(i)
            sensitivity_array = np.delete(sensitivity_array, nan_idx, axis=0)
            if sensitivity_array.size != 0:
                average = np.mean(sensitivity_array, axis=0)
                if sensitivity_array.shape[0] == 1:
                    stdev = np.zeros((sensitivity_array.shape[1], sensitivity_array.shape[2]))
                else:
                    stdev = np.std(sensitivity_array, axis=0, ddof=1)
                for l, condition in enumerate(self.model.problem.conditions):
                    plt.bar(
                        np.arange(len(reaction_indices)) + l * self._plotting.width,
                        average[:, l],
                        yerr=stdev[:, l],
                        ecolor=self._plotting.cmap[l],
                        capsize=2,
                        width=self._plotting.width,
                        color=self._plotting.cmap[l],
                        align="center",
                        label=condition,
                    )
                if show_indices:
                    self._write_reaction_indices(
                        reaction_indices, average, stdev, self._plotting.width
                    )
                plt.hlines([0], -self._plotting.width, len(reaction_indices), "k", lw=1)
                plt.xticks([])
                plt.ylabel(
                    "Control coefficients on\n" + metric + " (" + obs_name.replace("_", " ") + ")"
                )
                plt.xlim(-self._plotting.width, len(reaction_indices))
                if self._plotting.legend_kws is not None:
                    plt.legend(**self._plotting.legend_kws)
                plt.savefig(
                    os.path.join(
                        self._path_to_figs(metric),
                        "barplot",
                        f"{obs_name}",
                    ),
                )
                plt.close()

    def _heatmap_sensitivity(
        self,
        metric: str,
        sensitivity_coefficients: np.ndarray,
        reaction_indices: List[int],
    ) -> None:
        """
        Visualize sensitivity coefficients using heatmap.
        """
        os.makedirs(os.path.join(self._path_to_figs(metric), "heatmap"), exist_ok=True)
        # rcParams
        self.model.viz.set_sensitivity_rcParams()

        for k, obs_name in enumerate(self.model.observables):
            for l, condition in enumerate(self.model.problem.conditions):
                sensitivity_matrix = remove_nan(sensitivity_coefficients[:, :, k, l])
                if sensitivity_matrix.shape[0] > 1 and not np.all(sensitivity_matrix == 0.0):
                    g = sns.clustermap(
                        data=sensitivity_matrix,
                        center=0,
                        robust=True,
                        method="ward",
                        cmap="RdBu_r",
                        linewidth=0.5,
                        col_cluster=False,
                        figsize=self._plotting.figsize,
                        xticklabels=[str(j) for j in reaction_indices],
                        yticklabels=[],
                        # cbar_kws={"ticks": [-1, 0, 1]}
                    )
                    cbar = g.ax_heatmap.collections[0].colorbar
                    cbar.ax.tick_params(labelsize=8)
                    plt.savefig(
                        os.path.join(
                            self._path_to_figs(metric),
                            "heatmap",
                            f"{condition.replace('.', '')}_{obs_name}",
                        ),
                    )
                    plt.close()

    def analyze(self, *, metric: str, style: str, options: dict) -> None:
        """
        Perform sensitivity analysis.
        """
        if options["overwrite"] and os.path.isfile(self._coefficients(metric)):
            os.remove(self._coefficients(metric))
        if not self.model.rxn.reactions:
            raise ValueError("Define reaction indices (reactions) in reaction_network.py")
        biological_processes = self._group()
        reaction_indices = np.sum(biological_processes, axis=0)
        sensitivity_coefficients = self._load_sc(metric, reaction_indices)

        if style == "barplot":
            self._barplot_sensitivity(
                metric,
                sensitivity_coefficients,
                biological_processes,
                reaction_indices,
                options["show_indices"],
            )
        elif style == "heatmap":
            self._heatmap_sensitivity(
                metric,
                sensitivity_coefficients,
                reaction_indices,
            )
        else:
            raise ValueError("Available styles are: 'barplot', 'heatmap'")

    @staticmethod
    def _is_duplicate(
        reactions: Dict[str, List[int]],
        biological_processes: List[List[int]],
    ) -> bool:
        reaction_indices = (
            sum(biological_processes, []) if len(reactions) > 1 else biological_processes[0]
        )
        duplicate_reaction = [i for i in set(reaction_indices) if reaction_indices.count(i) > 1]
        if not duplicate_reaction:
            return False
        else:
            which_process = []
            for reaction_index in duplicate_reaction:
                for process, indices in reactions.items():
                    if reaction_index in indices:
                        which_process.append(process)
            raise ValueError(
                "Duplicate reaction: {} found in {}.".format(
                    ", ".join(map(str, duplicate_reaction)),
                    ", ".join(which_process),
                )
            )

    def _group(self) -> list:
        """
        Group reactions according to biological processes.
        """
        biological_processes = []
        for process, indices in self.model.rxn.reactions.items():
            if not isinstance(indices, list):
                raise TypeError("Use list for reaction indices in {}".format(process))
            biological_processes.append(indices)

        if not self._is_duplicate(self.model.rxn.reactions, biological_processes):
            return biological_processes
        assert False
