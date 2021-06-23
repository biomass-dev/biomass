import os
import sys
from dataclasses import dataclass
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from ..exec_model import ExecModel, ModelObject
from .util import dlnyi_dlnxj, get_signaling_metric


@dataclass
class ReactionSensitivity(ExecModel):
    """Sensitivity for rate equations"""

    model: ModelObject

    def _calc_sensitivity_coefficients(
        self,
        metric: str,
        reaction_indices: List[int],
        options: dict,
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
        rate = 1.01  # 1% change
        n_file = self.get_executable()
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
            optimized = self.load_param(nth_paramset)
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
                            signaling_metric[i, j, k, l] = get_signaling_metric(
                                metric, self.model.problem.simulations[k, :, l], options
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
                        signaling_metric[i, -1, k, l] = get_signaling_metric(
                            metric, self.model.problem.simulations[k, :, l], options
                        )
        sensitivity_coefficients = dlnyi_dlnxj(
            signaling_metric,
            n_file,
            reaction_indices,
            self.model.observables,
            self.model.problem.conditions,
            rate,
        )

        return sensitivity_coefficients

    def _load_sc(
        self,
        metric: str,
        reaction_indices: List[int],
        options: dict,
    ) -> np.ndarray:
        """
        Load (or calculate) sensitivity coefficients.
        """
        if not os.path.isfile(
            os.path.join(
                self.model.path,
                "sensitivity_coefficients",
                "reaction",
                f"{metric}.npy",
            )
        ):
            os.makedirs(
                os.path.join(
                    self.model.path,
                    "sensitivity_coefficients",
                    "reaction",
                ),
                exist_ok=True,
            )
            sensitivity_coefficients = self._calc_sensitivity_coefficients(
                metric,
                reaction_indices,
                options,
            )
            np.save(
                os.path.join(
                    self.model.path,
                    "sensitivity_coefficients",
                    "reaction",
                    f"{metric}",
                ),
                sensitivity_coefficients,
            )
        else:
            sensitivity_coefficients = np.load(
                os.path.join(
                    self.model.path,
                    "sensitivity_coefficients",
                    "reaction",
                    f"{metric}.npy",
                )
            )

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
        save_format: str,
        show_indices: bool,
    ) -> None:
        """
        Visualize sensitivity coefficients using barplot.
        """
        os.makedirs(
            os.path.join(
                self.model.path,
                "figure",
                "sensitivity",
                "reaction",
                f"{metric}",
                "barplot",
            ),
            exist_ok=True,
        )
        options = self.model.viz.sensitivity_options

        # rcParams
        self.model.viz.set_sensitivity_rcParams()

        if len(options["cmap"]) < len(self.model.problem.conditions):
            raise ValueError(
                "len(sensitivity_options['cmap']) must be equal to "
                "or greater than len(problem.conditions)."
            )
        for k, obs_name in enumerate(self.model.observables):
            plt.figure(figsize=options["figsize"])
            self._draw_vertical_span(biological_processes, options["width"])

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
                        np.arange(len(reaction_indices)) + l * options["width"],
                        average[:, l],
                        yerr=stdev[:, l],
                        ecolor=options["cmap"][l],
                        capsize=2,
                        width=options["width"],
                        color=options["cmap"][l],
                        align="center",
                        label=condition,
                    )
                if show_indices:
                    self._write_reaction_indices(
                        reaction_indices, average, stdev, options["width"]
                    )
                plt.hlines([0], -options["width"], len(reaction_indices), "k", lw=1)
                plt.xticks([])
                plt.ylabel(
                    "Control coefficients on\n" + metric + " (" + obs_name.replace("_", " ") + ")"
                )
                plt.xlim(-options["width"], len(reaction_indices))
                plt.legend(loc=options["legend_loc"], frameon=False)
                plt.savefig(
                    os.path.join(
                        self.model.path,
                        "figure",
                        "sensitivity",
                        "reaction",
                        f"{metric}",
                        "barplot",
                        f"{obs_name}.{save_format}",
                    ),
                    dpi=600 if save_format == "png" else None,
                    bbox_inches="tight",
                )
                plt.close()

    @staticmethod
    def _remove_nan(sensitivity_matrix: np.ndarray, normalize: bool) -> np.ndarray:
        """
        Remove NaN from sensitivity matrix.
        """
        nan_idx = []
        for i in range(sensitivity_matrix.shape[0]):
            if np.isnan(sensitivity_matrix[i, :]).any():
                nan_idx.append(i)
            else:
                pass
            if np.nanmax(np.abs(sensitivity_matrix[i, :])) == 0.0:
                sensitivity_matrix[i, :] = np.zeros(sensitivity_matrix.shape[1])
            else:
                sensitivity_matrix[i, :] = sensitivity_matrix[i, :] / (
                    np.nanmax(np.abs(sensitivity_matrix[i, :])) if normalize else 1
                )

        return np.delete(sensitivity_matrix, nan_idx, axis=0)

    def _heatmap_sensitivity(
        self,
        metric: str,
        sensitivity_coefficients: np.ndarray,
        reaction_indices: List[int],
        save_format: str,
    ) -> None:
        """
        Visualize sensitivity coefficients using heatmap.
        """
        os.makedirs(
            os.path.join(
                self.model.path,
                "figure",
                "sensitivity",
                "reaction",
                f"{metric}",
                "heatmap",
            ),
            exist_ok=True,
        )
        options = self.model.viz.sensitivity_options
        # rcParams
        self.model.viz.set_sensitivity_rcParams()

        for k, obs_name in enumerate(self.model.observables):
            for l, condition in enumerate(self.model.problem.conditions):
                sensitivity_matrix = self._remove_nan(
                    sensitivity_coefficients[:, :, k, l], normalize=False
                )
                if sensitivity_matrix.shape[0] > 1 and not np.all(sensitivity_matrix == 0.0):
                    g = sns.clustermap(
                        data=sensitivity_matrix,
                        center=0,
                        robust=True,
                        method="ward",
                        cmap="RdBu_r",
                        linewidth=0.5,
                        col_cluster=False,
                        figsize=options["figsize"],
                        xticklabels=[str(j) for j in reaction_indices],
                        yticklabels=[],
                        # cbar_kws={"ticks": [-1, 0, 1]}
                    )
                    cbar = g.ax_heatmap.collections[0].colorbar
                    cbar.ax.tick_params(labelsize=8)
                    plt.savefig(
                        os.path.join(
                            self.model.path,
                            "figure",
                            "sensitivity",
                            "reaction",
                            f"{metric}",
                            "heatmap",
                            f"{condition}_{obs_name}.{save_format}",
                        ),
                        dpi=600 if save_format == "png" else None,
                        bbox_inches="tight",
                    )
                    plt.close()

    def analyze(self, *, metric: str, style: str, save_format: str, options: dict) -> None:
        """
        Perform sensitivity analysis.
        """
        if not self.model.rxn.reactions:
            raise ValueError("Define reaction indices (reactions) in reaction_network.py")
        biological_processes = self.model.rxn.group()
        reaction_indices = np.sum(biological_processes, axis=0)
        sensitivity_coefficients = self._load_sc(metric, reaction_indices, options)

        if style == "barplot":
            self._barplot_sensitivity(
                metric,
                sensitivity_coefficients,
                biological_processes,
                reaction_indices,
                save_format,
                options["show_indices"],
            )
        elif style == "heatmap":
            self._heatmap_sensitivity(
                metric,
                sensitivity_coefficients,
                reaction_indices,
                save_format,
            )
        else:
            raise ValueError("Available styles are: 'barplot', 'heatmap'")


def is_duplicate(
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
