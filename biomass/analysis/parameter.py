import os
import sys
from dataclasses import dataclass
from typing import List

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from ..exec_model import ExecModel, ModelObject
from .util import dlnyi_dlnxj, get_signaling_metric


@dataclass
class ParameterSensitivity(ExecModel):
    """Sensitivity for parameters"""

    model: ModelObject

    def _get_param_indices(self, excluded_params: List[str]) -> List[int]:
        param_indices = []
        x = self.model.pval()
        for i, val in enumerate(x):
            if self.model.parameters[i] not in excluded_params and val != 0.0:
                param_indices.append(i)
        if not param_indices:
            raise ValueError("No nonzero parameters.")

        return param_indices

    def _calc_sensitivity_coefficients(
        self,
        metric: str,
        param_indices: List[int],
        options,
    ) -> np.ndarray:
        """Calculating Sensitivity Coefficients

        Parameters
        ----------
        metric : str
            The signaling metric used for sensitivity analysis.

        param_indices : list of int
            List of parameter indices for sensitivity analysis.

        Returns
        -------
        sensitivity_coefficients : numpy array

        """

        rate = 1.01  # 1% change
        n_file = self.get_executable()

        signaling_metric = np.full(
            (
                len(n_file),
                len(param_indices) + 1,
                len(self.model.observables),
                len(self.model.problem.conditions),
            ),
            np.nan,
        )
        for i, nth_paramset in enumerate(n_file):
            optimized = self.load_param(nth_paramset)
            for j, idx in enumerate(param_indices):
                x = optimized.params[:]
                x[idx] = optimized.params[idx] * rate
                if self.model.problem.simulate(x, optimized.initials) is None:
                    for k, _ in enumerate(self.model.observables):
                        for l, _ in enumerate(self.model.problem.conditions):
                            signaling_metric[i, j, k, l] = get_signaling_metric(
                                metric, self.model.problem.simulations[k, :, l], options
                            )
                sys.stdout.write(
                    "\r{:d} / {:d}".format(
                        i * len(param_indices) + j + 1, len(n_file) * len(param_indices)
                    )
                )
            # Signaling metric without perturbation (j=-1)
            x = optimized.params[:]
            if self.model.problem.simulate(x, optimized.initials) is None:
                for k, _ in enumerate(self.model.observables):
                    for l, _ in enumerate(self.model.problem.conditions):
                        signaling_metric[i, -1, k, l] = get_signaling_metric(
                            metric, self.model.problem.simulations[k, :, l], options
                        )
        sensitivity_coefficients = dlnyi_dlnxj(
            signaling_metric,
            n_file,
            param_indices,
            self.model.observables,
            self.model.problem.conditions,
            rate,
        )

        return sensitivity_coefficients

    def _load_sc(
        self,
        metric: str,
        param_indices: List[int],
        options: dict,
    ) -> np.ndarray:
        """
        Load (or calculate) sensitivity coefficients.
        """
        if not os.path.isfile(
            os.path.join(
                self.model.path,
                "sensitivity_coefficients",
                "parameter",
                f"{metric}.npy",
            )
        ):
            os.makedirs(
                os.path.join(
                    self.model.path,
                    "sensitivity_coefficients",
                    "parameter",
                ),
                exist_ok=True,
            )
            sensitivity_coefficients = self._calc_sensitivity_coefficients(
                metric,
                param_indices,
                options,
            )
            np.save(
                os.path.join(
                    self.model.path,
                    "sensitivity_coefficients",
                    "parameter",
                    f"{metric}",
                ),
                sensitivity_coefficients,
            )
        else:
            sensitivity_coefficients = np.load(
                os.path.join(
                    self.model.path,
                    "sensitivity_coefficients",
                    "parameter",
                    f"{metric}.npy",
                )
            )
            if len(param_indices) != sensitivity_coefficients.shape[1]:
                # User changed options['excluded_params'] after the last trial
                sensitivity_coefficients = self._calc_sensitivity_coefficients(
                    metric,
                    param_indices,
                    options,
                )
                np.save(
                    os.path.join(
                        self.model.path,
                        "sensitivity_coefficients",
                        "parameter",
                        f"{metric}",
                    ),
                    sensitivity_coefficients,
                )

        return sensitivity_coefficients

    def _barplot_sensitivity(
        self,
        metric: str,
        sensitivity_coefficients: np.ndarray,
        param_indices: List[int],
        save_format: str,
    ) -> None:
        """
        Visualize sensitivity coefficients using barplot.
        """
        os.makedirs(
            os.path.join(
                self.model.path,
                "figure",
                "sensitivity",
                "parameter",
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
            plt.hlines([0], -options["width"], len(param_indices), "k", lw=1)
            for l, condition in enumerate(self.model.problem.conditions):
                sensitivity_matrix = sensitivity_coefficients[:, :, k, l]
                nan_idx = []
                for m in range(sensitivity_matrix.shape[0]):
                    if np.isnan(sensitivity_matrix[m, :]).any():
                        nan_idx.append(m)
                sensitivity_matrix = np.delete(sensitivity_matrix, nan_idx, axis=0)
                if sensitivity_matrix.size != 0:
                    average = np.mean(sensitivity_matrix, axis=0)
                    if sensitivity_matrix.shape[0] == 1:
                        stdev = np.zeros(sensitivity_matrix.shape[1])
                    else:
                        stdev = np.std(sensitivity_matrix, axis=0, ddof=1)
                    plt.bar(
                        np.arange(len(param_indices)) + l * options["width"],
                        average,
                        yerr=stdev,
                        ecolor=options["cmap"][l],
                        capsize=2,
                        width=options["width"],
                        color=options["cmap"][l],
                        align="center",
                        label=condition,
                    )
            plt.xticks(
                np.arange(len(param_indices))
                + options["width"] * 0.5 * (len(self.model.problem.conditions) - 1),
                [self.model.parameters[i] for i in param_indices],
                fontsize=6,
                rotation=90,
            )
            plt.ylabel(
                "Control coefficients on\n" + metric + " (" + obs_name.replace("_", " ") + ")"
            )
            plt.xlim(-options["width"], len(param_indices))
            plt.legend(loc=options["legend_loc"], frameon=False)
            plt.savefig(
                os.path.join(
                    self.model.path,
                    "figure",
                    "sensitivity",
                    "parameter",
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
        param_indices: List[int],
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
                "parameter",
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
                        xticklabels=[self.model.parameters[i] for i in param_indices],
                        yticklabels=[],
                        # cbar_kws={"ticks": [-1, 0, 1]}
                    )
                    cbar = g.ax_heatmap.collections[0].colorbar
                    cbar.ax.tick_params(labelsize=8)
                    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)
                    plt.savefig(
                        os.path.join(
                            self.model.path,
                            "figure",
                            "sensitivity",
                            "parameter",
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
        param_indices = self._get_param_indices(options["excluded_params"])
        sensitivity_coefficients = self._load_sc(metric, param_indices, options)
        if style == "barplot":
            self._barplot_sensitivity(
                metric,
                sensitivity_coefficients,
                param_indices,
                save_format,
            )
        elif style == "heatmap":
            if len(param_indices) < 2:
                pass
            else:
                self._heatmap_sensitivity(
                    metric,
                    sensitivity_coefficients,
                    param_indices,
                    save_format,
                )
        else:
            raise ValueError("Available styles are: 'barplot', 'heatmap'")
