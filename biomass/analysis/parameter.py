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
class ParameterSensitivity(SignalingMetric):
    """Sensitivity for parameters"""

    model: ModelObject
    create_metrics: Optional[Dict[str, Callable[[np.ndarray], Union[int, float]]]]

    def __post_init__(self) -> None:
        self._plotting: SensitivityOptions = self.model.viz.get_sensitivity_options()
        self._coefficients: Callable[[str], str] = lambda metric: os.path.join(
            self.model.path,
            "sensitivity_coefficients",
            "parameter",
            f"{metric}.npy",
        )
        self._path_to_figs: Callable[[str], str] = lambda metric: os.path.join(
            self.model.path,
            "figure",
            "sensitivity",
            "parameter",
            f"{metric}",
        )
        if self.create_metrics is not None:
            for name, function in self.create_metrics.items():
                self.quantification[name] = function

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

        rate: Final[float] = 1.01  # 1% change
        n_file = self.model.get_executable()

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
            optimized = self.model.load_param(nth_paramset)
            for j, idx in enumerate(param_indices):
                x = optimized.params[:]
                x[idx] = optimized.params[idx] * rate
                if self.model.problem.simulate(x, optimized.initials) is None:
                    for k, _ in enumerate(self.model.observables):
                        for l, _ in enumerate(self.model.problem.conditions):
                            signaling_metric[i, j, k, l] = self.quantification[metric](
                                self.model.problem.simulations[k, l]
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
                        signaling_metric[i, -1, k, l] = self.quantification[metric](
                            self.model.problem.simulations[k, l]
                        )
        sensitivity_coefficients = dlnyi_dlnxj(
            signaling_metric,
            len(n_file),
            len(param_indices),
            len(self.model.observables),
            len(self.model.problem.conditions),
            rate,
        )

        return sensitivity_coefficients

    def _load_sc(
        self,
        metric: str,
        param_indices: List[int],
    ) -> np.ndarray:
        """
        Load (or calculate) sensitivity coefficients.
        """
        if not os.path.isfile(self._coefficients(metric)):
            os.makedirs(
                os.path.join(
                    self.model.path,
                    "sensitivity_coefficients",
                    "parameter",
                ),
                exist_ok=True,
            )
            sensitivity_coefficients = self._calc_sensitivity_coefficients(metric, param_indices)
            np.save(self._coefficients(metric), sensitivity_coefficients)
        else:
            sensitivity_coefficients = np.load(self._coefficients(metric))
            if len(param_indices) != sensitivity_coefficients.shape[1]:
                # User changed options['excluded_params'] after the last trial
                sensitivity_coefficients = self._calc_sensitivity_coefficients(
                    metric, param_indices
                )
                np.save(self._coefficients(metric), sensitivity_coefficients)

        return sensitivity_coefficients

    def _barplot_sensitivity(
        self,
        metric: str,
        sensitivity_coefficients: np.ndarray,
        param_indices: List[int],
    ) -> None:
        """
        Visualize sensitivity coefficients using barplot.
        """
        os.makedirs(os.path.join(self._path_to_figs(metric), "barplot"), exist_ok=True)
        # rcParams
        self.model.viz.set_sensitivity_rcParams()

        for k, obs_name in enumerate(self.model.observables):
            plt.figure(figsize=self._plotting.figsize)
            plt.hlines([0], -self._plotting.width, len(param_indices), "k", lw=1)
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
                        np.arange(len(param_indices)) + l * self._plotting.width,
                        average,
                        yerr=stdev,
                        ecolor=self._plotting.cmap[l],
                        capsize=2,
                        width=self._plotting.width,
                        color=self._plotting.cmap[l],
                        align="center",
                        label=condition,
                    )
            plt.xticks(
                np.arange(len(param_indices))
                + self._plotting.width * 0.5 * (len(self.model.problem.conditions) - 1),
                [self.model.parameters[i] for i in param_indices],
                fontsize=6,
                rotation=90,
            )
            plt.ylabel(
                "Control coefficients on\n" + metric + " (" + obs_name.replace("_", " ") + ")"
            )
            plt.xlim(-self._plotting.width, len(param_indices))
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
        param_indices: List[int],
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
                        xticklabels=[self.model.parameters[i] for i in param_indices],
                        yticklabels=[],
                        # cbar_kws={"ticks": [-1, 0, 1]}
                    )
                    cbar = g.ax_heatmap.collections[0].colorbar
                    cbar.ax.tick_params(labelsize=8)
                    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)
                    plt.savefig(
                        os.path.join(
                            self._path_to_figs(metric),
                            "heatmap",
                            f"{condition.replace(',', '')}_{obs_name}",
                        ),
                    )
                    plt.close()

    def analyze(self, *, metric: str, style: str, options: dict) -> None:
        """
        Perform sensitivity analysis.
        """
        if options["overwrite"] and os.path.isfile(self._coefficients(metric)):
            os.remove(self._coefficients(metric))
        param_indices = self._get_param_indices(options["excluded_params"])
        sensitivity_coefficients = self._load_sc(metric, param_indices)
        if style == "barplot":
            self._barplot_sensitivity(
                metric,
                sensitivity_coefficients,
                param_indices,
            )
        elif style == "heatmap":
            if len(param_indices) < 2:
                pass
            else:
                self._heatmap_sensitivity(
                    metric,
                    sensitivity_coefficients,
                    param_indices,
                )
        else:
            raise ValueError("Available styles are: 'barplot', 'heatmap'")
