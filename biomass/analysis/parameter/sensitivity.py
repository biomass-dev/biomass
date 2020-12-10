import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List

from biomass.exec_model import ExecModel
from biomass.analysis import get_signaling_metric, dlnyi_dlnxj


class ParameterSensitivity(ExecModel):
    """Sensitivity for parameters"""

    def __init__(self, model, excluded_params):
        super().__init__(model)
        self.excluded_params = excluded_params

    def _get_param_indices(self) -> List[int]:
        param_indices = []
        x = self.pval()
        for i, val in enumerate(x):
            if self.parameters[i] not in self.excluded_params and val != 0.0:
                param_indices.append(i)
        if not param_indices:
            raise ValueError("No nonzero parameters")

        return param_indices

    def _calc_sensitivity_coefficients(self, metric: str, param_indices: List[int]) -> np.ndarray:
        """Calculating Sensitivity Coefficients

        Parameters
        ----------
        metric : str
            - 'maximum': The maximum value.
            - 'minimum': The minimum value.
            - 'duration': The time it takes to decline below 10% of its maximum.
            - 'integral': The integral of concentration over the observation time.

        param_indices : list of int
            List of parameter indices for sensitivity analysis.

        Returns
        -------
        sensitivity_coefficients : numpy array

        """

        rate = 1.01  # 1% change
        param_indices = self._get_param_indices()
        n_file = self.get_executable()

        signaling_metric = np.full(
            (
                len(n_file),
                len(param_indices) + 1,
                len(self.obs),
                len(self.sim.conditions),
            ),
            np.nan,
        )
        for i, nth_paramset in enumerate(n_file):
            (x, y0) = self.load_param(nth_paramset)
            x_init = x[:]
            for j, idx in enumerate(param_indices):
                x = x_init[:]
                x[idx] = x_init[idx] * rate
                if self.sim.simulate(x, y0) is None:
                    for k, _ in enumerate(self.obs):
                        for l, _ in enumerate(self.sim.conditions):
                            signaling_metric[i, j, k, l] = get_signaling_metric(metric, self.sim.simulations[k, :, l])
                sys.stdout.write(
                    "\r{:d} / {:d}".format(i * len(param_indices) + j + 1, len(n_file) * len(param_indices))
                )
            # Signaling metric without perturbation (j=-1)
            x = x_init[:]
            if self.sim.simulate(x, y0) is None:
                for k, _ in enumerate(self.obs):
                    for l, _ in enumerate(self.sim.conditions):
                        signaling_metric[i, -1, k, l] = get_signaling_metric(metric, self.sim.simulations[k, :, l])
        sensitivity_coefficients = dlnyi_dlnxj(
            signaling_metric, n_file, param_indices, self.obs, self.sim.conditions, rate
        )

        return sensitivity_coefficients

    def _load_sc(self, metric: str, param_indices: List[int]):
        """
        Load (or calculate) sensitivity coefficients.
        """
        os.makedirs(
            self.model_path + "/figure/sensitivity/" f"parameter/{metric}/heatmap",
            exist_ok=True,
        )
        if not os.path.isfile(self.model_path + "/sensitivity_coefficients/" f"parameter/{metric}/sc.npy"):
            os.makedirs(
                self.model_path + "/sensitivity_coefficients/" f"parameter/{metric}",
                exist_ok=True,
            )
            sensitivity_coefficients = self._calc_sensitivity_coefficients(metric, param_indices)
            np.save(
                self.model_path + "/sensitivity_coefficients/" f"parameter/{metric}/sc",
                sensitivity_coefficients,
            )
        else:
            sensitivity_coefficients = np.load(
                self.model_path + "/sensitivity_coefficients/" f"parameter/{metric}/sc.npy"
            )

        return sensitivity_coefficients

    def _barplot_sensitivity(
        self,
        metric: str,
        sensitivity_coefficients: np.ndarray,
        param_indices: List[int],
    ):
        """
        Visualize sensitivity coefficients using barplot.
        """
        options = self.viz.sensitivity_options

        # rcParams
        self.viz.set_sensitivity_rcParams()

        if len(options["cmap"]) < len(self.sim.conditions):
            raise ValueError(
                "len(sensitivity_options['cmap']) must be equal to" " or greater than len(sim.conditions)."
            )
        for k, obs_name in enumerate(self.obs):
            plt.figure(figsize=options["figsize"])
            plt.hlines([0], -options["width"], len(param_indices), "k", lw=1)
            for l, condition in enumerate(self.sim.conditions):
                sensitivity_matrix = sensitivity_coefficients[:, :, k, l]
                nan_idx = []
                for m in range(sensitivity_matrix.shape[0]):
                    if any(np.isnan(sensitivity_matrix[m, :])):
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
                np.arange(len(param_indices)) + options["width"] * 0.5 * (len(self.sim.conditions) - 1),
                [self.parameters[i] for i in param_indices],
                fontsize=6,
                rotation=90,
            )
            plt.ylabel("Control coefficients on\n" + metric + " (" + obs_name.replace("_", " ") + ")")
            plt.xlim(-options["width"], len(param_indices))
            plt.legend(loc=options["legend_loc"], frameon=False)
            plt.savefig(
                self.model_path + "/figure/sensitivity/parameter/" f"{metric}/{obs_name}.pdf",
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
            if any(np.isnan(sensitivity_matrix[i, :])):
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
    ):
        """
        Visualize sensitivity coefficients using heatmap.
        """
        options = self.viz.sensitivity_options
        # rcParams
        self.viz.set_sensitivity_rcParams()

        for k, obs_name in enumerate(self.obs):
            for l, condition in enumerate(self.sim.conditions):
                sensitivity_matrix = self._remove_nan(sensitivity_coefficients[:, :, k, l], normalize=False)
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
                        xticklabels=[self.parameters[i] for i in param_indices],
                        yticklabels=[],
                        # cbar_kws={"ticks": [-1, 0, 1]}
                    )
                    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)
                    plt.savefig(
                        self.model_path + "/figure/sensitivity/parameter/"
                        f"{metric}/heatmap/{condition}_{obs_name}.pdf",
                        bbox_inches="tight",
                    )
                    plt.close()

    def analyze(self, metric: str, style: str):
        """
        Perform sensitivity analysis.
        """
        param_indices = self._get_param_indices()
        sensitivity_coefficients = self._load_sc(metric, param_indices)
        if style == "barplot":
            self._barplot_sensitivity(metric, sensitivity_coefficients, param_indices)
        elif style == "heatmap":
            if len(param_indices) < 2:
                pass
            else:
                self._heatmap_sensitivity(metric, sensitivity_coefficients, param_indices)
        else:
            raise ValueError("Available styles are: 'barplot', 'heatmap'")
