import os
import warnings
from typing import List, NoReturn, Optional

import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt

from ..exec_model import ExecModel
from ..template import BioMassModel


class VisualizeWarning(UserWarning):
    pass


def _check_unknown_options(unknown_options: dict) -> Optional[NoReturn]:
    if unknown_options:
        msg = ", ".join(map(str, unknown_options.keys()))
        warnings.warn(f"Unknown visualization options: {msg}", VisualizeWarning)
    return None


class TemporalDynamics(ExecModel):
    def __init__(self, model: BioMassModel) -> None:
        super().__init__(model)

    def plot_timecourse(
        self,
        n_file: List[int],
        viz_type: str,
        show_all: bool,
        stdev: bool,
        save_format: str,
        simulations_all: np.ndarray,
    ) -> None:
        """
        Plot time course of each observable.

        Parameters
        ----------
        n_file : list of integers
            Optimized parameter sets in out/.

        viz_type : str
            One of ['average', 'best', 'original', 'n(=1,2,...)', 'experiment'].

        show_all : bool
            Whether to show all simulation results.

        stdev : bool
            If True, the standard deviation of simulated values will be shown
            (only available for 'average' visualization type).

        save_format : str (default: "pdf")
            Either "png" or "pdf", indicating whether to save figures
            as png or pdf format.

        simulations_all : numpy array
            Array containing all simulated values.

        """
        os.makedirs(
            os.path.join(
                self.model.path,
                "figure",
                "simulation",
                f"{viz_type}",
            ),
            exist_ok=True,
        )
        self.model.exp.set_data()
        self.model.viz.set_timecourse_rcParams()
        singleplot = self.model.viz.get_timecourse_options()
        multiplot = self.model.viz.multiplot_observables()

        for mode in range(2):
            # mode 0 : timecourse_for_each_observable
            # mode 1 : multiplot_observables
            set_fig = False
            for i, obs_name in enumerate(self.model.obs):
                if len(singleplot[i]["cmap"]) < len(self.model.sim.conditions) or len(
                    singleplot[i]["shape"]
                ) < len(self.model.sim.conditions):
                    raise ValueError(
                        "len(cmap), len(shape) must be equal to"
                        " or greater than len(self.model.sim.conditions)."
                    )
                if mode == 1 and obs_name not in multiplot["observables"]:
                    continue
                if mode == 0:
                    plt.figure(figsize=(4, 3))
                elif mode == 1 and not set_fig:
                    plt.figure(figsize=(4, 3))
                    set_fig = True
                plt.gca().spines["right"].set_visible(False)
                plt.gca().spines["top"].set_visible(False)
                if viz_type != "experiment":
                    if show_all:
                        self._plot_show_all(
                            n_file, simulations_all, obs_name, mode, singleplot, multiplot
                        )
                    if viz_type == "average":
                        normalized = self._normalize_array(
                            n_file, simulations_all, obs_name, mode, singleplot, multiplot
                        )
                        if (
                            self.model.sim.normalization
                            and self.model.sim.normalization[obs_name]["timepoint"] is None
                        ):
                            normalized = self._divide_by_maximum(normalized, obs_name)
                        self._plot_average(normalized, obs_name, mode, singleplot, multiplot)
                        if stdev:
                            self._show_sd(normalized, obs_name, mode, singleplot, multiplot)
                    else:
                        self._plot_simulations(obs_name, mode, singleplot, multiplot)
                if singleplot[i]["exp_data"] and self.model.exp.experiments[i] is not None:
                    exp_t = self.model.exp.get_timepoint(obs_name)
                    if self.model.exp.error_bars[i] is not None:
                        self._plot_experimental_data_with_error_bars(
                            viz_type, exp_t, obs_name, mode, singleplot, multiplot
                        )
                    else:
                        self._plot_experimental_data_without_error_bars(
                            viz_type, exp_t, obs_name, mode, singleplot, multiplot
                        )
                if mode == 0:
                    self._save_mode_0(obs_name, singleplot, viz_type, save_format)
            if mode == 1 and multiplot["observables"]:
                self._save_mode_1(multiplot, viz_type, save_format)

    def plot_param_range(
        self,
        popt: np.ndarray,
        save_format: str,
        orientation: str,
        distribution: str,
        scatter: bool,
        **unknown_options,
    ) -> None:
        """
        Plot estimated parameter/initial values.

        Parameters
        ----------
        popt : numpy array
            Array containing all optimized parameter/initial values.

        save_format : str (default: "pdf")
            Either "png" or "pdf", indicating whether to save figures
            as png or pdf format.

        orientation : str (default: 'portrait')
            Either 'portrait' or 'landscape'.

        distribution : str (default: 'boxenplot')
            Either 'boxplot' or 'boxenplot'.

        scatter : bool (default: False)
            If True, draw a stripplot.

        """
        _check_unknown_options(unknown_options)
        os.makedirs(
            os.path.join(
                self.model.path,
                "figure",
                "param_range",
            ),
            exist_ok=True,
        )

        self.model.viz.set_param_range_rcParams()

        for val in ["parameter_value", "initial_value"]:
            if (val == "parameter_value" and len(self.model.sp.idx_params) == 0) or (
                val == "initial_value" and len(self.model.sp.idx_initials) == 0
            ):
                continue
            else:
                fig = plt.figure(
                    figsize=(
                        8,
                        (
                            len(self.model.sp.idx_params)
                            if val == "parameter_value"
                            else len(self.model.sp.idx_initials)
                        )
                        / 2.5,
                    )
                    if orientation == "portrait"
                    else (
                        (
                            len(self.model.sp.idx_params)
                            if val == "parameter_value"
                            else len(self.model.sp.idx_initials)
                        )
                        / 2.2,
                        6,
                    )
                )
                arguments_distribution = {
                    "data": popt[:, : len(self.model.sp.idx_params)]
                    if val == "parameter_value"
                    else popt[:, len(self.model.sp.idx_params) :],
                    "orient": "h" if orientation == "portrait" else "v",
                    "linewidth": 0.5,
                    "palette": "Set2",
                }
                arguments_scatter = {
                    "data": popt[:, : len(self.model.sp.idx_params)]
                    if val == "parameter_value"
                    else popt[:, len(self.model.sp.idx_params) :],
                    "orient": "h" if orientation == "portrait" else "v",
                    "size": 2.5,
                    "color": ".26",
                }
                if distribution == "boxplot":
                    ax = sns.boxplot(**arguments_distribution)
                else:
                    ax = sns.boxenplot(**arguments_distribution)
                if scatter:
                    ax = sns.stripplot(**arguments_scatter)
                sns.despine()
                if orientation == "portrait":
                    ax.set_xlabel(val.capitalize().replace("_", " "))
                    ax.set_xscale("log")
                    ax.set_ylabel("")
                    ax.set_yticklabels(
                        [self.model.parameters[i] for i in self.model.sp.idx_params]
                        if val == "parameter_value"
                        else [self.model.species[i] for i in self.model.sp.idx_initials]
                    )
                else:
                    ax.set_xlabel("")
                    ax.set_xticklabels(
                        [self.model.parameters[i] for i in self.model.sp.idx_params]
                        if val == "parameter_value"
                        else [self.model.species[i] for i in self.model.sp.idx_initials],
                        rotation=90,
                    )
                    ax.set_ylabel(val.capitalize().replace("_", " "))
                    ax.set_yscale("log")
                plt.savefig(
                    os.path.join(
                        self.model.path,
                        "figure",
                        "param_range",
                        f"estimated_{val}s." + save_format,
                    ),
                    dpi=600 if save_format == "png" else None,
                    bbox_inches="tight",
                )
                plt.close(fig)

    def _plot_show_all(
        self,
        n_file: List[int],
        simulations_all: np.ndarray,
        obs_name: str,
        mode: int,
        singleplot: List[dict],
        multiplot: dict,
    ) -> None:
        """
        Plot time course simulated values (show_all == True).
        """
        i = self.model.obs.index(obs_name)
        for j, _ in enumerate(n_file):
            for l, condition in enumerate(self.model.sim.conditions):
                if (mode == 0 and condition not in singleplot[i]["dont_show"]) or (
                    mode == 1 and condition == multiplot["condition"]
                ):
                    plt.plot(
                        np.array(self.model.sim.t) / singleplot[i]["divided_by"],
                        simulations_all[i, j, :, l]
                        / (
                            1
                            if not self.model.sim.normalization
                            or np.max(simulations_all[i, j, :, l]) == 0.0
                            else np.max(
                                simulations_all[
                                    i,
                                    j,
                                    self.model.sim.normalization[obs_name]["timepoint"],
                                    [
                                        self.model.sim.conditions.index(c)
                                        for c in self.model.sim.normalization[obs_name][
                                            "condition"
                                        ]
                                    ],
                                ]
                            )
                            if self.model.sim.normalization[obs_name]["timepoint"] is not None
                            else np.max(
                                simulations_all[
                                    i,
                                    j,
                                    :,
                                    [
                                        self.model.sim.conditions.index(c)
                                        for c in self.model.sim.normalization[obs_name][
                                            "condition"
                                        ]
                                    ],
                                ]
                            )
                        ),
                        color=singleplot[i]["cmap"][l]
                        if mode == 0
                        else multiplot["cmap"][multiplot["observables"].index(obs_name)],
                        alpha=0.05,
                    )

    def _normalize_array(
        self,
        n_file: List[int],
        simulations_all: np.ndarray,
        obs_name: str,
        mode: int,
        singleplot: List[dict],
        multiplot: dict,
    ) -> np.ndarray:
        """
        Normalize the array simulations_all using sim.normalization set in observable.py.
        """
        normalized: np.ndarray = np.empty_like(simulations_all)
        i = self.model.obs.index(obs_name)
        for j, _ in enumerate(n_file):
            for l, condition in enumerate(self.model.sim.conditions):
                if (mode == 0 and condition not in singleplot[i]["dont_show"]) or (
                    mode == 1 and condition == multiplot["condition"]
                ):
                    normalized[i, j, :, l] = simulations_all[i, j, :, l] / (
                        1
                        if not self.model.sim.normalization
                        or np.max(simulations_all[i, j, :, l]) == 0.0
                        else np.max(
                            simulations_all[
                                i,
                                j,
                                self.model.sim.normalization[obs_name]["timepoint"],
                                [
                                    self.model.sim.conditions.index(c)
                                    for c in self.model.sim.normalization[obs_name]["condition"]
                                ],
                            ]
                        )
                        if self.model.sim.normalization[obs_name]["timepoint"] is not None
                        else np.max(
                            simulations_all[
                                i,
                                j,
                                :,
                                [
                                    self.model.sim.conditions.index(c)
                                    for c in self.model.sim.normalization[obs_name]["condition"]
                                ],
                            ]
                        )
                    )
        return normalized

    def _divide_by_maximum(self, normalized: np.ndarray, obs_name: str) -> np.ndarray:
        """
        Divide the array by its maximum.
        """
        mean_vec = []
        for c in self.model.sim.normalization[obs_name]["condition"]:
            mean_vec.append(
                np.nanmean(
                    normalized[
                        self.model.obs.index(obs_name), :, :, self.model.sim.conditions.index(c)
                    ],
                    axis=0,
                )
            )
        norm_max = np.max(mean_vec)
        if not np.isnan(norm_max) and norm_max != 0.0:
            normalized[self.model.obs.index(obs_name), :, :, :] /= norm_max

        return normalized

    def _plot_average(
        self,
        normalized: np.ndarray,
        obs_name: str,
        mode: int,
        singleplot: List[dict],
        multiplot: dict,
    ) -> None:
        """
        Plot time course simulated values (viz_type == 'average').
        """
        i = self.model.obs.index(obs_name)
        for l, condition in enumerate(self.model.sim.conditions):
            if (mode == 0 and condition not in singleplot[i]["dont_show"]) or (
                mode == 1 and condition == multiplot["condition"]
            ):
                plt.plot(
                    np.array(self.model.sim.t) / singleplot[i]["divided_by"],
                    np.nanmean(normalized[i, :, :, l], axis=0),
                    color=singleplot[i]["cmap"][l]
                    if mode == 0
                    else multiplot["cmap"][multiplot["observables"].index(obs_name)],
                    label=condition if mode == 0 else singleplot[i]["ylabel"],
                )

    def _show_sd(
        self,
        normalized: np.ndarray,
        obs_name: str,
        mode: int,
        singleplot: List[dict],
        multiplot: dict,
    ) -> None:
        """
        Plot standard deviation (SD) as shaded area when stdev == True.
        """
        i = self.model.obs.index(obs_name)
        for l, condition in enumerate(self.model.sim.conditions):
            if (mode == 0 and condition not in singleplot[i]["dont_show"]) or (
                mode == 1 and condition == multiplot["condition"]
            ):
                y_mean = np.nanmean(normalized[i, :, :, l], axis=0)
                y_std = [
                    np.nanstd(normalized[i, :, k, l], ddof=1)
                    for k, _ in enumerate(self.model.sim.t)
                ]
                plt.fill_between(
                    np.array(self.model.sim.t) / singleplot[i]["divided_by"],
                    y_mean - y_std,
                    y_mean + y_std,
                    lw=0,
                    color=singleplot[i]["cmap"][l]
                    if mode == 0
                    else multiplot["cmap"][multiplot["observables"].index(obs_name)],
                    alpha=0.1,
                )

    def _plot_simulations(
        self,
        obs_name: str,
        mode: int,
        singleplot: List[dict],
        multiplot: dict,
    ) -> None:
        """
        Plot time course simulated values (viz_type not in ['average', 'experiment']).
        """
        i = self.model.obs.index(obs_name)
        for l, condition in enumerate(self.model.sim.conditions):
            if (mode == 0 and condition not in singleplot[i]["dont_show"]) or (
                mode == 1 and condition == multiplot["condition"]
            ):
                plt.plot(
                    np.array(self.model.sim.t) / singleplot[i]["divided_by"],
                    self.model.sim.simulations[i, :, l]
                    / (
                        1
                        if not self.model.sim.normalization
                        or np.max(self.model.sim.simulations[i, :, l]) == 0.0
                        else np.max(
                            self.model.sim.simulations[
                                i,
                                self.model.sim.normalization[obs_name]["timepoint"],
                                [
                                    self.model.sim.conditions.index(c)
                                    for c in self.model.sim.normalization[obs_name]["condition"]
                                ],
                            ]
                        )
                        if self.model.sim.normalization[obs_name]["timepoint"] is not None
                        else np.max(
                            self.model.sim.simulations[
                                i,
                                :,
                                [
                                    self.model.sim.conditions.index(c)
                                    for c in self.model.sim.normalization[obs_name]["condition"]
                                ],
                            ]
                        )
                    ),
                    color=singleplot[i]["cmap"][l]
                    if mode == 0
                    else multiplot["cmap"][multiplot["observables"].index(obs_name)],
                    label=condition if mode == 0 else singleplot[i]["ylabel"],
                )

    def _plot_experimental_data_with_error_bars(
        self,
        viz_type: str,
        exp_t: Optional[List[int]],
        obs_name: str,
        mode: int,
        singleplot: List[dict],
        multiplot: dict,
    ) -> None:
        """
        Plot experimental measurements with error bars.
        """
        i = self.model.obs.index(obs_name)
        for l, condition in enumerate(self.model.sim.conditions):
            if (
                condition in self.model.exp.experiments[i]
                and (mode == 0 and condition not in singleplot[i]["dont_show"])
                or (mode == 1 and condition == multiplot["condition"])
            ):
                exp_data = plt.errorbar(
                    np.array(exp_t) / singleplot[i]["divided_by"],
                    self.model.exp.experiments[i][condition],
                    yerr=self.model.exp.error_bars[i][condition],
                    color=singleplot[i]["cmap"][l]
                    if mode == 0
                    else multiplot["cmap"][multiplot["observables"].index(obs_name)],
                    ecolor=singleplot[i]["cmap"][l]
                    if mode == 0
                    else multiplot["cmap"][multiplot["observables"].index(obs_name)],
                    elinewidth=1,
                    capsize=8,
                    markerfacecolor="None",
                    markeredgecolor=singleplot[i]["cmap"][l]
                    if mode == 0
                    else multiplot["cmap"][multiplot["observables"].index(obs_name)],
                    fmt=singleplot[i]["shape"][l]
                    if mode == 0
                    else multiplot["shape"][multiplot["observables"].index(obs_name)],
                    clip_on=False,
                    label=singleplot[i]["ylabel"]
                    if mode == 1 and viz_type == "experiment"
                    else None,
                )
                for capline in exp_data[1]:
                    capline.set_clip_on(False)
                for barlinecol in exp_data[2]:
                    barlinecol.set_clip_on(False)

    def _plot_experimental_data_without_error_bars(
        self,
        viz_type: str,
        exp_t: Optional[List[int]],
        obs_name: str,
        mode: int,
        singleplot: List[dict],
        multiplot: dict,
    ) -> None:
        """
        Plot experimental measurements when model.exp.error_bars[i] is None.
        """
        i = self.model.obs.index(obs_name)
        for l, condition in enumerate(self.model.sim.conditions):
            if (
                condition in self.model.exp.experiments[i]
                and (mode == 0 and condition not in singleplot[i]["dont_show"])
                or (mode == 1 and condition == multiplot["condition"])
            ):
                plt.plot(
                    np.array(exp_t) / singleplot[i]["divided_by"],
                    self.model.exp.experiments[i][condition],
                    singleplot[i]["shape"][l]
                    if mode == 0
                    else multiplot["shape"][multiplot["observables"].index(obs_name)],
                    markerfacecolor="None",
                    markeredgecolor=singleplot[i]["cmap"][l]
                    if mode == 0
                    else multiplot["cmap"][multiplot["observables"].index(obs_name)],
                    color=singleplot[i]["cmap"][l]
                    if mode == 0
                    else multiplot["cmap"][multiplot["observables"].index(obs_name)],
                    clip_on=False,
                    label=singleplot[i]["ylabel"]
                    if mode == 1 and viz_type == "experiment"
                    else None,
                )

    def _save_mode_0(
        self,
        obs_name: str,
        singleplot: List[dict],
        viz_type: str,
        save_format: str,
    ) -> None:
        """
        Plot time course of each observable.
        """
        i = self.model.obs.index(obs_name)
        if singleplot[i]["xlim"]:
            plt.xlim(singleplot[i]["xlim"])
        if singleplot[i]["xticks"] is not None:
            plt.xticks(singleplot[i]["xticks"])
        plt.xlabel(singleplot[i]["xlabel"])
        if singleplot[i]["ylim"]:
            plt.ylim(singleplot[i]["ylim"])
        if singleplot[i]["yticks"] is not None:
            plt.yticks(singleplot[i]["yticks"])
        plt.ylabel(singleplot[i]["ylabel"])
        if singleplot[i]["legend_loc"] is not None:
            plt.legend(
                loc=singleplot[i]["legend_loc"],
                frameon=False,
                fontsize=12,
            )
        plt.savefig(
            os.path.join(
                self.model.path,
                "figure",
                "simulation",
                f"{viz_type}",
                f"{obs_name}." + save_format,
            ),
            dpi=600 if save_format == "png" else None,
            bbox_inches="tight",
        )
        plt.close()

    def _save_mode_1(
        self,
        multiplot: dict,
        viz_type: str,
        save_format: str,
    ) -> None:
        """
        Plot time course of multiple observables in one figure.
        """
        if multiplot["xlim"]:
            plt.xlim(multiplot["xlim"])
        if multiplot["xticks"] is not None:
            plt.xticks(multiplot["xticks"])
        plt.xlabel(multiplot["xlabel"])
        if multiplot["ylim"]:
            plt.ylim(multiplot["ylim"])
        if multiplot["yticks"] is not None:
            plt.yticks(multiplot["yticks"])
        plt.ylabel(multiplot["ylabel"])
        plt.legend(
            bbox_to_anchor=(1.05, 1),
            loc="upper left",
            borderaxespad=0,
            labelspacing=1.25,
            frameon=False,
            fontsize=12,
        )
        plt.savefig(
            os.path.join(
                self.model.path,
                "figure",
                "simulation",
                f"{viz_type}",
                f"{multiplot['fig_name']}." + save_format,
            ),
            dpi=600 if save_format == "png" else None,
            bbox_inches="tight",
        )
        plt.close()
