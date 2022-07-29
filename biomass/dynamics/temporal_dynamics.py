import os
from dataclasses import dataclass
from math import isnan
from typing import List, Optional

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes._axes import _log as matplotlib_axes_logger

from ..model_object import ModelObject
from ..plotting import MultipleObservables, SingleObservable

matplotlib_axes_logger.setLevel("ERROR")


@dataclass
class TemporalDynamics(object):
    model: ModelObject

    def plot_timecourse(
        self,
        n_file: List[int],
        viz_type: str,
        show_all: bool,
        stdev: bool,
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

        simulations_all : numpy array
            Array containing all simulated values.

        """
        os.makedirs(
            os.path.join(self.model.path, "figure", "simulation", f"{viz_type}"),
            exist_ok=True,
        )
        self.model.problem.set_data()
        self.model.viz.set_timecourse_rcParams()
        singleplotting = self.model.viz.get_single_observable_options()
        multiplotting = self.model.viz.get_multiple_observables_options()
        for mode in range(2):
            # mode 0 : timecourse_for_each_observable
            # mode 1 : multiple_observables
            set_fig = False
            for i, obs_name in enumerate(self.model.observables):
                if mode == 1 and obs_name not in multiplotting.observables:
                    continue
                if mode == 0:
                    plt.figure(figsize=singleplotting[i].figsize)
                elif mode == 1 and not set_fig:
                    plt.figure(figsize=multiplotting.figsize)
                    set_fig = True
                plt.gca().spines["right"].set_visible(False)
                plt.gca().spines["top"].set_visible(False)
                if viz_type != "experiment":
                    if show_all:
                        self._plot_show_all(
                            n_file, simulations_all, obs_name, mode, singleplotting, multiplotting
                        )
                    if viz_type == "average":
                        normalized = self._normalize_array(n_file, simulations_all, obs_name)
                        if (
                            self.model.problem.normalization
                            and self.model.problem.normalization[obs_name]["timepoint"] is None
                        ):
                            normalized = self._divide_by_maximum(normalized, obs_name)
                        self._plot_average(
                            normalized, obs_name, mode, singleplotting, multiplotting
                        )
                        if stdev:
                            self._show_sd(
                                normalized, obs_name, mode, singleplotting, multiplotting
                            )
                    else:
                        self._plot_simulations(obs_name, mode, singleplotting, multiplotting)
                if (
                    viz_type == "experiment" or singleplotting[i].exp_data
                ) and self.model.problem.experiments[i] is not None:
                    exp_t = self.model.problem.get_timepoint(obs_name)
                    if self.model.problem.error_bars[i] is not None:
                        self._plot_experimental_data_with_error_bars(
                            viz_type, exp_t, obs_name, mode, singleplotting, multiplotting
                        )
                    else:
                        self._plot_experimental_data_without_error_bars(
                            viz_type, exp_t, obs_name, mode, singleplotting, multiplotting
                        )
                if mode == 0:
                    self._save_mode_0(obs_name, singleplotting, viz_type)
            if mode == 1 and multiplotting.observables:
                self._save_mode_1(multiplotting, viz_type)

    def _plot_show_all(
        self,
        n_file: List[int],
        simulations_all: np.ndarray,
        obs_name: str,
        mode: int,
        singleplotting: List[SingleObservable],
        multiplotting: MultipleObservables,
    ) -> None:
        """
        Plot time course simulated values (show_all == True).
        """
        i = self.model.observables.index(obs_name)
        for j, _ in enumerate(n_file):
            for k, condition in enumerate(self.model.problem.conditions):
                if (mode == 0 and condition not in singleplotting[i].dont_show) or (
                    mode == 1 and condition == multiplotting.condition
                ):
                    plt.plot(
                        np.array(self.model.problem.t) / singleplotting[i].divided_by,
                        simulations_all[i, j, k]
                        / (
                            1
                            if not self.model.problem.normalization
                            or np.max(simulations_all[i, j, k]) == 0.0
                            else np.max(
                                simulations_all[
                                    i,
                                    j,
                                    [
                                        self.model.problem.conditions.index(c)
                                        for c in self.model.problem.normalization[obs_name][
                                            "condition"
                                        ]
                                    ],
                                    self.model.problem.normalization[obs_name]["timepoint"],
                                ]
                            )
                            if self.model.problem.normalization[obs_name]["timepoint"] is not None
                            else np.max(
                                simulations_all[
                                    i,
                                    j,
                                    [
                                        self.model.problem.conditions.index(c)
                                        for c in self.model.problem.normalization[obs_name][
                                            "condition"
                                        ]
                                    ],
                                ]
                            )
                        ),
                        color=singleplotting[i].cmap[k]
                        if mode == 0
                        else multiplotting.cmap[multiplotting.observables.index(obs_name)],
                        alpha=0.05,
                    )

    def _normalize_array(
        self,
        n_file: List[int],
        simulations_all: np.ndarray,
        obs_name: str,
    ) -> np.ndarray:
        """
        Normalize the array simulations_all using problem.normalization set in observable.py.
        """
        normalized: np.ndarray = np.empty_like(simulations_all)
        i = self.model.observables.index(obs_name)
        for j, _ in enumerate(n_file):
            for k, _ in enumerate(self.model.problem.conditions):
                normalized[i, j, k] = simulations_all[i, j, k] / (
                    1
                    if not self.model.problem.normalization
                    or np.max(simulations_all[i, j, k]) == 0.0
                    else np.max(
                        simulations_all[
                            i,
                            j,
                            [
                                self.model.problem.conditions.index(c)
                                for c in self.model.problem.normalization[obs_name]["condition"]
                            ],
                            self.model.problem.normalization[obs_name]["timepoint"],
                        ]
                    )
                    if self.model.problem.normalization[obs_name]["timepoint"] is not None
                    else np.max(
                        simulations_all[
                            i,
                            j,
                            [
                                self.model.problem.conditions.index(c)
                                for c in self.model.problem.normalization[obs_name]["condition"]
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
        for c in self.model.problem.normalization[obs_name]["condition"]:
            mean_vec.append(
                np.nanmean(
                    normalized[
                        self.model.observables.index(obs_name),
                        :,
                        self.model.problem.conditions.index(c),
                    ],
                    axis=0,
                )
            )
        norm_max = np.max(mean_vec)
        if not isnan(norm_max) and norm_max != 0.0:
            normalized[self.model.observables.index(obs_name)] /= norm_max

        return normalized

    def _plot_average(
        self,
        normalized: np.ndarray,
        obs_name: str,
        mode: int,
        singleplotting: List[SingleObservable],
        multiplotting: MultipleObservables,
    ) -> None:
        """
        Plot time course simulated values (viz_type == 'average').
        """
        i = self.model.observables.index(obs_name)
        for k, condition in enumerate(self.model.problem.conditions):
            if (mode == 0 and condition not in singleplotting[i].dont_show) or (
                mode == 1 and condition == multiplotting.condition
            ):
                plt.plot(
                    np.array(self.model.problem.t) / singleplotting[i].divided_by,
                    np.nanmean(normalized[i, :, k, :], axis=0),
                    color=singleplotting[i].cmap[k]
                    if mode == 0
                    else multiplotting.cmap[multiplotting.observables.index(obs_name)],
                    label=condition if mode == 0 else singleplotting[i].ylabel,
                )

    def _show_sd(
        self,
        normalized: np.ndarray,
        obs_name: str,
        mode: int,
        singleplotting: List[SingleObservable],
        multiplotting: MultipleObservables,
    ) -> None:
        """
        Plot standard deviation (SD) as shaded area when stdev == True.
        """
        i = self.model.observables.index(obs_name)
        for k, condition in enumerate(self.model.problem.conditions):
            if (mode == 0 and condition not in singleplotting[i].dont_show) or (
                mode == 1 and condition == multiplotting.condition
            ):
                y_mean = np.nanmean(normalized[i, :, k], axis=0)
                y_std = [
                    np.nanstd(normalized[i, :, k, l], ddof=1)
                    for l, _ in enumerate(self.model.problem.t)
                ]
                plt.fill_between(
                    np.array(self.model.problem.t) / singleplotting[i].divided_by,
                    y_mean - y_std,
                    y_mean + y_std,
                    lw=0,
                    color=singleplotting[i].cmap[k]
                    if mode == 0
                    else multiplotting.cmap[multiplotting.observables.index(obs_name)],
                    alpha=0.1,
                )

    def _plot_simulations(
        self,
        obs_name: str,
        mode: int,
        singleplotting: List[SingleObservable],
        multiplotting: MultipleObservables,
    ) -> None:
        """
        Plot time course simulated values (viz_type not in ['average', 'experiment']).
        """
        i = self.model.observables.index(obs_name)
        for k, condition in enumerate(self.model.problem.conditions):
            if (mode == 0 and condition not in singleplotting[i].dont_show) or (
                mode == 1 and condition == multiplotting.condition
            ):
                plt.plot(
                    np.array(self.model.problem.t) / singleplotting[i].divided_by,
                    self.model.problem.simulations[i, k, :]
                    / (
                        1
                        if not self.model.problem.normalization
                        or np.max(self.model.problem.simulations[i, k, :]) == 0.0
                        else np.max(
                            self.model.problem.simulations[
                                i,
                                [
                                    self.model.problem.conditions.index(c)
                                    for c in self.model.problem.normalization[obs_name][
                                        "condition"
                                    ]
                                ],
                                self.model.problem.normalization[obs_name]["timepoint"],
                            ]
                        )
                        if self.model.problem.normalization[obs_name]["timepoint"] is not None
                        else np.max(
                            self.model.problem.simulations[
                                i,
                                [
                                    self.model.problem.conditions.index(c)
                                    for c in self.model.problem.normalization[obs_name][
                                        "condition"
                                    ]
                                ],
                            ]
                        )
                    ),
                    color=singleplotting[i].cmap[k]
                    if mode == 0
                    else multiplotting.cmap[multiplotting.observables.index(obs_name)],
                    label=condition if mode == 0 else singleplotting[i].ylabel,
                )

    def _plot_experimental_data_with_error_bars(
        self,
        viz_type: str,
        exp_t: Optional[List[int]],
        obs_name: str,
        mode: int,
        singleplotting: List[SingleObservable],
        multiplotting: MultipleObservables,
    ) -> None:
        """
        Plot experimental measurements with error bars.
        """
        i = self.model.observables.index(obs_name)
        for k, condition in enumerate(self.model.problem.conditions):
            if (
                condition in self.model.problem.experiments[i]
                and (mode == 0 and condition not in singleplotting[i].dont_show)
                or (mode == 1 and condition == multiplotting.condition)
            ):
                try:
                    exp_data = plt.errorbar(
                        np.array(exp_t) / singleplotting[i].divided_by,
                        self.model.problem.experiments[i][condition],
                        yerr=self.model.problem.error_bars[i][condition],
                        color=singleplotting[i].cmap[k]
                        if mode == 0
                        else multiplotting.cmap[multiplotting.observables.index(obs_name)],
                        ecolor=singleplotting[i].cmap[k]
                        if mode == 0
                        else multiplotting.cmap[multiplotting.observables.index(obs_name)],
                        elinewidth=1,
                        capsize=8,
                        markerfacecolor="None",
                        markeredgecolor=singleplotting[i].cmap[k]
                        if mode == 0
                        else multiplotting.cmap[multiplotting.observables.index(obs_name)],
                        fmt=singleplotting[i].shape[k]
                        if mode == 0
                        else multiplotting.shape[multiplotting.observables.index(obs_name)],
                        clip_on=False,
                        label=singleplotting[i].ylabel
                        if mode == 1 and viz_type == "experiment"
                        else None,
                    )
                    for capline in exp_data[1]:
                        capline.set_clip_on(False)
                    for barlinecol in exp_data[2]:
                        barlinecol.set_clip_on(False)
                except ValueError as e:
                    print(e, f"in {obs_name}.")

    def _plot_experimental_data_without_error_bars(
        self,
        viz_type: str,
        exp_t: Optional[List[int]],
        obs_name: str,
        mode: int,
        singleplotting: List[SingleObservable],
        multiplotting: MultipleObservables,
    ) -> None:
        """
        Plot experimental measurements when model.problem.error_bars[i] is None.
        """
        i = self.model.observables.index(obs_name)
        for k, condition in enumerate(self.model.problem.conditions):
            if (
                condition in self.model.problem.experiments[i]
                and (mode == 0 and condition not in singleplotting[i].dont_show)
                or (mode == 1 and condition == multiplotting.condition)
            ):
                try:
                    plt.plot(
                        np.array(exp_t) / singleplotting[i].divided_by,
                        self.model.problem.experiments[i][condition],
                        singleplotting[i].shape[k]
                        if mode == 0
                        else multiplotting.shape[multiplotting.observables.index(obs_name)],
                        markerfacecolor="None",
                        markeredgecolor=singleplotting[i].cmap[k]
                        if mode == 0
                        else multiplotting.cmap[multiplotting.observables.index(obs_name)],
                        color=singleplotting[i].cmap[k]
                        if mode == 0
                        else multiplotting.cmap[multiplotting.observables.index(obs_name)],
                        clip_on=False,
                        label=singleplotting[i].ylabel
                        if mode == 1 and viz_type == "experiment"
                        else None,
                    )
                except ValueError as e:
                    print(e, f"in {obs_name}.")

    def _save_mode_0(
        self,
        obs_name: str,
        singleplotting: List[SingleObservable],
        viz_type: str,
    ) -> None:
        """
        Plot time course of each observable.
        """
        i = self.model.observables.index(obs_name)
        if singleplotting[i].xlim:
            plt.xlim(singleplotting[i].xlim)
        if singleplotting[i].xticks is not None:
            plt.xticks(singleplotting[i].xticks)
        plt.xlabel(singleplotting[i].xlabel)
        if singleplotting[i].ylim:
            plt.ylim(singleplotting[i].ylim)
        if singleplotting[i].yticks is not None:
            plt.yticks(singleplotting[i].yticks)
        plt.ylabel(singleplotting[i].ylabel)
        if singleplotting[i].legend_kws is not None:
            plt.legend(**singleplotting[i].legend_kws)
        plt.savefig(
            os.path.join(
                self.model.path,
                "figure",
                "simulation",
                f"{viz_type}",
                f"{obs_name}",
            ),
        )
        plt.close()

    def _save_mode_1(
        self,
        multiplotting: MultipleObservables,
        viz_type: str,
    ) -> None:
        """
        Plot time course of multiple observables in one figure.
        """
        if multiplotting.xlim:
            plt.xlim(multiplotting.xlim)
        if multiplotting.xticks is not None:
            plt.xticks(multiplotting.xticks)
        plt.xlabel(multiplotting.xlabel)
        if multiplotting.ylim:
            plt.ylim(multiplotting.ylim)
        if multiplotting.yticks is not None:
            plt.yticks(multiplotting.yticks)
        plt.ylabel(multiplotting.ylabel)
        plt.legend(**multiplotting.legend_kws)
        plt.savefig(
            os.path.join(
                self.model.path,
                "figure",
                "simulation",
                f"{viz_type}",
                f"{multiplotting.fname}",
            ),
        )
        plt.close()
