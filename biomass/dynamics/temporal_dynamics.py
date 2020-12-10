import os
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

from biomass.exec_model import ExecModel


class TemporalDynamics(ExecModel):
    def __init__(self, model):
        super().__init__(model)

    def plot_timecourse(self, sim, n_file, viz_type, show_all, stdev, simulations_all):
        """
        Plot time course of each observable.

        Parameters
        ----------
        sim : NumericalSimulation
            Information about simulation.

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
        os.makedirs(self.model_path + "/figure/simulation/{}".format(viz_type), exist_ok=True)
        self.exp.set_data()
        self.viz.set_timecourse_rcParams()
        timecourse = self.viz.get_timecourse_options()
        multiplot = self.viz.multiplot_observables()

        for mode in range(2):
            # mode 0 : timecourse_for_each_observable
            # mode 1 : multiplot_observables
            set_fig = False
            for i, obs_name in enumerate(self.obs):
                if len(timecourse[i]["cmap"]) < len(sim.conditions) or len(timecourse[i]["shape"]) < len(
                    sim.conditions
                ):
                    raise ValueError("len(cmap), len(shape) must be equal to" " or greater than len(sim.conditions).")
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
                        for j, _ in enumerate(n_file):
                            for l, condition in enumerate(sim.conditions):
                                if (mode == 0 and condition not in timecourse[i]["dont_show"]) or (
                                    mode == 1 and condition == multiplot["condition"]
                                ):
                                    plt.plot(
                                        np.array(sim.t) / timecourse[i]["divided_by"],
                                        simulations_all[i, j, :, l]
                                        / (
                                            1
                                            if not sim.normalization or np.max(simulations_all[i, j, :, l]) == 0.0
                                            else np.max(
                                                simulations_all[
                                                    i,
                                                    j,
                                                    sim.normalization[obs_name]["timepoint"],
                                                    [
                                                        sim.conditions.index(c)
                                                        for c in sim.normalization[obs_name]["condition"]
                                                    ],
                                                ]
                                            )
                                            if sim.normalization[obs_name]["timepoint"] is not None
                                            else np.max(
                                                simulations_all[
                                                    i,
                                                    j,
                                                    :,
                                                    [
                                                        sim.conditions.index(c)
                                                        for c in sim.normalization[obs_name]["condition"]
                                                    ],
                                                ]
                                            )
                                        ),
                                        color=timecourse[i]["cmap"][l]
                                        if mode == 0
                                        else multiplot["cmap"][multiplot["observables"].index(obs_name)],
                                        alpha=0.05,
                                    )
                    if viz_type == "average":
                        normalized = np.empty_like(simulations_all)
                        for j, _ in enumerate(n_file):
                            for l, condition in enumerate(sim.conditions):
                                if (mode == 0 and condition not in timecourse[i]["dont_show"]) or (
                                    mode == 1 and condition == multiplot["condition"]
                                ):
                                    normalized[i, j, :, l] = simulations_all[i, j, :, l] / (
                                        1
                                        if not sim.normalization or np.max(simulations_all[i, j, :, l]) == 0.0
                                        else np.max(
                                            simulations_all[
                                                i,
                                                j,
                                                sim.normalization[obs_name]["timepoint"],
                                                [
                                                    sim.conditions.index(c)
                                                    for c in sim.normalization[obs_name]["condition"]
                                                ],
                                            ]
                                        )
                                        if sim.normalization[obs_name]["timepoint"] is not None
                                        else np.max(
                                            simulations_all[
                                                i,
                                                j,
                                                :,
                                                [
                                                    sim.conditions.index(c)
                                                    for c in sim.normalization[obs_name]["condition"]
                                                ],
                                            ]
                                        )
                                    )
                        if sim.normalization and sim.normalization[obs_name]["timepoint"] is None:
                            mean_vec = []
                            for c in sim.normalization[obs_name]["condition"]:
                                mean_vec.append(np.nanmean(normalized[i, :, :, sim.conditions.index(c)], axis=0))
                            norm_max = np.max(mean_vec)
                            if not np.isnan(norm_max) and norm_max != 0.0:
                                normalized[i, :, :, :] /= norm_max
                        for l, condition in enumerate(sim.conditions):
                            if (mode == 0 and condition not in timecourse[i]["dont_show"]) or (
                                mode == 1 and condition == multiplot["condition"]
                            ):
                                plt.plot(
                                    np.array(sim.t) / timecourse[i]["divided_by"],
                                    np.nanmean(normalized[i, :, :, l], axis=0),
                                    color=timecourse[i]["cmap"][l]
                                    if mode == 0
                                    else multiplot["cmap"][multiplot["observables"].index(obs_name)],
                                    label=condition if mode == 0 else timecourse[i]["ylabel"],
                                )
                        if stdev:
                            for l, condition in enumerate(sim.conditions):
                                if (mode == 0 and condition not in timecourse[i]["dont_show"]) or (
                                    mode == 1 and condition == multiplot["condition"]
                                ):
                                    y_mean = np.nanmean(normalized[i, :, :, l], axis=0)
                                    y_std = [np.nanstd(normalized[i, :, k, l], ddof=1) for k, _ in enumerate(sim.t)]
                                    plt.fill_between(
                                        np.array(sim.t) / timecourse[i]["divided_by"],
                                        y_mean - y_std,
                                        y_mean + y_std,
                                        lw=0,
                                        color=timecourse[i]["cmap"][l]
                                        if mode == 0
                                        else multiplot["cmap"][multiplot["observables"].index(obs_name)],
                                        alpha=0.1,
                                    )
                    else:
                        for l, condition in enumerate(sim.conditions):
                            if (mode == 0 and condition not in timecourse[i]["dont_show"]) or (
                                mode == 1 and condition == multiplot["condition"]
                            ):
                                plt.plot(
                                    np.array(sim.t) / timecourse[i]["divided_by"],
                                    sim.simulations[i, :, l]
                                    / (
                                        1
                                        if not sim.normalization or np.max(sim.simulations[i, :, l]) == 0.0
                                        else np.max(
                                            sim.simulations[
                                                i,
                                                sim.normalization[obs_name]["timepoint"],
                                                [
                                                    sim.conditions.index(c)
                                                    for c in sim.normalization[obs_name]["condition"]
                                                ],
                                            ]
                                        )
                                        if sim.normalization[obs_name]["timepoint"] is not None
                                        else np.max(
                                            sim.simulations[
                                                i,
                                                :,
                                                [
                                                    sim.conditions.index(c)
                                                    for c in sim.normalization[obs_name]["condition"]
                                                ],
                                            ]
                                        )
                                    ),
                                    color=timecourse[i]["cmap"][l]
                                    if mode == 0
                                    else multiplot["cmap"][multiplot["observables"].index(obs_name)],
                                    label=condition if mode == 0 else timecourse[i]["ylabel"],
                                )
                if timecourse[i]["exp_data"] and self.exp.experiments[i] is not None:
                    exp_t = self.exp.get_timepoint(obs_name)
                    if self.exp.error_bars[i] is not None:
                        for l, condition in enumerate(sim.conditions):
                            if (
                                condition in self.exp.experiments[i]
                                and (mode == 0 and condition not in timecourse[i]["dont_show"])
                                or (mode == 1 and condition == multiplot["condition"])
                            ):
                                exp_data = plt.errorbar(
                                    np.array(exp_t) / timecourse[i]["divided_by"],
                                    self.exp.experiments[i][condition],
                                    yerr=self.exp.error_bars[i][condition],
                                    color=timecourse[i]["cmap"][l]
                                    if mode == 0
                                    else multiplot["cmap"][multiplot["observables"].index(obs_name)],
                                    ecolor=timecourse[i]["cmap"][l]
                                    if mode == 0
                                    else multiplot["cmap"][multiplot["observables"].index(obs_name)],
                                    elinewidth=1,
                                    capsize=8,
                                    markerfacecolor="None",
                                    markeredgecolor=timecourse[i]["cmap"][l]
                                    if mode == 0
                                    else multiplot["cmap"][multiplot["observables"].index(obs_name)],
                                    fmt=timecourse[i]["shape"][l]
                                    if mode == 0
                                    else multiplot["shape"][multiplot["observables"].index(obs_name)],
                                    clip_on=False,
                                    label=timecourse[i]["ylabel"] if mode == 1 and viz_type == "experiment" else None,
                                )
                                for capline in exp_data[1]:
                                    capline.set_clip_on(False)
                                for barlinecol in exp_data[2]:
                                    barlinecol.set_clip_on(False)
                    else:
                        for l, condition in enumerate(sim.conditions):
                            if (
                                condition in self.exp.experiments[i]
                                and (mode == 0 and condition not in timecourse[i]["dont_show"])
                                or (mode == 1 and condition == multiplot["condition"])
                            ):
                                plt.plot(
                                    np.array(exp_t) / timecourse[i]["divided_by"],
                                    self.exp.experiments[i][condition],
                                    timecourse[i]["shape"][l]
                                    if mode == 0
                                    else multiplot["shape"][multiplot["observables"].index(obs_name)],
                                    markerfacecolor="None",
                                    markeredgecolor=timecourse[i]["cmap"][l]
                                    if mode == 0
                                    else multiplot["cmap"][multiplot["observables"].index(obs_name)],
                                    color=timecourse[i]["cmap"][l]
                                    if mode == 0
                                    else multiplot["cmap"][multiplot["observables"].index(obs_name)],
                                    clip_on=False,
                                    label=timecourse[i]["ylabel"] if mode == 1 and viz_type == "experiment" else None,
                                )
                if mode == 0:
                    if timecourse[i]["xlim"]:
                        plt.xlim(timecourse[i]["xlim"])
                    if timecourse[i]["xticks"] is not None:
                        plt.xticks(timecourse[i]["xticks"])
                    plt.xlabel(timecourse[i]["xlabel"])
                    if timecourse[i]["ylim"]:
                        plt.ylim(timecourse[i]["ylim"])
                    if timecourse[i]["yticks"] is not None:
                        plt.yticks(timecourse[i]["yticks"])
                    plt.ylabel(timecourse[i]["ylabel"])
                    if timecourse[i]["legend_loc"] is not None:
                        plt.legend(loc=timecourse[i]["legend_loc"], frameon=False, fontsize=12)
                    plt.savefig(self.model_path + f"/figure/simulation/{viz_type}/{obs_name}.pdf", bbox_inches="tight")
                    plt.close()
            if mode == 1 and multiplot["observables"]:
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
                    self.model_path + f'/figure/simulation/{viz_type}/{multiplot["fig_name"]}.pdf', bbox_inches="tight"
                )
                plt.close()

    def plot_param_range(self, popt: np.ndarray, portrait: bool):
        """
        Plot estimated parameter/initial values.

        Parameters
        ----------
        popt : numpy array
            Array containing all optimized parameter/initial values.

        portrait : bool
        """
        os.makedirs(self.model_path + "/figure/param_range", exist_ok=True)

        self.viz.set_param_range_rcParams()

        for val in ["parameter_value", "initial_value"]:
            if (val == "parameter_value" and len(self.sp.idx_params) == 0) or (
                val == "initial_value" and len(self.sp.idx_initials) == 0
            ):
                continue
            else:
                fig = plt.figure(
                    figsize=(
                        8,
                        (len(self.sp.idx_params) if val == "parameter_value" else len(self.sp.idx_initials)) / 2.5,
                    )
                    if portrait
                    else ((len(self.sp.idx_params) if val == "parameter_value" else len(self.sp.idx_initials)) / 2.2, 6)
                )
                ax = sns.boxenplot(
                    data=popt[:, : len(self.sp.idx_params)]
                    if val == "parameter_value"
                    else popt[:, len(self.sp.idx_params) :],
                    orient="h" if portrait else "v",
                    linewidth=0.5,
                    palette="Set2",
                )
                sns.despine()
                if portrait:
                    ax.set_xlabel(val.capitalize().replace("_", " "))
                    ax.set_xscale("log")
                    ax.set_ylabel("")
                    ax.set_yticklabels(
                        [self.parameters[i] for i in self.sp.idx_params]
                        if val == "parameter_value"
                        else [self.species[i] for i in self.sp.idx_initials]
                    )
                else:
                    ax.set_xlabel("")
                    ax.set_xticklabels(
                        [self.parameters[i] for i in self.sp.idx_params]
                        if val == "parameter_value"
                        else [self.species[i] for i in self.sp.idx_initials],
                        rotation=90,
                    )
                    ax.set_ylabel(val.capitalize().replace("_", " "))
                    ax.set_yscale("log")
                plt.savefig(
                    self.model_path + "/figure/param_range/" + "estimated_{}s.pdf".format(val), bbox_inches="tight"
                )
                plt.close(fig)
