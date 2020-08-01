import os
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

from biomass.exec_model import ExecModel


class TemporalDynamics(ExecModel):
    def __init__(self, model):
        super().__init__(model)

    def plot_timecourse(self, sim, n_file, viz_type, show_all, stdev, simulations_all):
        os.makedirs(
            self.model_path
            + '/figure/simulation/{}'.format(viz_type), exist_ok=True
        )
        self.viz.set_timecourse_rcParams()
        options = self.viz.get_timecourse_options()
        multiple = self.viz.multiplot_observables()

        for rule in ['default', 'custom']:
            for i, obs_name in enumerate(self.obs):
                if len(options[i]['cmap']) < len(sim.conditions) or \
                        len(options[i]['shape']) < len(sim.conditions):
                    raise ValueError(
                        'len(cmap), len(shape) must be equal to'
                        ' or greater than len(sim.conditions).'
                    )
                if rule == 'custom' and \
                        obs_name not in multiple['observables']:
                    continue
                if rule == 'default':
                    plt.figure(figsize=(4, 3))
                plt.gca().spines['right'].set_visible(False)
                plt.gca().spines['top'].set_visible(False)
                if viz_type != 'experiment':
                    if show_all:
                        for j, _ in enumerate(n_file):
                            for l, condition in enumerate(sim.conditions):
                                if (rule == 'default' and condition not in 
                                        options[i]['dont_show']) or (
                                        rule == 'custom' and
                                        condition == multiple['condition']):
                                    plt.plot(
                                        np.array(sim.t) / 
                                            options[i]['divided_by'],
                                        simulations_all[i, j, :, l] / (
                                            1 if not sim.normalization else
                                            np.max(
                                                simulations_all[i, j, :, :]
                                            )
                                        ),
                                        color=options[i]['cmap'][l] \
                                            if rule == 'default' \
                                            else multiple['cmap'][
                                                multiple['observables'].index(
                                                    obs_name
                                                )
                                            ],
                                        alpha=0.05
                                    )
                    if viz_type == 'average':
                        normalized = np.empty_like(simulations_all)
                        for j, _ in enumerate(n_file):
                            for l, condition in enumerate(sim.conditions):
                                if (rule == 'default' and condition not in 
                                        options[i]['dont_show']) or (
                                            rule == 'custom' and
                                            condition == multiple['condition']):
                                    normalized[i, j, :, l] = (
                                        simulations_all[i, j, :, l] / (
                                            1 if not sim.normalization else
                                            np.max(
                                                simulations_all[i, j, :, :]
                                            )
                                        )
                                    )
                        normalized[i, :, :, :] /= \
                            1 if not sim.normalization else np.max(
                                np.nanmean(normalized[i, :, :, :], axis=0)
                            )
                        for l, condition in enumerate(sim.conditions):
                            if (rule == 'default' and condition not in
                                    options[i]['dont_show']) or (
                                        rule == 'custom' and
                                        condition == multiple['condition']):
                                plt.plot(
                                    np.array(sim.t) / options[i]['divided_by'], 
                                    np.nanmean(normalized[i, :, :, l], axis=0),
                                    color=options[i]['cmap'][l] \
                                        if rule == 'default' \
                                        else multiple['cmap'][
                                            multiple['observables'].index(
                                                obs_name
                                            )
                                        ],
                                    label = options[i]['ylabel']
                                )
                        if stdev:
                            for l, condition in enumerate(sim.conditions):
                                if (rule == 'default' and condition not in 
                                        options[i]['dont_show']) or (
                                            rule == 'custom' and
                                            condition == multiple['condition']):
                                    y_mean = np.nanmean(
                                        normalized[i, :, :, l], axis=0
                                    )
                                    y_std = [
                                        np.nanstd(
                                            normalized[i, :, k, l], ddof=1
                                        ) for k, _ in enumerate(sim.t)
                                    ]
                                    plt.fill_between(
                                        np.array(sim.t) / 
                                            options[i]['divided_by'], 
                                        y_mean - y_std, y_mean + y_std,
                                        lw=0, color=options[i]['cmap'][l] \
                                            if rule == 'default' \
                                            else multiple['cmap'][
                                                multiple['observables'].index(
                                                    obs_name
                                                )
                                            ],
                                        alpha=0.1
                                    )
                    else:
                        for l, condition in enumerate(sim.conditions):
                            if (rule == 'default' and condition not in
                                    options[i]['dont_show']) or (
                                        rule == 'custom' and
                                        condition == multiple['condition']):
                                plt.plot(
                                    np.array(sim.t) / options[i]['divided_by'], 
                                    sim.simulations[i, :, l] / (
                                        1 if not sim.normalization else
                                        np.max(
                                            sim.simulations[i]
                                        )
                                    ),
                                    color=options[i]['cmap'][l] \
                                        if rule == 'default' \
                                        else multiple['cmap'][
                                            multiple['observables'].index(
                                                obs_name
                                            )
                                        ],
                                    label = options[i]['ylabel']
                                )
                if options[i]['exp_data'] and self.exp.experiments[i] is not None:
                    exp_t = self.exp.get_timepoint(i)
                    if self.exp.error_bar[i] is not None:
                        for l, condition in enumerate(sim.conditions):
                            if condition in self.exp.experiments[i] and \
                                    (rule == 'default' and condition not in
                                    options[i]['dont_show']) or (
                                        rule == 'custom' and
                                        condition == multiple['condition']):
                                exp_data = plt.errorbar(
                                    np.array(exp_t) / options[i]['divided_by'],
                                    self.exp.experiments[i][condition],
                                    yerr=self.exp.error_bar[i][condition], 
                                    color=options[i]['cmap'][l] \
                                        if rule == 'default' \
                                        else multiple['cmap'][
                                            multiple['observables'].index(
                                                obs_name
                                            )
                                        ],
                                    ecolor=options[i]['cmap'][l] \
                                        if rule == 'default' \
                                        else multiple['cmap'][
                                            multiple['observables'].index(
                                                obs_name
                                            )
                                        ],
                                    elinewidth=1, capsize=8,
                                    markerfacecolor='None',
                                    markeredgecolor=options[i]['cmap'][l] \
                                        if rule == 'default' \
                                        else multiple['cmap'][
                                            multiple['observables'].index(
                                                obs_name
                                            )
                                        ],
                                    fmt=options[i]['shape'][l] \
                                        if rule == 'default' \
                                        else multiple['shape'][
                                            multiple['observables'].index(
                                                obs_name
                                            )
                                        ],
                                    clip_on=False
                                )
                                for capline in exp_data[1]:
                                    capline.set_clip_on(False)
                                for barlinecol in exp_data[2]:
                                    barlinecol.set_clip_on(False)
                    else:
                        for l, condition in enumerate(sim.conditions):
                            if condition in self.exp.experiments[i] and \
                                    (rule == 'default' and condition not in
                                    options[i]['dont_show']) or (
                                        rule == 'custom' and
                                        condition == multiple['condition']):
                                plt.plot(
                                    np.array(exp_t) / options[i]['divided_by'], 
                                    self.exp.experiments[i][condition],
                                    fmt=options[i]['shape'][l] \
                                        if rule == 'default' \
                                        else multiple['shape'][
                                            multiple['observables'].index(
                                                obs_name
                                            )
                                        ], 
                                    markerfacecolor='None', 
                                    markeredgecolor=options[i]['cmap'][l] \
                                        if rule == 'default' \
                                        else multiple['cmap'][
                                            multiple['observables'].index(
                                                obs_name
                                            )
                                        ],
                                    color=options[i]['cmap'][l] \
                                        if rule == 'default' \
                                        else multiple['cmap'][
                                            multiple['observables'].index(
                                                obs_name
                                            )
                                        ], 
                                    clip_on=False
                                )
                if rule == 'default':
                    if options[i]['xlim']:
                        plt.xlim(options[i]['xlim'])
                    if options[i]['xticks'] is not None:
                        plt.xticks(options[i]['xticks'])
                    plt.xlabel(options[i]['xlabel'])
                    if options[i]['ylim']:
                        plt.ylim(options[i]['ylim'])
                    if options[i]['yticks'] is not None:
                        plt.yticks(options[i]['yticks'])
                    plt.ylabel(options[i]['ylabel'])
                    plt.savefig(
                        self.model_path 
                        + '/figure/simulation/{}/{}.pdf'.format(
                            viz_type, obs_name
                        ), bbox_inches='tight'
                    )
                    plt.close()
            if rule == 'custom' and multiple['observables']:
                if multiple['xlim']:
                    plt.xlim(multiple['xlim'])
                if multiple['xticks'] is not None:
                    plt.xticks(multiple['xticks'])
                plt.xlabel(multiple['xlabel'])
                if multiple['ylim']:
                    plt.ylim(multiple['ylim'])
                if multiple['yticks'] is not None:
                    plt.yticks(multiple['yticks'])
                plt.ylabel(multiple['ylabel'])
                plt.legend(
                    bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0,
                    labelspacing=1.25, frameon = False, fontsize=12
                )
                plt.savefig(
                    self.model_path 
                    + '/figure/simulation/{}/multiple_observables.pdf'.format(
                        viz_type
                    ), bbox_inches='tight'
                )
                plt.close()


    def plot_param_range(self, popt, portrait):
        os.makedirs(self.model_path + '/figure/param_range', exist_ok=True)

        self.viz.set_param_range_rcParams()

        if len(self.sp.idx_params) > 0:
            fig = plt.figure(
                figsize=(8, len(self.sp.idx_params) / 2.5) if portrait else
                        (len(self.sp.idx_params) / 2.2, 6)
            )
            ax = sns.boxenplot(
                data=popt[:, :len(self.sp.idx_params)],
                orient='h' if portrait else 'v',
                linewidth=0.5,
                palette='Set2'
            )
            sns.despine()
            if portrait:
                ax.set_xlabel('Parameter value')
                ax.set_xscale('log')
                ax.set_ylabel('')
                ax.set_yticklabels(
                    [self.parameters[i] for i in self.sp.idx_params]
                )
            else:
                ax.set_xlabel('')
                ax.set_xticklabels(
                    [self.parameters[i] for i in self.sp.idx_params],
                    rotation=90
                )
                ax.set_ylabel('Parameter value')
                ax.set_yscale('log')
            plt.savefig(
                self.model_path 
                + '/figure/param_range/estimated_param_values_range.pdf',
                bbox_inches='tight'
            )
            plt.close(fig)
        if len(self.sp.idx_initials) > 0:
            fig = plt.figure(
                figsize=(8, len(self.sp.idx_initials) / 2.5) if portrait else
                        (len(self.sp.idx_initials) / 2.2, 6)
            )
            ax = sns.boxenplot(
                data=popt[:, len(self.sp.idx_params):],
                orient='h' if portrait else 'v',
                linewidth=0.5,
                palette='Set2'
            )
            sns.despine()
            if portrait:
                ax.set_xlabel('Initial value')
                ax.set_xscale('log')
                ax.set_ylabel('')
                ax.set_yticklabels(
                    [self.species[i] for i in self.sp.idx_initials]
                )
            else:
                ax.set_xlabel('')
                ax.set_yticklabels(
                    [self.species[i] for i in self.sp.idx_initials],
                    rotation=90
                )
                ax.set_ylabel('Initial value')
                ax.set_yscale('log')
            plt.savefig(
                self.model_path 
                + '/figure/param_range/estimated_initial_values_range.pdf',
                bbox_inches='tight'
            )
            plt.close(fig)
