import os
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns


class PlotFunc(object):
    def __init__(self, model):
        self.model_path = model.__path__[0]
        self.parameters = model.C.NAMES
        self.species = model.V.NAMES
        self.obs = model.observables
        self.viz = model.Visualization()
        self.exp = model.ExperimentalData()
        self.sp = model.SearchParam()

    def timecourse(self, sim, n_file, viz_type, show_all, stdev, simulations_all):
        os.makedirs(
            self.model_path
            + '/figure/simulation/{}'.format(viz_type), exist_ok=True
        )
        self.viz.set_timecourse_rcParams()
        options = self.viz.get_timecourse_options()

        for i, obs_name in enumerate(self.obs):
            if len(options[i]['cmap']) < len(sim.conditions) or \
                    len(options[i]['shape']) < len(sim.conditions):
                raise ValueError(
                    'len(cmap), len(shape) must be equal to'
                    ' or greater than len(sim.conditions).'
                )
            plt.figure(figsize=(4, 3))
            plt.gca().spines['right'].set_visible(False)
            plt.gca().spines['top'].set_visible(False)

            if viz_type != 'experiment':
                if show_all:
                    for j, _ in enumerate(n_file):
                        for l, _ in enumerate(sim.conditions):
                            plt.plot(
                                np.array(sim.t) / options[i]['divided_by'],
                                simulations_all[i, j, :, l] /
                                    np.max(simulations_all[i, j, :, :]),
                                color=options[i]['cmap'][l],
                                alpha=0.05
                            )
                if viz_type == 'average':
                    normalized = np.empty_like(simulations_all)
                    for j, _ in enumerate(n_file):
                        for l, _ in enumerate(sim.conditions):
                            normalized[i, j, :, l] = (
                                simulations_all[i, j, :, l] /
                                np.max(simulations_all[i, j, :, :])
                            )
                    normalized[i, :, :, :] = normalized[i, :, :, :] / np.max(
                        np.nanmean(normalized[i, :, :, :], axis=0)
                    )
                    for l, _ in enumerate(sim.conditions):
                        plt.plot(
                            np.array(sim.t) / options[i]['divided_by'], 
                            np.nanmean(normalized[i, :, :, l], axis=0),
                            color=options[i]['cmap'][l]
                        )
                    if stdev:
                        for l, _ in enumerate(sim.conditions):
                            y_mean = np.nanmean(normalized[i, :, :, l], axis=0)
                            y_std = [
                                np.nanstd(normalized[i, :, k, l], ddof=1)
                                for k, _ in enumerate(sim.t)
                            ]
                            plt.fill_between(
                                np.array(sim.t) / options[i]['divided_by'], 
                                y_mean - y_std, y_mean + y_std,
                                lw=0, color=options[i]['cmap'][l],
                                alpha=0.1
                            )
                else:
                    for l, _ in enumerate(sim.conditions):
                        plt.plot(
                            np.array(sim.t) / options[i]['divided_by'], 
                            sim.simulations[i, :, l]/np.max(sim.simulations[i]),
                            color=options[i]['cmap'][l]
                        )
            if self.exp.experiments[i] is not None:
                exp_t = self.exp.get_timepoint(i)
                if self.exp.standard_error[i] is not None:
                    for l, condition in enumerate(sim.conditions):
                        if condition in self.exp.experiments[i]:
                            exp_data = plt.errorbar(
                                np.array(exp_t) / options[i]['divided_by'],
                                self.exp.experiments[i][condition],
                                yerr=self.exp.standard_error[i][condition], 
                                color=options[i]['cmap'][l],
                                ecolor=options[i]['cmap'][l],
                                elinewidth=1, capsize=8, markerfacecolor='None',
                                markeredgecolor=options[i]['cmap'][l],
                                fmt=options[i]['shape'][l],
                                clip_on=False
                            )
                            for capline in exp_data[1]:
                                capline.set_clip_on(False)
                            for barlinecol in exp_data[2]:
                                barlinecol.set_clip_on(False)
                else:
                    for l, condition in enumerate(sim.conditions):
                        if condition in self.exp.experiments[i]:
                            plt.plot(
                                np.array(exp_t) / options[i]['divided_by'], 
                                self.exp.experiments[i][condition],
                                options[i]['shape'][l], 
                                markerfacecolor='None', 
                                markeredgecolor=options[i]['cmap'][l],
                                color=options[i]['cmap'][l], 
                                clip_on=False
                            )
            if options[i]['xlim']:
                plt.xlim(*options[i]['xlim'])
            if options[i]['xticks']:
                plt.xticks(options[i]['xticks'])
            if options[i]['xlabel'] is not None:
                plt.xlabel(options[i]['xlabel'])
            if options[i]['ylim']:
                plt.ylim(*options[i]['ylim'])
            if options[i]['yticks']:
                plt.yticks(options[i]['yticks'])
            plt.ylabel(options[i]['ylabel'])
            plt.savefig(
                self.model_path + '/figure/simulation/{}/{}.pdf'.format(
                    viz_type, obs_name
                ), bbox_inches='tight'
            )
            plt.close()


    def param_range(self, popt, portrait):
        os.makedirs(self.model_path + '/figure/param_range', exist_ok=True)

        self.viz.set_param_range_rcParams()

        if portrait:
            if len(self.sp.idx_params) > 0:
                fig = plt.figure(
                    figsize=(8, len(self.sp.idx_params) / 2.5)
                )
                ax = sns.boxenplot(
                    data=popt[:, :len(self.sp.idx_params)],
                    orient='h',
                    linewidth=0.5,
                    palette='Set2'
                )
                sns.despine()
                ax.set_xlabel('Parameter value')
                ax.set_ylabel('')
                ax.set_yticklabels(
                    [self.parameters[i] for i in self.sp.idx_params]
                )
                ax.set_xscale('log')

                plt.savefig(
                    self.model_path 
                    + '/figure/param_range/param_range.pdf',
                    bbox_inches='tight'
                )
                plt.close(fig)
            if len(self.sp.idx_initials) > 0:
                fig = plt.figure(
                    figsize=(8, len(self.sp.idx_initials) / 2.5)
                )
                ax = sns.boxenplot(
                    data=popt[:, len(self.sp.idx_params):],
                    orient='h',
                    linewidth=0.5,
                    palette='Set2'
                )
                sns.despine()
                ax.set_xlabel('Initial value')
                ax.set_ylabel('')
                ax.set_yticklabels(
                    [self.species[i] for i in self.sp.idx_initials]
                )
                ax.set_xscale('log')

                plt.savefig(
                    self.model_path 
                    + '/figure/param_range/initial_value_range.pdf',
                    bbox_inches='tight'
                )
                plt.close(fig)
        else:
            if len(self.sp.idx_params) > 0:
                fig = plt.figure(
                    figsize=(len(self.sp.idx_params) / 2.2, 6)
                )
                ax = sns.boxenplot(
                    data=popt[:, :len(self.sp.idx_params)],
                    linewidth=0.5,
                    palette='Set2'
                )
                sns.despine()
                ax.set_xlabel('')
                ax.set_xticklabels(
                    [self.parameters[i] for i in self.sp.idx_params], rotation=45
                )
                ax.set_ylabel('Parameter value')
                ax.set_yscale('log')

                plt.savefig(
                    self.model_path 
                    + '/figure/param_range/param_range_h.pdf',
                    bbox_inches='tight'
                )
                plt.close(fig)
            if len(self.sp.idx_initials) > 0:
                fig = plt.figure(
                    figsize=(len(self.sp.idx_initials) / 2.2, 6)
                )
                ax = sns.boxenplot(
                    data=popt[:, len(self.sp.idx_params):],
                    linewidth=0.5,
                    palette='Set2'
                )
                ax.set_xlabel('')
                ax.set_xticklabels(
                    [self.species[i] for i in self.sp.idx_initials], rotation=45
                )
                sns.despine()
                ax.set_ylabel('Initial value')
                ax.set_yscale('log')

                plt.savefig(
                    self.model_path 
                    + '/figure/param_range/initail_value_range_h.pdf',
                    bbox_inches='tight'
                )
                plt.close(fig)
