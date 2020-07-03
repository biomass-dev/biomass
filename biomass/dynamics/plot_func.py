import os
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns


def timecourse(sim, n_file, viz_type, show_all, stdev,
               simulations_all, model_path, obs, exp):
    os.makedirs(
        model_path + '/figure/simulation/{}'.format(viz_type), exist_ok=True
    )

    plt.rcParams['font.size'] = 20
    plt.rcParams['axes.linewidth'] = 1.5
    plt.rcParams['xtick.major.width'] = 1.5
    plt.rcParams['ytick.major.width'] = 1.5
    plt.rcParams['lines.linewidth'] = 1.8
    plt.rcParams['lines.markersize'] = 12
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.it'] = 'Arial:italic'

    cmap = ['goldenrod', 'seagreen']
    shape = ['^', 'o']
    if len(cmap) < len(sim.conditions) or len(shape) < len(sim.conditions):
        raise ValueError(
            'len(cmap), len(shape) must be equal to or greater than len(sim.conditions).'
        )

    for i, obs_name in enumerate(obs):

        plt.figure(figsize=(4, 3))
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['top'].set_visible(False)

        if viz_type != 'experiment':
            if show_all:
                for j, _ in enumerate(n_file):
                    for l, _ in enumerate(sim.conditions):
                        plt.plot(
                            sim.t, simulations_all[i, j, :, l] /
                            np.max(simulations_all[i, j, :, :]),
                            color=cmap[l], alpha=0.05
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
                        sim.t, np.nanmean(normalized[i, :, :, l], axis=0),
                        color=cmap[l]
                    )
                if stdev:
                    for l, _ in enumerate(sim.conditions):
                        y_mean = np.nanmean(normalized[i, :, :, l], axis=0)
                        y_std = [
                            np.nanstd(normalized[i, :, k, l], ddof=1)
                            for k, _ in enumerate(sim.t)
                        ]
                        plt.fill_between(
                            sim.t, y_mean - y_std, y_mean + y_std,
                            lw=0, color=cmap[l], alpha=0.1
                        )
            else:
                for l, _ in enumerate(sim.conditions):
                    plt.plot(
                        sim.t, sim.simulations[i, :, l]/np.max(sim.simulations[i]),
                        color=cmap[l]
                    )
        if exp.experiments[i] is not None:
            exp_t = exp.get_timepoint(i)
            if exp.standard_error[i] is not None:
                for l, condition in enumerate(sim.conditions):
                    if condition in exp.experiments[i]:
                        exp_data = plt.errorbar(
                            np.array(exp_t) / 60., exp.experiments[i][condition],
                            yerr=exp.standard_error[i][condition], color=cmap[l],
                            ecolor=cmap[l], elinewidth=1, capsize=8,
                            markerfacecolor='None', markeredgecolor=cmap[l],
                            fmt=shape[l], clip_on=False
                        )
                        for capline in exp_data[1]:
                            capline.set_clip_on(False)
                        for barlinecol in exp_data[2]:
                            barlinecol.set_clip_on(False)
            else:
                for l, condition in enumerate(sim.conditions):
                    if condition in exp.experiments[i]:
                        plt.plot(
                            np.array(exp_t) / 60., exp.experiments[i][condition],
                            shape[l], markerfacecolor='None', markeredgecolor=cmap[l],
                            color=cmap[l], clip_on=False
                        )
        plt.xlim(-5, 95)
        plt.xticks([0, 30, 60, 90])
        plt.yticks([0, 0.3, 0.6, 0.9, 1.2])
        plt.ylim(-0.1, 1.3)
        plt.xlabel('Time (min)')
        plt.ylabel(obs_name.replace('__', '\n').replace('_', ' '))

        plt.savefig(
            model_path + '/figure/simulation/{}/{}.pdf'.format(
                viz_type, obs_name
            ), bbox_inches='tight'
        )
        plt.close()


def param_range(search_idx, popt, model_path, parameters, species, sp, portrait):
    os.makedirs(model_path + '/figure/param_range', exist_ok=True)

    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.linewidth'] = 1.2
    plt.rcParams['xtick.major.width'] = 1.2
    plt.rcParams['ytick.major.width'] = 1.2
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.it'] = 'Arial:italic'

    if portrait:
        if len(sp.idx_params) > 0:
            fig = plt.figure(
                figsize=(8, len(sp.idx_params) / 2.5)
            )
            ax = sns.boxenplot(
                data=popt[:, :len(sp.idx_params)],
                orient='h',
                linewidth=0.5,
                palette='Set2'
            )
            sns.despine()
            ax.set_xlabel('Parameter value')
            ax.set_ylabel('')
            ax.set_yticklabels([parameters[i] for i in sp.idx_params])
            ax.set_xscale('log')

            plt.savefig(
                model_path + '/figure/param_range/param_range.pdf',
                bbox_inches='tight'
            )
            plt.close(fig)
        if len(sp.idx_initials) > 0:
            fig = plt.figure(
                figsize=(8, len(sp.idx_initials) / 2.5)
            )
            ax = sns.boxenplot(
                data=popt[:, len(sp.idx_params):],
                orient='h',
                linewidth=0.5,
                palette='Set2'
            )
            sns.despine()
            ax.set_xlabel('Initial value')
            ax.set_ylabel('')
            ax.set_yticklabels([species[i] for i in sp.idx_initials])
            ax.set_xscale('log')

            plt.savefig(
                model_path + '/figure/param_range/initial_value_range.pdf',
                bbox_inches='tight'
            )
            plt.close(fig)
    else:
        if len(sp.idx_params) > 0:
            fig = plt.figure(
                figsize=(len(sp.idx_params) / 2.2, 6)
            )
            ax = sns.boxenplot(
                data=popt[:, :len(sp.idx_params)],
                linewidth=0.5,
                palette='Set2'
            )
            sns.despine()
            ax.set_xlabel('')
            ax.set_xticklabels(
                [parameters[i] for i in sp.idx_params], rotation=45
            )
            ax.set_ylabel('Parameter value')
            ax.set_yscale('log')

            plt.savefig(
                model_path + '/figure/param_range/param_range_h.pdf',
                bbox_inches='tight'
            )
            plt.close(fig)
        if len(sp.idx_initials) > 0:
            fig = plt.figure(
                figsize=(len(sp.idx_initials) / 2.2, 6)
            )
            ax = sns.boxenplot(
                data=popt[:, len(sp.idx_params):],
                linewidth=0.5,
                palette='Set2'
            )
            ax.set_xlabel('')
            ax.set_xticklabels(
                [species[i] for i in sp.idx_initials], rotation=45
            )
            sns.despine()
            ax.set_ylabel('Initial value')
            ax.set_yscale('log')

            plt.savefig(
                model_path + '/figure/param_range/initail_value_range_h.pdf',
                bbox_inches='tight'
            )
            plt.close(fig)
