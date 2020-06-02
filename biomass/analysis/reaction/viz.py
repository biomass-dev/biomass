import os
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

from biomass.current_model import observables, NumericalSimulation


def barplot_sensitivity(metric, sensitivity_coefficients, biological_processes, 
                        n_reaction, sort_idx, reaction_indices):
    sim = NumericalSimulation()
    width = 0.3

    # rcParams
    plt.rcParams['font.size'] = 15
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.it'] = 'Arial:italic'
    plt.rcParams['axes.linewidth'] = 1.2
    plt.rcParams['xtick.major.width'] = 1.2
    plt.rcParams['ytick.major.width'] = 1.2

    colors = ['mediumblue', 'red']
    if len(colors) < len(sim.conditions):
        raise ValueError(
            'len(colors) must be equal to or greater than len(sim.conditions).'
        )
    for k, obs_name in enumerate(observables):
        plt.figure(figsize=(12, 5))
        # draw_vertical_span
        if len(biological_processes) > 1:
            left_end = 0
            for i, ith_module in enumerate(biological_processes):
                if i % 2 == 0:
                    plt.axvspan(
                        left_end - width,
                        left_end - width + len(ith_module),
                        facecolor='k', alpha=0.1
                    )
                left_end += len(ith_module)
        plt.hlines(
            [0], -width, n_reaction-1-width, 'k', lw=1
        )
        sensitivity_array = sensitivity_coefficients[:, :, k, :]
        nan_idx = []
        for i in range(sensitivity_array.shape[0]):
            for j in range(sensitivity_array.shape[1]):
                if any(np.isnan(sensitivity_array[i, j, :])):
                    nan_idx.append(i)
        sensitivity_array = np.delete(
            sensitivity_array, nan_idx, axis=0
        )
        if sensitivity_array.size != 0:
            average = np.mean(sensitivity_array, axis=0)
            stdev = np.std(sensitivity_array, axis=0, ddof=1)
            for l, condition in enumerate(sim.conditions):
                plt.bar(
                    np.arange(n_reaction) + l * width,
                    average[sort_idx, l], yerr=stdev[sort_idx, l],
                    ecolor=colors[l], capsize=2, width=width, color=colors[l],
                    align='center', label=condition
                )
            distance = np.max(average) * 0.05
            for i, j in enumerate(sort_idx):
                if j != 0:
                    xp = i + width/2
                    yp = average[j, np.argmax(np.abs(average[j, :]))]
                    yerr = stdev[j, np.argmax(stdev[j, :])]
                    if yp > 0:
                        plt.text(
                            xp, yp + yerr + distance, reaction_indices[i],
                            ha='center', va='bottom', fontsize=10, rotation=90
                        )
                    else:
                        plt.text(
                            xp, yp - yerr - distance, reaction_indices[i],
                            ha='center', va='top', fontsize=10, rotation=90
                        )
            plt.xticks([])
            plt.ylabel(
                'Control coefficients on\n'+metric +
                ' (' + obs_name.replace('_', ' ') + ')'
            )
            plt.xlim(-width, n_reaction-1-width)
            # plt.ylim(-1.2,0.6)
            # plt.yticks([-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6])
            plt.legend(loc='lower right', frameon=False)
            plt.savefig(
                'figure/sensitivity/reaction/{}/{}.pdf'.format(
                    metric, obs_name
                ), bbox_inches='tight'
            )
            plt.close()


def heatmap_sensitivity(metric, sensitivity_coefficients, biological_processes,
                        n_reaction, sort_idx, reaction_indices):
        sim = NumericalSimulation()

        # rcParams
        plt.rcParams['font.size'] = 8
        plt.rcParams['font.family'] = 'Arial'
        plt.rcParams['mathtext.fontset'] = 'custom'
        plt.rcParams['mathtext.it'] = 'Arial:italic'
        plt.rcParams['axes.linewidth'] = 1.2
        plt.rcParams['xtick.major.width'] = 1.2
        plt.rcParams['ytick.major.width'] = 1.2

        for k, obs_name in enumerate(observables):
            for l, condition in enumerate(sim.conditions):
                sensitivity_matrix = \
                    sensitivity_coefficients[:, sort_idx[:-1], k, l]
                # Normalize from -1 to 1
                nan_idx = []
                for i in range(sensitivity_matrix.shape[0]):
                    if any(np.isnan(sensitivity_matrix[i, :])):
                        nan_idx.append(i)
                    if np.nanmax(np.abs(sensitivity_matrix[i, :])) == 0.0:
                        sensitivity_matrix[i, :] = np.zeros(
                            sensitivity_matrix.shape[1]
                        )
                    else:
                        sensitivity_matrix[i, :] = (
                            sensitivity_matrix[i, :] / 
                            np.nanmax(
                                np.abs(
                                    sensitivity_matrix[i, :]
                                )
                            )
                        )
                sensitivity_matrix = np.delete(
                    sensitivity_matrix, nan_idx, axis=0
                )
                if sensitivity_matrix.size != 0 and not np.all(sensitivity_matrix == 0.0):
                    sns.clustermap(
                        sensitivity_matrix,
                        center=0,
                        method='ward',
                        cmap='RdBu_r',
                        linewidth=.5,
                        col_cluster=False,
                        figsize=(16, 8),
                        xticklabels=[
                            reaction_indices[i] for i in range(n_reaction-1)
                        ],
                        yticklabels=[],
                        cbar_kws={"ticks": [-1, 0, 1]}
                    )
                    plt.savefig(
                        'figure/sensitivity/reaction/{}/heatmap/{}_{}.pdf'.format(
                            metric, condition, obs_name
                        ), bbox_inches='tight'
                    )
                    plt.close()
