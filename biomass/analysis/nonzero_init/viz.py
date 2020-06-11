import os
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns


def barplot_sensitivity(metric, sensitivity_coefficients, nonzero_idx,
                        species, obs, sim):
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
    for k, obs_name in enumerate(obs):
        plt.figure(figsize=(9, 5))
        plt.hlines(
            [0], -width, len(nonzero_idx)-width, 'k', lw=1
        )
        for l, condition in enumerate(sim.conditions):
            sensitivity_matrix = sensitivity_coefficients[:, :, k, l]
            nan_idx = []
            for m in range(sensitivity_matrix.shape[0]):
                if any(np.isnan(sensitivity_matrix[m, :])):
                    nan_idx.append(m)
            sensitivity_matrix = np.delete(
                sensitivity_matrix, nan_idx, axis=0
            )
            if sensitivity_matrix.size != 0:
                average = np.mean(sensitivity_matrix, axis=0)
                if sensitivity_matrix.shape[0] == 1:
                    stdev = np.zeros(sensitivity_matrix.shape[1])
                else:
                    stdev = np.std(sensitivity_matrix, axis=0, ddof=1)
                plt.bar(
                    np.arange(len(nonzero_idx))+l*width, average, yerr=stdev,
                    ecolor=colors[l], capsize=2, width=width, color=colors[l],
                    align='center', label=condition
                )
        plt.xticks(
            np.arange(len(nonzero_idx)) + width/2,
            [species[i] for i in nonzero_idx],
            rotation=30
        )
        plt.ylabel(
            'Control coefficients on\n' + metric +
            ' (' + obs_name.replace('_', ' ') + ')'
        )
        plt.xlim(-width, len(nonzero_idx)-width)
        plt.legend(loc='upper left', frameon=False)
        plt.savefig(
            'figure/sensitivity/nonzero_init/{}/{}.pdf'.format(
                metric, obs_name
            ), bbox_inches='tight'
        )
        plt.close()


def _remove_nan(sensitivity_matrix, normalize):
    nan_idx = []
    for i in range(sensitivity_matrix.shape[0]):
        if any(np.isnan(sensitivity_matrix[i, :])):
            nan_idx.append(i)
        else:
            pass
        if np.nanmax(np.abs(sensitivity_matrix[i, :])) == 0.0:
            sensitivity_matrix[i, :] = np.zeros(
                sensitivity_matrix.shape[1]
            )
        else:
            sensitivity_matrix[i, :] = sensitivity_matrix[i, :] / (
                np.nanmax(np.abs(sensitivity_matrix[i, :])) if normalize else 1
            )

    return np.delete(sensitivity_matrix, nan_idx, axis=0)


def heatmap_sensitivity(metric, sensitivity_coefficients, nonzero_idx,
                        species, obs, sim):
    width = 0.3

    # rcParams
    plt.rcParams['font.size'] = 8
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.it'] = 'Arial:italic'

    for k, obs_name in enumerate(obs):
        for l, condition in enumerate(sim.conditions):
            sensitivity_matrix = _remove_nan(
                sensitivity_coefficients[:, :, k, l], normalize=False
            )
            if sensitivity_matrix.shape[0] > 1 and \
                    not np.all(sensitivity_matrix == 0.0):
                sns.clustermap(
                    data=sensitivity_matrix,
                    center=0,
                    robust=True,
                    method='ward',
                    cmap='RdBu_r',
                    linewidth=.5,
                    col_cluster=False,
                    figsize=(16, 8),
                    xticklabels=[species[i] for i in nonzero_idx],
                    yticklabels=[],
                    #cbar_kws={"ticks": [-1, 0, 1]}
                )
                plt.savefig(
                    'figure/sensitivity/nonzero_init/{}/heatmap/{}_{}.pdf'.format(
                        metric, condition, obs_name
                    ), bbox_inches='tight'
                )
                plt.close()