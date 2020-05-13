import os
import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

from biomass.model import V, initial_values
from biomass.observable import observables, NumericalSimulation
from .sensitivity import calc_sensitivity_coefficients


def get_nonzero_idx():
    nonzero_idx = []
    y0 = initial_values()
    for i, val in enumerate(y0):
        if val != 0.0:
            nonzero_idx.append(i)
    if len(nonzero_idx) == 0:
        print('No nonzero initial values')
        sys.exit()
    
    return nonzero_idx


def _load_sc(metric):
    nonzero_idx = get_nonzero_idx()
    os.makedirs(
        './figure/sensitivity/nonzero_init/{}/heatmap'.format(
            metric
        ), exist_ok=True
    )
    if not os.path.isfile('sc_npy/nonzero_init/{}/sc.npy'.format(metric)):
        os.makedirs(
            './sc_npy/nonzero_init/{}'.format(
                metric
            ), exist_ok=True
        )
        sensitivity_coefficients = calc_sensitivity_coefficients(
            metric, nonzero_idx
        )
        np.save(
            'sc_npy/nonzero_init/{}/sc'.format(
                metric
            ), sensitivity_coefficients
        )
    else:
        sensitivity_coefficients = np.load(
            'sc_npy/nonzero_init/{}/sc.npy'.format(
                metric
            )
        )
        
    return sensitivity_coefficients


def analyze(metric, style):
    sim = NumericalSimulation()
    width = 0.3

    nonzero_idx = get_nonzero_idx()
    if style == 'barplot':
        sensitivity_coefficients = _load_sc(metric)
        # rcParams
        plt.rcParams['font.size'] = 15
        plt.rcParams['font.family'] = 'Arial'
        plt.rcParams['mathtext.fontset'] = 'custom'
        plt.rcParams['mathtext.it'] = 'Arial:italic'
        plt.rcParams['axes.linewidth'] = 1.2
        plt.rcParams['xtick.major.width'] = 1.2
        plt.rcParams['ytick.major.width'] = 1.2

        colors = ['mediumblue', 'red']
        for k, obs_name in enumerate(observables):
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
                    stdev = np.std(sensitivity_matrix, axis=0, ddof=1)
                    plt.bar(
                        np.arange(len(nonzero_idx))+l*width, average, yerr=stdev,
                        ecolor=colors[l], capsize=2, width=width, color=colors[l],
                        align='center', label=condition
                    )
            plt.xticks(
                np.arange(len(nonzero_idx)) + width/2,
                [V.var_names[i] for i in nonzero_idx],
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
    elif style == 'heatmap':
        if len(nonzero_idx) < 2:
            pass
        else:
            sensitivity_coefficients = _load_sc(metric)
            # rcParams
            plt.rcParams['font.size'] = 8
            plt.rcParams['font.family'] = 'Arial'
            plt.rcParams['mathtext.fontset'] = 'custom'
            plt.rcParams['mathtext.it'] = 'Arial:italic'

            for k, obs_name in enumerate(observables):
                for l, condition in enumerate(sim.conditions):
                    sensitivity_matrix = sensitivity_coefficients[:, :, k, l]
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
                            xticklabels=[V.var_names[i] for i in nonzero_idx],
                            yticklabels=[],
                            cbar_kws={"ticks": [-1, 0, 1]}
                        )
                        plt.savefig(
                            'figure/sensitivity/nonzero_init/{}/heatmap/{}_{}.pdf'.format(
                                metric, condition, obs_name
                            ), bbox_inches='tight'
                        )
                        plt.close()
    else:
        raise ValueError(
            "Available styles are: 'barplot', 'heatmap'"
    )