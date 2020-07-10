import os
import sys
import re
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from biomass import ExecModel
from biomass.dynamics import load_param, get_executable
from biomass.analysis import get_signaling_metric, dlnyi_dlnxj


class ReactionSensitivity(ExecModel):
    def __init__(self, model):
        super().__init__(model)
        self.model_path = model.__path__[0]
        self.reaction_system = model.set_model
        self.obs = model.observables
        self.sim = model.NumericalSimulation()
        self.viz = model.Visualization()
        self.sp = model.SearchParam()
        self.rxn = model.ReactionNetwork()
    
    def _calc_sensitivity_coefficients(self, metric, n_reaction):
        """ Calculating Sensitivity Coefficients

        Parameters
        ----------
        metric: str
            - 'amplitude':
                The maximum value.
            - 'duration':
                The time it takes to decline below 10% of its maximum.
            - 'integral': 
                The integral of concentration over the observation time.
        n_reaction: int
            len(v) in set_model.py/diffeq

        Returns
        -------
        sensitivity_coefficients: numpy array
        
        """
        rate = 1.01  # 1% change
        n_file = get_executable(self.model_path)
        signaling_metric = np.full(
            (len(n_file), n_reaction, len(self.obs), len(self.sim.conditions)),
            np.nan
        )
        for i, nth_paramset in enumerate(n_file):
            (x, y0) = load_param(self.model_path, nth_paramset, self.sp.update)
            for j in range(n_reaction):
                self.reaction_system.perturbation = [1] * n_reaction
                self.reaction_system.perturbation[j] = rate
                if self.sim.simulate(x, y0) is None:
                    for k, _ in enumerate(self.obs):
                        for l, _ in enumerate(self.sim.conditions):
                            signaling_metric[i, j, k, l] = \
                                get_signaling_metric(
                                    metric, self.sim.simulations[k, :, l]
                                )
                sys.stdout.write(
                    '\r{:d} / {:d}'.format(
                        i*n_reaction+j+1, len(n_file)*n_reaction
                    )
                )
        sensitivity_coefficients = dlnyi_dlnxj(
            signaling_metric, n_file, range(n_reaction),
            self.obs, self.sim.conditions, rate, metric_idx=0
        )

        return sensitivity_coefficients

    def _load_sc(self, metric, n_reaction):
        os.makedirs(
            self.model_path + '/figure/sensitivity/' \
            'reaction/{}/heatmap'.format(metric), exist_ok=True
        )
        if not os.path.isfile(
                self.model_path + '/sensitivity_coefficients/' \
                'reaction/{}/sc.npy'.format(metric)):
            os.makedirs(
                self.model_path + '/sensitivity_coefficients/' \
                'reaction/{}'.format(metric), exist_ok=True
            )
            sensitivity_coefficients = \
                self._calc_sensitivity_coefficients(metric, n_reaction)
            np.save(
                self.model_path + '/sensitivity_coefficients/' \
                'reaction/{}/sc'.format(metric), sensitivity_coefficients
            )
        else:
            sensitivity_coefficients = np.load(
                self.model_path + '/sensitivity_coefficients/' \
                'reaction/{}/sc.npy'.format(metric)
            )
            
        return sensitivity_coefficients

    @staticmethod
    def _draw_vertical_span(biological_processes, width):
        if len(biological_processes) > 1:
            left_end = 0
            for i, proc in enumerate(biological_processes):
                if i % 2 == 0:
                    plt.axvspan(
                        left_end - width,
                        left_end - width + len(proc),
                        facecolor='k', alpha=0.1
                    )
                left_end += len(proc)

    @staticmethod
    def _write_reaction_indices(sort_idx, reaction_indices, average, stdev, width):
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


    def _barplot_sensitivity(
            self,
            metric,
            sensitivity_coefficients,
            biological_processes, 
            n_reaction,
            sort_idx,
            reaction_indices
    ):
        options = self.viz.sensitivity_options

        # rcParams
        plt.rcParams['font.size'] = 15
        plt.rcParams['font.family'] = 'Arial'
        plt.rcParams['mathtext.fontset'] = 'custom'
        plt.rcParams['mathtext.it'] = 'Arial:italic'
        plt.rcParams['axes.linewidth'] = 1.2
        plt.rcParams['xtick.major.width'] = 1.2
        plt.rcParams['ytick.major.width'] = 1.2

        if len(options['cmap']) < len(self.sim.conditions):
            raise ValueError(
                "len(sensitivity_options['cmap']) must be equal to"
                " or greater than len(sim.conditions)."
            )
        for k, obs_name in enumerate(self.obs):
            plt.figure(figsize=(12, 5))
            self._draw_vertical_span(biological_processes, options['width'])
            sensitivity_array = sensitivity_coefficients[:, :, k, :]
            # Remove NaN
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
                if sensitivity_array.shape[0] == 1:
                    stdev = np.zeros(
                        (sensitivity_array.shape[1], sensitivity_array.shape[2])
                    )
                else:
                    stdev = np.std(sensitivity_array, axis=0, ddof=1)
                for l, condition in enumerate(self.sim.conditions):
                    plt.bar(
                        np.arange(n_reaction) + l * options['width'],
                        average[sort_idx, l], yerr=stdev[sort_idx, l],
                        ecolor=options['cmap'][l], capsize=2, width=options['width'],
                        color=options['cmap'][l], align='center', label=condition
                    )
                self._write_reaction_indices(
                    sort_idx, reaction_indices, average, stdev, options['width']
                )
                plt.hlines(
                    [0], -options['width'], n_reaction-1-options['width'],
                    'k', lw=1
                )
                plt.xticks([])
                plt.ylabel(
                    'Control coefficients on\n'+metric +
                    ' (' + obs_name.replace('_', ' ') + ')'
                )
                plt.xlim(-options['width'], n_reaction-1-options['width'])
                # plt.ylim(-1.2,0.6)
                # plt.yticks([-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6])
                plt.legend(loc='lower right', frameon=False)
                plt.savefig(
                    model_path + '/figure/sensitivity/reaction/'\
                    '{}/{}.pdf'.format(
                        metric, obs_name
                    ), bbox_inches='tight'
                )
                plt.close()

    @staticmethod
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
                    np.nanmax(
                        np.abs(sensitivity_matrix[i, :])
                    ) if normalize else 1
                )
        
        return np.delete(sensitivity_matrix, nan_idx, axis=0)

    def _heatmap_sensitivity(
            self,
            metric,
            sensitivity_coefficients,
            biological_processes,
            n_reaction,
            sort_idx,
            reaction_indices
    ):
        # rcParams
        plt.rcParams['font.size'] = 8
        plt.rcParams['font.family'] = 'Arial'
        plt.rcParams['mathtext.fontset'] = 'custom'
        plt.rcParams['mathtext.it'] = 'Arial:italic'
        plt.rcParams['axes.linewidth'] = 1.2
        plt.rcParams['xtick.major.width'] = 1.2
        plt.rcParams['ytick.major.width'] = 1.2

        for k, obs_name in enumerate(self.obs):
            for l, condition in enumerate(self.sim.conditions):
                sensitivity_matrix = self._remove_nan(
                    sensitivity_coefficients[:, sort_idx[:-1], k, l],
                    normalize=False
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
                        xticklabels=[
                            reaction_indices[i] for i in range(n_reaction-1)
                        ],
                        yticklabels=[],
                        #cbar_kws={"ticks": [-1, 0, 1]}
                    )
                    plt.savefig(
                        self.model_path + '/figure/sensitivity/reaction/'\
                        '{}/heatmap/{}_{}.pdf'.format(
                            metric, condition, obs_name
                        ), bbox_inches='tight'
                    )
                    plt.close()

    def analyze(self, metric, style):
        biological_processes = self.rxn.group()
        n_reaction = 1
        for reactions_in_process in biological_processes:
            n_reaction += len(reactions_in_process)
        sort_idx = [0] * n_reaction
        sort_idx[:-1] = np.sum(biological_processes, axis=0)
        reaction_indices = list(
            map(
                lambda x: str(x), sort_idx
            )
        )
        sensitivity_coefficients = self._load_sc(metric, n_reaction)

        if style == 'barplot':
            self._barplot_sensitivity(
                metric, sensitivity_coefficients, biological_processes,
                n_reaction, sort_idx, reaction_indices
            )
        elif style == 'heatmap':
            self._heatmap_sensitivity(
                metric, sensitivity_coefficients, biological_processes,
                n_reaction, sort_idx, reaction_indices
            )
        else:
            raise ValueError("Available styles are: 'barplot', 'heatmap'")
