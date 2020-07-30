import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from biomass.exec_model import ExecModel
from biomass.analysis import get_signaling_metric, dlnyi_dlnxj

class InitialConditionSensitivity(ExecModel):
    """ sensitivity for species with nonzero initial conditions
    """
    def __init__(self, model):
        super().__init__(model)

    def _get_nonzero_indices(self):
        nonzero_indices = []
        y0 = self.ival()
        for i, val in enumerate(y0):
            if val != 0.0:
                nonzero_indices.append(i)
        if not nonzero_indices:
            raise ValueError('No nonzero initial conditions')
        
        return nonzero_indices

    def _calc_sensitivity_coefficients(self, metric, nonzero_indices):
        """ Calculating Sensitivity Coefficients

        Parameters
        ----------
        metric: str
            - 'amplitude': The maximum value.
            - 'duration': The time it takes to decline below 10% of its maximum.
            - 'integral': The integral of concentration over the observation time.
        nonzero_indices: list
            for i in nonzero_indices:
                y0[i] != 0.0

        Returns
        -------
        sensitivity_coefficients: numpy array
        
        """

        rate = 1.01  # 1% change
        y0 = self.ival()
        nonzero_indices = self._get_nonzero_indices()
        n_file = self.get_executable()

        signaling_metric = np.full(
            (
                len(n_file),
                len(nonzero_indices)+1,
                len(self.obs),
                len(self.sim.conditions)
            ), np.nan
        )
        for i, nth_paramset in enumerate(n_file):
            (x, y0) = self.load_param(nth_paramset)
            y_init = y0[:]
            for j, idx in enumerate(nonzero_indices):
                y0 = y_init[:]
                y0[idx] = y_init[idx] * rate
                if self.sim.simulate(x, y0) is None:
                    for k, _ in enumerate(self.obs):
                        for l, _ in enumerate(self.sim.conditions):
                            signaling_metric[i, j, k, l] = \
                                get_signaling_metric(
                                    metric, self.sim.simulations[k, :, l]
                                )
                sys.stdout.write(
                    '\r{:d} / {:d}'.format(
                        i*len(nonzero_indices)+j+1,
                        len(n_file)*len(nonzero_indices)
                    )
                )
            # Signaling metric without perturbation (j=-1)
            y0 = y_init[:]
            if self.sim.simulate(x, y0) is None:
                for k, _ in enumerate(self.obs):
                    for l, _ in enumerate(self.sim.conditions):
                        signaling_metric[i, -1, k, l] = get_signaling_metric(
                            metric, self.sim.simulations[k, :, l]
                        )
        sensitivity_coefficients = dlnyi_dlnxj(
            signaling_metric, n_file, nonzero_indices,
            self.obs, self.sim.conditions, rate
        )

        return sensitivity_coefficients

    def _load_sc(self, metric, nonzero_indices):
        os.makedirs(
            self.model_path + '/figure/sensitivity/' \
            'initial_condition/{}/heatmap'.format(metric), exist_ok=True
        )
        if not os.path.isfile(
                self.model_path + '/sensitivity_coefficients/' \
                'initial_condition/{}/sc.npy'.format(metric)):
            os.makedirs(
                self.model_path + '/sensitivity_coefficients/' \
                'initial_condition/{}'.format(metric), exist_ok=True
            )
            sensitivity_coefficients = \
                self._calc_sensitivity_coefficients(metric, nonzero_indices)
            np.save(
                self.model_path + '/sensitivity_coefficients/' \
                'initial_condition/{}/sc'.format(metric), sensitivity_coefficients
            )
        else:
            sensitivity_coefficients = np.load(
                self.model_path + '/sensitivity_coefficients/' \
                'initial_condition/{}/sc.npy'.format(metric)
            )
            
        return sensitivity_coefficients

    def _barplot_sensitivity(
            self,
            metric,
            sensitivity_coefficients,
            nonzero_indices
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
            plt.figure(figsize=(9, 5))
            plt.hlines(
                [0], -options['width'], len(nonzero_indices) - options['width'],
                'k', lw=1
            )
            for l, condition in enumerate(self.sim.conditions):
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
                        np.arange(len(nonzero_indices)) + l * options['width'],
                        average, yerr=stdev, ecolor=options['cmap'][l], 
                        capsize=2, width=options['width'],
                        color=options['cmap'][l], align='center', label=condition
                    )
            plt.xticks(
                np.arange(len(nonzero_indices)) + options['width'] / 2,
                [self.viz.convert_species_name(self.species[i]) for i in nonzero_indices],
                rotation=90
            )
            plt.ylabel(
                'Control coefficients on\n' + metric +
                ' (' + obs_name.replace('_', ' ') + ')'
            )
            plt.xlim(-options['width'], len(nonzero_indices)-options['width'])
            plt.legend(loc='upper left', frameon=False)
            plt.savefig(
                self.model_path
                + '/figure/sensitivity/initial_condition/'\
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
                    np.nanmax(np.abs(sensitivity_matrix[i, :])) if normalize else 1
                )

        return np.delete(sensitivity_matrix, nan_idx, axis=0)


    def _heatmap_sensitivity(
            self,
            metric,
            sensitivity_coefficients,
            nonzero_indices
    ):
        # rcParams
        plt.rcParams['font.size'] = 12
        plt.rcParams['font.family'] = 'Arial'
        plt.rcParams['mathtext.fontset'] = 'custom'
        plt.rcParams['mathtext.it'] = 'Arial:italic'

        for k, obs_name in enumerate(self.obs):
            for l, condition in enumerate(self.sim.conditions):
                sensitivity_matrix = self._remove_nan(
                    sensitivity_coefficients[:, :, k, l], normalize=False
                )
                if sensitivity_matrix.shape[0] > 1 and \
                        not np.all(sensitivity_matrix == 0.0):
                    g=sns.clustermap(
                        data=sensitivity_matrix,
                        center=0,
                        robust=True,
                        method='ward',
                        cmap='RdBu_r',
                        linewidth=.5,
                        col_cluster=False,
                        figsize=(sensitivity_matrix.shape[0]*0.5,
                                 sensitivity_matrix.shape[0]*0.35),
                        xticklabels=[
                            self.viz.convert_species_name(self.species[i])
                            for i in nonzero_indices],
                        yticklabels=[],
                        #cbar_kws={"ticks": [-1, 0, 1]}
                    )
                    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)
                    plt.savefig(
                        self.model_path
                        + '/figure/sensitivity/initial_condition/'\
                        '{}/heatmap/{}_{}.pdf'.format(
                            metric, condition, obs_name
                        ), bbox_inches='tight'
                    )
                    plt.close()

    def analyze(self, metric, style):
        nonzero_indices = self._get_nonzero_indices()
        sensitivity_coefficients = self._load_sc(metric, nonzero_indices)
        if style == 'barplot':
            self._barplot_sensitivity(
                metric, sensitivity_coefficients, nonzero_indices
            )
        elif style == 'heatmap':
            if len(nonzero_indices) < 2:
                pass
            else:
                self._heatmap_sensitivity(
                    metric, sensitivity_coefficients, nonzero_indices
                )
        else:
            raise ValueError("Available styles are: 'barplot', 'heatmap'")