import os
import sys
import re
import numpy as np

from biomass.dynamics import load_param, get_executable
from biomass.analysis import get_signaling_metric, dlnyi_dlnxj
from .viz import barplot_sensitivity, heatmap_sensitivity


class ReactionSensitivity(object):
    def __init__(self, reaction_system, obs, sim, sp, rxn):
        self.reaction_system = reaction_system
        self.obs = obs
        self.sim = sim
        self.sp = sp
        self.rxn = rxn
    
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
        n_file = get_executable()
        signaling_metric = np.full(
            (len(n_file), n_reaction, len(self.obs), len(self.sim.conditions)),
            np.nan
        )
        for i, nth_paramset in enumerate(n_file):
            (x, y0) = load_param(nth_paramset, self.sp.update)
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
            './figure/sensitivity/' \
            'reaction/{}/heatmap'.format(metric), exist_ok=True
        )
        if not os.path.isfile(
                './sensitivity_coefficients/' \
                'reaction/{}/sc.npy'.format(metric)):
            os.makedirs(
                './sensitivity_coefficients/' \
                'reaction/{}'.format(metric), exist_ok=True
            )
            sensitivity_coefficients = \
                self._calc_sensitivity_coefficients(metric, n_reaction)
            np.save(
                './sensitivity_coefficients/' \
                'reaction/{}/sc'.format(metric), sensitivity_coefficients
            )
        else:
            sensitivity_coefficients = np.load(
                './sensitivity_coefficients/' \
                'reaction/{}/sc.npy'.format(metric)
            )
            
        return sensitivity_coefficients

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
            barplot_sensitivity(
                metric, sensitivity_coefficients, biological_processes,
                n_reaction, sort_idx, reaction_indices, self.obs, self.sim
            )
        elif style == 'heatmap':
            heatmap_sensitivity(
                metric, sensitivity_coefficients, biological_processes,
                n_reaction, sort_idx, reaction_indices, self.obs, self.sim
            )
        else:
            raise ValueError("Available styles are: 'barplot', 'heatmap'")
