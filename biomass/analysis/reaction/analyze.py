import os
import numpy as np

from .sensitivity import calc_sensitivity_coefficients
from .viz import *


class ReactionSensitivity(object):
    def __init__(self, reaction_system, obs, sim, sp, rxn):
        self.reaction_system = reaction_system
        self.obs = obs
        self.sim = sim
        self.sp = sp
        self.rxn = rxn
    
    def _load_sc(self, metric, n_reaction):
        os.makedirs(
            './figure/sensitivity/reaction/{}/heatmap'.format(
                metric
            ), exist_ok=True
        )
        if not os.path.isfile('sc_npy/reaction/{}/sc.npy'.format(metric)):
            os.makedirs(
                './sc_npy/reaction/{}'.format(
                    metric
                ), exist_ok=True
            )
            sensitivity_coefficients = calc_sensitivity_coefficients(
                metric, n_reaction, self.reaction_system, self.obs, self.sim, self.sp
            )
            np.save(
                'sc_npy/reaction/{}/sc'.format(
                    metric
                ), sensitivity_coefficients
            )
        else:
            sensitivity_coefficients = np.load(
                'sc_npy/reaction/{}/sc.npy'.format(
                    metric
                )
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
            raise ValueError(
                "Available styles are: 'barplot', 'heatmap'"
        )
