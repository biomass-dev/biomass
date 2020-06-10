import os
import sys
import numpy as np

from .sensitivity import calc_sensitivity_coefficients
from .viz import *

class NonZeroInitSensitivity(object):
    def __init__(self, species, ival, obs, sim, sp):
        self.species = species
        self.ival = ival
        self.obs = obs
        self.sim = sim
        self.sp = sp

    def _get_nonzero_idx(self):
        nonzero_idx = []
        y0 = self.ival()
        for i, val in enumerate(y0):
            if val != 0.0:
                nonzero_idx.append(i)
        if len(nonzero_idx) == 0:
            print('No nonzero initial values')
            sys.exit()
        
        return nonzero_idx


    def _load_sc(self, metric, nonzero_idx):
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
                metric, nonzero_idx, self.ival, self.obs, self.sim, self.sp
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


    def analyze(self, metric, style):
        nonzero_idx = self._get_nonzero_idx()
        sensitivity_coefficients = self._load_sc(metric, nonzero_idx)
        if style == 'barplot':
            barplot_sensitivity(
                metric, sensitivity_coefficients, nonzero_idx,
                self.species, self.obs, self.sim
            )
        elif style == 'heatmap':
            if len(nonzero_idx) < 2:
                pass
            else:
                heatmap_sensitivity(
                    metric, sensitivity_coefficients, nonzero_idx,
                    self.species, self.obs, self.sim
                )
        else:
            raise ValueError(
                "Available styles are: 'barplot', 'heatmap'"
        )