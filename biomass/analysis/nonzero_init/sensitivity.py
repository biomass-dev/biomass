import os
import sys
import re
import numpy as np

from biomass.dynamics import load_param, get_executable
from biomass.analysis import get_signaling_metric, dlnyi_dlnxj
from .viz import barplot_sensitivity, heatmap_sensitivity

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
        if not nonzero_idx:
            raise ValueError('No nonzero initial conditions')
        
        return nonzero_idx

    def _calc_sensitivity_coefficients(self, metric, nonzero_idx):
        """ Calculating Sensitivity Coefficients

        Parameters
        ----------
        metric: str
            - 'amplitude': The maximum value.
            - 'duration': The time it takes to decline below 10% of its maximum.
            - 'integral': The integral of concentration over the observation time.
        nonzero_idx: list
            for i in nonzero_idx:
                y0[i] != 0.0

        Returns
        -------
        sensitivity_coefficients: numpy array
        
        """

        rate = 1.01  # 1% change
        y0 = self.ival()
        nonzero_idx = self._get_nonzero_idx()
        n_file = get_executable()

        signaling_metric = np.full(
            (len(n_file), len(nonzero_idx)+1, len(self.obs), len(self.sim.conditions)),
            np.nan
        )
        for i, nth_paramset in enumerate(n_file):
            (x, y0) = load_param(nth_paramset, self.sp.update)
            copy_y0 = y0[:]
            for j, idx in enumerate(nonzero_idx):
                y0 = copy_y0[:]
                y0[idx] = copy_y0[idx] * rate
                if self.sim.simulate(x, y0) is None:
                    for k, _ in enumerate(self.obs):
                        for l, _ in enumerate(self.sim.conditions):
                            signaling_metric[i, j, k, l] = \
                                get_signaling_metric(
                                    metric, self.sim.simulations[k, :, l]
                                )
                sys.stdout.write(
                    '\r{:d} / {:d}'.format(
                        i*len(nonzero_idx)+j+1, len(n_file)*len(nonzero_idx)
                    )
                )
            # Signaling metric without perturbation (j=-1)
            y0 = copy_y0[:]
            if self.sim.simulate(x, y0) is None:
                for k, _ in enumerate(self.obs):
                    for l, _ in enumerate(self.sim.conditions):
                        signaling_metric[i, -1, k, l] = get_signaling_metric(
                            metric, self.sim.simulations[k, :, l]
                        )
        sensitivity_coefficients = dlnyi_dlnxj(
            signaling_metric, n_file, nonzero_idx,
            self.obs, self.sim.conditions, rate, metric_idx=-1
        )

        return sensitivity_coefficients

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
            sensitivity_coefficients = \
                self._calc_sensitivity_coefficients(metric, nonzero_idx)
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