import os
import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

from biomass.model import initial_values
from .sensitivity import calc_sensitivity_coefficients
from .viz import *


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
    sensitivity_coefficients = _load_sc(metric)
    nonzero_idx = get_nonzero_idx()
    if style == 'barplot':
        barplot_sensitivity(
            metric, sensitivity_coefficients, nonzero_idx
        )
    elif style == 'heatmap':
        if len(nonzero_idx) < 2:
            pass
        else:
            heatmap_sensitivity(
                metric, sensitivity_coefficients, nonzero_idx
            )
    else:
        raise ValueError(
            "Available styles are: 'barplot', 'heatmap'"
    )