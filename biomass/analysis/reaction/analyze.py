import os
import numpy as np

from .sensitivity import calc_sensitivity_coefficients
from .reaction import *
from .viz import *


def _load_sc(metric):
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
            metric, num_reaction
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


def analyze(metric, style):
    sensitivity_coefficients = _load_sc(metric)
    reaction_module = get_reaction_module()
    sort_idx = get_sort_idx()
    reaction_number = list(
        map(
            lambda x: str(x), sort_idx
        )
    )
    
    if style == 'barplot':
        barplot_sensitivity(
            metric, sensitivity_coefficients, num_reaction,
            reaction_module, sort_idx, reaction_number
        )
    elif style == 'heatmap':
        heatmap_sensitivity(
            metric, sensitivity_coefficients, num_reaction,
            reaction_module, sort_idx, reaction_number
        )
    else:
        raise ValueError(
            "Available styles are: 'barplot', 'heatmap'"
    )
