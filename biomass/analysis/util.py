import sys
from dataclasses import dataclass, field
from math import fabs, isnan, log
from typing import Callable, Dict, List, Union

import numpy as np
from scipy.integrate import simpson


@dataclass
class SignalingMetric(object):
    """
    Signaling metric used for sensitivity analysis.

    Attributes
    ----------
    quantification : Dict[str, Callable[[np.ndarray], Union[int, float]]]
        Functions to quantify signaling metrics.
    """

    quantification: Dict[str, Callable[[np.ndarray], Union[int, float]]] = field(
        default_factory=lambda: dict(
            maximum=np.max,
            minimum=np.min,
            integral=simpson,
        ),
        init=False,
    )


def dlnyi_dlnxj(
    signaling_metric: np.ndarray,
    n_file: List[int],
    perturbed_idx: List[int],
    observables: List[str],
    conditions: List[str],
    rate: float,
) -> np.ndarray:
    """
    Numerical computation of sensitivities using finite difference approximations
    with 1% changes in the reaction rates, parameter values or non-zero initial values.

    Parameters
    ----------
    signaling_metric : numpy array
        Signaling metric

    n_file : list of integers
        Optimized parameter sets in out/

    perturbed_idx : list of integers
        Indices of rate equations or non-zero initial values.

    observables : list of strings
        observables in observable.py

    conditions : list of strings
        Experimental conditions.

    rate : float ~ 1
        1.01 for 1% change.

    Returns
    -------
    sensitivity_coefficients: numpy array

    """
    sensitivity_coefficients = np.empty(
        (len(n_file), len(perturbed_idx), len(observables), len(conditions))
    )
    for i, _ in enumerate(n_file):
        for j, _ in enumerate(perturbed_idx):
            for k, _ in enumerate(observables):
                for l, _ in enumerate(conditions):
                    if isnan(signaling_metric[i, j, k, l]):
                        sensitivity_coefficients[i, j, k, l] = np.nan
                    elif (
                        fabs(signaling_metric[i, -1, k, l]) < sys.float_info.epsilon
                        or fabs(signaling_metric[i, j, k, l] - signaling_metric[i, -1, k, l])
                        < sys.float_info.epsilon
                        or (signaling_metric[i, j, k, l] / signaling_metric[i, -1, k, l]) <= 0
                    ):
                        # 1. Signaling metric before adding perturbation is zero
                        # 2. Absolute change caused by perturbation is too small
                        # 3. Antilogarithm <= 0
                        sensitivity_coefficients[i, j, k, l] = 0.0
                    else:
                        sensitivity_coefficients[i, j, k, l] = log(
                            signaling_metric[i, j, k, l] / signaling_metric[i, -1, k, l]
                        ) / log(rate)

    return sensitivity_coefficients
